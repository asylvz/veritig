#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <set>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include "project.h"
#include "edlib.h"


void Project::run_mapping(parameters& params, const std::string& svtig_fasta, const std::string& paf_out)
{
	std::string threads_str = std::to_string(params.threads);
	std::cerr << "  Mapping " << svtig_fasta << " to reference...";
	std::string cmd = "minimap2 -cx " + shell_quote(effective_preset(params)) + " --cs -r2k -t " + threads_str + " "
		+ shell_quote(params.reference_path) + " " + shell_quote(svtig_fasta)
		+ " --secondary=no > " + shell_quote(paf_out) + " 2>/dev/null";
	run_command(cmd, paf_out);
	std::cerr << " done\n";
}


std::map<std::string, std::string> Project::load_fasta(const std::string& fasta_path)
{
	std::map<std::string, std::string> out;
	std::ifstream fp(fasta_path);
	if (!fp.good())
	{
		std::cerr << "[veritig] Cannot open FASTA: " << fasta_path << "\n";
		return out;
	}

	std::string line;
	std::string cur_name;
	std::string cur_seq;
	while (getline(fp, line))
	{
		if (line.empty()) continue;
		if (line[0] == '>')
		{
			if (!cur_name.empty())
				out[cur_name] = cur_seq;
			std::string name = line.substr(1);
			size_t sp = name.find_first_of(" \t");
			if (sp != std::string::npos) name = name.substr(0, sp);
			cur_name = name;
			cur_seq.clear();
		}
		else
		{
			for (char& c : line) c = std::toupper(c);
			cur_seq += line;
		}
	}
	if (!cur_name.empty())
		out[cur_name] = cur_seq;

	return out;
}


// Parse the cs:Z: tag and emit >=min_svlen INS/DEL events.
// cs short-form grammar:
//   :N     -> N identical bases (advances ref and query by N)
//   *xy    -> substitution of x to y (advances ref and query by 1)
//   +seq   -> insertion in query (advances query by len(seq))
//   -seq   -> deletion from ref (advances ref by len(seq))
static void extract_svs_from_cs(const std::string& cs, int tgt_start,
                                const std::string& chrom, const std::string& svtig_name,
                                int mapq, int haplo, int min_svlen,
                                std::vector<ProjectedSV>& out)
{
	int ref_offset = 0;
	size_t i = 0;
	while (i < cs.size())
	{
		char op = cs[i++];
		if (op == ':')
		{
			int len = 0;
			while (i < cs.size() && std::isdigit(cs[i]))
				len = len * 10 + (cs[i++] - '0');
			ref_offset += len;
		}
		else if (op == '*')
		{
			// substitution: two base letters
			if (i + 1 < cs.size()) i += 2;
			ref_offset += 1;
		}
		else if (op == '+')
		{
			// insertion in query (ref unchanged)
			std::string seq;
			while (i < cs.size() && std::isalpha(cs[i]))
				seq += std::toupper(cs[i++]);
			int len = (int)seq.size();
			if (len >= min_svlen)
			{
				ProjectedSV sv;
				sv.chrom = chrom;
				sv.pos = tgt_start + ref_offset;   // 1-based anchor = last ref base consumed
				if (sv.pos < 1) sv.pos = 1;
				sv.svtype = 'I';
				sv.svlen = len;
				sv.ref_base = "N";
				sv.alt_seq = "N" + seq;
				sv.svtig_name = svtig_name;
				sv.mapq = mapq;
				sv.haplo = haplo;
				out.push_back(sv);
			}
		}
		else if (op == '-')
		{
			// deletion from ref (query unchanged)
			std::string seq;
			while (i < cs.size() && std::isalpha(cs[i]))
				seq += std::toupper(cs[i++]);
			int len = (int)seq.size();
			if (len >= min_svlen)
			{
				ProjectedSV sv;
				sv.chrom = chrom;
				sv.pos = tgt_start + ref_offset;
				if (sv.pos < 1) sv.pos = 1;
				sv.svtype = 'D';
				sv.svlen = len;
				sv.ref_base = "N" + seq;
				sv.alt_seq = "N";
				sv.svtig_name = svtig_name;
				sv.mapq = mapq;
				sv.haplo = haplo;
				out.push_back(sv);
			}
			ref_offset += len;
		}
		else
		{
			// unknown op; skip until next op char
			while (i < cs.size() && cs[i] != ':' && cs[i] != '*' && cs[i] != '+' && cs[i] != '-')
				i++;
		}
	}
}


std::vector<ProjectedSV> Project::parse_paf(const std::string& paf_path, int haplo, parameters& params,
                                            const std::map<std::string, std::string>& svtig_seqs,
                                            std::vector<std::string>* unmapped_out)
{
	std::vector<ProjectedSV> svs;
	std::ifstream fp(paf_path);
	if (!fp.good())
	{
		std::cerr << "[veritig] Cannot open PAF: " << paf_path << "\n";
		return svs;
	}

	std::string line;
	int n_alignments = 0;
	std::map<std::string, std::vector<AlignmentRecord>> svtig_alns;
	std::set<std::string> aligned_svtig_names;

	while (getline(fp, line))
	{
		if (line.empty()) continue;

		std::vector<std::string> tokens;
		std::stringstream ss(line);
		std::string tok;
		while (getline(ss, tok, '\t'))
			tokens.push_back(tok);

		if (tokens.size() < 12) continue;

		std::string svtig_name = tokens[0];
		int q_len = std::atoi(tokens[1].c_str());
		int q_start = std::atoi(tokens[2].c_str());
		int q_end = std::atoi(tokens[3].c_str());
		char strand = tokens[4][0];
		std::string chrom = tokens[5];
		int tgt_start_0b = std::atoi(tokens[7].c_str());
		int tgt_end_0b = std::atoi(tokens[8].c_str());
		int mapq = std::atoi(tokens[11].c_str());

		if (q_len <= 0) continue;   // skip malformed PAF lines

		aligned_svtig_names.insert(svtig_name);

		if (mapq < params.min_mapq) continue;

		// Find cs:Z: tag; also check tp:A: for primary/supplementary
		bool is_primary = true;
		std::string cs_tag;
		for (size_t t = 12; t < tokens.size(); t++)
		{
			if (tokens[t].rfind("tp:A:", 0) == 0)
				is_primary = (tokens[t][5] == 'P');
			else if (tokens[t].rfind("cs:Z:", 0) == 0)
				cs_tag = tokens[t].substr(5);
		}

		// Track this alignment for split-aware analysis
		AlignmentRecord ar;
		ar.q_start = q_start; ar.q_end = q_end; ar.q_len = q_len;
		ar.t_chrom = chrom; ar.t_start = tgt_start_0b; ar.t_end = tgt_end_0b;
		ar.strand = strand; ar.mapq = mapq; ar.is_primary = is_primary;
		svtig_alns[svtig_name].push_back(ar);

		if (cs_tag.empty()) continue;

		n_alignments++;
		extract_svs_from_cs(cs_tag, tgt_start_0b, chrom, svtig_name, mapq, haplo, params.min_svlen, svs);
	}

	std::cerr << "  " << paf_path << ": " << n_alignments << " alignments, "
		<< svs.size() << " SVs >=" << params.min_svlen << "bp from CIGAR\n";

	// Split-alignment / soft-clip INS detection from query coverage gaps
	std::vector<ProjectedSV> split_svs = extract_split_ins(svtig_alns, svtig_seqs, haplo, params);
	if (!split_svs.empty())
	{
		std::cerr << "  Split/soft-clip analysis: " << split_svs.size()
			<< " additional INS candidates (>=" << params.min_split_svlen << "bp)\n";
		svs.insert(svs.end(), split_svs.begin(), split_svs.end());
	}

	// Identify unmapped svtigs (in FASTA but never aligned)
	if (unmapped_out)
	{
		for (const auto& kv : svtig_seqs)
			if (aligned_svtig_names.find(kv.first) == aligned_svtig_names.end())
				unmapped_out->push_back(kv.first);
	}

	// Within-svtig clustering: merge fragmented events from same svtig (e.g. tandem-repeat artefacts)
	svs = cluster_within_svtig(svs, params.cluster_window, params.cluster_size_ratio);

	// Tolerance-based dedup across svtigs (within same haplotype):
	// Two events are merged if same chrom + svtype, |pos diff| <= pos_tol, |size diff| <= size_tol * larger.
	// When merged, kept record's svtig_name accumulates the svtig list (comma-separated) and support_count increments.
	std::sort(svs.begin(), svs.end(), [](const ProjectedSV& a, const ProjectedSV& b) {
		if (a.chrom != b.chrom) return a.chrom < b.chrom;
		if (a.svtype != b.svtype) return a.svtype < b.svtype;
		if (a.pos != b.pos) return a.pos < b.pos;
		return a.svlen < b.svlen;
	});

	std::vector<ProjectedSV> dedup;
	int dedup_dropped = 0;
	int pos_tol = params.dedup_pos_tol;
	double size_tol = params.dedup_size_tol;

	for (auto& sv : svs)
	{
		bool merged = false;
		// Compare against all recent dedup records on the same chrom+svtype within pos_tol.
		// Since input is sorted by (chrom, svtype, pos, svlen), we walk back through tail.
		for (int i = (int)dedup.size() - 1; i >= 0; i--)
		{
			ProjectedSV& p = dedup[i];
			if (p.chrom != sv.chrom || p.svtype != sv.svtype) break;
			if (sv.pos - p.pos > pos_tol) break;  // sorted: further back records are even smaller pos
			// size tolerance: |size_a - size_b| / max(a, b) <= size_tol
			int size_diff = std::abs(p.svlen - sv.svlen);
			int size_max = std::max(p.svlen, sv.svlen);
			if (size_max > 0 && (double)size_diff / size_max <= size_tol)
			{
				// Same SV, different svtigs -> merge svtig list using pipe separator
				// (VCF Number=1 String fields cannot contain commas — they would be parsed as list)
				const size_t SVTIG_LIST_MAX = 200;
				if (p.svtig_name.find(sv.svtig_name) == std::string::npos)
				{
					p.support_count++;
					if (p.svtig_name.size() < SVTIG_LIST_MAX)
						p.svtig_name += "|" + sv.svtig_name;
					else if (p.svtig_name.size() < SVTIG_LIST_MAX + 10
					         && p.svtig_name.find("...") == std::string::npos)
						p.svtig_name += "|...";
				}
				if (sv.mapq > p.mapq) p.mapq = sv.mapq;
				merged = true;
				dedup_dropped++;
				break;
			}
		}
		if (!merged)
			dedup.push_back(sv);
	}

	if (dedup_dropped > 0)
		std::cerr << "  Tolerance dedup: dropped " << dedup_dropped
			<< " duplicates (pos_tol=" << pos_tol << ", size_tol=" << size_tol << ")\n";

	return dedup;
}


// Extract INS from unaligned query gaps (soft-clips and inter-alignment gaps).
// Each svtig's alignments are sorted by query position; a gap >= min_split_svlen is
// emitted as an INS anchored at the adjacent ref position. Internal gaps must be same
// chrom/strand within max_split_distance; strand flips (INV) and cross-chrom (TRA) are
// left for other passes.
std::vector<ProjectedSV> Project::extract_split_ins(
    const std::map<std::string, std::vector<AlignmentRecord>>& svtig_alns,
    const std::map<std::string, std::string>& svtig_seqs,
    int haplo, parameters& params)
{
	std::vector<ProjectedSV> out;
	int min_size = params.min_split_svlen;
	int max_dist = params.max_split_distance;
	int n_inv_skipped = 0;
	int n_tra_skipped = 0;
	int n_distant_skipped = 0;

	for (const auto& kv : svtig_alns)
	{
		const std::string& svtig_name = kv.first;
		auto alns = kv.second;
		if (alns.empty()) continue;

		// Sort by query start
		std::sort(alns.begin(), alns.end(), [](const AlignmentRecord& a, const AlignmentRecord& b) {
			return a.q_start < b.q_start;
		});

		int q_len = alns[0].q_len;
		auto seq_it = svtig_seqs.find(svtig_name);
		bool have_seq = (seq_it != svtig_seqs.end() && (int)seq_it->second.size() == q_len);

		// Helper to emit one INS candidate for a query gap [q_lo, q_hi) anchored at ref_pos (1-based)
		auto emit_gap = [&](int q_lo, int q_hi, const std::string& chrom, int ref_pos_1b, int mapq) {
			int gap_size = q_hi - q_lo;
			if (gap_size < min_size) return;
			ProjectedSV sv;
			sv.chrom = chrom;
			sv.pos = (ref_pos_1b < 1) ? 1 : ref_pos_1b;
			sv.svtype = 'I';
			sv.svlen = gap_size;
			sv.ref_base = "N";
			if (have_seq)
				sv.alt_seq = "N" + seq_it->second.substr(q_lo, gap_size);
			else
				sv.alt_seq = "N" + std::string(gap_size, 'N');
			sv.svtig_name = svtig_name;
			sv.mapq = mapq;
			sv.haplo = haplo;
			out.push_back(sv);
		};

		// Soft-clip on a single alignment is intentionally NOT emitted as INS:
		// the unaligned svtig portion has no anchoring evidence on the opposite side,
		// so the true INS position on the reference cannot be localised. We only
		// emit INS from internal gaps between two alignments (split-alignment case),
		// where both flanks anchor the gap to a specific ref locus.

		// First pass: detect 3-segment INV patterns (+, -, +) or (-, +, -).
		// In a true inversion, minimap2 produces three adjacent alignments where the middle
		// one is on the opposite strand. The INV reference span is bounded by the outer two
		// alignments, NOT the reverse alignment's full span (which often extends past the
		// real inversion when query content matches both strands locally).
		for (size_t i = 0; i + 2 < alns.size(); i++)
		{
			const AlignmentRecord& A = alns[i];
			const AlignmentRecord& B = alns[i + 1];
			const AlignmentRecord& C = alns[i + 2];

			if (A.t_chrom != B.t_chrom || B.t_chrom != C.t_chrom) continue;
			// require A,C same strand and B opposite
			if (A.strand != C.strand || A.strand == B.strand) continue;

			// INV bounded by (A.t_end .. C.t_start) when A,C are forward (the typical case)
			// or by the equivalent when A,C are reverse.
			int inv_start, inv_end;
			if (A.strand == '+')
			{
				inv_start = A.t_end;     // 0-based last aligned ref base of A; 1-based start of inversion = A.t_end + 1; we use A.t_end as 1-based anchor
				inv_end   = C.t_start;   // 1-based last inverted ref base = C.t_start
			}
			else
			{
				inv_start = C.t_end;
				inv_end   = A.t_start;
			}
			if (inv_end <= inv_start) continue;
			int inv_svlen = inv_end - inv_start;
			// INV uses --min-svlen (default 50), NOT --min-split-svlen which is for soft-clip
			// artefacts. INV from 3-segment pattern has unambiguous breakpoints.
			if (inv_svlen < params.min_svlen) continue;

			ProjectedSV sv;
			sv.chrom = A.t_chrom;
			sv.pos = inv_start + 1;     // 1-based VCF POS
			if (sv.pos < 1) sv.pos = 1;
			sv.svtype = 'V';
			sv.svlen = inv_svlen;
			sv.ref_base = "N";
			sv.alt_seq = "<INV>";
			sv.svtig_name = svtig_name;
			sv.mapq = std::max({A.mapq, B.mapq, C.mapq});
			sv.haplo = haplo;
			out.push_back(sv);
		}


		// Second pass: pairwise transitions for INS (same strand, query gap)
		for (size_t i = 1; i < alns.size(); i++)
		{
			const AlignmentRecord& L = alns[i - 1];
			const AlignmentRecord& R = alns[i];

			if (L.t_chrom != R.t_chrom)
			{
				n_tra_skipped++;
				continue;
			}

			// Strand-flip in pairwise: INV emitted only from 3-segment pattern above.
			// 2-segment INVs (without flanking forward alignment on the other side) have
			// ambiguous breakpoints that produce many FPs in benchmarks.
			if (L.strand != R.strand)
			{
				n_inv_skipped++;
				continue;
			}

			// Same strand: check query gap for INS detection
			int q_gap = R.q_start - L.q_end;
			if (q_gap < min_size) continue;

			int t_distance = std::abs(R.t_start - L.t_end);
			if (t_distance > max_dist)
			{
				n_distant_skipped++;
				continue;
			}

			// Anchor at L.t_end (last aligned ref base of left alignment, 1-based)
			emit_gap(L.q_end, R.q_start, L.t_chrom, L.t_end, std::max(L.mapq, R.mapq));
		}
	}

	if (n_inv_skipped + n_tra_skipped + n_distant_skipped > 0)
		std::cerr << "  Split skipped: " << n_inv_skipped << " strand-flip (INV candidate), "
			<< n_tra_skipped << " cross-chrom (TRA), " << n_distant_skipped << " distant (>max-split-distance)\n";

	return out;
}


// Merge fragmented events within the same svtig: same SVTYPE, gap <= window, similar size.
// Conservative: only events with size_min/size_max >= ratio are merged.
// Skipped if window <= 0.
std::vector<ProjectedSV> Project::cluster_within_svtig(std::vector<ProjectedSV>& events, int window, double size_ratio)
{
	if (window <= 0 || events.empty())
		return events;

	// Group by (svtig_name, chrom, svtype). Same type only — DEL+INS pairs stay separate.
	std::map<std::tuple<std::string, std::string, char>, std::vector<ProjectedSV>> groups;
	for (auto& ev : events)
		groups[{ev.svtig_name, ev.chrom, ev.svtype}].push_back(ev);

	std::vector<ProjectedSV> result;
	int merge_events = 0;
	int merge_clusters = 0;

	for (auto& kv : groups)
	{
		auto& group = kv.second;
		std::sort(group.begin(), group.end(), [](const ProjectedSV& a, const ProjectedSV& b) {
			return a.pos < b.pos;
		});

		ProjectedSV current = group[0];
		bool in_cluster = false;

		for (size_t i = 1; i < group.size(); i++)
		{
			const ProjectedSV& nxt = group[i];
			// Reference-space gap between events:
			// DEL consumes ref bases (span = svlen), INS consumes 0 ref bases
			int end_curr = current.pos + (current.svtype == 'D' ? current.svlen : 0);
			int gap = nxt.pos - end_curr;
			double smin = std::min((double)current.svlen, (double)nxt.svlen);
			double smax = std::max((double)current.svlen, (double)nxt.svlen);
			bool size_ok = (smax > 0) && (smin / smax >= size_ratio);

			if (gap <= window && size_ok)
			{
				// Merge nxt into current
				int end_curr = current.pos + current.svlen;
				int end_nxt = nxt.pos + nxt.svlen;
				int new_end = std::max(end_curr, end_nxt);

				if (current.svtype == 'D')
				{
					// DEL: span = leftmost_pos to rightmost_end
					current.svlen = new_end - current.pos;
					// REF sequence content unknown for the gap region; placeholder: anchor + N's
					current.ref_base = "N" + std::string(current.svlen, 'N');
				}
				else
				{
					// INS: combine inserted sequences (drop anchor of nxt)
					current.svlen += nxt.svlen;
					if (nxt.alt_seq.size() > 1)
						current.alt_seq += nxt.alt_seq.substr(1);
				}
				current.mapq = std::max(current.mapq, nxt.mapq);
				in_cluster = true;
				merge_events++;
			}
			else
			{
				if (in_cluster) merge_clusters++;
				result.push_back(current);
				current = nxt;
				in_cluster = false;
			}
		}
		if (in_cluster) merge_clusters++;
		result.push_back(current);
	}

	if (merge_events > 0)
		std::cerr << "  Clustering: merged " << merge_events << " fragmented events into "
			<< merge_clusters << " composite SVs (window=" << window
			<< ", size_ratio=" << size_ratio << ")\n";

	return result;
}


// Relabel an INS as DUP when its inserted sequence S (length L) closely matches a
// nearby reference copy: edlib HW-aligns S to a +/- factor*L window and relabels if
// identity (1 - editDistance/L) >= dup_similarity. Post-hoc; only the svtype changes,
// not pos/svlen/alt_seq.
void Project::classify_dups(std::vector<ProjectedSV>& events, parameters& params)
{
	if (ref_seqs.empty())
	{
		std::cerr << "  DUP detection: reference not loaded, skipping\n";
		return;
	}

	int n_checked = 0;
	int n_dup = 0;
	double sim_threshold = params.dup_similarity;
	double window_factor = params.dup_window_factor;
	int min_size = params.dup_min_svlen;

	for (auto& sv : events)
	{
		if (sv.svtype != 'I') continue;
		if (sv.svlen < min_size) continue;
		// alt_seq starts with anchor 'N'; the actual inserted sequence is alt_seq[1..]
		if ((int)sv.alt_seq.size() < sv.svlen + 1) continue;

		auto ref_it = ref_seqs.find(sv.chrom);
		if (ref_it == ref_seqs.end()) continue;
		const std::string& ref = ref_it->second;

		int half_window = (int)(window_factor * sv.svlen);
		int win_start = std::max(0, sv.pos - 1 - half_window);
		int win_end = std::min((int)ref.size(), sv.pos + half_window);
		if (win_end - win_start < sv.svlen) continue;

		std::string query = sv.alt_seq.substr(1);   // strip anchor
		const char* target = ref.data() + win_start;
		int target_len = win_end - win_start;

		// HW = infix mode: find best alignment of query inside target
		EdlibAlignResult result = edlibAlign(
			query.data(), (int)query.size(),
			target, target_len,
			edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));

		n_checked++;
		if (result.status == EDLIB_STATUS_OK && result.editDistance >= 0)
		{
			double identity = 1.0 - (double)result.editDistance / (double)query.size();
			if (identity >= sim_threshold)
			{
				sv.svtype = 'U';   // reclassify as DUP
				n_dup++;
			}
		}
		edlibFreeAlignResult(result);
	}

	if (n_checked > 0)
		std::cerr << "  DUP detection: " << n_dup << "/" << n_checked
			<< " INS events reclassified as DUP (similarity >= " << sim_threshold << ")\n";
}


// Flag an INS/DUP as ambiguous when its inserted sequence sits in a repeat: such calls
// are often minimap2 artefacts in tandem-repeat regions. Count k-mer abundances of S
// against a +/- max(window_factor*|S|, ambig_min_window) reference window; if the median
// abundance >= ambig_threshold, set sv.ambig (written as FILTER=Ambig).
void Project::flag_ambiguous_repeats(std::vector<ProjectedSV>& events, parameters& params)
{
	if (ref_seqs.empty()) return;

	int k = params.ambig_kmer_size;
	int threshold = params.ambig_threshold;
	int min_half_window = params.ambig_min_window;
	double window_factor = params.dup_window_factor;

	int n_checked = 0;
	int n_flagged = 0;

	for (auto& sv : events)
	{
		// Only check insertion-like events (INS, DUP). DEL events do not insert sequence.
		if (sv.svtype != 'I' && sv.svtype != 'U') continue;
		if ((int)sv.alt_seq.size() < sv.svlen + 1) continue;

		auto ref_it = ref_seqs.find(sv.chrom);
		if (ref_it == ref_seqs.end()) continue;
		const std::string& ref = ref_it->second;

		int half_window = std::max((int)(window_factor * sv.svlen), min_half_window);
		int win_start = std::max(0, sv.pos - 1 - half_window);
		int win_end = std::min((int)ref.size(), sv.pos + half_window);
		if (win_end - win_start < k) continue;

		// Count k-mer abundance in the ref window
		std::map<std::string, int> ref_kmers;
		for (int i = win_start; i + k <= win_end; i++)
			ref_kmers[ref.substr(i, k)]++;

		// Look up abundance for each k-mer in the inserted sequence
		std::string inserted = sv.alt_seq.substr(1);
		std::vector<int> abundances;
		for (int i = 0; i + k <= (int)inserted.size(); i++)
		{
			auto it = ref_kmers.find(inserted.substr(i, k));
			abundances.push_back(it == ref_kmers.end() ? 0 : it->second);
		}
		if (abundances.empty()) continue;

		std::sort(abundances.begin(), abundances.end());
		int median = abundances[abundances.size() / 2];

		n_checked++;
		if (median >= threshold)
		{
			sv.ambig = true;
			n_flagged++;
		}
	}

	if (n_checked > 0)
		std::cerr << "  Ambig flag: " << n_flagged << "/" << n_checked
			<< " INS/DUP events flagged as repeat-ambiguous (k=" << k
			<< ", threshold=" << threshold << ", min_window=" << min_half_window << ")\n";
}


std::vector<ProjectedSV> Project::merge_haplotypes(std::vector<ProjectedSV>& h1_svs, std::vector<ProjectedSV>& h2_svs)
{
	return merge_haplotypes(h1_svs, h2_svs, 50, 0.1);
}

std::vector<ProjectedSV> Project::merge_haplotypes(std::vector<ProjectedSV>& h1_svs, std::vector<ProjectedSV>& h2_svs,
		int pos_tol, double size_tol)
{
	std::vector<ProjectedSV> merged;

	// Group H1 indices by (chrom, svtype) and sort within each group by pos for binary search
	std::map<std::pair<std::string,char>, std::vector<size_t>> h1_groups;
	for (size_t i = 0; i < h1_svs.size(); i++)
	{
		auto key = std::make_pair(h1_svs[i].chrom, h1_svs[i].svtype);
		h1_groups[key].push_back(i);
	}
	for (auto& kv : h1_groups)
	{
		std::sort(kv.second.begin(), kv.second.end(),
			[&h1_svs](size_t a, size_t b) { return h1_svs[a].pos < h1_svs[b].pos; });
	}

	std::set<size_t> h1_matched;

	// For each H2 sv, search candidates in same chrom+svtype within pos_tol; match if size also within tolerance
	for (auto& sv2 : h2_svs)
	{
		auto key = std::make_pair(sv2.chrom, sv2.svtype);
		auto git = h1_groups.find(key);
		size_t best_idx = SIZE_MAX;
		int best_pos_diff = pos_tol + 1;

		if (git != h1_groups.end())
		{
			const auto& sorted_idx = git->second;
			// Binary search for first index with pos >= sv2.pos - pos_tol
			int target_lo = sv2.pos - pos_tol;
			int target_hi = sv2.pos + pos_tol;
			auto lower = std::lower_bound(sorted_idx.begin(), sorted_idx.end(), target_lo,
				[&h1_svs](size_t i, int v) { return h1_svs[i].pos < v; });

			for (auto it = lower; it != sorted_idx.end() && h1_svs[*it].pos <= target_hi; ++it)
			{
				if (h1_matched.count(*it)) continue;
				const auto& sv1 = h1_svs[*it];
				int pos_diff = std::abs(sv1.pos - sv2.pos);
				int max_len = std::max(sv1.svlen, sv2.svlen);
				int min_len = std::min(sv1.svlen, sv2.svlen);
				double size_diff = max_len > 0 ? (double)(max_len - min_len) / (double)max_len : 0.0;
				if (size_diff > size_tol) continue;
				if (pos_diff < best_pos_diff)
				{
					best_pos_diff = pos_diff;
					best_idx = *it;
				}
			}
		}

		if (best_idx != SIZE_MAX)
		{
			// homozygous: take H1 record (sequences usually nearly identical; pick one)
			ProjectedSV sv = h1_svs[best_idx];
			sv.haplo = 3;
			merged.push_back(sv);
			h1_matched.insert(best_idx);
		}
		else
		{
			merged.push_back(sv2);  // H2-only
		}
	}

	// Add remaining H1-only
	for (size_t i = 0; i < h1_svs.size(); i++)
	{
		if (h1_matched.count(i) == 0)
			merged.push_back(h1_svs[i]);
	}

	// Sort final list by chrom, pos
	std::sort(merged.begin(), merged.end(), [](const ProjectedSV& a, const ProjectedSV& b) {
		if (a.chrom != b.chrom) return a.chrom < b.chrom;
		return a.pos < b.pos;
	});

	return merged;
}


// Merge untagged (unphased) SVarp svtigs into an already-merged H1+H2 set.
//   - If an untagged SV matches an existing call (chrom, svtype, ±pos_tol, size within tol):
//       increment support_count, append svtig name -- the existing call already covers it.
//   - If an untagged SV is novel (no match): emit as a homozygous call (haplo=3) with
//       untagged_only flag set (FILTER=Untagged in VCF). This is a default assumption;
//       sequence-level verification against the assembly should refine the genotype.
void Project::merge_untagged(std::vector<ProjectedSV>& merged, std::vector<ProjectedSV>& untagged,
	int pos_tol, double size_tol, int& n_matched, int& n_novel)
{
	// Index merged calls by (chrom, svtype) sorted by pos for fast tolerance lookup
	std::map<std::pair<std::string,char>, std::vector<size_t>> idx;
	for (size_t i = 0; i < merged.size(); i++)
	{
		auto key = std::make_pair(merged[i].chrom, merged[i].svtype);
		idx[key].push_back(i);
	}
	for (auto& kv : idx)
	{
		std::sort(kv.second.begin(), kv.second.end(),
			[&merged](size_t a, size_t b) { return merged[a].pos < merged[b].pos; });
	}

	n_matched = 0;
	n_novel = 0;
	for (auto& uv : untagged)
	{
		auto key = std::make_pair(uv.chrom, uv.svtype);
		auto it = idx.find(key);
		size_t best_idx = SIZE_MAX;
		int best_pos_diff = pos_tol + 1;

		if (it != idx.end())
		{
			int target_lo = uv.pos - pos_tol;
			int target_hi = uv.pos + pos_tol;
			auto lower = std::lower_bound(it->second.begin(), it->second.end(), target_lo,
				[&merged](size_t i, int v) { return merged[i].pos < v; });
			for (auto k = lower; k != it->second.end() && merged[*k].pos <= target_hi; ++k)
			{
				const auto& mv = merged[*k];
				int max_len = std::max(mv.svlen, uv.svlen);
				int min_len = std::min(mv.svlen, uv.svlen);
				double size_diff = max_len > 0 ? (double)(max_len - min_len) / (double)max_len : 0.0;
				if (size_diff > size_tol) continue;
				int pos_diff = std::abs(mv.pos - uv.pos);
				if (pos_diff < best_pos_diff)
				{
					best_pos_diff = pos_diff;
					best_idx = *k;
				}
			}
		}

		if (best_idx != SIZE_MAX)
		{
			merged[best_idx].support_count += uv.support_count;
			if (!uv.svtig_name.empty())
			{
				if (!merged[best_idx].svtig_name.empty())
					merged[best_idx].svtig_name += "|";
				merged[best_idx].svtig_name += uv.svtig_name;
			}
			n_matched++;
		}
		else
		{
			ProjectedSV sv = uv;
			sv.haplo = 3;             // default homozygous; verify can correct
			sv.untagged_only = true;  // mark for FILTER=Untagged
			merged.push_back(sv);
			n_novel++;
		}
	}

	// Re-sort to keep merged in genome order
	std::sort(merged.begin(), merged.end(), [](const ProjectedSV& a, const ProjectedSV& b) {
		if (a.chrom != b.chrom) return a.chrom < b.chrom;
		return a.pos < b.pos;
	});
}


void Project::write_vcf(const std::vector<ProjectedSV>& svs, parameters& params)
{
	std::string vcf_path = params.log_path + params.sample_name + ".projected.vcf";
	std::ofstream fp(vcf_path);

	fp << "##fileformat=VCFv4.2\n";
	fp << "##source=veritig-project\n";
	fp << "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
	fp << "##FILTER=<ID=Ambig,Description=\"Inserted sequence has high k-mer abundance in flanking ref window; likely alignment artefact in repeat-rich region\">\n";
	fp << "##FILTER=<ID=Untagged,Description=\"Call derived only from SVarp untagged (unphased) svtigs; default GT=1|1 should be confirmed via sequence-level verification\">\n";

	// Emit ##contig=<ID=...,length=...> from reference .fai (if present)
	std::string fai_path = params.reference_path + ".fai";
	std::ifstream fai(fai_path);
	if (fai.good())
	{
		std::string line;
		while (getline(fai, line))
		{
			std::stringstream ss(line);
			std::string name, length_str;
			getline(ss, name, '\t');
			getline(ss, length_str, '\t');
			if (!name.empty() && !length_str.empty())
				fp << "##contig=<ID=" << name << ",length=" << length_str << ">\n";
		}
	}
	else
	{
		std::cerr << "  Note: no .fai found at " << fai_path << "; ##contig headers omitted\n";
	}

	fp << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
	fp << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variant (positive for INS/DUP/INV, negative for DEL)\">\n";
	fp << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant (1-based, inclusive)\">\n";
	fp << "##INFO=<ID=SVTIG,Number=1,Type=String,Description=\"Source svtig name(s); pipe-separated when multiple svtigs support this SV\">\n";
	fp << "##INFO=<ID=SVTIG_COUNT,Number=1,Type=Integer,Description=\"Number of svtigs supporting this SV\">\n";
	fp << "##INFO=<ID=ALN_MAPQ,Number=1,Type=Integer,Description=\"Maximum svtig alignment mapping quality across supporting svtigs\">\n";
	fp << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	fp << "##ALT=<ID=INS,Description=\"Insertion\">\n";
	fp << "##ALT=<ID=DEL,Description=\"Deletion\">\n";
	fp << "##ALT=<ID=DUP,Description=\"Tandem duplication (INS whose sequence matches a flanking ref window)\">\n";
	fp << "##ALT=<ID=INV,Description=\"Inversion (svtig segment aligned to reference on opposite strand)\">\n";
	fp << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << params.sample_name << "\n";

	int idx = 0;
	for (const auto& sv : svs)
	{
		std::string id = "veritig_" + std::to_string(++idx);
		std::string svtype_str;
		switch (sv.svtype) {
			case 'I': svtype_str = "INS"; break;
			case 'D': svtype_str = "DEL"; break;
			case 'U': svtype_str = "DUP"; break;
			case 'V': svtype_str = "INV"; break;
			default:  svtype_str = "INS"; break;
		}
		// SVLEN sign: negative for DEL only; INS, DUP, INV positive
		int signed_len = (sv.svtype == 'D') ? -sv.svlen : sv.svlen;

		std::string gt;
		if (sv.haplo == 0)      gt = "1/.";  // unphased single FASTA
		else if (sv.haplo == 3) gt = "1|1";
		else if (sv.haplo == 1) gt = "1|0";
		else                    gt = "0|1";  // haplo == 2

		// END for symbolic ALTs (INV, DEL with span). For point insertions (INS/DUP), END = POS.
		int end_pos;
		if (sv.svtype == 'D' || sv.svtype == 'V')
			end_pos = sv.pos + sv.svlen;       // span: last affected ref base
		else
			end_pos = sv.pos;                  // point insertion: anchor only

		fp << sv.chrom << "\t"
			<< sv.pos << "\t"
			<< id << "\t"
			<< sv.ref_base << "\t"
			<< sv.alt_seq << "\t"
			<< "." << "\t"
			<< (sv.ambig && sv.untagged_only ? "Ambig;Untagged" :
			    sv.ambig ? "Ambig" :
			    sv.untagged_only ? "Untagged" : "PASS") << "\t"
			<< "SVTYPE=" << svtype_str
			<< ";END=" << end_pos
			<< ";SVLEN=" << signed_len
			<< ";SVTIG=" << sv.svtig_name
			<< ";SVTIG_COUNT=" << sv.support_count
			<< ";ALN_MAPQ=" << sv.mapq << "\t"
			<< "GT" << "\t"
			<< gt << "\n";
	}
	std::cerr << "  VCF written: " << vcf_path << " (" << svs.size() << " records)\n";
}


void Project::write_bed(const std::vector<ProjectedSV>& svs, parameters& params)
{
	std::string bed_path = params.log_path + params.sample_name + ".projected.bed";
	std::ofstream fp(bed_path);

	fp << "#chrom\tstart\tend\tsvtig_id\tsvtype\tsvlen\tmapq\thaplo\talt_seq\n";
	for (const auto& sv : svs)
	{
		int start = sv.pos - 1;  // BED is 0-based
		// DEL and INV both span ref[pos..pos+svlen]; INS/DUP are point insertions
		int end = (sv.svtype == 'D' || sv.svtype == 'V') ? (start + sv.svlen) : start + 1;
		std::string svtype_str;
		switch (sv.svtype) {
			case 'I': svtype_str = "INS"; break;
			case 'D': svtype_str = "DEL"; break;
			case 'U': svtype_str = "DUP"; break;
			case 'V': svtype_str = "INV"; break;
			default:  svtype_str = "INS"; break;
		}
		std::string haplo_str;
		if (sv.haplo == 0)      haplo_str = "unphased";
		else if (sv.haplo == 3) haplo_str = "homo";
		else if (sv.haplo == 1) haplo_str = "H1";
		else                    haplo_str = "H2";

		// For BED, write only the inserted sequence for INS/DUP (without anchor) to save space,
		// "." for DEL/INV (no payload sequence; INV is symbolic)
		std::string alt_for_bed;
		if (sv.svtype == 'I' || sv.svtype == 'U')
			alt_for_bed = (sv.alt_seq.size() > 1) ? sv.alt_seq.substr(1) : ".";
		else
			alt_for_bed = ".";

		fp << sv.chrom << "\t"
			<< start << "\t"
			<< end << "\t"
			<< sv.svtig_name << "\t"
			<< svtype_str << "\t"
			<< sv.svlen << "\t"
			<< sv.mapq << "\t"
			<< haplo_str << "\t"
			<< alt_for_bed << "\n";
	}
	std::cerr << "  BED written: " << bed_path << "\n";
}


// Write a list of svtigs that did not align to the reference. These may represent
// novel insertions absent from the reference. Sequence is included for downstream use.
void Project::write_unmapped(const std::vector<std::string>& unmapped,
                             const std::map<std::string, std::string>& svtig_seqs,
                             parameters& params, const std::string& suffix)
{
	if (unmapped.empty()) return;
	std::string path = params.log_path + params.sample_name + suffix;
	std::ofstream fp(path);
	fp << "#svtig_name\tsize\tsequence\n";
	for (const auto& name : unmapped)
	{
		auto it = svtig_seqs.find(name);
		std::string seq = (it != svtig_seqs.end()) ? it->second : std::string();
		fp << name << "\t" << seq.size() << "\t" << seq << "\n";
	}
	std::cerr << "  Unmapped svtigs: " << unmapped.size() << " -> " << path << "\n";
}


void Project::run(parameters& params)
{
	bool phased = !params.svtig2_path.empty();
	std::cerr << "  Projecting svtigs to linear reference ("
		<< (phased ? "phased mode" : "unphased mode")
		<< ", min_svlen=" << params.min_svlen
		<< ", min_mapq=" << params.min_mapq << ")...\n";

	std::string paf_dir = params.log_path + "paf/";

	std::cerr << "  Loading reference FASTA (this may take a while)...\n";
	ref_seqs = load_fasta(params.reference_path);
	size_t ref_total = 0;
	for (const auto& kv : ref_seqs) ref_total += kv.second.size();
	std::cerr << "  Reference loaded: " << ref_seqs.size() << " contigs, "
		<< (ref_total / 1024 / 1024) << " MB total sequence\n";

	if (phased)
	{
		this->h1_paf_path = paf_dir + params.sample_name + ".project.H1.paf";
		this->h2_paf_path = paf_dir + params.sample_name + ".project.H2.paf";

		if (!params.skip_mapping)
		{
			run_mapping(params, params.svtig1_path, this->h1_paf_path);
			run_mapping(params, params.svtig2_path, this->h2_paf_path);
		}

		std::cerr << "  Loading H1 svtig FASTA...\n";
		auto h1_seqs = load_fasta(params.svtig1_path);
		std::cerr << "  Loading H2 svtig FASTA...\n";
		auto h2_seqs = load_fasta(params.svtig2_path);

		std::cerr << "  Parsing H1 PAF...\n";
		std::vector<std::string> h1_unmapped, h2_unmapped;
		std::vector<ProjectedSV> h1_svs = parse_paf(this->h1_paf_path, 1, params, h1_seqs, &h1_unmapped);
		std::cerr << "  Parsing H2 PAF...\n";
		std::vector<ProjectedSV> h2_svs = parse_paf(this->h2_paf_path, 2, params, h2_seqs, &h2_unmapped);

		write_unmapped(h1_unmapped, h1_seqs, params, ".unmapped_svtigs.H1.tsv");
		write_unmapped(h2_unmapped, h2_seqs, params, ".unmapped_svtigs.H2.tsv");

		std::cerr << "  Merging haplotypes (pos_tol=" << params.merge_pos_tol
			<< ", size_tol=" << params.merge_size_tol << ")...\n";
		size_t n_h1_before = h1_svs.size();
		size_t n_h2_before = h2_svs.size();
		std::vector<ProjectedSV> merged = merge_haplotypes(h1_svs, h2_svs,
			params.merge_pos_tol, params.merge_size_tol);
		size_t n_homo_merged = n_h1_before + n_h2_before - merged.size();
		std::cerr << "  Tolerance merge: " << n_homo_merged << " homozygous matches\n";

		// Optional: process untagged (unphased) svtigs and merge into the H1/H2 result
		if (!params.svtig3_path.empty())
		{
			this->untagged_paf_path = paf_dir + params.sample_name + ".project.untagged.paf";
			if (!params.skip_mapping)
				run_mapping(params, params.svtig3_path, this->untagged_paf_path);

			std::cerr << "  Loading untagged svtig FASTA...\n";
			auto un_seqs = load_fasta(params.svtig3_path);
			std::cerr << "  Parsing untagged PAF...\n";
			std::vector<std::string> un_unmapped;
			std::vector<ProjectedSV> un_svs = parse_paf(this->untagged_paf_path, 0, params,
				un_seqs, &un_unmapped);
			write_unmapped(un_unmapped, un_seqs, params, ".unmapped_svtigs.untagged.tsv");

			std::cerr << "  Merging untagged...\n";
			int n_matched = 0, n_novel = 0;
			merge_untagged(merged, un_svs, params.merge_pos_tol, params.merge_size_tol,
				n_matched, n_novel);
			std::cerr << "  Untagged: " << n_matched << " matched existing calls, "
				<< n_novel << " novel calls (default GT=1|1, FILTER=Untagged)\n";
		}

		std::cerr << "  Classifying tandem duplications...\n";
		classify_dups(merged, params);
		std::cerr << "  Flagging repeat-ambiguous events...\n";
		flag_ambiguous_repeats(merged, params);

		int n_homo = 0, n_h1 = 0, n_h2 = 0, n_untagged = 0;
		for (const auto& sv : merged)
		{
			if (sv.untagged_only) n_untagged++;
			if (sv.haplo == 3) n_homo++;
			else if (sv.haplo == 1) n_h1++;
			else n_h2++;
		}
		std::cerr << "  Final: " << merged.size() << " SVs ("
			<< n_homo << " homozygous, " << n_h1 << " H1-only, " << n_h2 << " H2-only";
		if (n_untagged > 0) std::cerr << ", " << n_untagged << " untagged-only";
		std::cerr << ")\n";

		write_vcf(merged, params);
		write_bed(merged, params);
	}
	else
	{
		// Unphased: single svtig FASTA, no haplotype assignment
		this->h1_paf_path = paf_dir + params.sample_name + ".project.paf";

		if (!params.skip_mapping)
			run_mapping(params, params.svtig1_path, this->h1_paf_path);

		std::cerr << "  Loading svtig FASTA...\n";
		auto svtig_seqs = load_fasta(params.svtig1_path);

		std::cerr << "  Parsing PAF...\n";
		std::vector<std::string> unmapped;
		std::vector<ProjectedSV> svs = parse_paf(this->h1_paf_path, 0, params, svtig_seqs, &unmapped);

		write_unmapped(unmapped, svtig_seqs, params, ".unmapped_svtigs.tsv");

		std::cerr << "  Classifying tandem duplications...\n";
		classify_dups(svs, params);
		std::cerr << "  Flagging repeat-ambiguous events...\n";
		flag_ambiguous_repeats(svs, params);

		std::cerr << "  Final: " << svs.size() << " SVs (unphased)\n";

		write_vcf(svs, params);
		write_bed(svs, params);
	}
}
