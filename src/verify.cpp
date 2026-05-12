#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <set>
#include <map>
#include "verify.h"

extern std::map<std::string, std::string> load_fasta_external(const std::string& path);
void run_command(const std::string& cmd, const std::string& output_file);


// --- helpers ---

std::string Verify::reverse_complement(const std::string& s)
{
	std::string out;
	out.resize(s.size());
	auto comp = [](char c) {
		switch (c) {
			case 'A': return 'T'; case 'a': return 't';
			case 'T': return 'A'; case 't': return 'a';
			case 'C': return 'G'; case 'c': return 'g';
			case 'G': return 'C'; case 'g': return 'c';
			case 'N': return 'N'; case 'n': return 'n';
			default: return 'N';
		}
	};
	for (size_t i = 0; i < s.size(); i++)
		out[s.size() - 1 - i] = comp(s[i]);
	return out;
}


int Verify::parse_int_info(const std::string& info_field)
{
	try { return std::stoi(info_field); } catch (...) { return 0; }
}


// Walk a CIGAR like "200M50I300M" and record the largest I and D op sizes.
void Verify::parse_cigar_max_indel(const std::string& cigar, int& max_ins, int& max_del)
{
	max_ins = 0; max_del = 0;
	if (cigar.empty()) return;
	int n = 0;
	for (char c : cigar)
	{
		if (c >= '0' && c <= '9')
		{
			n = n * 10 + (c - '0');
		}
		else
		{
			if (c == 'I') { if (n > max_ins) max_ins = n; }
			else if (c == 'D') { if (n > max_del) max_del = n; }
			n = 0;
		}
	}
}


// --- VCF reader (line-based, VCF v4.2 dialect produced by veritig --project) ---

std::vector<VerifySV> Verify::read_vcf(const std::string& vcf_path)
{
	std::vector<VerifySV> out;
	std::ifstream fp(vcf_path);
	if (!fp.good())
	{
		std::cerr << "[veritig] Cannot open input VCF: " << vcf_path << "\n";
		return out;
	}

	std::string line;
	while (getline(fp, line))
	{
		if (line.empty() || line[0] == '#') continue;
		std::vector<std::string> cols;
		{
			std::stringstream ss(line);
			std::string tok;
			while (getline(ss, tok, '\t')) cols.push_back(tok);
		}
		if (cols.size() < 8) continue;

		VerifySV sv;
		sv.chrom = cols[0];
		sv.pos = std::atoi(cols[1].c_str());
		sv.id = cols[2];
		sv.ref = cols[3];
		// Take first ALT only (split on comma)
		std::string alt_field = cols[4];
		size_t comma = alt_field.find(',');
		sv.alt = (comma == std::string::npos) ? alt_field : alt_field.substr(0, comma);
		sv.filter_str = cols.size() > 6 ? cols[6] : "PASS";

		// Parse INFO (cols[7])
		std::string info_str = cols[7];
		std::stringstream is(info_str);
		std::string kv;
		while (getline(is, kv, ';'))
		{
			size_t eq = kv.find('=');
			std::string key = (eq == std::string::npos) ? kv : kv.substr(0, eq);
			std::string val = (eq == std::string::npos) ? "" : kv.substr(eq + 1);
			if (key == "SVTYPE") sv.svtype = val;
			else if (key == "SVLEN")
			{
				int v = std::atoi(val.c_str());
				sv.svlen = std::abs(v);
			}
			else if (key == "END") sv.end_pos = std::atoi(val.c_str());
			else if (key == "SVTIG") sv.svtig = val;
			else if (key == "SVTIG_COUNT") sv.svtig_count = std::atoi(val.c_str());
			else if (key == "ALN_MAPQ") sv.aln_mapq = std::atoi(val.c_str());
		}

		// Parse FORMAT/sample (cols[8] and cols[9..]). Find GT index.
		if (cols.size() >= 10)
		{
			std::string format = cols[8];
			std::string sample = cols[9];
			std::stringstream fs(format), ss(sample);
			std::vector<std::string> fkeys, svals;
			std::string t;
			while (getline(fs, t, ':')) fkeys.push_back(t);
			while (getline(ss, t, ':')) svals.push_back(t);
			for (size_t i = 0; i < fkeys.size() && i < svals.size(); i++)
			{
				if (fkeys[i] == "GT") { sv.gt_str = svals[i]; break; }
			}
		}

		out.push_back(sv);
	}
	return out;
}


// --- Expected sequence builder ---

std::string Verify::build_expected_seq(const VerifySV& sv,
	const std::map<std::string, std::string>& ref_seqs, int flank, std::string& reason)
{
	reason = "ok";
	auto it = ref_seqs.find(sv.chrom);
	if (it == ref_seqs.end()) { reason = "ref_chrom_missing"; return ""; }
	const std::string& chrom_seq = it->second;
	int ctg_len = (int)chrom_seq.size();
	int start = sv.pos - 1;  // 0-based
	if (sv.svlen <= 0) { reason = "svlen_zero"; return ""; }

	int up_start = std::max(0, start - flank);
	int up_end = start;
	if (up_end <= up_start) { reason = "no_upstream_flank"; return ""; }
	std::string upstream = chrom_seq.substr(up_start, up_end - up_start);

	std::string expected;

	if (sv.svtype == "DEL")
	{
		// anchor + skip deleted region + downstream
		std::string anchor = sv.ref.empty() ? "" : sv.ref.substr(0, 1);
		int ds_start = start + (int)sv.ref.size();
		int ds_end = std::min(ctg_len, ds_start + flank);
		if (ds_end <= ds_start) { reason = "no_downstream_flank"; return ""; }
		std::string downstream = chrom_seq.substr(ds_start, ds_end - ds_start);
		expected = upstream + anchor + downstream;
	}
	else if (sv.svtype == "INS" || sv.svtype == "DUP")
	{
		if (sv.alt.empty() || sv.alt[0] == '<') { reason = "symbolic_ins_dup_unsupported"; return ""; }
		if (sv.alt.size() <= 1) { reason = "alt_too_short"; return ""; }
		int ds_start = start + (int)sv.ref.size();
		int ds_end = std::min(ctg_len, ds_start + flank);
		if (ds_end <= ds_start) { reason = "no_downstream_flank"; return ""; }
		std::string downstream = chrom_seq.substr(ds_start, ds_end - ds_start);
		// Uppercase ALT
		std::string alt_up = sv.alt;
		for (char& c : alt_up) c = std::toupper(c);
		expected = upstream + alt_up + downstream;
	}
	else if (sv.svtype == "INV")
	{
		int inv_start = start;
		int inv_end = std::min(ctg_len, start + sv.svlen);
		if (inv_end <= inv_start) { reason = "no_inv_segment"; return ""; }
		std::string inv_seq = chrom_seq.substr(inv_start, inv_end - inv_start);
		int ds_start = inv_end;
		int ds_end = std::min(ctg_len, ds_start + flank);
		std::string downstream = ds_end > ds_start ? chrom_seq.substr(ds_start, ds_end - ds_start) : "";
		expected = upstream + reverse_complement(inv_seq) + downstream;
	}
	else
	{
		reason = "unsupported_svtype";
		return "";
	}

	if ((int)expected.size() < 100) { reason = "expected_too_short"; return ""; }
	return expected;
}


void Verify::write_expected_fasta(const std::vector<VerifySV>& svs,
	const std::map<std::string, std::string>& ref_seqs, int flank,
	const std::string& out_path, std::vector<std::string>& skip_reasons)
{
	std::ofstream out(out_path);
	skip_reasons.clear();
	skip_reasons.resize(svs.size(), "");
	int n_emitted = 0;
	std::map<std::string, int> reason_counts;
	for (size_t i = 0; i < svs.size(); i++)
	{
		std::string reason;
		std::string seq = build_expected_seq(svs[i], ref_seqs, flank, reason);
		if (seq.empty())
		{
			skip_reasons[i] = reason;
			reason_counts[reason]++;
			continue;
		}
		// Use stable per-record name: prefer VCF ID if non-empty/non-".", else "sv_<index>"
		std::string sv_id = (svs[i].id.empty() || svs[i].id == ".") ?
			("sv_" + std::to_string(i + 1)) : svs[i].id;
		out << ">" << sv_id << "\n" << seq << "\n";
		n_emitted++;
	}
	std::cerr << "  Built " << n_emitted << " expected sequences from " << svs.size() << " VCF records\n";
	if (!reason_counts.empty())
	{
		std::cerr << "  Skipped (";
		bool first = true;
		for (auto& kv : reason_counts)
		{
			if (!first) std::cerr << ", ";
			std::cerr << kv.first << "=" << kv.second;
			first = false;
		}
		std::cerr << ")\n";
	}
}


// --- minimap2 wrapper (target = assembly, query = expected.fa) ---

void Verify::run_minimap(parameters& params, const std::string& query_fa,
	const std::string& target_fa, const std::string& paf_out)
{
	std::string threads_str = std::to_string(params.threads);
	std::cerr << "  Mapping expected_seq to " << target_fa << "...";
	std::string cmd = "minimap2 -cx asm5 -t " + threads_str + " "
		+ "--secondary=no " + target_fa + " " + query_fa
		+ " > " + paf_out + " 2>/dev/null";
	run_command(cmd, paf_out);
	std::cerr << " done\n";
}


// --- PAF parser with bounding-box span + CIGAR gap tracking ---

std::map<std::string, VerifyResult> Verify::parse_paf(const std::string& paf_path)
{
	// Per-(qname, tname, strand) groups
	std::map<std::tuple<std::string,std::string,char>, VerifyAlnGroup> groups;
	std::ifstream fp(paf_path);
	if (!fp.good())
	{
		std::cerr << "[veritig] Cannot open PAF: " << paf_path << "\n";
		return {};
	}
	std::string line;
	while (getline(fp, line))
	{
		if (line.empty()) continue;
		std::vector<std::string> cols;
		std::stringstream ss(line);
		std::string tok;
		while (getline(ss, tok, '\t')) cols.push_back(tok);
		if (cols.size() < 12) continue;

		std::string qname = cols[0];
		int qlen = std::atoi(cols[1].c_str());
		int qstart = std::atoi(cols[2].c_str());
		int qend = std::atoi(cols[3].c_str());
		char strand = cols[4][0];
		std::string tname = cols[5];
		int tstart = std::atoi(cols[7].c_str());
		int tend = std::atoi(cols[8].c_str());
		int matches = std::atoi(cols[9].c_str());
		int mapq = std::atoi(cols[11].c_str());

		std::string cigar;
		for (size_t i = 12; i < cols.size(); i++)
		{
			if (cols[i].compare(0, 5, "cg:Z:") == 0)
			{
				cigar = cols[i].substr(5);
				break;
			}
		}
		int ins_gap = 0, del_gap = 0;
		parse_cigar_max_indel(cigar, ins_gap, del_gap);

		auto key = std::make_tuple(qname, tname, strand);
		auto& g = groups[key];
		g.qlen = qlen;
		g.qstarts.push_back(qstart); g.qends.push_back(qend);
		g.tstarts.push_back(tstart); g.tends.push_back(tend);
		g.matches += matches;
		if (ins_gap > g.max_ins_gap) g.max_ins_gap = ins_gap;
		if (del_gap > g.max_del_gap) g.max_del_gap = del_gap;
		if (mapq > g.best_mapq) g.best_mapq = mapq;
		g.tname = tname;
		g.strand = strand;
	}

	// For each query, pick the best group (largest aligned q span) and compute final metrics
	std::map<std::string, std::vector<std::tuple<std::string,char,VerifyAlnGroup>>> per_query;
	for (auto& kv : groups)
	{
		const std::string& qname = std::get<0>(kv.first);
		per_query[qname].push_back(std::make_tuple(std::get<1>(kv.first), std::get<2>(kv.first), kv.second));
	}

	std::map<std::string, VerifyResult> out;
	for (auto& kv : per_query)
	{
		const std::string& qname = kv.first;
		size_t best_idx = 0;
		int best_aligned = -1;
		for (size_t i = 0; i < kv.second.size(); i++)
		{
			const VerifyAlnGroup& g = std::get<2>(kv.second[i]);
			int aligned_sum = 0;
			for (size_t j = 0; j < g.qstarts.size(); j++)
				aligned_sum += g.qends[j] - g.qstarts[j];
			if (aligned_sum > best_aligned) { best_aligned = aligned_sum; best_idx = i; }
		}
		const std::string& tname = std::get<0>(kv.second[best_idx]);
		char strand = std::get<1>(kv.second[best_idx]);
		const VerifyAlnGroup& g = std::get<2>(kv.second[best_idx]);

		int qlen = g.qlen;
		int aligned = 0;
		for (size_t j = 0; j < g.qstarts.size(); j++)
			aligned += g.qends[j] - g.qstarts[j];
		if (aligned > qlen) aligned = qlen;

		int qmin = *std::min_element(g.qstarts.begin(), g.qstarts.end());
		int qmax = *std::max_element(g.qends.begin(), g.qends.end());
		int tmin = *std::min_element(g.tstarts.begin(), g.tstarts.end());
		int tmax = *std::max_element(g.tends.begin(), g.tends.end());
		int qspan = qmax - qmin;
		int tspan = tmax - tmin;

		VerifyResult r;
		r.map_ratio = qlen ? (double)aligned / (double)qlen : 0.0;
		r.identity = aligned ? (double)g.matches / (double)aligned : 0.0;
		r.mapq = g.best_mapq;
		r.qspan = qspan;
		r.tspan = tspan;
		r.span_delta = tspan - qspan;
		r.max_ins_gap = g.max_ins_gap;
		r.max_del_gap = g.max_del_gap;
		r.present = true;
		r.tname = tname;
		(void)strand;
		out[qname] = r;
	}
	return out;
}


// --- Verification logic ---

bool Verify::passes(const VerifyResult& r, int svlen, double gap_frac, int min_indel_gap,
	double min_map_ratio, double min_identity)
{
	if (!r.present) return false;
	if (r.map_ratio < min_map_ratio) return false;
	if (r.identity < min_identity) return false;
	int gap_threshold = std::max(min_indel_gap, (int)(svlen * gap_frac));
	if (r.max_ins_gap >= gap_threshold) return false;
	if (r.max_del_gap >= gap_threshold) return false;
	if (std::abs(r.span_delta) >= gap_threshold) return false;
	return true;
}


// --- Output VCF writer (passthrough header from input VCF + new INFO/FILTER lines) ---

void Verify::write_vcf(const std::string& in_vcf_path, const std::string& out_vcf_path,
	const std::vector<VerifySV>& svs,
	const std::map<std::string, VerifyResult>& h1_res,
	const std::map<std::string, VerifyResult>& h2_res,
	parameters& params,
	int& n_total, int& n_verified,
	std::map<std::string, std::map<std::string,int>>& by_type)
{
	n_total = 0; n_verified = 0;
	std::ifstream in(in_vcf_path);
	std::ofstream out(out_vcf_path);
	if (!in.good() || !out.good())
	{
		std::cerr << "[veritig] verify write_vcf: failed to open files\n";
		return;
	}

	bool header_done = false;
	bool injected = false;
	std::string line;

	// We re-use the original VCF header, append our INFO lines just before #CHROM.
	auto inject_headers = [&out]() {
		out << "##INFO=<ID=VERIFIED,Number=1,Type=String,Description=\"Sequence-level verification result vs assembly: yes/no\">\n";
		out << "##INFO=<ID=VAL_HAP,Number=1,Type=String,Description=\"Haplotype where verified (H1/H2/Both/None)\">\n";
		out << "##INFO=<ID=VAL_MAP_RATIO,Number=1,Type=Float,Description=\"Best aligned fraction of expected_seq\">\n";
		out << "##INFO=<ID=VAL_IDENTITY,Number=1,Type=Float,Description=\"Best identity of expected_seq alignment\">\n";
		out << "##INFO=<ID=VAL_REASON,Number=1,Type=String,Description=\"Reason for verification result\">\n";
	};

	size_t sv_idx = 0;

	while (getline(in, line))
	{
		if (line.empty()) { out << "\n"; continue; }
		if (line[0] == '#')
		{
			if (line.compare(0, 6, "#CHROM") == 0)
			{
				if (!injected) { inject_headers(); injected = true; }
			}
			out << line << "\n";
			if (line.compare(0, 6, "#CHROM") == 0) header_done = true;
			continue;
		}
		if (!header_done) continue;

		// Parse this record minimally (we already have parsed list aligned by index)
		if (sv_idx >= svs.size()) { out << line << "\n"; continue; }
		const VerifySV& sv = svs[sv_idx];

		// Lookup verify results
		std::string sv_id = (sv.id.empty() || sv.id == ".") ?
			("sv_" + std::to_string(sv_idx + 1)) : sv.id;
		auto it1 = h1_res.find(sv_id);
		auto it2 = h2_res.find(sv_id);
		const VerifyResult* r1 = (it1 != h1_res.end()) ? &it1->second : nullptr;
		const VerifyResult* r2 = (it2 != h2_res.end()) ? &it2->second : nullptr;

		bool h1_ok = r1 && passes(*r1, sv.svlen, params.verify_gap_frac, params.verify_min_indel_gap,
		                          params.min_map_ratio, params.min_aln_identity);
		bool h2_ok = r2 && passes(*r2, sv.svlen, params.verify_gap_frac, params.verify_min_indel_gap,
		                          params.min_map_ratio, params.min_aln_identity);

		// Pick best hap for reporting metrics
		double s1 = r1 ? r1->map_ratio * r1->identity : -1.0;
		double s2 = r2 ? r2->map_ratio * r2->identity : -1.0;
		const VerifyResult* best = (s1 >= s2) ? r1 : r2;
		double best_R = best ? best->map_ratio : 0.0;
		double best_I = best ? best->identity : 0.0;

		std::string val_hap = "None";
		std::string reason;
		bool verified = false;
		std::string new_gt = sv.gt_str;
		if (h1_ok && h2_ok)      { val_hap = "Both"; reason = "verified_in_both"; verified = true; new_gt = "1|1"; }
		else if (h1_ok)          { val_hap = "H1"; reason = "verified_in_H1"; verified = true; new_gt = "1|0"; }
		else if (h2_ok)          { val_hap = "H2"; reason = "verified_in_H2"; verified = true; new_gt = "0|1"; }
		else
		{
			val_hap = "None";
			if (!best) reason = "no_alignment";
			else if (best->map_ratio < params.min_map_ratio) reason = "low_map_ratio";
			else if (best->identity < params.min_aln_identity) reason = "low_identity";
			else reason = "gap_inconsistent";
		}

		// Re-split the original line and rewrite INFO + GT
		std::vector<std::string> cols;
		{
			std::stringstream ss(line);
			std::string tok;
			while (getline(ss, tok, '\t')) cols.push_back(tok);
		}
		if (cols.size() < 8) { out << line << "\n"; sv_idx++; continue; }

		// Append new INFO fields
		std::ostringstream ann;
		ann << ";VERIFIED=" << (verified ? "yes" : "no")
		    << ";VAL_HAP=" << val_hap
		    << ";VAL_MAP_RATIO=" << best_R
		    << ";VAL_IDENTITY=" << best_I
		    << ";VAL_REASON=" << reason;
		cols[7] += ann.str();

		// Update GT in sample column if present (pos 0 of FORMAT)
		if (cols.size() >= 10)
		{
			std::string format = cols[8];
			std::string sample = cols[9];
			std::vector<std::string> fkeys, svals;
			std::stringstream fs(format), ss(sample);
			std::string t;
			while (getline(fs, t, ':')) fkeys.push_back(t);
			while (getline(ss, t, ':')) svals.push_back(t);
			for (size_t i = 0; i < fkeys.size() && i < svals.size(); i++)
			{
				if (fkeys[i] == "GT")
				{
					svals[i] = new_gt;
					break;
				}
			}
			std::ostringstream so;
			for (size_t i = 0; i < svals.size(); i++)
			{
				if (i) so << ":";
				so << svals[i];
			}
			cols[9] = so.str();
		}

		// Emit
		for (size_t i = 0; i < cols.size(); i++)
		{
			if (i) out << "\t";
			out << cols[i];
		}
		out << "\n";

		// Stats
		n_total++;
		if (verified) n_verified++;
		auto& bt = by_type[sv.svtype];
		bt["total"]++;
		if (verified) bt["verified"]++;

		sv_idx++;
	}
}


void Verify::write_summary(const std::string& path, int n_total, int n_verified,
	const std::map<std::string, std::map<std::string,int>>& by_type)
{
	std::ofstream out(path);
	out << "svtype\ttotal\tverified\tverified_pct\n";
	for (const auto& kv : by_type)
	{
		int t = 0, v = 0;
		auto it_t = kv.second.find("total");
		auto it_v = kv.second.find("verified");
		if (it_t != kv.second.end()) t = it_t->second;
		if (it_v != kv.second.end()) v = it_v->second;
		double pct = t > 0 ? 100.0 * v / t : 0.0;
		out << kv.first << "\t" << t << "\t" << v << "\t" << pct << "\n";
	}
	double pct = n_total > 0 ? 100.0 * n_verified / n_total : 0.0;
	out << "TOTAL\t" << n_total << "\t" << n_verified << "\t" << pct << "\n";
}


// --- main entry ---

void Verify::run(parameters& params)
{
	if (params.verify_vcf_path.empty())
	{
		std::cerr << "[veritig] --verify requires --vcf <input.vcf>\n";
		return;
	}
	if (params.haplo1_assembly_path.empty() || params.reference_path.empty())
	{
		std::cerr << "[veritig] --verify requires --ref and --h1 (and --h2 for diploid)\n";
		return;
	}

	std::cerr << "  Reading input VCF " << params.verify_vcf_path << "...\n";
	std::vector<VerifySV> svs = read_vcf(params.verify_vcf_path);
	std::cerr << "  Loaded " << svs.size() << " SV records\n";

	std::cerr << "  Loading reference FASTA...\n";
	std::map<std::string, std::string> ref_seqs = load_fasta_external(params.reference_path);
	std::cerr << "  Reference: " << ref_seqs.size() << " contigs\n";

	std::string expected_fa = params.log_path + params.sample_name + ".verify.expected.fa";
	std::string h1_paf = params.log_path + params.sample_name + ".verify.h1.paf";
	std::string h2_paf = params.log_path + params.sample_name + ".verify.h2.paf";

	std::cerr << "  Building expected sequences (flank=" << params.verify_flank << ")...\n";
	std::vector<std::string> skip_reasons;
	write_expected_fasta(svs, ref_seqs, params.verify_flank, expected_fa, skip_reasons);

	// Reference is only needed for expected_seq construction. Free it before launching
	// minimap2 (which itself loads the reference into ~15 GB hash tables). This avoids
	// holding a redundant copy of the reference sequence in veritig memory.
	std::map<std::string, std::string>().swap(ref_seqs);
	std::cerr << "  Reference released from memory after expected_seq construction\n";

	if (!params.skip_mapping)
	{
		run_minimap(params, expected_fa, params.haplo1_assembly_path, h1_paf);
		if (!params.haplo2_assembly_path.empty())
			run_minimap(params, expected_fa, params.haplo2_assembly_path, h2_paf);
	}

	std::cerr << "  Parsing H1 PAF...\n";
	auto h1_results = parse_paf(h1_paf);
	std::cerr << "  H1: " << h1_results.size() << "/" << svs.size() << " SVs have any alignment\n";

	std::map<std::string, VerifyResult> h2_results;
	if (!params.haplo2_assembly_path.empty())
	{
		std::cerr << "  Parsing H2 PAF...\n";
		h2_results = parse_paf(h2_paf);
		std::cerr << "  H2: " << h2_results.size() << "/" << svs.size() << " SVs have any alignment\n";
	}

	std::string out_vcf = params.log_path + params.sample_name + ".verified.vcf";
	std::string summary_path = params.log_path + params.sample_name + ".verify.summary.tsv";
	std::cerr << "  Writing verified VCF " << out_vcf << "...\n";
	int n_total = 0, n_verified = 0;
	std::map<std::string, std::map<std::string,int>> by_type;
	write_vcf(params.verify_vcf_path, out_vcf, svs, h1_results, h2_results, params,
	          n_total, n_verified, by_type);
	write_summary(summary_path, n_total, n_verified, by_type);

	std::cerr << "  Final: " << n_verified << "/" << n_total
	          << " verified (" << (n_total ? (100.0 * n_verified / n_total) : 0.0) << "%)\n";
	std::cerr << "  Summary: " << summary_path << "\n";
}


// Bridge to project.cpp's load_fasta. To avoid duplicating the FASTA parser, we expose
// project's via this thin helper. Implemented here to keep linkage simple.
std::map<std::string, std::string> load_fasta_external(const std::string& path)
{
	std::map<std::string, std::string> out;
	std::ifstream fp(path);
	if (!fp.good())
	{
		std::cerr << "[veritig] Cannot open FASTA: " << path << "\n";
		return out;
	}
	std::string line, cur_name, cur_seq;
	while (getline(fp, line))
	{
		if (line.empty()) continue;
		if (line[0] == '>')
		{
			if (!cur_name.empty()) out[cur_name] = cur_seq;
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
	if (!cur_name.empty()) out[cur_name] = cur_seq;
	return out;
}
