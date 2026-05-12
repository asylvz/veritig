#ifndef __PROJECT
#define __PROJECT

#include <string>
#include <vector>
#include <map>
#include "common.h"


struct ProjectedSV
{
	std::string chrom;
	int pos;               // 1-based reference position (VCF convention)
	char svtype;           // 'I' for INS, 'D' for DEL, 'U' for DUP (tandem duplication)
	int svlen;             // positive for INS, positive for DEL (sign handled at VCF write)
	std::string ref_base;  // anchor base on reference
	std::string alt_seq;   // full ALT sequence for INS (ref_base + inserted), or ref_base for DEL
	std::string svtig_name;   // comma-separated list when multiple svtigs support this SV
	int support_count = 1;    // number of svtigs supporting this SV
	int mapq;
	int haplo;             // 0 = unphased, 1 = H1, 2 = H2, 3 = both (homozygous)
	bool ambig = false;    // true if event is in a repeat-rich ref window (Ambig FILTER tag)
	bool untagged_only = false;  // true if event came only from SVarp untagged svtigs (FILTER=Untagged)
};


// One PAF alignment record retained per (svtig, alignment) for split-aware analysis
struct AlignmentRecord
{
	int q_start = 0, q_end = 0, q_len = 0;
	std::string t_chrom;
	int t_start = 0, t_end = 0;     // 0-based ref coords from PAF
	char strand = '+';
	int mapq = 0;
	bool is_primary = true;
};


class Project
{
public:
	void run(parameters& params);

private:
	void run_mapping(parameters& params, const std::string& svtig_fasta, const std::string& paf_out);
	std::vector<ProjectedSV> parse_paf(const std::string& paf_path, int haplo, parameters& params,
	                                   const std::map<std::string, std::string>& svtig_seqs,
	                                   std::vector<std::string>* unmapped_out);
	std::vector<ProjectedSV> extract_split_ins(
	    const std::map<std::string, std::vector<AlignmentRecord>>& svtig_alns,
	    const std::map<std::string, std::string>& svtig_seqs,
	    int haplo, parameters& params);
	std::vector<ProjectedSV> cluster_within_svtig(std::vector<ProjectedSV>& events, int window, double size_ratio);
	void classify_dups(std::vector<ProjectedSV>& events, parameters& params);
	void flag_ambiguous_repeats(std::vector<ProjectedSV>& events, parameters& params);
	std::vector<ProjectedSV> merge_haplotypes(std::vector<ProjectedSV>& h1_svs, std::vector<ProjectedSV>& h2_svs);
	std::vector<ProjectedSV> merge_haplotypes(std::vector<ProjectedSV>& h1_svs, std::vector<ProjectedSV>& h2_svs,
		int pos_tol, double size_tol);
	void merge_untagged(std::vector<ProjectedSV>& merged, std::vector<ProjectedSV>& untagged,
		int pos_tol, double size_tol, int& n_matched, int& n_novel);
	void write_vcf(const std::vector<ProjectedSV>& svs, parameters& params);
	void write_bed(const std::vector<ProjectedSV>& svs, parameters& params);
	void write_unmapped(const std::vector<std::string>& unmapped,
	                    const std::map<std::string, std::string>& svtig_seqs,
	                    parameters& params, const std::string& suffix);
	std::map<std::string, std::string> load_fasta(const std::string& fasta_path);

	std::string h1_paf_path;
	std::string h2_paf_path;
	std::string untagged_paf_path;
	std::map<std::string, std::string> ref_seqs;   // reference contig name -> sequence (uppercase)
};

#endif
