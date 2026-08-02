#ifndef __VERIFY
#define __VERIFY

#include <string>
#include <vector>
#include <map>
#include "common.h"


// One SV record parsed from the input VCF (typically a veritig --project output)
struct VerifySV
{
	// VCF fields
	std::string chrom;
	int pos = 0;             // 1-based VCF POS
	std::string id;
	std::string ref;
	std::string alt;         // first ALT only; symbolic forms (<INV>) handled
	std::string filter_str;  // raw FILTER column ("PASS" or "Ambig" or "Ambig;Untagged" etc.)
	// INFO subset we care about
	std::string svtype;      // INS / DEL / DUP / INV
	int svlen = 0;           // absolute value
	std::string svtig;       // SVTIG INFO (pipe-separated svtig names)
	int svtig_count = 0;
	int aln_mapq = 0;
	int end_pos = 0;         // INFO/END (for symbolic ALT)
	// FORMAT/sample
	std::string gt_str;      // original GT like "1|0"
	// Raw rest-of-line for passthrough write (everything after CHROM\tPOS\tID\tREF\tALT)
	std::string raw_qual_filter_info_format_sample;
};


// Per-(query, target_contig, strand) aggregated alignment stats
struct VerifyAlnGroup
{
	int qlen = 0;
	int matches = 0;
	int best_mapq = 0;
	int max_ins_gap = 0;     // largest insertion gap inside any CIGAR (relative to query)
	int max_del_gap = 0;     // largest deletion gap inside any CIGAR
	std::vector<int> qstarts, qends;
	std::vector<int> tstarts, tends;
	std::string tname;
	char strand = '+';
};


// Final per-query verification metrics on a single haplotype
struct VerifyResult
{
	double map_ratio = 0.0;
	double identity = 0.0;
	int mapq = 0;
	int qspan = 0;
	int tspan = 0;
	int span_delta = 0;      // tspan - qspan (positive = extra in target)
	int max_ins_gap = 0;
	int max_del_gap = 0;
	bool present = false;    // does this query have any alignment to this haplotype?
	std::string tname;
};


class Verify
{
public:
	void run(parameters& params);

	// Sequence utility, also used when projecting minus-strand alignments.
	static std::string reverse_complement(const std::string& s);

private:
	// Pipeline steps
	std::vector<VerifySV> read_vcf(const std::string& vcf_path);
	std::string build_expected_seq(const VerifySV& sv,
	                               const std::map<std::string, std::string>& ref_seqs,
	                               int flank, std::string& reason);
	void write_expected_fasta(const std::vector<VerifySV>& svs,
	                          const std::map<std::string, std::string>& ref_seqs,
	                          int flank, const std::string& out_path,
	                          std::vector<std::string>& skip_reasons);
	void run_minimap(parameters& params, const std::string& query_fa,
	                 const std::string& target_fa, const std::string& paf_out);
	std::map<std::string, VerifyResult> parse_paf(const std::string& paf_path);
	bool passes(const VerifyResult& r, int svlen, double gap_frac, int min_indel_gap,
	            double min_map_ratio, double min_identity);
	void write_vcf(const std::string& in_vcf_path, const std::string& out_vcf_path,
	               const std::vector<VerifySV>& svs,
	               const std::map<std::string, VerifyResult>& h1_res,
	               const std::map<std::string, VerifyResult>& h2_res,
	               parameters& params,
	               int& n_total, int& n_verified,
	               std::map<std::string, std::map<std::string,int>>& by_type);
	void write_summary(const std::string& path, int n_total, int n_verified,
	                   const std::map<std::string, std::map<std::string,int>>& by_type);

	// Helpers
	static int parse_int_info(const std::string& info_field);
	static void parse_cigar_max_indel(const std::string& cigar, int& max_ins, int& max_del);
};

#endif
