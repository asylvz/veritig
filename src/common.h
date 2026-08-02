#ifndef __COMMON
#define __COMMON

#include <fstream>
#include <vector>
#include <string>

#define RETURN_SUCCESS 0
#define RETURN_ERROR -1

#define MINSVSIZE 50

#define DELETION 'D'
#define INSERTION 'I'


typedef struct _parameters
{
	std::string fasta;
	std::string output_path;
	std::string svtig1_path;
	std::string svtig2_path;
	std::string svtig3_path;   // optional: SVarp untagged svtigs (cannot be phased)
	std::string log_path;
	std::string haplo1_assembly_path;
	std::string haplo2_assembly_path;
	std::string reference_path;
	std::string sample_name;
	bool skip_mapping = false;
	bool concordance = false;
	bool validate = false;
	bool stats = false;
	bool filter = false;
	bool compare = false;
	bool project = false;
	bool verify = false;
	bool phase = false;
	bool detailed = false;
	int threads = 16;
	int min_svlen = 50;
	int min_mapq = 20;
	int cluster_window = 200;
	double cluster_size_ratio = 0.7;
	int dedup_pos_tol = 50;
	double dedup_size_tol = 0.1;
	bool split_ins = false;           // recover INS from split alignments; off by default (see --split-ins)
	int min_split_svlen = 2000;       // matches minimap2 -r2k chain gap threshold; smaller gaps may be chain heuristics, not real SVs
	int max_split_distance = 100000;
	double dup_similarity = 0.85;     // identity threshold to reclassify INS as DUP (tandem duplication)
	double dup_window_factor = 2.0;   // search inserted seq in ref window of size +/- factor*svlen around breakpoint
	int dup_min_svlen = 50;           // skip DUP detection for very small INS (noise risk)
	int ambig_kmer_size = 15;         // k-mer size for ref-window abundance check
	int ambig_threshold = 5;          // median k-mer abundance >= threshold => mark as Ambig (FILTER)
	int ambig_min_window = 5000;      // minimum half-window (bp) for k-mer abundance counting
	int merge_pos_tol = 50;           // H1/H2 homozygous match: positional tolerance (bp)
	double merge_size_tol = 0.1;      // H1/H2 homozygous match: relative size tolerance (10%)
	// --- verify mode parameters ---
	std::string verify_vcf_path;      // input VCF (e.g. veritig --project output)
	int verify_flank = 500;           // flanking ref bp on each side of expected_seq
	double verify_gap_frac = 0.5;     // reject if max_indel or |span_delta| >= svlen * this
	int verify_min_indel_gap = 30;    // absolute floor for gap rejection (bp)
	double min_map_ratio = 0.85;
	double min_aln_identity = 0.85;
	std::string minimap_preset = "asm10";
	bool preset_set = false;          // true once --preset/-P is given explicitly
} parameters;

// Preset per mode: --project/--verify default to asm5, all other modes to asm10;
// --preset overrides both. Resolved here so one mode cannot alter another's preset.
inline std::string effective_preset(const parameters& params)
{
	if (!params.preset_set && (params.project || params.verify))
		return "asm5";
	return params.minimap_preset;
}

typedef struct _read
{
	std::string name;
	int freq = 0;
	int sv_count = 0;
	int ins_count = 0;
	int del_count = 0;
	int mapq = 0;
	int edit_dist = 0;
	int aln_score = 0;
	int svtig_size = 0;
	int total_aligned_bases = 0;
	double highest_map_ratio = 0.0;
	double aln_identity = 0.0;
	int haplo = 0;
	bool homo = false;
	double h1_mr = 0.0, h1_id = 0.0;  // best H1 alignment metrics (for concordant-on-both check)
	double h2_mr = 0.0, h2_id = 0.0;  // best H2 alignment metrics (for concordant-on-both check)
	std::string contig;
} Read;

struct PafRecord
{
	std::string read_name;
	std::string contig;
	int mapq = 0;
	int edit_dist = 0;
	int aln_score = 0;
	int svtig_size = 0;
	int sv_count = 0;
	int ins_count = 0;
	int del_count = 0;
	int aligned_bases = 0;
	double map_ratio = 0.0;
	double aln_identity = 0.0;
	bool is_primary = true;
};

int parse_paf_line(std::string& line, PafRecord& rec);
std::string shell_quote(const std::string& s);
double compute_veriscore(double map_ratio, double aln_identity, int mapq);
int decompose_cigar(std::string cigar, std::vector<int>& cigarLen, std::vector<char>& cigarOp);
void run_command(const std::string& cmd, const std::string& output_file);

#endif
