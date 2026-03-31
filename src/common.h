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
	std::string log_path;
	std::string haplo1_assembly_path;
	std::string haplo2_assembly_path;
	std::string sample_name;
	bool skip_mapping = false;
	bool concordance = false;
	bool validate = false;
	bool stats = false;
	bool filter = false;
	bool compare = false;
	bool phase = false;
	bool detailed = false;
	int threads = 16;
	double min_map_ratio = 0.85;
	double min_aln_identity = 0.85;
} parameters;

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
double compute_veriscore(double map_ratio, double aln_identity, int mapq);
int decompose_cigar(std::string cigar, std::vector<int>& cigarLen, std::vector<char>& cigarOp);
void run_command(const std::string& cmd, const std::string& output_file);

#endif
