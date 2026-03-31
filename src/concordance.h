#ifndef __CONCORDANCE
#define __CONCORDANCE

#include <map>
#include <vector>
#include "common.h"

struct SvtigResult
{
	std::string name;
	int svtig_size = 0;
	double highest_map_ratio = 0.0;
	double combined_map_ratio = 0.0;
	double aln_identity = 0.0;
	double mapq = 0.0;
	int edit_dist = 0;
	double aln_score = 0.0;
	int ins_count = 0;
	int del_count = 0;
	int sv_count = 0;
	std::string sv_type;
	bool concordant = false;
	double veriscore = 0.0;
	std::string haplotype;
};

class Concordance
{
private:
	std::string paf_H1_svtig1_path;
	std::string paf_H1_svtig2_path;
	std::string paf_H2_svtig1_path;
	std::string paf_H2_svtig2_path;

	int h1_count = 0;
	int h2_count = 0;
	std::map<std::string, int> svtig_sizes;
	std::map<std::string, Read*> reads;
	std::vector<SvtigResult> results;

	void count_svtigs(parameters& params);
	void run_mapping(parameters& params);
	int mapping_stats(std::string filename, int total_svtig_count, std::string haplotype, parameters& params);
	void write_report(parameters& params);
	void write_detailed(std::string haplotype_label, int total_svtig_count, int total_concordant, parameters& params);

public:
	void run(parameters& params);
};

#endif
