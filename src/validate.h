#ifndef __VALIDATE
#define __VALIDATE

#include <map>
#include <vector>
#include "common.h"

struct SvResult
{
	std::string name;
	std::string contig;
	std::string haplotype;
	int sv_size = 0;
	int ins_count = 0;
	int del_count = 0;
	int sv_count = 0;
	double map_ratio = 0.0;
	double aln_identity = 0.0;
	int mapq = 0;
	int edit_dist = 0;
	int aln_score = 0;
	std::string sv_type;
	bool concordant = false;
	double veriscore = 0.0;
};

class Validate
{
private:
	std::string paf_H1_sv_path;
	std::string paf_H2_sv_path;

	int sv_count = 0;
	std::map<std::string, int> sv_sizes;
	std::map<std::string, Read*> reads;
	std::vector<SvResult> results;

	void count_svs(parameters& params);
	void run_mapping(parameters& params);
	void compute_stats(parameters& params);
	void write_report(parameters& params);
	void write_detailed(std::string filename, std::string haplotype, parameters& params);

public:
	void run(parameters& params);
};

#endif
