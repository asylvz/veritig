#ifndef __SVTIG_COMPARE
#define __SVTIG_COMPARE

#include <map>
#include <vector>
#include "common.h"

class SvtigCompare
{
private:
	struct SetStats
	{
		int total_input = 0;
		int total_mapped = 0;
		int total_unmapped = 0;
		int total_concordant = 0;
		double concordance_rate = 0.0;
		double avg_map_ratio = 0.0;
		double avg_aln_identity = 0.0;
		double avg_mapq = 0.0;
		double avg_svtig_size = 0.0;
		int total_ins = 0;
		int total_del = 0;
		int size_50_100 = 0;
		int size_100_500 = 0;
		int size_500_1000 = 0;
		int size_1000_10000 = 0;
		int size_10000_plus = 0;
		int conc_50_100 = 0;
		int conc_100_500 = 0;
		int conc_500_1000 = 0;
		int conc_1000_10000 = 0;
		int conc_10000_plus = 0;
	};

	std::string run_mapping(std::string fasta, std::string label, parameters& params);
	SetStats compute_stats(std::string paf_path, int total_input, parameters& params);
	int count_seqs(std::string fasta_path);
	void write_report(SetStats& s1, SetStats& s2, parameters& params);

public:
	void run(parameters& params);
};

#endif
