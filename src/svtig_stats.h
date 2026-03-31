#ifndef __SVTIG_STATS
#define __SVTIG_STATS

#include <vector>
#include "common.h"

class SvtigStats
{
private:
	struct SeqStats
	{
		int count = 0;
		long total_bases = 0;
		int min_size = 0;
		int max_size = 0;
		double mean_size = 0.0;
		int median_size = 0;
		int n50 = 0;
		double gc_content = 0.0;
		double n_content = 0.0;
		int size_50_100 = 0;
		int size_100_500 = 0;
		int size_500_1000 = 0;
		int size_1000_10000 = 0;
		int size_10000_plus = 0;
	};

	SeqStats compute(std::string fasta_path);
	void write_report(SeqStats& s1, SeqStats& s2, SeqStats& combined, bool has_svtig2, parameters& params);

public:
	void run(parameters& params);
};

#endif
