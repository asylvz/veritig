#ifndef __SVTIG_FILTER
#define __SVTIG_FILTER

#include <map>
#include <set>
#include "common.h"

class SvtigFilter
{
private:
	std::string paf_path;
	std::map<std::string, Read*> reads;

	void run_mapping(parameters& params);
	std::set<std::string> compute_passing(parameters& params);
	void write_filtered_fasta(std::set<std::string>& passing, parameters& params);

public:
	void run(parameters& params);
};

#endif
