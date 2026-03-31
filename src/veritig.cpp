#include <iostream>
#include "cmdline.h"
#include "concordance.h"
#include "validate.h"
#include "svtig_stats.h"
#include "svtig_filter.h"
#include "svtig_compare.h"


int main(int argc, char** argv)
{
	parameters params = parameters();
	if (parse_command_line(argc, argv, params) != RETURN_SUCCESS)
		return RETURN_ERROR;

	init_logs(params);

	if (params.stats)
	{
		std::cerr<<"\nRunning svtig statistics...\n";
		SvtigStats s;
		s.run(params);
	}

	if (params.compare)
	{
		std::cerr<<"\nComparing svtig sets...\n";
		SvtigCompare c;
		c.run(params);
	}

	if (params.filter)
	{
		std::cerr<<"\nFiltering svtigs...\n";
		SvtigFilter f;
		f.run(params);
	}

	if (params.concordance)
	{
		std::cerr<<"\nRunning concordance analysis...\n";
		Concordance c;
		c.run(params);
	}

	if (params.validate)
	{
		std::cerr<<"\nRunning SV validation...\n";
		Validate v;
		v.run(params);
	}

	std::cerr<<"\nDone.\n";

	return RETURN_SUCCESS;
}
