#include <iostream>
#include <string>
#include <filesystem>
#include <fstream>
#include <getopt.h>
#include "cmdline.h"


int parse_command_line(int argc, char** argv, parameters& params)
{
	int index, o;

	static struct option long_options[] =
	{
		{"h1" , required_argument, NULL, 'i'},
		{"h2" , required_argument, NULL, 'j'},
		{"fasta" , required_argument, NULL, 'f'},
		{"help"   , no_argument, 0, 'h'},
		{"skip-mapping" , no_argument, NULL, 'm'},
		{"out" , required_argument, NULL, 'o'},
		{"phase" , no_argument, NULL, 'p'},
		{"sample" , required_argument, NULL, 's'},
		{"svtig1" , required_argument, NULL, 't'},
		{"svtig2" , required_argument, NULL, 'u'},
		{"concordance"   , no_argument, 0, 'v'},
		{"validate"   , no_argument, 0, 'y'},
		{"stats"   , no_argument, 0, 'q'},
		{"filter"  , no_argument, 0, 'r'},
		{"compare" , no_argument, 0, 'c'},
		{"threads" , required_argument, NULL, 'T'},
		{"min-map-ratio" , required_argument, NULL, 'M'},
		{"min-identity" , required_argument, NULL, 'I'},
		{"detailed" , no_argument, NULL, 'D'},
		{"preset" , required_argument, NULL, 'P'},

		{NULL, 0, NULL, 0}
	};

	while((o = getopt_long( argc, argv, "cf:hi:j:mo:pqrs:t:u:vyT:M:I:DP:", long_options, &index)) != -1)
	{
		switch(o)
		{
			case 'f':
				params.fasta = optarg;
				break;
			case 'i':
				params.haplo1_assembly_path = optarg;
				break;
			case 'j':
				params.haplo2_assembly_path = optarg;
				break;
			case 'm':
				params.skip_mapping = true;
				break;
			case 'o':
				params.output_path = optarg;
				break;
			case 'p':
				params.phase = true;
				break;
			case 's':
				params.sample_name = optarg;
				break;
			case 't':
				params.svtig1_path = optarg;
				break;
			case 'u':
				params.svtig2_path = optarg;
				break;
			case 'v':
				params.concordance = true;
				break;
			case 'y':
				params.validate = true;
				break;
			case 'q':
				params.stats = true;
				break;
			case 'r':
				params.filter = true;
				break;
			case 'c':
				params.compare = true;
				break;
			case 'T':
				params.threads = atoi(optarg);
				break;
			case 'M':
				params.min_map_ratio = atof(optarg);
				break;
			case 'I':
				params.min_aln_identity = atof(optarg);
				break;
			case 'D':
				params.detailed = true;
				break;
			case 'P':
				params.minimap_preset = optarg;
				break;
			case 'h':
				print_help();
				exit(0);
		}
	}
	if (params.concordance)
	{
		if((params.haplo1_assembly_path).empty())
		{
			std::cerr<<"[veritig] Please enter assembly file path for haplotype 1 using \"--h1\".\n";
			return RETURN_ERROR;
		}

		if((params.svtig1_path).empty())
		{
			std::cerr<<"[veritig] Please enter svtig path for the first haplotype using \"--svtig1\".\n";
			return RETURN_ERROR;
		}

		if (params.phase)
		{
			if((params.haplo2_assembly_path).empty())
			{
				std::cerr<<"[veritig] Please enter assembly file path for haplotype 2 using \"--h2\".\n";
				return RETURN_ERROR;
			}
			if((params.svtig2_path).empty())
			{
				std::cerr<<"[veritig] Please enter svtig path for the second haplotype using \"--svtig2\".\n";
				return RETURN_ERROR;
			}
		}
	}

	if (params.validate)
	{
		if((params.haplo1_assembly_path).empty())
		{
			std::cerr<<"[veritig] Please enter assembly file path for haplotype 1 using \"--h1\".\n";
			return RETURN_ERROR;
		}

		if((params.haplo2_assembly_path).empty())
		{
			std::cerr<<"[veritig] Please enter assembly file path for haplotype 2 using \"--h2\".\n";
			return RETURN_ERROR;
		}

		if((params.fasta).empty())
		{
			std::cerr<<"[veritig] Please enter SV sequences FASTA file using \"--fasta\".\n";
			return RETURN_ERROR;
		}
	}

	if (params.stats)
	{
		if ((params.svtig1_path).empty() && (params.fasta).empty())
		{
			std::cerr<<"[veritig] Please enter svtig FASTA file using \"--svtig1\" or \"--fasta\".\n";
			return RETURN_ERROR;
		}
	}

	if (params.compare)
	{
		if ((params.haplo1_assembly_path).empty())
		{
			std::cerr<<"[veritig] Please enter assembly file using \"--h1\".\n";
			return RETURN_ERROR;
		}
		if ((params.svtig1_path).empty())
		{
			std::cerr<<"[veritig] Please enter first svtig FASTA file using \"--svtig1\".\n";
			return RETURN_ERROR;
		}
		if ((params.svtig2_path).empty())
		{
			std::cerr<<"[veritig] Please enter second svtig FASTA file using \"--svtig2\".\n";
			return RETURN_ERROR;
		}
	}

	if (params.filter)
	{
		if ((params.haplo1_assembly_path).empty())
		{
			std::cerr<<"[veritig] Please enter assembly file using \"--h1\".\n";
			return RETURN_ERROR;
		}
		if ((params.svtig1_path).empty() && (params.fasta).empty())
		{
			std::cerr<<"[veritig] Please enter svtig FASTA file using \"--svtig1\" or \"--fasta\".\n";
			return RETURN_ERROR;
		}
	}

	if((params.output_path).empty())
	{
		std::string cwd = std::filesystem::current_path().string();
		params.log_path = cwd + "/veritig_results/";
	}
	else
		params.log_path = params.output_path + "veritig_results/";

	if((params.sample_name).empty())
		params.sample_name = "sample";

	return RETURN_SUCCESS;
}


void init_logs(parameters& params)
{
	std::cerr << "\n...veritig is running...\n";

	std::filesystem::create_directories(params.log_path);

	bool needs_mapping = (params.concordance || params.validate || params.filter || params.compare);
	if (needs_mapping)
		std::filesystem::create_directories(params.log_path + "paf/");

	if (params.detailed)
		std::filesystem::create_directories(params.log_path + "detailed/");

	std::cerr << "Output folder: " << params.log_path << "\n";
	if (needs_mapping)
	{
		std::cerr << "Concordance thresholds: map_ratio >= " << params.min_map_ratio
			<< ", identity >= " << params.min_aln_identity << "\n";
		std::cerr << "minimap2 preset: " << params.minimap_preset << "\n";
		std::cerr << "Threads: " << params.threads << "\n";
	}
}


void print_help()
{
	std::cerr << std::endl;
	std::cout << "veritig - svtig Verification Tool" << std::endl;
	std::cout << "\tVersion " << VERITIG_VERSION << ", Last update: " << VERITIG_UPDATE << "\n";
	std::cerr << std::endl;
	std::cerr << "Concordance analysis (--concordance)" << std::endl;
	std::cerr << "\t--h1                : Haplotype 1 assembly file" << std::endl;
	std::cerr << "\t--svtig1            : svtig FASTA file for haplotype 1" << std::endl;
	std::cerr << "\t--phase             : Enable phased mode (requires --h2, --svtig2)" << std::endl;
	std::cerr << "\t--h2                : Haplotype 2 assembly file" << std::endl;
	std::cerr << "\t--svtig2            : svtig FASTA file for haplotype 2" << std::endl;
	std::cerr << std::endl;
	std::cerr << "SV validation (--validate)" << std::endl;
	std::cerr << "\t--h1                : Haplotype 1 assembly file" << std::endl;
	std::cerr << "\t--h2                : Haplotype 2 assembly file" << std::endl;
	std::cerr << "\t--fasta             : SV sequences FASTA file" << std::endl;
	std::cerr << std::endl;
	std::cerr << "svtig statistics (--stats)" << std::endl;
	std::cerr << "\t--svtig1            : svtig FASTA file (or --fasta)" << std::endl;
	std::cerr << "\t--svtig2            : Second svtig FASTA file (optional)" << std::endl;
	std::cerr << std::endl;
	std::cerr << "svtig comparison (--compare)" << std::endl;
	std::cerr << "\t--h1                : Haplotype assembly file" << std::endl;
	std::cerr << "\t--svtig1            : First svtig FASTA file" << std::endl;
	std::cerr << "\t--svtig2            : Second svtig FASTA file" << std::endl;
	std::cerr << std::endl;
	std::cerr << "svtig filtering (--filter)" << std::endl;
	std::cerr << "\t--h1                : Haplotype assembly file" << std::endl;
	std::cerr << "\t--svtig1            : svtig FASTA file (or --fasta)" << std::endl;
	std::cerr << std::endl;
	std::cerr << "General options" << std::endl;
	std::cerr << "\t--preset (-P)       : minimap2 preset [asm10] (e.g., asm5, asm10, asm20)" << std::endl;
	std::cerr << "\t--threads (-T)      : Number of threads for minimap2 [16]" << std::endl;
	std::cerr << "\t--min-map-ratio (-M): Min mapping ratio for concordance [0.85]" << std::endl;
	std::cerr << "\t--min-identity (-I) : Min alignment identity for concordance [0.85]" << std::endl;
	std::cerr << "\t--detailed (-D)     : Write additional analysis files (summary, size bins, SV types)" << std::endl;
	std::cerr << "\t--skip-mapping (-m) : Skip minimap2 mapping, use existing PAF files" << std::endl;
	std::cerr << "\t--out (-o)          : Output folder path" << std::endl;
	std::cerr << "\t--sample (-s)       : Sample name" << std::endl;
	std::cerr << "\t--help              : Print this help menu" << std::endl;
	std::cerr << std::endl;
}
