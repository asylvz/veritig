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
		{"project" , no_argument, NULL, 'w'},
		{"ref" , required_argument, NULL, 'R'},
		{"min-svlen" , required_argument, NULL, 'L'},
		{"min-mapq" , required_argument, NULL, 'Q'},
		{"cluster-window" , required_argument, NULL, 'W'},
		{"cluster-size-ratio" , required_argument, NULL, 'Z'},
		{"dedup-pos-tol" , required_argument, NULL, 'X'},
		{"dedup-size-tol" , required_argument, NULL, 'Y'},
		{"min-split-svlen" , required_argument, NULL, 'A'},
		{"max-split-distance" , required_argument, NULL, 'B'},
		{"dup-similarity" , required_argument, NULL, 'C'},
		{"dup-window-factor" , required_argument, NULL, 'E'},
		{"ambig-kmer-size" , required_argument, NULL, 'F'},
		{"ambig-threshold" , required_argument, NULL, 'G'},
		{"ambig-min-window" , required_argument, NULL, 'K'},
		{"merge-pos-tol" , required_argument, NULL, 'N'},
		{"merge-size-tol" , required_argument, NULL, 'O'},
		{"svtig3" , required_argument, NULL, 'U'},
		{"verify" , no_argument, NULL, 'V'},
		{"vcf" , required_argument, NULL, 'a'},
		{"gap-frac" , required_argument, NULL, 'b'},
		{"min-indel-gap" , required_argument, NULL, 'd'},
		{"flank" , required_argument, NULL, 'e'},

		{NULL, 0, NULL, 0}
	};

	while((o = getopt_long( argc, argv, "a:b:cd:e:f:hi:j:mo:pqrs:t:u:vyT:M:I:DP:wR:L:Q:W:Z:X:Y:A:B:C:E:F:G:K:N:O:U:V", long_options, &index)) != -1)
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
			case 'w':
				params.project = true;
				break;
			case 'R':
				params.reference_path = optarg;
				break;
			case 'L':
				params.min_svlen = atoi(optarg);
				break;
			case 'Q':
				params.min_mapq = atoi(optarg);
				break;
			case 'W':
				params.cluster_window = atoi(optarg);
				break;
			case 'Z':
				params.cluster_size_ratio = atof(optarg);
				break;
			case 'X':
				params.dedup_pos_tol = atoi(optarg);
				break;
			case 'Y':
				params.dedup_size_tol = atof(optarg);
				break;
			case 'A':
				params.min_split_svlen = atoi(optarg);
				break;
			case 'B':
				params.max_split_distance = atoi(optarg);
				break;
			case 'C':
				params.dup_similarity = atof(optarg);
				break;
			case 'E':
				params.dup_window_factor = atof(optarg);
				break;
			case 'F':
				params.ambig_kmer_size = atoi(optarg);
				break;
			case 'G':
				params.ambig_threshold = atoi(optarg);
				break;
			case 'K':
				params.ambig_min_window = atoi(optarg);
				break;
			case 'N':
				params.merge_pos_tol = atoi(optarg);
				break;
			case 'O':
				params.merge_size_tol = atof(optarg);
				break;
			case 'U':
				params.svtig3_path = optarg;
				break;
			case 'V':
				params.verify = true;
				break;
			case 'a':
				params.verify_vcf_path = optarg;
				break;
			case 'b':
				params.verify_gap_frac = atof(optarg);
				break;
			case 'd':
				params.verify_min_indel_gap = atoi(optarg);
				break;
			case 'e':
				params.verify_flank = atoi(optarg);
				break;
			case 'h':
				print_help();
				exit(0);
		}
	}

	if (!params.concordance && !params.validate && !params.stats
		&& !params.filter && !params.compare && !params.project && !params.verify)
	{
		std::cerr << "[veritig] No mode specified. Use one of:\n"
			<< "    --stats, --concordance, --validate, --filter, --compare, --project, --verify\n"
			<< "Run 'veritig --help' for full usage.\n";
		return RETURN_ERROR;
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

	if (params.project)
	{
		if ((params.reference_path).empty())
		{
			std::cerr<<"[veritig] Please enter linear reference FASTA file using \"--ref\".\n";
			return RETURN_ERROR;
		}
		if ((params.svtig1_path).empty())
		{
			std::cerr<<"[veritig] Please enter svtig FASTA file using \"--svtig1\" (and optionally --svtig2 for phased input).\n";
			return RETURN_ERROR;
		}
		// --svtig2 is optional: omit for unphased single-FASTA SVarp output (no --phase)
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

	bool needs_mapping = (params.concordance || params.validate || params.filter || params.compare || params.project);
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
	std::cerr << "svtig projection (--project)" << std::endl;
	std::cerr << "\t--ref (-R)          : Linear reference FASTA file" << std::endl;
	std::cerr << "\t--svtig1            : svtig FASTA file (H1 if phased, or unified output if unphased)" << std::endl;
	std::cerr << "\t--svtig2            : H2 svtig FASTA file (optional; phased mode if provided)" << std::endl;
	std::cerr << "\t--svtig3 (-U)       : SVarp untagged (unphased) svtigs (FASTA, optional)" << std::endl;
	std::cerr << "\t--min-svlen (-L)    : Minimum SV length [50]" << std::endl;
	std::cerr << "\t--min-mapq (-Q)     : Minimum svtig alignment MAPQ [0]" << std::endl;
	std::cerr << "\t--dup-similarity (-C) : Identity threshold to reclassify INS as DUP [0.85]" << std::endl;
	std::cerr << std::endl;
	std::cerr << "svtig VCF verification (--verify)" << std::endl;
	std::cerr << "\t--vcf (-a)          : Input VCF for verification" << std::endl;
	std::cerr << "\t--ref (-R)          : Linear reference FASTA file (matching the VCF coordinates)" << std::endl;
	std::cerr << "\t--h1                : Haplotype 1 assembly file" << std::endl;
	std::cerr << "\t--h2                : Haplotype 2 assembly file (optional)" << std::endl;
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
