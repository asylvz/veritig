#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <cstdlib>
#include "common.h"


void run_command(const std::string& cmd, const std::string& output_file)
{
	int ret = system(cmd.c_str());
	if (ret != 0)
	{
		std::cerr << "[VERITIG ERROR] Command failed (exit code " << ret << "): " << cmd << "\n";
		exit(EXIT_FAILURE);
	}
	if (!output_file.empty() && (!std::filesystem::exists(output_file) || std::filesystem::file_size(output_file) == 0))
	{
		std::cerr << "[VERITIG ERROR] Expected output file missing or empty: " << output_file << "\n";
		exit(EXIT_FAILURE);
	}
}


int decompose_cigar(std::string cigar, std::vector<int>& cigarLen, std::vector<char>& cigarOp)
{
	int count = 0;
	int num = 0;

	for (char c : cigar)
	{
		if (std::isdigit(c))
		{
			num = num * 10 + (c - '0');
		}
		else
		{
			cigarLen.push_back(num);
			cigarOp.push_back(c);
			num = 0;
			count++;
		}
	}

	return count;
}


double compute_veriscore(double map_ratio, double aln_identity, int mapq)
{
	// Floor=3 (uninformative prior, P≈0.5), cap=60 (numerical stability)
	int effective_mapq = std::min(std::max(mapq, 3), 60);

	// Phred-scale to probability: P_correct = 1 - 10^(-MAPQ/10)
	double p_correct = 1.0 - std::pow(10.0, -effective_mapq / 10.0);

	// Geometric mean of three quality components
	if (map_ratio <= 0.0 || aln_identity <= 0.0)
		return 0.0;

	return std::cbrt(map_ratio * aln_identity * p_correct);
}


int parse_paf_line(std::string& line, PafRecord& rec)
{
	std::vector<std::string> tokens;
	std::string tmp_str;
	std::stringstream s(line);

	while (getline(s, tmp_str, '\t'))
		tokens.push_back(tmp_str);

	if (tokens.size() < 12)
		return RETURN_ERROR;

	rec.read_name = tokens[0];
	rec.svtig_size = stoi(tokens[1]);
	rec.contig = tokens[5];
	rec.mapq = stoi(tokens[11]);

	int query_start = stoi(tokens[2]);
	int query_end = stoi(tokens[3]);

	if (query_start == 0 && query_end == 0)
		return RETURN_ERROR;

	rec.aligned_bases = query_end - query_start;
	rec.map_ratio = (double)rec.aligned_bases / rec.svtig_size;
	rec.aln_identity = (double)stoi(tokens[9]) / stoi(tokens[10]);

	rec.sv_count = 0;
	rec.ins_count = 0;
	rec.del_count = 0;
	rec.is_primary = true;

	for (auto& tok : tokens)
	{
		if (tok.substr(0, 5) == "tp:A:")
		{
			rec.is_primary = (tok[5] == 'P');
		}
		else if (tok.substr(0, 5) == "cg:Z:")
		{
			std::vector<int> cigarLen;
			std::vector<char> cigarOp;
			std::string cigar = tok.substr(5);
			int cigar_cnt = decompose_cigar(cigar, cigarLen, cigarOp);
			for (int c = 0; c < cigar_cnt; c++)
			{
				if (cigarLen[c] > MINSVSIZE)
				{
					if (cigarOp[c] == INSERTION)
					{
						rec.ins_count++;
						rec.sv_count++;
					}
					else if (cigarOp[c] == DELETION)
					{
						rec.del_count++;
						rec.sv_count++;
					}
				}
			}
		}
		else if (tok.substr(0, 5) == "NM:i:")
			rec.edit_dist = stoi(tok.substr(5));
		else if (tok.substr(0, 5) == "AS:i:")
			rec.aln_score = stoi(tok.substr(5));
	}

	return rec.sv_count;
}
