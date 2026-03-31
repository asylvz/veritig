#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include "svtig_filter.h"


void SvtigFilter::run_mapping(parameters& params)
{
	std::string paf_dir = params.log_path + "paf/";
	this->paf_path = paf_dir + params.sample_name + ".filter.paf";
	std::string threads_str = std::to_string(params.threads);

	std::string fasta = params.svtig1_path.empty() ? params.fasta : params.svtig1_path;

	std::cerr << "  Mapping svtigs to assembly...";
	std::string minimap_cmd = "minimap2 -cx asm10 -t " + threads_str + " "
		+ params.haplo1_assembly_path + " " + fasta
		+ " --secondary=no > " + this->paf_path + " 2>/dev/null";
	run_command(minimap_cmd, this->paf_path);
	std::cerr << " done\n";
}


std::set<std::string> SvtigFilter::compute_passing(parameters& params)
{
	std::set<std::string> passing;
	std::string line;
	std::ifstream fp(this->paf_path);

	if (!fp.good())
	{
		std::cerr << "  Error opening '" << this->paf_path << "'\n";
		return passing;
	}

	while (getline(fp, line))
	{
		PafRecord rec;
		int sv_count = parse_paf_line(line, rec);
		if (sv_count == RETURN_ERROR)
			continue;

		auto it = reads.find(rec.read_name);
		if (it != reads.end())
		{
			if (!rec.is_primary)
			{
				it->second->total_aligned_bases += rec.aligned_bases;
				continue;
			}
			it->second->freq++;
			it->second->aln_identity += rec.aln_identity;
			it->second->total_aligned_bases += rec.aligned_bases;
			if (rec.map_ratio > it->second->highest_map_ratio)
				it->second->highest_map_ratio = rec.map_ratio;
		}
		else
		{
			if (!rec.is_primary)
				continue;

			Read *r = new Read();
			r->name = rec.read_name;
			r->freq = 1;
			r->svtig_size = rec.svtig_size;
			r->aln_identity = rec.aln_identity;
			r->highest_map_ratio = rec.map_ratio;
			r->total_aligned_bases = rec.aligned_bases;
			reads.insert(std::pair<std::string, Read*>(rec.read_name, r));
		}
	}

	for (auto it = reads.begin(); it != reads.end(); ++it)
	{
		Read* r = it->second;
		double combined_map_ratio = std::min(1.0, (double)r->total_aligned_bases / r->svtig_size);
		double avg_aln_identity = r->aln_identity / r->freq;

		if (combined_map_ratio >= params.min_map_ratio && avg_aln_identity >= params.min_aln_identity)
			passing.insert(it->first);
	}

	return passing;
}


void SvtigFilter::write_filtered_fasta(std::set<std::string>& passing, parameters& params)
{
	std::string fasta = params.svtig1_path.empty() ? params.fasta : params.svtig1_path;
	std::string output_path = params.log_path + params.sample_name + ".filtered.fa";

	std::ifstream fp(fasta);
	std::ofstream out(output_path);
	std::string line;
	bool writing = false;
	int total = 0, passed = 0;

	while (getline(fp, line))
	{
		if (line.empty())
			continue;

		if (line[0] == '>')
		{
			total++;
			std::string name = line.substr(1);
			// Trim whitespace and everything after first space
			size_t sp = name.find(' ');
			if (sp != std::string::npos)
				name = name.substr(0, sp);
			sp = name.find('\t');
			if (sp != std::string::npos)
				name = name.substr(0, sp);

			if (passing.count(name))
			{
				writing = true;
				passed++;
				out << line << "\n";
			}
			else
			{
				writing = false;
			}
		}
		else if (writing)
		{
			out << line << "\n";
		}
	}

	std::cerr << "  " << passed << "/" << total << " svtigs passed ("
		<< ((total > 0) ? (double)passed / total * 100 : 0) << "%)\n";
	std::cerr << "  Filtered FASTA written to " << output_path << "\n";
}


void SvtigFilter::run(parameters& params)
{
	std::cerr << "  Filtering svtigs (map_ratio >= " << params.min_map_ratio
		<< ", identity >= " << params.min_aln_identity << ")...\n";

	if (!params.skip_mapping)
		run_mapping(params);
	else
		this->paf_path = params.log_path + "paf/" + params.sample_name + ".filter.paf";

	std::set<std::string> passing = compute_passing(params);
	write_filtered_fasta(passing, params);

	for (auto it = reads.begin(); it != reads.end(); ++it)
		delete it->second;
	reads.clear();
}
