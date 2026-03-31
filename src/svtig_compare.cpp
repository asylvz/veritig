#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <map>
#include "svtig_compare.h"


int SvtigCompare::count_seqs(std::string fasta_path)
{
	int count = 0;
	std::string line;
	std::ifstream fp(fasta_path);
	while (getline(fp, line))
	{
		if (!line.empty() && line[0] == '>')
			count++;
	}
	return count;
}


std::string SvtigCompare::run_mapping(std::string fasta, std::string label, parameters& params)
{
	std::string paf_dir = params.log_path + "paf/";
	std::string paf_path = paf_dir + params.sample_name + ".compare." + label + ".paf";
	std::string threads_str = std::to_string(params.threads);

	std::cerr << "  Mapping " << label << " to assembly...";
	std::string cmd = "minimap2 -cx asm10 -t " + threads_str + " "
		+ params.haplo1_assembly_path + " " + fasta
		+ " --secondary=no > " + paf_path + " 2>/dev/null";
	run_command(cmd, paf_path);
	std::cerr << " done\n";

	return paf_path;
}


SvtigCompare::SetStats SvtigCompare::compute_stats(std::string paf_path, int total_input, parameters& params)
{
	SetStats s;
	s.total_input = total_input;

	std::map<std::string, Read*> reads;
	std::string line;
	std::ifstream fp(paf_path);

	if (!fp.good())
	{
		std::cerr << "  Error opening '" << paf_path << "'\n";
		return s;
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
				it->second->ins_count += rec.ins_count;
				it->second->del_count += rec.del_count;
				continue;
			}
			it->second->freq++;
			it->second->sv_count += sv_count;
			it->second->ins_count += rec.ins_count;
			it->second->del_count += rec.del_count;
			it->second->mapq += rec.mapq;
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
			r->sv_count = sv_count;
			r->ins_count = rec.ins_count;
			r->del_count = rec.del_count;
			r->mapq = rec.mapq;
			r->highest_map_ratio = rec.map_ratio;
			r->aln_identity = rec.aln_identity;
			r->total_aligned_bases = rec.aligned_bases;
			reads.insert(std::pair<std::string, Read*>(rec.read_name, r));
		}
	}

	s.total_mapped = (int)reads.size();
	s.total_unmapped = s.total_input - s.total_mapped;

	double sum_map_ratio = 0, sum_identity = 0, sum_mapq = 0;
	long sum_size = 0;

	for (auto it = reads.begin(); it != reads.end(); ++it)
	{
		Read* r = it->second;
		double combined_map_ratio = std::min(1.0, (double)r->total_aligned_bases / r->svtig_size);
		double avg_identity = r->aln_identity / r->freq;

		bool concordant = (combined_map_ratio >= params.min_map_ratio && avg_identity >= params.min_aln_identity);
		if (concordant)
			s.total_concordant++;

		sum_map_ratio += combined_map_ratio;
		sum_identity += avg_identity;
		sum_mapq += (double)r->mapq / r->freq;
		sum_size += r->svtig_size;

		if (r->ins_count > 0) s.total_ins++;
		if (r->del_count > 0) s.total_del++;

		int sz = r->svtig_size;
		if (sz >= 50 && sz < 100) { s.size_50_100++; if (concordant) s.conc_50_100++; }
		else if (sz >= 100 && sz < 500) { s.size_100_500++; if (concordant) s.conc_100_500++; }
		else if (sz >= 500 && sz < 1000) { s.size_500_1000++; if (concordant) s.conc_500_1000++; }
		else if (sz >= 1000 && sz < 10000) { s.size_1000_10000++; if (concordant) s.conc_1000_10000++; }
		else if (sz >= 10000) { s.size_10000_plus++; if (concordant) s.conc_10000_plus++; }
	}

	if (s.total_mapped > 0)
	{
		s.concordance_rate = (double)s.total_concordant / s.total_mapped;
		s.avg_map_ratio = sum_map_ratio / s.total_mapped;
		s.avg_aln_identity = sum_identity / s.total_mapped;
		s.avg_mapq = sum_mapq / s.total_mapped;
		s.avg_svtig_size = (double)sum_size / s.total_mapped;
	}

	for (auto it = reads.begin(); it != reads.end(); ++it)
		delete it->second;

	return s;
}


void SvtigCompare::write_report(SetStats& s1, SetStats& s2, parameters& params)
{
	std::string report_path = params.log_path + params.sample_name + ".compare.tsv";
	std::ofstream fp(report_path);

	fp << "metric\tsvtig1\tsvtig2\n";
	fp << "total_input\t" << s1.total_input << "\t" << s2.total_input << "\n";
	fp << "total_mapped\t" << s1.total_mapped << "\t" << s2.total_mapped << "\n";
	fp << "total_unmapped\t" << s1.total_unmapped << "\t" << s2.total_unmapped << "\n";
	fp << "total_concordant\t" << s1.total_concordant << "\t" << s2.total_concordant << "\n";
	fp << "concordance_rate\t" << s1.concordance_rate << "\t" << s2.concordance_rate << "\n";
	fp << "avg_map_ratio\t" << s1.avg_map_ratio << "\t" << s2.avg_map_ratio << "\n";
	fp << "avg_aln_identity\t" << s1.avg_aln_identity << "\t" << s2.avg_aln_identity << "\n";
	fp << "avg_mapq\t" << s1.avg_mapq << "\t" << s2.avg_mapq << "\n";
	fp << "avg_svtig_size\t" << s1.avg_svtig_size << "\t" << s2.avg_svtig_size << "\n";
	fp << "total_with_ins\t" << s1.total_ins << "\t" << s2.total_ins << "\n";
	fp << "total_with_del\t" << s1.total_del << "\t" << s2.total_del << "\n";

	auto pct = [](int conc, int total) -> double {
		return (total > 0) ? (double)conc / total * 100 : 0;
	};

	fp << "size_50_100\t" << s1.size_50_100 << "\t" << s2.size_50_100 << "\n";
	fp << "conc_50_100_pct\t" << pct(s1.conc_50_100, s1.size_50_100) << "\t" << pct(s2.conc_50_100, s2.size_50_100) << "\n";
	fp << "size_100_500\t" << s1.size_100_500 << "\t" << s2.size_100_500 << "\n";
	fp << "conc_100_500_pct\t" << pct(s1.conc_100_500, s1.size_100_500) << "\t" << pct(s2.conc_100_500, s2.size_100_500) << "\n";
	fp << "size_500_1000\t" << s1.size_500_1000 << "\t" << s2.size_500_1000 << "\n";
	fp << "conc_500_1000_pct\t" << pct(s1.conc_500_1000, s1.size_500_1000) << "\t" << pct(s2.conc_500_1000, s2.size_500_1000) << "\n";
	fp << "size_1000_10000\t" << s1.size_1000_10000 << "\t" << s2.size_1000_10000 << "\n";
	fp << "conc_1000_10000_pct\t" << pct(s1.conc_1000_10000, s1.size_1000_10000) << "\t" << pct(s2.conc_1000_10000, s2.size_1000_10000) << "\n";
	fp << "size_10000_plus\t" << s1.size_10000_plus << "\t" << s2.size_10000_plus << "\n";
	fp << "conc_10000_plus_pct\t" << pct(s1.conc_10000_plus, s1.size_10000_plus) << "\t" << pct(s2.conc_10000_plus, s2.size_10000_plus) << "\n";
}


void SvtigCompare::run(parameters& params)
{
	int count1 = count_seqs(params.svtig1_path);
	int count2 = count_seqs(params.svtig2_path);
	std::cerr << "  SVtig counts: svtig1=" << count1 << ", svtig2=" << count2 << "\n";

	std::string paf1, paf2;

	if (!params.skip_mapping)
	{
		paf1 = run_mapping(params.svtig1_path, "svtig1", params);
		paf2 = run_mapping(params.svtig2_path, "svtig2", params);
	}
	else
	{
		std::string paf_dir = params.log_path + "paf/";
		paf1 = paf_dir + params.sample_name + ".compare.svtig1.paf";
		paf2 = paf_dir + params.sample_name + ".compare.svtig2.paf";
	}

	std::cerr << "  Computing statistics...\n";
	SetStats s1 = compute_stats(paf1, count1, params);
	SetStats s2 = compute_stats(paf2, count2, params);

	write_report(s1, s2, params);

	std::cerr << "  svtig1: " << s1.total_concordant << "/" << s1.total_mapped
		<< " concordant (" << s1.concordance_rate * 100 << "%)\n";
	std::cerr << "  svtig2: " << s2.total_concordant << "/" << s2.total_mapped
		<< " concordant (" << s2.concordance_rate * 100 << "%)\n";
	std::cerr << "  Results written to " << params.log_path << "\n";
}
