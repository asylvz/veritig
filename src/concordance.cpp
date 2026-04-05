#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <climits>
#include <set>
#include "concordance.h"


void Concordance::count_svtigs(parameters& params)
{
	std::string line;
	std::string current_name;
	int current_bases = 0;

	std::ifstream fp(params.svtig1_path);
	while (getline(fp, line))
	{
		if (line[0] == '>')
		{
			if (!current_name.empty())
				svtig_sizes[current_name] = current_bases;
			current_name = line.substr(1);
			size_t space_pos = current_name.find(' ');
			if (space_pos != std::string::npos)
				current_name = current_name.substr(0, space_pos);
			current_bases = 0;
			this->h1_count++;
		}
		else
		{
			current_bases += (int)line.size();
		}
	}
	if (!current_name.empty())
		svtig_sizes[current_name] = current_bases;

	if (params.phase)
	{
		current_name.clear();
		current_bases = 0;

		std::ifstream fp2(params.svtig2_path);
		while (getline(fp2, line))
		{
			if (line[0] == '>')
			{
				if (!current_name.empty())
					svtig_sizes[current_name] = current_bases;
				current_name = line.substr(1);
				size_t space_pos = current_name.find(' ');
				if (space_pos != std::string::npos)
					current_name = current_name.substr(0, space_pos);
				current_bases = 0;
				this->h2_count++;
			}
			else
			{
				current_bases += (int)line.size();
			}
		}
		if (!current_name.empty())
			svtig_sizes[current_name] = current_bases;
	}
}

void Concordance::run_mapping(parameters& params)
{
	std::string minimap_cmd;
	std::string threads_str = std::to_string(params.threads);

	std::string paf_dir = params.log_path + "paf/";
	this->paf_H1_svtig1_path = paf_dir + params.sample_name + ".svtig.H1.paf";
	this->paf_H2_svtig1_path = paf_dir + params.sample_name + ".svtig.H2.paf";

	std::cerr << "  Mapping svtig1 to H1...";
	minimap_cmd = "minimap2 -cx " + params.minimap_preset + " -t " + threads_str + " " + params.haplo1_assembly_path + " " + params.svtig1_path + " --secondary=no > " + this->paf_H1_svtig1_path + " 2>/dev/null";
	run_command(minimap_cmd, this->paf_H1_svtig1_path);
	std::cerr << " done\n";

	std::cerr << "  Mapping svtig1 to H2...";
	minimap_cmd = "minimap2 -cx " + params.minimap_preset + " -t " + threads_str + " " + params.haplo2_assembly_path + " " + params.svtig1_path + " --secondary=no > " + this->paf_H2_svtig1_path + " 2>/dev/null";
	run_command(minimap_cmd, this->paf_H2_svtig1_path);
	std::cerr << " done\n";

	if (params.phase)
	{
		this->paf_H1_svtig2_path = paf_dir + params.sample_name + ".svtig2.H1.paf";
		this->paf_H2_svtig2_path = paf_dir + params.sample_name + ".svtig2.H2.paf";

		std::cerr << "  Mapping svtig2 to H1...";
		minimap_cmd = "minimap2 -cx " + params.minimap_preset + " -t " + threads_str + " " + params.haplo1_assembly_path + " " + params.svtig2_path + " --secondary=no > " + this->paf_H1_svtig2_path + " 2>/dev/null";
		run_command(minimap_cmd, this->paf_H1_svtig2_path);
		std::cerr << " done\n";

		std::cerr << "  Mapping svtig2 to H2...";
		minimap_cmd = "minimap2 -cx " + params.minimap_preset + " -t " + threads_str + " " + params.haplo2_assembly_path + " " + params.svtig2_path + " --secondary=no > " + this->paf_H2_svtig2_path + " 2>/dev/null";
		run_command(minimap_cmd, this->paf_H2_svtig2_path);
		std::cerr << " done\n";
	}
}


int Concordance::mapping_stats(std::string filename, int total_svtig_count, std::string haplotype, parameters& params)
{
	std::string line;
	int total_supplementary = 0, total_uniq = 0;
	std::ifstream fp(filename);

	if (!fp.good())
	{
		std::cerr << "Error opening '" << filename << std::endl;
		return RETURN_ERROR;
	}

	while (getline(fp, line))
	{
		PafRecord rec;
		int sv_count = parse_paf_line(line, rec);
		if (sv_count == RETURN_ERROR)
			continue;

		if (!rec.is_primary)
		{
			total_supplementary++;
			auto it = reads.find(rec.read_name);
			if (it != reads.end())
			{
				it->second->total_aligned_bases += rec.aligned_bases;
				it->second->sv_count += sv_count;
				it->second->ins_count += rec.ins_count;
				it->second->del_count += rec.del_count;
			}
			continue;
		}

		auto it = reads.find(rec.read_name);
		if (it != reads.end())
		{
			it->second->freq++;
			it->second->sv_count += sv_count;
			it->second->ins_count += rec.ins_count;
			it->second->del_count += rec.del_count;
			it->second->mapq += rec.mapq;
			it->second->edit_dist += rec.edit_dist;
			it->second->aln_score += rec.aln_score;
			it->second->aln_identity += rec.aln_identity;
			it->second->total_aligned_bases += rec.aligned_bases;
			if (rec.map_ratio > it->second->highest_map_ratio)
				it->second->highest_map_ratio = rec.map_ratio;
		}
		else
		{
			total_uniq++;

			Read *r = new Read();
			r->name = rec.read_name;
			r->freq = 1;
			r->sv_count = sv_count;
			r->ins_count = rec.ins_count;
			r->del_count = rec.del_count;
			r->svtig_size = rec.svtig_size;
			r->edit_dist = rec.edit_dist;
			r->aln_score = rec.aln_score;
			r->mapq = rec.mapq;
			r->highest_map_ratio = rec.map_ratio;
			r->aln_identity = rec.aln_identity;
			r->total_aligned_bases = rec.aligned_bases;
			reads.insert(std::pair<std::string, Read*>(rec.read_name, r));
		}

	}

	// Compute per-read stats and collect results
	int total_concordant = 0;

	for (auto it = reads.begin(); it != reads.end(); ++it)
	{
		Read* r = it->second;
		double combined_map_ratio = std::min(1.0, (double)r->total_aligned_bases / r->svtig_size);
		double avg_mapq = (double)r->mapq / r->freq;
		double avg_aln_score = (double)r->aln_score / r->freq;
		double avg_aln_identity = r->aln_identity / r->freq;

		bool is_concordant = (combined_map_ratio >= params.min_map_ratio && avg_aln_identity >= params.min_aln_identity);
		if (is_concordant)
			total_concordant++;

		std::string sv_type;
		if (r->ins_count > 0 && r->del_count == 0) sv_type = "INS";
		else if (r->del_count > 0 && r->ins_count == 0) sv_type = "DEL";
		else if (r->ins_count > 0 && r->del_count > 0) sv_type = "MIXED";
		else sv_type = "NONE";

		SvtigResult res;
		res.name = it->first;
		res.svtig_size = r->svtig_size;
		res.highest_map_ratio = r->highest_map_ratio;
		res.combined_map_ratio = combined_map_ratio;
		res.aln_identity = avg_aln_identity;
		res.mapq = avg_mapq;
		res.edit_dist = r->edit_dist;
		res.aln_score = avg_aln_score;
		res.ins_count = r->ins_count;
		res.del_count = r->del_count;
		res.sv_count = r->sv_count;
		res.sv_type = sv_type;
		res.concordant = is_concordant;
		res.veriscore = compute_veriscore(combined_map_ratio, avg_aln_identity, (int)avg_mapq);
		res.haplotype = haplotype;
		results.push_back(res);
	}

	double concordance_rate = (total_uniq > 0) ? (double)total_concordant / total_uniq : 0;

	std::cerr << "  " << filename << ": "
		<< total_uniq << "/" << total_svtig_count << " mapped, "
		<< total_supplementary << " supplementary, "
		<< total_concordant << " concordant (" << concordance_rate * 100 << "%)\n";

	if (params.detailed)
		write_detailed(haplotype, total_svtig_count, total_concordant, params);

	for (auto it = reads.begin(); it != reads.end(); ++it)
		delete it->second;
	reads.clear();

	return total_concordant;
}


void Concordance::write_detailed(std::string haplotype_label, int total_svtig_count, int total_concordant, parameters& params)
{
	std::string detailed_dir = params.log_path + "detailed/";
	std::string base = detailed_dir + params.sample_name + ".svtig." + haplotype_label;

	std::string details_path = base + ".details.tsv";
	std::ofstream fp_details(details_path);
	fp_details << "read_name\tsvtig_size\tmap_ratio\tcombined_map_ratio\taln_identity\tmapq\tedit_dist\taln_score\tins_count\tdel_count\tsv_count\tsv_type\tconcordant\tveriscore\n";

	int total_read_count = 0, total_sv_count = 0;
	int total_ins = 0, total_del = 0;
	long total_svtig_size = 0;
	double total_map_ratio = 0.0, total_mapq = 0, total_aln_score = 0, total_edit_dist = 0, total_aln_identity = 0.0;
	int total_unmapped = 0;

	struct SizeBin { int min_size; int max_size; std::string label; int total; int concordant; };
	std::vector<SizeBin> size_bins = {
		{50, 100, "50-100", 0, 0},
		{100, 500, "100-500", 0, 0},
		{500, 1000, "500-1000", 0, 0},
		{1000, 10000, "1000-10000", 0, 0},
		{10000, INT_MAX, "10000+", 0, 0}
	};

	struct TypeBin { std::string label; int total; int concordant; };
	std::vector<TypeBin> type_bins = {
		{"INS", 0, 0},
		{"DEL", 0, 0},
		{"MIXED", 0, 0},
		{"NONE", 0, 0}
	};

	// Find the results batch for this PAF (contiguous at the end, same haplotype)
	int batch_size = 0;
	for (int i = (int)results.size() - 1; i >= 0; i--)
	{
		if (results[i].haplotype == results.back().haplotype)
			batch_size++;
		else
			break;
	}
	size_t batch_start = results.size() - batch_size;

	for (size_t i = batch_start; i < results.size(); i++)
	{
		auto& res = results[i];
		total_read_count++;
		total_sv_count += res.sv_count;
		total_ins += res.ins_count;
		total_del += res.del_count;
		total_map_ratio += res.combined_map_ratio;
		total_mapq += res.mapq;
		total_aln_score += res.aln_score;
		total_aln_identity += res.aln_identity;
		total_edit_dist += res.edit_dist;
		total_svtig_size += res.svtig_size;

		for (auto& bin : size_bins)
		{
			if (res.svtig_size >= bin.min_size && res.svtig_size < bin.max_size)
			{
				bin.total++;
				if (res.concordant) bin.concordant++;
				break;
			}
		}

		for (auto& tb : type_bins)
		{
			if (tb.label == res.sv_type)
			{
				tb.total++;
				if (res.concordant) tb.concordant++;
				break;
			}
		}

		fp_details << res.name
			<< "\t" << res.svtig_size
			<< "\t" << res.highest_map_ratio
			<< "\t" << res.combined_map_ratio
			<< "\t" << res.aln_identity
			<< "\t" << res.mapq
			<< "\t" << res.edit_dist
			<< "\t" << res.aln_score
			<< "\t" << res.ins_count
			<< "\t" << res.del_count
			<< "\t" << res.sv_count
			<< "\t" << res.sv_type
			<< "\t" << (res.concordant ? "yes" : "no")
			<< "\t" << res.veriscore
			<< "\n";
	}

	total_unmapped = total_svtig_count - total_read_count;

	double avg_sv_count = (total_read_count > 0) ? (double)total_sv_count / total_read_count : 0;
	double avg_map_ratio = (total_read_count > 0) ? total_map_ratio / total_read_count : 0;
	double avg_mapq_val = (total_read_count > 0) ? total_mapq / total_read_count : 0;
	double avg_edit_dist = (total_read_count > 0) ? total_edit_dist / total_read_count : 0;
	double avg_aln_score = (total_read_count > 0) ? total_aln_score / total_read_count : 0;
	double avg_aln_identity = (total_read_count > 0) ? total_aln_identity / total_read_count : 0;
	double avg_svtig_size = (total_read_count > 0) ? (double)total_svtig_size / total_read_count : 0;
	double concordance_rate = (total_read_count > 0) ? (double)total_concordant / total_read_count : 0;

	std::string summary_path = base + ".summary.tsv";
	std::ofstream fp_summary(summary_path);

	fp_summary << "metric\tvalue\n";
	fp_summary << "total_svtigs_input\t" << total_svtig_count << "\n";
	fp_summary << "total_mapped\t" << total_read_count << "\n";
	fp_summary << "total_unmapped\t" << total_unmapped << "\n";
	fp_summary << "total_concordant\t" << total_concordant << "\n";
	fp_summary << "concordance_rate\t" << concordance_rate << "\n";
	fp_summary << "total_ins\t" << total_ins << "\n";
	fp_summary << "total_del\t" << total_del << "\n";
	fp_summary << "avg_sv_count\t" << avg_sv_count << "\n";
	fp_summary << "avg_svtig_size\t" << avg_svtig_size << "\n";
	fp_summary << "avg_map_ratio\t" << avg_map_ratio << "\n";
	fp_summary << "avg_mapq\t" << avg_mapq_val << "\n";
	fp_summary << "avg_edit_dist\t" << avg_edit_dist << "\n";
	fp_summary << "avg_aln_score\t" << avg_aln_score << "\n";
	fp_summary << "avg_aln_identity\t" << avg_aln_identity << "\n";
	fp_summary << "concordance_min_map_ratio\t" << params.min_map_ratio << "\n";
	fp_summary << "concordance_min_aln_identity\t" << params.min_aln_identity << "\n";

	// Aggregate VeriScore for this batch
	{
		std::vector<double> scores;
		double sum = 0.0;
		for (size_t i = batch_start; i < results.size(); i++)
		{
			scores.push_back(results[i].veriscore);
			sum += results[i].veriscore;
		}
		if (!scores.empty())
		{
			std::sort(scores.begin(), scores.end());
			double mean_vs = sum / scores.size();
			double median_vs = (scores.size() % 2 == 0)
				? (scores[scores.size()/2 - 1] + scores[scores.size()/2]) / 2.0
				: scores[scores.size()/2];
			fp_summary << "mean_veriscore\t" << mean_vs << "\n";
			fp_summary << "median_veriscore\t" << median_vs << "\n";
		}
	}

	std::string sizebin_path = base + ".size_bins.tsv";
	std::ofstream fp_sizebin(sizebin_path);
	fp_sizebin << "size_bin\ttotal\tconcordant\tconcordance_pct\n";
	for (auto& bin : size_bins)
	{
		double pct = (bin.total > 0) ? (double)bin.concordant / bin.total * 100 : 0;
		fp_sizebin << bin.label << "\t" << bin.total << "\t" << bin.concordant << "\t" << pct << "\n";
	}

	std::string svtype_path = base + ".svtype.tsv";
	std::ofstream fp_svtype(svtype_path);
	fp_svtype << "sv_type\ttotal\tconcordant\tconcordance_pct\n";
	for (auto& tb : type_bins)
	{
		double pct = (tb.total > 0) ? (double)tb.concordant / tb.total * 100 : 0;
		fp_svtype << tb.label << "\t" << tb.total << "\t" << tb.concordant << "\t" << pct << "\n";
	}
}


void Concordance::write_report(parameters& params)
{
	std::string report_path = params.log_path + params.sample_name + ".svtig.report.tsv";
	std::ofstream fp(report_path);

	fp << "svtig_name\tsvtig_size\thaplotype\tcombined_map_ratio\taln_identity\tmapq\tedit_dist\taln_score\tins_count\tdel_count\tsv_count\tsv_type\tconcordant\tveriscore\n";

	for (auto& res : results)
	{
		fp << res.name
			<< "\t" << res.svtig_size
			<< "\t" << res.haplotype
			<< "\t" << res.combined_map_ratio
			<< "\t" << res.aln_identity
			<< "\t" << res.mapq
			<< "\t" << res.edit_dist
			<< "\t" << res.aln_score
			<< "\t" << res.ins_count
			<< "\t" << res.del_count
			<< "\t" << res.sv_count
			<< "\t" << res.sv_type
			<< "\t" << (res.concordant ? "yes" : "no")
			<< "\t" << res.veriscore
			<< "\n";
	}
}


void Concordance::run(parameters& params)
{
	double conc_H1_svtig1 = 0, conc_H1_svtig2 = 0, conc_H2_svtig1 = 0, conc_H2_svtig2 = 0;

	count_svtigs(params);
	std::cerr << "  svtig count: " << this->h1_count << " (H1)";
	if (params.phase)
		std::cerr << ", " << this->h2_count << " (H2)";
	std::cerr << "\n";

	std::string paf_dir = params.log_path + "paf/";
	this->paf_H1_svtig1_path = paf_dir + params.sample_name + ".svtig.H1.paf";
	this->paf_H2_svtig1_path = paf_dir + params.sample_name + ".svtig.H2.paf";
	if (params.phase)
	{
		this->paf_H1_svtig2_path = paf_dir + params.sample_name + ".svtig2.H1.paf";
		this->paf_H2_svtig2_path = paf_dir + params.sample_name + ".svtig2.H2.paf";
	}

	if (!params.skip_mapping)
		run_mapping(params);

	std::cerr << "  Computing statistics...\n";
	conc_H1_svtig1 = mapping_stats(this->paf_H1_svtig1_path, h1_count, "svtig1_H1", params);
	conc_H2_svtig1 = mapping_stats(this->paf_H2_svtig1_path, h1_count, "svtig1_H2", params);

	if (params.phase)
	{
		conc_H1_svtig2 = mapping_stats(this->paf_H1_svtig2_path, h2_count, "svtig2_H1", params);
		conc_H2_svtig2 = mapping_stats(this->paf_H2_svtig2_path, h2_count, "svtig2_H2", params);
	}

	// Add unmapped svtigs to results
	{
		std::set<std::string> mapped_names;
		for (auto& r : results)
			mapped_names.insert(r.name);
		for (auto& kv : svtig_sizes)
		{
			if (mapped_names.find(kv.first) == mapped_names.end())
			{
				SvtigResult res;
				res.name = kv.first;
				res.svtig_size = kv.second;
				res.haplotype = "unmapped";
				res.sv_type = "NA";
				results.push_back(res);
			}
		}
	}

	write_report(params);

	std::string concordance_path = params.log_path + params.sample_name + ".svtig.concordance.tsv";
	std::ofstream fp_conc(concordance_path);
	fp_conc << "svtig_set\thaplotype\tconcordant\ttotal\tconcordance_pct\n";

	if (params.phase)
	{
		fp_conc << "svtig1\tH1\t" << conc_H1_svtig1 << "\t" << h1_count << "\t" << (conc_H1_svtig1 / (double)h1_count) * 100 << "\n";
		fp_conc << "svtig1\tH2\t" << conc_H2_svtig1 << "\t" << h1_count << "\t" << (conc_H2_svtig1 / (double)h1_count) * 100 << "\n";
		fp_conc << "svtig2\tH1\t" << conc_H1_svtig2 << "\t" << h2_count << "\t" << (conc_H1_svtig2 / (double)h2_count) * 100 << "\n";
		fp_conc << "svtig2\tH2\t" << conc_H2_svtig2 << "\t" << h2_count << "\t" << (conc_H2_svtig2 / (double)h2_count) * 100 << "\n";

		double conc_H1, conc_H2;
		if (conc_H1_svtig1 > conc_H2_svtig1)
		{
			conc_H1 = (conc_H1_svtig1 / (double)h1_count) * 100;
			conc_H2 = (conc_H2_svtig2 / (double)h2_count) * 100;
		}
		else
		{
			conc_H1 = (conc_H1_svtig2 / (double)h2_count) * 100;
			conc_H2 = (conc_H2_svtig1 / (double)h1_count) * 100;
		}

		std::cerr << "  Concordance: H1=" << conc_H1 << "%, H2=" << conc_H2
			<< "%, Average=" << (conc_H1 + conc_H2) / 2.0 << "%\n";
	}
	else
	{
		fp_conc << "svtig1\tH1\t" << conc_H1_svtig1 << "\t" << h1_count << "\t" << (conc_H1_svtig1 / (double)h1_count) * 100 << "\n";
		fp_conc << "svtig1\tH2\t" << conc_H2_svtig1 << "\t" << h1_count << "\t" << (conc_H2_svtig1 / (double)h1_count) * 100 << "\n";

		double conc_H1 = (conc_H1_svtig1 / (double)h1_count) * 100;
		double conc_H2 = (conc_H2_svtig1 / (double)h1_count) * 100;
		std::cerr << "  Concordance: H1=" << conc_H1 << "%, H2=" << conc_H2 << "%\n";
	}

	// Aggregate VeriScore
	if (!results.empty())
	{
		std::vector<double> scores;
		double sum = 0.0;
		for (auto& r : results)
		{
			scores.push_back(r.veriscore);
			sum += r.veriscore;
		}
		std::sort(scores.begin(), scores.end());
		double mean_vs = sum / scores.size();
		double median_vs = (scores.size() % 2 == 0)
			? (scores[scores.size()/2 - 1] + scores[scores.size()/2]) / 2.0
			: scores[scores.size()/2];

		std::cerr << "  VeriScore: mean=" << mean_vs << ", median=" << median_vs << "\n";
		fp_conc << "mean_veriscore\t" << mean_vs << "\n";
		fp_conc << "median_veriscore\t" << median_vs << "\n";
	}

	std::cerr << "  Results written to " << params.log_path << "\n";
}
