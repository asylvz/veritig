#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <climits>
#include "validate.h"


void Validate::count_svs(parameters& params)
{
	std::string line;
	std::string current_name;
	int current_bases = 0;

	std::ifstream fp(params.fasta);
	while (getline(fp, line))
	{
		if (line[0] == '>')
		{
			if (!current_name.empty())
				sv_sizes[current_name] = current_bases;
			current_name = line.substr(1);
			size_t space_pos = current_name.find(' ');
			if (space_pos != std::string::npos)
				current_name = current_name.substr(0, space_pos);
			current_bases = 0;
			this->sv_count++;
		}
		else
		{
			current_bases += (int)line.size();
		}
	}
	if (!current_name.empty())
		sv_sizes[current_name] = current_bases;
}

void Validate::run_mapping(parameters& params)
{
	std::string minimap_cmd;
	std::string threads_str = std::to_string(params.threads);

	std::string paf_dir = params.log_path + "paf/";
	this->paf_H1_sv_path = paf_dir + params.sample_name + ".sv.H1.paf";
	this->paf_H2_sv_path = paf_dir + params.sample_name + ".sv.H2.paf";

	std::cerr << "  Mapping SVs to H1...";
	minimap_cmd = "minimap2 -cx " + params.minimap_preset + " -t " + threads_str + " --paf-no-hit " + params.haplo1_assembly_path + " " + params.fasta + " --secondary=no > " + this->paf_H1_sv_path + " 2>/dev/null";
	run_command(minimap_cmd, this->paf_H1_sv_path);
	std::cerr << " done\n";

	std::cerr << "  Mapping SVs to H2...";
	minimap_cmd = "minimap2 -cx " + params.minimap_preset + " -t " + threads_str + " --paf-no-hit " + params.haplo2_assembly_path + " " + params.fasta + " --secondary=no > " + this->paf_H2_sv_path + " 2>/dev/null";
	run_command(minimap_cmd, this->paf_H2_sv_path);
	std::cerr << " done\n";
}


void Validate::compute_stats(parameters& params)
{
	// --- Pass 1: Read H1 PAF, populate reads map ---
	{
		std::string line;
		std::ifstream fp(this->paf_H1_sv_path);
		if (!fp.good())
		{
			std::cerr << "Error opening '" << this->paf_H1_sv_path << std::endl;
			return;
		}

		while (getline(fp, line))
		{
			PafRecord rec;
			int sv_count = parse_paf_line(line, rec);
			if (sv_count == RETURN_ERROR)
				continue;
			if (!rec.is_primary)
				continue;

			auto it = reads.find(rec.read_name);
			if (it != reads.end())
			{
				// Same SV, multiple H1 mappings: keep best
				if (rec.map_ratio > it->second->highest_map_ratio)
				{
					it->second->svtig_size = rec.svtig_size;
					it->second->edit_dist = rec.edit_dist;
					it->second->aln_score = rec.aln_score;
					it->second->mapq = rec.mapq;
					it->second->highest_map_ratio = rec.map_ratio;
					it->second->aln_identity = rec.aln_identity;
					it->second->contig = rec.contig;
					it->second->ins_count = rec.ins_count;
					it->second->del_count = rec.del_count;
					it->second->sv_count = sv_count;
				}
			}
			else
			{
				Read *r = new Read();
				r->haplo = 1;
				r->name = rec.read_name;
				r->svtig_size = rec.svtig_size;
				r->edit_dist = rec.edit_dist;
				r->aln_score = rec.aln_score;
				r->mapq = rec.mapq;
				r->highest_map_ratio = rec.map_ratio;
				r->aln_identity = rec.aln_identity;
				r->contig = rec.contig;
				r->ins_count = rec.ins_count;
				r->del_count = rec.del_count;
				r->sv_count = sv_count;
				reads.insert(std::pair<std::string, Read*>(rec.read_name, r));
			}
		}
	}

	// --- Pass 2: Read H2 PAF, update/add to reads map ---
	{
		std::string line;
		std::ifstream fp(this->paf_H2_sv_path);
		if (!fp.good())
		{
			std::cerr << "Error opening '" << this->paf_H2_sv_path << std::endl;
			return;
		}

		while (getline(fp, line))
		{
			PafRecord rec;
			int sv_count = parse_paf_line(line, rec);
			if (sv_count == RETURN_ERROR)
				continue;
			if (!rec.is_primary)
				continue;

			auto it = reads.find(rec.read_name);
			if (it != reads.end())
			{
				// Already seen in H1 → mark as homozygous
				it->second->homo = true;

				// If H2 mapping is better, update metrics and mark H2
				if (rec.map_ratio > it->second->highest_map_ratio)
				{
					it->second->haplo = 2;
					it->second->svtig_size = rec.svtig_size;
					it->second->edit_dist = rec.edit_dist;
					it->second->aln_score = rec.aln_score;
					it->second->mapq = rec.mapq;
					it->second->highest_map_ratio = rec.map_ratio;
					it->second->aln_identity = rec.aln_identity;
					it->second->contig = rec.contig;
					it->second->ins_count = rec.ins_count;
					it->second->del_count = rec.del_count;
					it->second->sv_count = sv_count;
				}
			}
			else
			{
				Read *r = new Read();
				r->haplo = 2;
				r->name = rec.read_name;
				r->svtig_size = rec.svtig_size;
				r->edit_dist = rec.edit_dist;
				r->aln_score = rec.aln_score;
				r->mapq = rec.mapq;
				r->highest_map_ratio = rec.map_ratio;
				r->aln_identity = rec.aln_identity;
				r->contig = rec.contig;
				r->ins_count = rec.ins_count;
				r->del_count = rec.del_count;
				r->sv_count = sv_count;
				reads.insert(std::pair<std::string, Read*>(rec.read_name, r));
			}
		}
	}

	// --- Build results ---
	for (auto it = reads.begin(); it != reads.end(); ++it)
	{
		Read* r = it->second;

		bool is_concordant = (r->highest_map_ratio >= params.min_map_ratio && r->aln_identity >= params.min_aln_identity);

		std::string sv_type;
		if (r->ins_count > 0 && r->del_count == 0) sv_type = "INS";
		else if (r->del_count > 0 && r->ins_count == 0) sv_type = "DEL";
		else if (r->ins_count > 0 && r->del_count > 0) sv_type = "MIXED";
		else sv_type = "NONE";

		std::string hap_label;
		if (r->homo)
			hap_label = "Homo";
		else if (r->haplo == 1)
			hap_label = "H1";
		else
			hap_label = "H2";

		SvResult res;
		res.name = it->first;
		res.sv_size = r->svtig_size;
		res.contig = r->contig;
		res.haplotype = hap_label;
		res.map_ratio = r->highest_map_ratio;
		res.aln_identity = r->aln_identity;
		res.mapq = r->mapq;
		res.edit_dist = r->edit_dist;
		res.aln_score = r->aln_score;
		res.ins_count = r->ins_count;
		res.del_count = r->del_count;
		res.sv_count = r->sv_count;
		res.sv_type = sv_type;
		res.concordant = is_concordant;
		res.veriscore = compute_veriscore(r->highest_map_ratio, r->aln_identity, r->mapq);
		results.push_back(res);
	}

	// Add unmapped SVs
	for (auto& kv : sv_sizes)
	{
		if (reads.find(kv.first) == reads.end())
		{
			SvResult res;
			res.name = kv.first;
			res.sv_size = kv.second;
			res.haplotype = "unmapped";
			res.sv_type = "NA";
			results.push_back(res);
		}
	}
}


void Validate::write_report(parameters& params)
{
	// --- Main report (always) ---
	std::string report_path = params.log_path + params.sample_name + ".sv.report.tsv";
	std::ofstream fp(report_path);

	fp << "sv_name\tsv_size\tcontig\thaplotype\tmap_ratio\taln_identity\tmapq\tedit_dist\taln_score\tins_count\tdel_count\tsv_count\tsv_type\tconcordant\tveriscore\n";

	int haplo1 = 0, haplo2 = 0, homo = 0, unmapped = 0;
	int concordant_h1 = 0, concordant_h2 = 0, concordant_homo = 0;

	for (auto& res : results)
	{
		fp << res.name
			<< "\t" << res.sv_size
			<< "\t" << res.contig
			<< "\t" << res.haplotype
			<< "\t" << res.map_ratio
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

		if (res.haplotype == "Homo") { homo++; if (res.concordant) concordant_homo++; }
		else if (res.haplotype == "H1") { haplo1++; if (res.concordant) concordant_h1++; }
		else if (res.haplotype == "H2") { haplo2++; if (res.concordant) concordant_h2++; }
	}

	unmapped = sv_count - (int)results.size();

	// --- Concordance summary (always) ---
	std::string concordance_path = params.log_path + params.sample_name + ".sv.concordance.tsv";
	std::ofstream fp_conc(concordance_path);
	fp_conc << "haplotype\ttotal\tconcordant\tconcordance_pct\n";

	int total_assigned = haplo1 + haplo2 + homo;
	int total_concordant = concordant_h1 + concordant_h2 + concordant_homo;

	if (haplo1 > 0) fp_conc << "H1\t" << haplo1 << "\t" << concordant_h1 << "\t" << (double)concordant_h1 / haplo1 * 100 << "\n";
	if (haplo2 > 0) fp_conc << "H2\t" << haplo2 << "\t" << concordant_h2 << "\t" << (double)concordant_h2 / haplo2 * 100 << "\n";
	if (homo > 0) fp_conc << "Homo\t" << homo << "\t" << concordant_homo << "\t" << (double)concordant_homo / homo * 100 << "\n";
	fp_conc << "Total\t" << total_assigned << "\t" << total_concordant << "\t" << ((total_assigned > 0) ? (double)total_concordant / total_assigned * 100 : 0) << "\n";

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

		fp_conc << "mean_veriscore\t" << mean_vs << "\n";
		fp_conc << "median_veriscore\t" << median_vs << "\n";

		std::cerr << "  VeriScore: mean=" << mean_vs << ", median=" << median_vs << "\n";
	}

	// --- Terminal output ---
	std::cerr << "  Mapped: " << results.size() << "/" << sv_count;
	if (unmapped > 0) std::cerr << " (" << unmapped << " unmapped)";
	std::cerr << "\n";
	std::cerr << "  Haplotype assignment: H1=" << haplo1 << ", H2=" << haplo2 << ", Homo=" << homo << "\n";
	std::cerr << "  Concordant: " << total_concordant << "/" << total_assigned
		<< " (" << ((total_assigned > 0) ? (double)total_concordant / total_assigned * 100 : 0) << "%)\n";
	std::cerr << "  Results written to " << params.log_path << "\n";
}


void Validate::write_detailed(std::string filename, std::string haplotype, parameters& params)
{
	std::string detailed_dir = params.log_path + "detailed/";
	std::string base = detailed_dir + params.sample_name + ".sv." + haplotype;

	// Filter results for this haplotype category
	std::vector<SvResult*> filtered;
	for (auto& res : results)
	{
		if (haplotype == "H1" && (res.haplotype == "H1" || res.haplotype == "Homo"))
			filtered.push_back(&res);
		else if (haplotype == "H2" && (res.haplotype == "H2" || res.haplotype == "Homo"))
			filtered.push_back(&res);
	}

	std::string details_path = base + ".details.tsv";
	std::ofstream fp_details(details_path);
	fp_details << "sv_name\tsv_size\tcontig\thaplotype\tmap_ratio\taln_identity\tmapq\tedit_dist\taln_score\tins_count\tdel_count\tsv_count\tsv_type\tconcordant\tveriscore\n";

	int total_concordant = 0, total_ins = 0, total_del = 0;
	long total_sv_size = 0;
	double total_map_ratio = 0, total_aln_identity = 0, total_mapq = 0, total_edit_dist = 0, total_aln_score = 0;

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
		{"INS", 0, 0}, {"DEL", 0, 0}, {"MIXED", 0, 0}, {"NONE", 0, 0}
	};

	for (auto* res : filtered)
	{
		if (res->concordant) total_concordant++;
		total_ins += res->ins_count;
		total_del += res->del_count;
		total_sv_size += res->sv_size;
		total_map_ratio += res->map_ratio;
		total_aln_identity += res->aln_identity;
		total_mapq += res->mapq;
		total_edit_dist += res->edit_dist;
		total_aln_score += res->aln_score;

		for (auto& bin : size_bins)
		{
			if (res->sv_size >= bin.min_size && res->sv_size < bin.max_size)
			{
				bin.total++;
				if (res->concordant) bin.concordant++;
				break;
			}
		}
		for (auto& tb : type_bins)
		{
			if (tb.label == res->sv_type)
			{
				tb.total++;
				if (res->concordant) tb.concordant++;
				break;
			}
		}

		fp_details << res->name << "\t" << res->sv_size << "\t" << res->contig
			<< "\t" << res->haplotype << "\t" << res->map_ratio << "\t" << res->aln_identity
			<< "\t" << res->mapq << "\t" << res->edit_dist << "\t" << res->aln_score
			<< "\t" << res->ins_count << "\t" << res->del_count << "\t" << res->sv_count
			<< "\t" << res->sv_type << "\t" << (res->concordant ? "yes" : "no")
			<< "\t" << res->veriscore << "\n";
	}

	int n = (int)filtered.size();

	// Summary
	std::string summary_path = base + ".summary.tsv";
	std::ofstream fp_summary(summary_path);
	fp_summary << "metric\tvalue\n";
	fp_summary << "total_svs_input\t" << sv_count << "\n";
	fp_summary << "total_mapped\t" << n << "\n";
	fp_summary << "total_concordant\t" << total_concordant << "\n";
	fp_summary << "concordance_rate\t" << ((n > 0) ? (double)total_concordant / n : 0) << "\n";
	fp_summary << "total_ins\t" << total_ins << "\n";
	fp_summary << "total_del\t" << total_del << "\n";
	fp_summary << "avg_sv_size\t" << ((n > 0) ? (double)total_sv_size / n : 0) << "\n";
	fp_summary << "avg_map_ratio\t" << ((n > 0) ? total_map_ratio / n : 0) << "\n";
	fp_summary << "avg_mapq\t" << ((n > 0) ? total_mapq / n : 0) << "\n";
	fp_summary << "avg_edit_dist\t" << ((n > 0) ? total_edit_dist / n : 0) << "\n";
	fp_summary << "avg_aln_score\t" << ((n > 0) ? total_aln_score / n : 0) << "\n";
	fp_summary << "avg_aln_identity\t" << ((n > 0) ? total_aln_identity / n : 0) << "\n";
	fp_summary << "concordance_min_map_ratio\t" << params.min_map_ratio << "\n";
	fp_summary << "concordance_min_aln_identity\t" << params.min_aln_identity << "\n";

	// Aggregate VeriScore
	if (!filtered.empty())
	{
		std::vector<double> scores;
		double sum = 0.0;
		for (auto* r : filtered)
		{
			scores.push_back(r->veriscore);
			sum += r->veriscore;
		}
		std::sort(scores.begin(), scores.end());
		double mean_vs = sum / scores.size();
		double median_vs = (scores.size() % 2 == 0)
			? (scores[scores.size()/2 - 1] + scores[scores.size()/2]) / 2.0
			: scores[scores.size()/2];
		fp_summary << "mean_veriscore\t" << mean_vs << "\n";
		fp_summary << "median_veriscore\t" << median_vs << "\n";
	}

	// Size bins
	std::string sizebin_path = base + ".size_bins.tsv";
	std::ofstream fp_sizebin(sizebin_path);
	fp_sizebin << "size_bin\ttotal\tconcordant\tconcordance_pct\n";
	for (auto& bin : size_bins)
	{
		double pct = (bin.total > 0) ? (double)bin.concordant / bin.total * 100 : 0;
		fp_sizebin << bin.label << "\t" << bin.total << "\t" << bin.concordant << "\t" << pct << "\n";
	}

	// SV types
	std::string svtype_path = base + ".svtype.tsv";
	std::ofstream fp_svtype(svtype_path);
	fp_svtype << "sv_type\ttotal\tconcordant\tconcordance_pct\n";
	for (auto& tb : type_bins)
	{
		double pct = (tb.total > 0) ? (double)tb.concordant / tb.total * 100 : 0;
		fp_svtype << tb.label << "\t" << tb.total << "\t" << tb.concordant << "\t" << pct << "\n";
	}
}


void Validate::run(parameters& params)
{
	count_svs(params);
	std::cerr << "  SV count: " << this->sv_count << "\n";

	std::string paf_dir = params.log_path + "paf/";
	this->paf_H1_sv_path = paf_dir + params.sample_name + ".sv.H1.paf";
	this->paf_H2_sv_path = paf_dir + params.sample_name + ".sv.H2.paf";

	if (!params.skip_mapping)
		run_mapping(params);

	std::cerr << "  Computing statistics...\n";
	compute_stats(params);
	write_report(params);

	if (params.detailed)
	{
		write_detailed(this->paf_H1_sv_path, "H1", params);
		write_detailed(this->paf_H2_sv_path, "H2", params);
	}

	for (auto it = reads.begin(); it != reads.end(); ++it)
		delete it->second;
	reads.clear();
}
