#include <iostream>
#include <algorithm>
#include "svtig_stats.h"


SvtigStats::SeqStats SvtigStats::compute(std::string fasta_path)
{
	SeqStats s;
	std::vector<int> sizes;
	std::string line;
	std::string current_seq;
	long gc_count = 0, n_count = 0;

	std::ifstream fp(fasta_path);
	if (!fp.good())
	{
		std::cerr << "  Error opening '" << fasta_path << "'\n";
		return s;
	}

	while (getline(fp, line))
	{
		if (line.empty())
			continue;

		if (line[0] == '>')
		{
			if (!current_seq.empty())
			{
				int len = (int)current_seq.size();
				sizes.push_back(len);
				s.total_bases += len;

				for (char c : current_seq)
				{
					if (c == 'G' || c == 'C' || c == 'g' || c == 'c')
						gc_count++;
					else if (c == 'N' || c == 'n')
						n_count++;
				}

				if (len >= 50 && len < 100) s.size_50_100++;
				else if (len >= 100 && len < 500) s.size_100_500++;
				else if (len >= 500 && len < 1000) s.size_500_1000++;
				else if (len >= 1000 && len < 10000) s.size_1000_10000++;
				else if (len >= 10000) s.size_10000_plus++;

				current_seq.clear();
			}
		}
		else
		{
			current_seq += line;
		}
	}

	// Last sequence
	if (!current_seq.empty())
	{
		int len = (int)current_seq.size();
		sizes.push_back(len);
		s.total_bases += len;

		for (char c : current_seq)
		{
			if (c == 'G' || c == 'C' || c == 'g' || c == 'c')
				gc_count++;
			else if (c == 'N' || c == 'n')
				n_count++;
		}

		if (len >= 50 && len < 100) s.size_50_100++;
		else if (len >= 100 && len < 500) s.size_100_500++;
		else if (len >= 500 && len < 1000) s.size_500_1000++;
		else if (len >= 1000 && len < 10000) s.size_1000_10000++;
		else if (len >= 10000) s.size_10000_plus++;
	}

	s.count = (int)sizes.size();
	if (s.count == 0)
		return s;

	std::sort(sizes.begin(), sizes.end());

	s.min_size = sizes.front();
	s.max_size = sizes.back();
	s.mean_size = (double)s.total_bases / s.count;
	s.gc_content = (s.total_bases > 0) ? (double)gc_count / s.total_bases : 0;
	s.n_content = (s.total_bases > 0) ? (double)n_count / s.total_bases : 0;

	// Median
	if (s.count % 2 == 0)
		s.median_size = (sizes[s.count / 2 - 1] + sizes[s.count / 2]) / 2;
	else
		s.median_size = sizes[s.count / 2];

	// N50: sort descending, cumulative sum until >= total/2
	std::sort(sizes.begin(), sizes.end(), std::greater<int>());
	long cumulative = 0;
	long half = s.total_bases / 2;
	for (int sz : sizes)
	{
		cumulative += sz;
		if (cumulative >= half)
		{
			s.n50 = sz;
			break;
		}
	}

	return s;
}


void SvtigStats::write_report(SeqStats& s1, SeqStats& s2, SeqStats& combined, bool has_svtig2, parameters& params)
{
	std::string report_path = params.log_path + params.sample_name + ".stats.tsv";
	std::ofstream fp(report_path);

	if (has_svtig2)
	{
		fp << "metric\tsvtig1\tsvtig2\tcombined\n";
		fp << "total_svtigs\t" << s1.count << "\t" << s2.count << "\t" << combined.count << "\n";
		fp << "total_bases\t" << s1.total_bases << "\t" << s2.total_bases << "\t" << combined.total_bases << "\n";
		fp << "min_size\t" << s1.min_size << "\t" << s2.min_size << "\t" << combined.min_size << "\n";
		fp << "max_size\t" << s1.max_size << "\t" << s2.max_size << "\t" << combined.max_size << "\n";
		fp << "mean_size\t" << s1.mean_size << "\t" << s2.mean_size << "\t" << combined.mean_size << "\n";
		fp << "median_size\t" << s1.median_size << "\t" << s2.median_size << "\t" << combined.median_size << "\n";
		fp << "N50\t" << s1.n50 << "\t" << s2.n50 << "\t" << combined.n50 << "\n";
		fp << "GC_content\t" << s1.gc_content << "\t" << s2.gc_content << "\t" << combined.gc_content << "\n";
		fp << "N_content\t" << s1.n_content << "\t" << s2.n_content << "\t" << combined.n_content << "\n";
		fp << "size_50_100\t" << s1.size_50_100 << "\t" << s2.size_50_100 << "\t" << combined.size_50_100 << "\n";
		fp << "size_100_500\t" << s1.size_100_500 << "\t" << s2.size_100_500 << "\t" << combined.size_100_500 << "\n";
		fp << "size_500_1000\t" << s1.size_500_1000 << "\t" << s2.size_500_1000 << "\t" << combined.size_500_1000 << "\n";
		fp << "size_1000_10000\t" << s1.size_1000_10000 << "\t" << s2.size_1000_10000 << "\t" << combined.size_1000_10000 << "\n";
		fp << "size_10000_plus\t" << s1.size_10000_plus << "\t" << s2.size_10000_plus << "\t" << combined.size_10000_plus << "\n";
	}
	else
	{
		fp << "metric\tvalue\n";
		fp << "total_svtigs\t" << s1.count << "\n";
		fp << "total_bases\t" << s1.total_bases << "\n";
		fp << "min_size\t" << s1.min_size << "\n";
		fp << "max_size\t" << s1.max_size << "\n";
		fp << "mean_size\t" << s1.mean_size << "\n";
		fp << "median_size\t" << s1.median_size << "\n";
		fp << "N50\t" << s1.n50 << "\n";
		fp << "GC_content\t" << s1.gc_content << "\n";
		fp << "N_content\t" << s1.n_content << "\n";
		fp << "size_50_100\t" << s1.size_50_100 << "\n";
		fp << "size_100_500\t" << s1.size_100_500 << "\n";
		fp << "size_500_1000\t" << s1.size_500_1000 << "\n";
		fp << "size_1000_10000\t" << s1.size_1000_10000 << "\n";
		fp << "size_10000_plus\t" << s1.size_10000_plus << "\n";
	}

	// Terminal summary
	std::cerr << "  SVtig statistics:\n";
	if (has_svtig2)
	{
		std::cerr << "    svtig1: " << s1.count << " SVtigs, " << s1.total_bases << " bases, N50=" << s1.n50 << "\n";
		std::cerr << "    svtig2: " << s2.count << " SVtigs, " << s2.total_bases << " bases, N50=" << s2.n50 << "\n";
		std::cerr << "    combined: " << combined.count << " SVtigs, " << combined.total_bases << " bases, N50=" << combined.n50 << "\n";
	}
	else
	{
		std::cerr << "    " << s1.count << " SVtigs, " << s1.total_bases << " bases, N50=" << s1.n50 << "\n";
	}
	std::cerr << "  Results written to " << params.log_path << "\n";
}


void SvtigStats::run(parameters& params)
{
	std::cerr << "  Computing SVtig statistics...\n";

	bool has_svtig2 = !params.svtig2_path.empty();
	std::string fasta = params.svtig1_path.empty() ? params.fasta : params.svtig1_path;

	SeqStats s1 = compute(fasta);
	SeqStats s2;
	SeqStats combined;

	if (has_svtig2)
	{
		s2 = compute(params.svtig2_path);

		// Recompute N50/median for combined set
		std::vector<int> all_sizes;
		long gc_count = 0, n_count = 0;

		for (const auto& path : {fasta, params.svtig2_path})
		{
			std::ifstream fp(path);
			std::string line, seq;
			while (getline(fp, line))
			{
				if (line.empty()) continue;
				if (line[0] == '>')
				{
					if (!seq.empty())
					{
						all_sizes.push_back((int)seq.size());
						for (char c : seq)
						{
							if (c == 'G' || c == 'C' || c == 'g' || c == 'c') gc_count++;
							else if (c == 'N' || c == 'n') n_count++;
						}
						seq.clear();
					}
				}
				else
					seq += line;
			}
			if (!seq.empty())
			{
				all_sizes.push_back((int)seq.size());
				for (char c : seq)
				{
					if (c == 'G' || c == 'C' || c == 'g' || c == 'c') gc_count++;
					else if (c == 'N' || c == 'n') n_count++;
				}
			}
		}

		combined.count = (int)all_sizes.size();
		combined.total_bases = s1.total_bases + s2.total_bases;
		combined.gc_content = (combined.total_bases > 0) ? (double)gc_count / combined.total_bases : 0;
		combined.n_content = (combined.total_bases > 0) ? (double)n_count / combined.total_bases : 0;

		if (combined.count > 0)
		{
			std::sort(all_sizes.begin(), all_sizes.end());
			combined.min_size = all_sizes.front();
			combined.max_size = all_sizes.back();
			combined.mean_size = (double)combined.total_bases / combined.count;

			if (combined.count % 2 == 0)
				combined.median_size = (all_sizes[combined.count / 2 - 1] + all_sizes[combined.count / 2]) / 2;
			else
				combined.median_size = all_sizes[combined.count / 2];

			std::sort(all_sizes.begin(), all_sizes.end(), std::greater<int>());
			long cumulative = 0, half = combined.total_bases / 2;
			for (int sz : all_sizes)
			{
				cumulative += sz;
				if (cumulative >= half) { combined.n50 = sz; break; }
			}

			for (int sz : all_sizes)
			{
				if (sz >= 50 && sz < 100) combined.size_50_100++;
				else if (sz >= 100 && sz < 500) combined.size_100_500++;
				else if (sz >= 500 && sz < 1000) combined.size_500_1000++;
				else if (sz >= 1000 && sz < 10000) combined.size_1000_10000++;
				else if (sz >= 10000) combined.size_10000_plus++;
			}
		}
	}

	write_report(s1, s2, combined, has_svtig2, params);
}
