# veritig

A toolkit for analyzing and verifying structural variant sequences (svtigs) using haplotype-resolved assemblies.

veritig maps svtig sequences to haplotype assemblies via minimap2, computes concordance metrics and VeriScore (a composite quality metric), and provides tools for quality assessment and filtering.

## Quick Start

```bash
git clone <repo-url>
cd veritig
make
```

## Requirements

- C++17 compatible compiler (g++)
- [minimap2](https://github.com/lh3/minimap2) (must be in PATH; not required for `--stats`)

## Modes

veritig provides five modes for different stages of svtig analysis:

### `--stats` — svtig Statistics

Computes basic statistics of svtig sequences directly from FASTA files. No mapping required.

**Use when:** You want a quick quality overview of your svtig set before running any analysis.

```bash
# Single set
veritig --stats --svtig1 svtigs.fa -o output/ -s sample_name

# Two sets (e.g., phased svtigs)
veritig --stats --svtig1 svtigs_h1.fa --svtig2 svtigs_h2.fa -o output/ -s sample_name
```

**Output:** `{sample}.stats.tsv` with total count, total bases, min/max/mean/median size, N50, GC content, N content, and size distribution. When two sets are provided, each set is reported separately along with combined statistics.

### `--concordance` — Concordance Analysis

Measures how well svtig sequences map to each haplotype assembly. All mappings for each svtig are **accumulated**: mapping metrics are averaged across multiple alignments, and supplementary alignments contribute to a combined mapping ratio. The result is a concordance score and VeriScore per svtig.

**Use when:** You want to measure overall mapping quality of svtigs against assemblies.

```bash
# Single haplotype
veritig --concordance --h1 hap1.fa --svtig1 svtigs.fa -o output/ -s sample_name

# Phased mode (both haplotypes)
veritig --concordance --h1 hap1.fa --svtig1 svtigs_h1.fa \
                      --h2 hap2.fa --svtig2 svtigs_h2.fa \
                      --phase -o output/ -s sample_name
```

**Output:** `{sample}.svtig.report.tsv` (per-svtig metrics including VeriScore) and `{sample}.svtig.concordance.tsv` (summary with mean/median VeriScore). Unmapped svtigs are included in the report with `haplotype=unmapped` and `veriscore=0`.

### `--validate` — SV Validation

Evaluates each SV sequence individually by mapping it to both haplotype assemblies. Only the **best mapping** (highest mapping ratio) is kept per SV. Each SV is assigned a haplotype (H1, H2, or Homo) based on which assemblies it maps to.

**Use when:** You want per-SV mapping metrics (sequence identity, mapping ratio) and haplotype assignment.

```bash
veritig --validate --h1 hap1.fa --h2 hap2.fa --fasta svs.fa -o output/ -s sample_name
```

**Output:** `{sample}.sv.report.tsv` (per-SV metrics with haplotype and VeriScore) and `{sample}.sv.concordance.tsv` (summary with mean/median VeriScore). Unmapped SVs are included with `haplotype=unmapped` and `veriscore=0`.

### `--compare` — svtig Set Comparison

Compares two svtig sets by mapping both to the same assembly and reporting side-by-side concordance metrics. Useful for comparing svtigs from different callers, parameters, or samples.

**Use when:** You want to evaluate which svtig set has better mapping quality, or how two sets differ across size ranges.

```bash
veritig --compare --h1 hap1.fa --svtig1 setA.fa --svtig2 setB.fa -o output/ -s sample_name
```

**Output:** `{sample}.compare.tsv` with per-set metrics including concordance rate, average mapping ratio, average identity, and size-binned concordance percentages.

### `--filter` — svtig Filtering

Maps svtig sequences to a single haplotype assembly and outputs only those that pass concordance thresholds as a new FASTA file.

**Use when:** You want to produce a cleaned svtig set for downstream analysis.

```bash
# Default thresholds (0.85/0.85)
veritig --filter --h1 hap1.fa --svtig1 svtigs.fa -o output/ -s sample_name

# Stricter thresholds
veritig --filter --h1 hap1.fa --svtig1 svtigs.fa -M 0.95 -I 0.95 -o output/ -s sample_name
```

**Output:** `{sample}.filtered.fa` containing only concordant svtigs.

### Key differences between `--concordance` and `--validate`

| | `--concordance` | `--validate` |
|---|---|---|
| Purpose | Measure concordance | Per-SV metrics and haplotype assignment |
| Multiple mappings | Accumulate and average | Keep only the best |
| Supplementary alignments | Included in combined mapping ratio | Not accumulated |
| Haplotype assignment | No | Yes (H1 / H2 / Homo) |
| Contig info | Not reported | Reported per SV |
| Phased mode | Yes (`--phase`) | No |

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `--threads (-T)` | Number of threads for minimap2 | 16 |
| `--min-map-ratio (-M)` | Minimum mapping ratio for concordance | 0.85 |
| `--min-identity (-I)` | Minimum alignment identity for concordance | 0.85 |
| `--detailed (-D)` | Write additional analysis files | off |
| `--skip-mapping (-m)` | Skip minimap2 mapping, use existing PAF files | off |
| `--out (-o)` | Output folder path | current directory |
| `--sample (-s)` | Sample name | sample |

## Output

Results are written to `veritig_results/` under the output directory.

```
veritig_results/
├── {sample}.stats.tsv                 # svtig statistics (--stats)
├── {sample}.svtig.report.tsv          # Per-svtig concordance report (--concordance)
├── {sample}.svtig.concordance.tsv     # svtig concordance summary (--concordance)
├── {sample}.sv.report.tsv             # Per-SV report with haplotype (--validate)
├── {sample}.sv.concordance.tsv        # SV concordance summary (--validate)
├── {sample}.compare.tsv               # Side-by-side comparison (--compare)
├── {sample}.filtered.fa               # Filtered svtigs (--filter)
├── paf/                               # Minimap2 alignment files
│   ├── {sample}.svtig.H1.paf
│   ├── {sample}.svtig.H2.paf
│   ├── {sample}.sv.H1.paf
│   ├── {sample}.sv.H2.paf
│   ├── {sample}.filter.paf
│   ├── {sample}.compare.svtig1.paf
│   └── {sample}.compare.svtig2.paf
└── detailed/                          # Only with --detailed
    ├── {sample}.svtig.{hap}.details.tsv
    ├── {sample}.svtig.{hap}.summary.tsv
    ├── {sample}.svtig.{hap}.size_bins.tsv
    ├── {sample}.svtig.{hap}.svtype.tsv
    └── ...
```

### Report columns

**svtig report** (`--concordance`):
`svtig_name`, `svtig_size`, `haplotype`, `combined_map_ratio`, `aln_identity`, `mapq`, `edit_dist`, `aln_score`, `ins_count`, `del_count`, `sv_count`, `sv_type`, `concordant`, `veriscore`

**SV report** (`--validate`):
`sv_name`, `sv_size`, `contig`, `haplotype`, `map_ratio`, `aln_identity`, `mapq`, `edit_dist`, `aln_score`, `ins_count`, `del_count`, `sv_count`, `sv_type`, `concordant`, `veriscore`

An svtig/SV is **concordant** if `map_ratio >= min_map_ratio` AND `aln_identity >= min_identity`.

Unmapped svtigs/SVs are listed with `haplotype=unmapped`, `sv_type=NA`, all metrics at zero, and `veriscore=0`.

## VeriScore

VeriScore is a composite quality metric that summarizes mapping quality as a single continuous value in [0, 1]:

```
VeriScore = (R × I × C)^(1/3)
```

where:
- **R** = mapping ratio (alignment completeness)
- **I** = alignment identity (alignment accuracy)
- **C** = 1 - 10^(-Q/10) (Phred-scaled mapping confidence)
- **Q** = max(MAPQ, 3), capped at 60

VeriScore is the geometric mean of three quality dimensions:

- **Geometric mean** is the statistically appropriate average for [0,1] ratio-scale variables. Unlike the arithmetic mean, it does not mask a weak component; unlike the simple product, it does not over-penalize.
- **Phred normalization** converts MAPQ back to its natural probability scale, preserving the non-linear diminishing returns (MAPQ 10 vs 20 is a bigger quality jump than 40 vs 60).
- **MAPQ floor = 3** corresponds to the maximum-entropy (uninformative) prior P_correct ≈ 0.5. When MAPQ = 0 (confidence unavailable), this avoids discarding otherwise well-mapped sequences while still penalizing ambiguous placements.
- **MAPQ cap = 60** avoids numerical issues as P_correct approaches 1.0.

Key properties:
- **Range**: [0, 1], interpretable as overall mapping quality
- **Monotonic**: increasing in all three inputs
- **Zero-dominant**: if R = 0 or I = 0, VeriScore = 0
- **Balanced**: no implicit weighting; a 10% drop in any metric reduces the score by ~3.4%
- **Complementary to concordance**: concordance is binary (pass/fail); VeriScore provides continuous ranking within both concordant and discordant groups

**Aggregate VeriScore** (mean and median) is reported in the concordance summary and terminal output as a single-number quality indicator for the entire svtig/SV set.

### With `--detailed`

Per-haplotype analysis files in the `detailed/` subdirectory:

| File | Description |
|------|-------------|
| `*.details.tsv` | Per-svtig/SV metrics for each haplotype mapping |
| `*.summary.tsv` | Aggregate statistics (counts, averages, mean/median VeriScore) |
| `*.size_bins.tsv` | Concordance by SV size range (50-100, 100-500, 500-1k, 1k-10k, 10k+) |
| `*.svtype.tsv` | Concordance by SV type (INS, DEL, MIXED, NONE) |

## How It Works

1. svtig/SV sequences are aligned to haplotype assemblies using minimap2 (`-cx asm10`)
2. PAF alignments are parsed to extract mapping ratio, alignment identity, MAPQ, edit distance, alignment score, and SV counts from CIGAR strings
3. Supplementary alignments are accumulated into a combined mapping ratio (`--concordance`) or ignored (`--validate`)
4. Structural variants (>50bp insertions/deletions) are detected from CIGAR and classified by type (INS, DEL, MIXED)
5. Each svtig/SV is marked concordant or discordant based on configurable thresholds
6. VeriScore is computed as the geometric mean of mapping ratio, alignment identity, and Phred-transformed MAPQ
7. Unmapped sequences are included in the report with VeriScore = 0

## Testing

```bash
cd test
bash run_test.sh
```

VeriScore mathematical validation:
```bash
python3 test/veriscore_analysis.py
```

## License

TBD
