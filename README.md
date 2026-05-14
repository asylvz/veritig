# veritig

[![tests](https://github.com/asylvz/veritig/actions/workflows/test.yml/badge.svg)](https://github.com/asylvz/veritig/actions/workflows/test.yml)

Sequence-level verification of structural variant assemblies against haplotype-resolved references.

veritig aligns assembled SV sequences to haplotype-resolved references and computes concordance metrics and VeriScore, complementing VCF-based benchmarking with sequence-level evaluation. Companion modes provide filtering and comparison, assembly-based SV calling against a linear reference (`--project`), and gap-aware verification of VCF SV calls (`--verify`).

An svtig (structural variant contig) is a sequence-resolved local assembly of an SV-bearing region together with flanking context, typically a few kb to tens of kb of contiguous sequence that captures the full variant allele. Unlike a VCF record, which only encodes coordinates and length, an svtig carries the actual variant sequence.

svtigs are currently produced by [SVarp](https://github.com/asylvz/svarp), a pangenome-based SV caller that clusters reads at candidate SV loci on a pangenome graph and locally assembles them. veritig accepts FASTA input: one file for unphased data, or two files (`--svtig1` for H1, `--svtig2` for H2) for haplotype-phased data. Although designed for svtigs, any FASTA of SV-bearing or local-assembly sequences is accepted, and the same alignment-based metrics apply.

## Installation

```bash
git clone https://github.com/asylvz/veritig.git
cd veritig
make
```

**Requirements:** C++17 compiler (g++), [minimap2](https://github.com/lh3/minimap2) in PATH (not needed for `--stats`)

## Modes

| Mode | Description |
|------|-------------|
| `--concordance` | Concordance analysis against haplotype assemblies |
| `--validate` | De novo haplotype assignment from a single unphased FASTA |
| `--filter` | Quality-based filtering, outputs concordant svtigs as FASTA |
| `--compare` | Side-by-side comparison of two svtig sets |
| `--project` | Assembly-based SV caller (svtigs or FASTA → single-sample VCF) |
| `--verify` | Gap-aware sequence-level verification of VCF SV calls |
| `--stats` | Descriptive statistics of svtig sets (no assembly needed) |

## Quick Examples

```bash
# Statistics (no assembly needed)
veritig --stats --svtig1 svtigs_h1.fa --svtig2 svtigs_h2.fa -o output/ -s sample

# Concordance with phasing
veritig --concordance --h1 hap1.fa --svtig1 svtigs_h1.fa \
        --h2 hap2.fa --svtig2 svtigs_h2.fa --phase -o output/ -s sample

# Validate (single unphased FASTA)
veritig --validate --h1 hap1.fa --h2 hap2.fa --fasta svs.fa -o output/ -s sample

# Filter low-quality svtigs
veritig --filter --h1 hap1.fa --svtig1 svtigs.fa -M 0.95 -I 0.95 -o output/ -s sample

# Compare two svtig sets
veritig --compare --h1 hap1.fa --svtig1 setA.fa --svtig2 setB.fa -o output/ -s sample

# Project svtigs to a single-sample VCF on a linear reference
# (accepts 1, 2, or 3 svtig FASTAs; phased GT when --svtig1 + --svtig2 are given)
veritig --project --ref ref.fa --svtig1 svtigs_h1.fa --svtig2 svtigs_h2.fa -o output/ -s sample

# Verify VCF SV calls against haplotype assemblies (gap-aware sequence-level check)
veritig --verify --vcf calls.vcf --ref ref.fa --h1 hap1.fa --h2 hap2.fa -o output/ -s sample
```

## Documentation

See the [Wiki](https://github.com/asylvz/veritig/wiki) for detailed documentation:

- [Modes and Usage](https://github.com/asylvz/veritig/wiki/Modes-and-Usage) — all modes, parameters, and examples
- [VeriScore](https://github.com/asylvz/veritig/wiki/VeriScore) — formulation, properties, and interpretation
- [Output Format](https://github.com/asylvz/veritig/wiki/Output-Format) — file structure, report columns, detailed output

## License

MIT
