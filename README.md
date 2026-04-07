# veritig

Sequence-level verification of structural variant assemblies against haplotype-resolved references.

veritig aligns svtig sequences to haplotype assemblies via minimap2, computes concordance metrics and VeriScore, and provides tools for quality assessment and filtering.

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
| `--stats` | Descriptive statistics of svtig sets (no assembly needed) |
| `--concordance` | Concordance analysis against haplotype assemblies |
| `--validate` | De novo haplotype assignment from a single unphased FASTA |
| `--filter` | Quality-based filtering, outputs concordant svtigs as FASTA |
| `--compare` | Side-by-side comparison of two svtig sets |

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
```

## Documentation

See the [Wiki](https://github.com/asylvz/veritig/wiki) for detailed documentation:

- [Modes and Usage](https://github.com/asylvz/veritig/wiki/Modes-and-Usage) — all modes, parameters, and examples
- [VeriScore](https://github.com/asylvz/veritig/wiki/VeriScore) — formulation, properties, and interpretation
- [Output Format](https://github.com/asylvz/veritig/wiki/Output-Format) — file structure, report columns, detailed output

## License

MIT
