#!/usr/bin/env python3
"""Generate synthetic test data for --project and --verify modes.

Produces:
  test/ref_test.fa        : 1 chromosome (chr_test, 2000 bp)
  test/svtigs_test.fa     : 2 svtigs (one with 120 bp INS, one with 100 bp DEL)
  test/h1_test.fa         : haplotype 1 assembly (carries the INS at the same locus)
  test/h2_test.fa         : haplotype 2 assembly (carries the DEL at the same locus)
  test/expected_pv.vcf    : VCF (INS+DEL records that --project should output)

Sequences are pseudo-random with a fixed seed so the data is reproducible.
"""

import random

random.seed(42)
BASES = "ACGT"

def random_seq(n):
    return "".join(random.choice(BASES) for _ in range(n))

def write_fasta(path, name, seq, line_width=80):
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), line_width):
            f.write(seq[i:i + line_width] + "\n")

# Reference: 2000 bp single chromosome
REF_LEN = 2000
ref = random_seq(REF_LEN)

# SV definitions
INS_POS = 600       # 0-based; insertion anchored just after this base (1-based pos = 600)
INS_SEQ = random_seq(120)
DEL_POS = 1200      # 0-based start of the deleted span
DEL_LEN = 100

# Svtig with INS: ref[0:INS_POS] + INS_SEQ + ref[INS_POS:1100]
# (cover ~1100 bp of reference, includes 120 bp insertion)
svtig_ins = ref[0:INS_POS] + INS_SEQ + ref[INS_POS:1100]

# Svtig with DEL: ref[800:DEL_POS] + ref[DEL_POS+DEL_LEN:1700]
# (cover ~900 bp of reference, includes 100 bp deletion)
svtig_del = ref[800:DEL_POS] + ref[DEL_POS + DEL_LEN:1700]

# H1 carries the INS, H2 carries the DEL (simulating heterozygous variants)
# H1: full reference but with INS inserted
h1 = ref[0:INS_POS] + INS_SEQ + ref[INS_POS:]
# H2: full reference but with DEL removed
h2 = ref[0:DEL_POS] + ref[DEL_POS + DEL_LEN:]

write_fasta("test/ref_test.fa", "chr_test", ref)

# Two svtigs in one FASTA
with open("test/svtigs_test.fa", "w") as f:
    f.write(">svtig_ins\n")
    for i in range(0, len(svtig_ins), 80):
        f.write(svtig_ins[i:i + 80] + "\n")
    f.write(">svtig_del\n")
    for i in range(0, len(svtig_del), 80):
        f.write(svtig_del[i:i + 80] + "\n")

write_fasta("test/h1_test.fa", "chr_test", h1)
write_fasta("test/h2_test.fa", "chr_test", h2)

# Expected VCF (for reference; --project should produce equivalent records)
# Note: --project uses 'N' for REF/ALT when constructing records (see project.cpp).
# We emit explicit sequences here just as a human-readable reference document.
with open("test/expected_pv.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    f.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">\n")
    f.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length\">\n")
    f.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate\">\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    # INS at 1-based position INS_POS (anchor base)
    f.write(f"chr_test\t{INS_POS}\t.\tN\tN{INS_SEQ}\t.\tPASS\tSVTYPE=INS;SVLEN={len(INS_SEQ)};END={INS_POS}\n")
    # DEL anchored at DEL_POS (1-based), spans DEL_LEN
    deleted_seq = ref[DEL_POS:DEL_POS + DEL_LEN]
    f.write(f"chr_test\t{DEL_POS}\t.\tN{deleted_seq}\tN\t.\tPASS\tSVTYPE=DEL;SVLEN={DEL_LEN};END={DEL_POS + DEL_LEN}\n")

print(f"ref:    {REF_LEN} bp")
print(f"svtig_ins: {len(svtig_ins)} bp (120 bp INS)")
print(f"svtig_del: {len(svtig_del)} bp (100 bp DEL)")
print(f"h1:     {len(h1)} bp (carries INS)")
print(f"h2:     {len(h2)} bp (carries DEL)")
print("Written: test/ref_test.fa, svtigs_test.fa, h1_test.fa, h2_test.fa, expected_pv.vcf")
