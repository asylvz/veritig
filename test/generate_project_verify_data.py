#!/usr/bin/env python3
"""Generate synthetic test data for --project and --verify modes.

Produces:
  test/ref_test.fa        : 1 chromosome (chr_test, 5000 bp)
  test/svtigs_test.fa     : 4 svtigs (INS, DEL, INV, DUP)
  test/h1_test.fa         : haplotype 1 assembly (carries INS + INV)
  test/h2_test.fa         : haplotype 2 assembly (carries DEL + DUP)
  test/expected_pv.vcf    : VCF (4 records that --project should output)

SV layout on the reference (anchor positions, 1-based):
   600  INS  (120 bp insertion)
  1200  DEL  (100 bp deletion)
  2400  INV  (300 bp inverted segment)
  3600  DUP  (150 bp tandem duplication)

Sequences are pseudo-random with a fixed seed so the data is reproducible.
"""

import random

random.seed(42)
BASES = "ACGT"
COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def random_seq(n):
    return "".join(random.choice(BASES) for _ in range(n))


def revcomp(s):
    return "".join(COMP[c] for c in reversed(s))


def write_fasta(path, name, seq, line_width=80):
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), line_width):
            f.write(seq[i:i + line_width] + "\n")


def write_multi_fasta(path, entries, line_width=80):
    with open(path, "w") as f:
        for name, seq in entries:
            f.write(f">{name}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i + line_width] + "\n")


# Reference: 5000 bp single chromosome
REF_LEN = 5000
ref = random_seq(REF_LEN)

# SV definitions (0-based reference coordinates)
INS_POS = 600        # 1-based anchor = 600
INS_SEQ = random_seq(120)

DEL_POS = 1200
DEL_LEN = 100

INV_START = 2400
INV_LEN = 300        # inverted segment 2400..2700

DUP_START = 3600     # source copy of the DUP unit lives at 3600..3750
DUP_LEN = 150        # tandem duplication: an extra copy is appended in svtig/H2

# Svtigs (each spans ~1000-1300 bp of reference with the variant applied)
svtig_ins = ref[0:INS_POS] + INS_SEQ + ref[INS_POS:1200]                       # ~1320 bp
svtig_del = ref[800:DEL_POS] + ref[DEL_POS + DEL_LEN:1900]                     # ~800 bp
svtig_inv = (
    ref[1900:INV_START]                                # left flank (~500 bp)
    + revcomp(ref[INV_START:INV_START + INV_LEN])      # inverted segment (300 bp)
    + ref[INV_START + INV_LEN:2900]                    # right flank (~200 bp)
)                                                                              # ~1000 bp
svtig_dup = (
    ref[3100:DUP_START + DUP_LEN]                       # left flank + first copy
    + ref[DUP_START:DUP_START + DUP_LEN]                # extra copy (tandem)
    + ref[DUP_START + DUP_LEN:4400]                     # right flank
)                                                                              # ~1450 bp

# Haplotype 1 carries INS + INV (positions shifted by INS_SEQ insertion length)
h1 = (
    ref[0:INS_POS]
    + INS_SEQ
    + ref[INS_POS:INV_START]
    + revcomp(ref[INV_START:INV_START + INV_LEN])
    + ref[INV_START + INV_LEN:]
)

# Haplotype 2 carries DEL + DUP
h2 = (
    ref[0:DEL_POS]
    + ref[DEL_POS + DEL_LEN:DUP_START + DUP_LEN]
    + ref[DUP_START:DUP_START + DUP_LEN]            # extra copy
    + ref[DUP_START + DUP_LEN:]
)

write_fasta("test/ref_test.fa", "chr_test", ref)
write_multi_fasta(
    "test/svtigs_test.fa",
    [
        ("svtig_ins", svtig_ins),
        ("svtig_del", svtig_del),
        ("svtig_inv", svtig_inv),
        ("svtig_dup", svtig_dup),
    ],
)
write_fasta("test/h1_test.fa", "chr_test", h1)
write_fasta("test/h2_test.fa", "chr_test", h2)

# Reference VCF (human-readable; --project writes equivalent records with 'N' anchors)
with open("test/expected_pv.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    f.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">\n")
    f.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length\">\n")
    f.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate\">\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    f.write(f"chr_test\t{INS_POS}\t.\tN\tN{INS_SEQ}\t.\tPASS\tSVTYPE=INS;SVLEN={len(INS_SEQ)};END={INS_POS}\n")
    deleted_seq = ref[DEL_POS:DEL_POS + DEL_LEN]
    f.write(f"chr_test\t{DEL_POS}\t.\tN{deleted_seq}\tN\t.\tPASS\tSVTYPE=DEL;SVLEN={DEL_LEN};END={DEL_POS + DEL_LEN}\n")
    f.write(f"chr_test\t{INV_START}\t.\tN\t<INV>\t.\tPASS\tSVTYPE=INV;SVLEN={INV_LEN};END={INV_START + INV_LEN}\n")
    f.write(f"chr_test\t{DUP_START}\t.\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;SVLEN={DUP_LEN};END={DUP_START + DUP_LEN}\n")

print(f"ref:          {REF_LEN} bp")
print(f"svtig_ins:    {len(svtig_ins)} bp (120 bp INS)")
print(f"svtig_del:    {len(svtig_del)} bp (100 bp DEL)")
print(f"svtig_inv:    {len(svtig_inv)} bp (300 bp INV)")
print(f"svtig_dup:    {len(svtig_dup)} bp (150 bp DUP)")
print(f"h1:           {len(h1)} bp (carries INS + INV)")
print(f"h2:           {len(h2)} bp (carries DEL + DUP)")
print("Written: test/ref_test.fa, svtigs_test.fa, h1_test.fa, h2_test.fa, expected_pv.vcf")
