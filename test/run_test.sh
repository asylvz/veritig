#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
VERITIG="$SCRIPT_DIR/../build/veritig"
TEST_DIR="$SCRIPT_DIR"
OUT_DIR="$TEST_DIR/output"

rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR/veritig_results/paf"

# Copy pre-made PAF files
cp "$TEST_DIR/H1_svtig.paf" "$OUT_DIR/veritig_results/paf/test.svtig.H1.paf"
cp "$TEST_DIR/H2_svtig.paf" "$OUT_DIR/veritig_results/paf/test.svtig.H2.paf"

echo "=== Test 1: concordance (default output) ==="
$VERITIG --concordance \
    --h1 dummy_h1.fa \
    --h2 dummy_h2.fa \
    --svtig1 "$TEST_DIR/svtigs.fa" \
    --skip-mapping \
    --out "$OUT_DIR/" \
    --sample test \
    --threads 4 \
    --min-map-ratio 0.85 \
    --min-identity 0.85

echo ""
echo "=== Output files (default) ==="
find "$OUT_DIR/veritig_results" -type f | sort

echo ""
echo "=== svtig Report ==="
cat "$OUT_DIR/veritig_results/test.svtig.report.tsv"

echo ""
echo "=== Concordance ==="
cat "$OUT_DIR/veritig_results/test.svtig.concordance.tsv"

echo ""
echo "=== Test 2: concordance (--detailed output) ==="
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR/veritig_results/paf"
cp "$TEST_DIR/H1_svtig.paf" "$OUT_DIR/veritig_results/paf/test.svtig.H1.paf"
cp "$TEST_DIR/H2_svtig.paf" "$OUT_DIR/veritig_results/paf/test.svtig.H2.paf"

$VERITIG --concordance \
    --h1 dummy_h1.fa \
    --h2 dummy_h2.fa \
    --svtig1 "$TEST_DIR/svtigs.fa" \
    --skip-mapping \
    --out "$OUT_DIR/" \
    --sample test \
    --threads 4 \
    --detailed

echo ""
echo "=== Output files (detailed) ==="
find "$OUT_DIR/veritig_results" -type f | sort

echo ""
echo "=== Details (H1) ==="
cat "$OUT_DIR/veritig_results/detailed/test.svtig.svtig1_H1.details.tsv"

echo ""
echo "=== Summary (H1) ==="
cat "$OUT_DIR/veritig_results/detailed/test.svtig.svtig1_H1.summary.tsv"

echo ""
echo "=== Size bins (H1) ==="
cat "$OUT_DIR/veritig_results/detailed/test.svtig.svtig1_H1.size_bins.tsv"

echo ""
echo "=== Test 3: validate (default output) ==="
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR/veritig_results/paf"
cp "$TEST_DIR/H1_sv.paf" "$OUT_DIR/veritig_results/paf/test.sv.H1.paf"
cp "$TEST_DIR/H2_sv.paf" "$OUT_DIR/veritig_results/paf/test.sv.H2.paf"

$VERITIG --validate \
    --h1 dummy_h1.fa \
    --h2 dummy_h2.fa \
    --fasta "$TEST_DIR/svs.fa" \
    --skip-mapping \
    --out "$OUT_DIR/" \
    --sample test \
    --threads 4

echo ""
echo "=== SV Report ==="
cat "$OUT_DIR/veritig_results/test.sv.report.tsv"

echo ""
echo "=== SV Concordance ==="
cat "$OUT_DIR/veritig_results/test.sv.concordance.tsv"

echo ""
echo "=== Test 4: validate (--detailed output) ==="
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR/veritig_results/paf"
cp "$TEST_DIR/H1_sv.paf" "$OUT_DIR/veritig_results/paf/test.sv.H1.paf"
cp "$TEST_DIR/H2_sv.paf" "$OUT_DIR/veritig_results/paf/test.sv.H2.paf"

$VERITIG --validate \
    --h1 dummy_h1.fa \
    --h2 dummy_h2.fa \
    --fasta "$TEST_DIR/svs.fa" \
    --skip-mapping \
    --out "$OUT_DIR/" \
    --sample test \
    --threads 4 \
    --detailed

echo ""
echo "=== SV Detailed files ==="
find "$OUT_DIR/veritig_results" -type f | sort

echo ""
echo "=== SV Details (H1) ==="
cat "$OUT_DIR/veritig_results/detailed/test.sv.H1.details.tsv"

echo ""
echo "=== Test 5: stats mode ==="
rm -rf "$OUT_DIR"

$VERITIG --stats \
    --svtig1 "$TEST_DIR/svtigs.fa" \
    --out "$OUT_DIR/" \
    --sample test

echo ""
echo "=== Stats output ==="
find "$OUT_DIR/veritig_results" -type f | sort
cat "$OUT_DIR/veritig_results/test.stats.tsv"

echo ""
echo "=== Test 6: filter mode ==="
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR/veritig_results/paf"
cp "$TEST_DIR/H1_svtig.paf" "$OUT_DIR/veritig_results/paf/test.filter.paf"

$VERITIG --filter \
    --h1 dummy_h1.fa \
    --svtig1 "$TEST_DIR/svtigs.fa" \
    --skip-mapping \
    --out "$OUT_DIR/" \
    --sample test \
    --threads 4

echo ""
echo "=== Filter output ==="
find "$OUT_DIR/veritig_results" -type f | sort
echo ""
echo "=== Filtered FASTA ==="
cat "$OUT_DIR/veritig_results/test.filtered.fa"

echo ""
echo "=== Test 7: compare mode ==="
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR/veritig_results/paf"
cp "$TEST_DIR/H1_svtig.paf" "$OUT_DIR/veritig_results/paf/test.compare.svtig1.paf"
cp "$TEST_DIR/H2_svtig.paf" "$OUT_DIR/veritig_results/paf/test.compare.svtig2.paf"

$VERITIG --compare \
    --h1 dummy_h1.fa \
    --svtig1 "$TEST_DIR/svtigs.fa" \
    --svtig2 "$TEST_DIR/svtigs.fa" \
    --skip-mapping \
    --out "$OUT_DIR/" \
    --sample test \
    --threads 4

echo ""
echo "=== Compare output ==="
cat "$OUT_DIR/veritig_results/test.compare.tsv"

echo ""
echo "=== DONE ==="
