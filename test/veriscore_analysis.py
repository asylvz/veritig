#!/usr/bin/env python3
"""
VeriScore Mathematical Analysis
================================
Tests the consistency and behavior of the VeriScore composite metric.

VeriScore = (map_ratio * aln_identity * P_correct)^(1/3)

where P_correct = 1 - 10^(-MAPQ/10)   [Phred-scale probability of correct mapping]

This is the geometric mean of three quality dimensions:
  1. map_ratio    : fraction of query covered by alignment (completeness)
  2. aln_identity : fraction of matching bases in alignment (accuracy)
  3. P_correct    : probability of correct mapping position (confidence)

Why geometric mean?
  - All three are [0,1] ratio-scale measurements
  - Geometric mean of ratios is statistically appropriate (unlike arithmetic mean)
  - A low value in ANY dimension dominates: you can't compensate bad identity with good coverage
  - Produces a single [0,1] score that's interpretable as "overall mapping quality"

Why Phred normalization for MAPQ?
  - MAPQ is defined as -10*log10(P_error), so P_correct = 1 - 10^(-MAPQ/10)
  - This converts the log-scale MAPQ back to its natural probability scale
  - Matches the [0,1] range of the other two components
  - Preserves the non-linear nature: MAPQ 0->0.0, 10->0.9, 20->0.99, 60->~1.0
"""

import math

def phred_to_prob(mapq):
    """Convert MAPQ (Phred-scale) to probability of correct mapping."""
    if mapq <= 0:
        return 0.0
    return 1.0 - math.pow(10.0, -mapq / 10.0)

def veriscore(map_ratio, aln_identity, mapq):
    """Compute VeriScore as geometric mean of three quality components."""
    p_correct = phred_to_prob(mapq)
    if map_ratio <= 0 or aln_identity <= 0 or p_correct <= 0:
        return 0.0
    return math.pow(map_ratio * aln_identity * p_correct, 1.0/3.0)


print("=" * 90)
print("TEST 1: Monotonicity — Score increases with each metric (others fixed)")
print("=" * 90)

print("\n--- Varying map_ratio (identity=0.95, MAPQ=60) ---")
print(f"{'map_ratio':>12} {'P_correct':>12} {'VeriScore':>12}")
for mr in [0.50, 0.60, 0.70, 0.80, 0.85, 0.90, 0.95, 1.00]:
    vs = veriscore(mr, 0.95, 60)
    print(f"{mr:12.2f} {phred_to_prob(60):12.6f} {vs:12.4f}")

print("\n--- Varying aln_identity (map_ratio=0.95, MAPQ=60) ---")
print(f"{'aln_identity':>12} {'P_correct':>12} {'VeriScore':>12}")
for ai in [0.80, 0.85, 0.90, 0.92, 0.95, 0.97, 0.99, 1.00]:
    vs = veriscore(0.95, ai, 60)
    print(f"{ai:12.2f} {phred_to_prob(60):12.6f} {vs:12.4f}")

print("\n--- Varying MAPQ (map_ratio=0.95, identity=0.95) ---")
print(f"{'MAPQ':>12} {'P_correct':>12} {'VeriScore':>12}")
for mq in [0, 1, 3, 5, 10, 20, 30, 40, 50, 60]:
    vs = veriscore(0.95, 0.95, mq)
    pc = phred_to_prob(mq)
    print(f"{mq:12d} {pc:12.6f} {vs:12.4f}")


print("\n" + "=" * 90)
print("TEST 2: Dominance — A single bad metric pulls the score down")
print("=" * 90)
print(f"\n{'Scenario':>40} {'map_r':>7} {'ident':>7} {'MAPQ':>6} {'VeriScore':>10} {'Concordant':>12}")

scenarios = [
    ("Perfect mapping",                   1.00, 1.00, 60),
    ("Excellent all-around",              0.98, 0.99, 60),
    ("Good mapping",                      0.95, 0.95, 40),
    ("Typical concordant",                0.90, 0.92, 30),
    ("Borderline concordant",             0.85, 0.85, 20),
    ("High coverage, low identity",       0.98, 0.75, 60),
    ("Low coverage, high identity",       0.70, 0.99, 60),
    ("Good metrics, ambiguous mapping",   0.95, 0.95,  0),
    ("Good metrics, low confidence",      0.95, 0.95,  5),
    ("Mediocre all-around",               0.80, 0.80, 10),
    ("Poor mapping",                      0.50, 0.70, 10),
    ("Very poor",                         0.30, 0.60,  5),
]

for name, mr, ai, mq in scenarios:
    vs = veriscore(mr, ai, mq)
    conc = "yes" if (mr >= 0.85 and ai >= 0.85) else "no"
    print(f"{name:>40} {mr:7.2f} {ai:7.2f} {mq:6d} {vs:10.4f} {conc:>12}")


print("\n" + "=" * 90)
print("TEST 3: Score vs Concordance threshold — Does VeriScore refine binary classification?")
print("=" * 90)
print("\nSVtigs near the concordance boundary (map_ratio=0.85, identity=0.85):")
print(f"{'map_r':>7} {'ident':>7} {'MAPQ':>6} {'VeriScore':>10} {'Concordant':>12} {'Interpretation':>30}")

boundary_cases = [
    (0.86, 0.86, 60, "Just concordant, confident"),
    (0.86, 0.86, 10, "Just concordant, low confidence"),
    (0.86, 0.86,  3, "Just concordant, very low confidence"),
    (0.84, 0.86, 60, "Just discordant (map_ratio)"),
    (0.86, 0.84, 60, "Just discordant (identity)"),
    (0.90, 0.90, 60, "Concordant, high quality"),
    (0.90, 0.90,  5, "Concordant, poor confidence"),
]

for mr, ai, mq, interp in boundary_cases:
    vs = veriscore(mr, ai, mq)
    conc = "yes" if (mr >= 0.85 and ai >= 0.85) else "no"
    print(f"{mr:7.2f} {ai:7.2f} {mq:6d} {vs:10.4f} {conc:>12}   {interp}")


print("\n" + "=" * 90)
print("TEST 4: Comparison with alternative formulas")
print("=" * 90)

def arithmetic_mean(mr, ai, mq):
    return (mr + ai + phred_to_prob(mq)) / 3.0

def simple_product(mr, ai, mq):
    return mr * ai * phred_to_prob(mq)

def harmonic_mean(mr, ai, mq):
    pc = phred_to_prob(mq)
    if mr <= 0 or ai <= 0 or pc <= 0:
        return 0.0
    return 3.0 / (1.0/mr + 1.0/ai + 1.0/pc)

test_cases = [
    ("Excellent",    0.98, 0.99, 60),
    ("Good",         0.90, 0.92, 40),
    ("Mixed quality",0.95, 0.70, 60),
    ("Low MAPQ",     0.95, 0.95,  3),
    ("Poor",         0.50, 0.60,  5),
]

print(f"\n{'Case':>15} {'Arith':>8} {'Geometric':>10} {'Harmonic':>10} {'Product':>10}")
for name, mr, ai, mq in test_cases:
    am = arithmetic_mean(mr, ai, mq)
    gm = veriscore(mr, ai, mq)
    hm = harmonic_mean(mr, ai, mq)
    sp = simple_product(mr, ai, mq)
    print(f"{name:>15} {am:8.4f} {gm:10.4f} {hm:10.4f} {sp:10.4f}")

print("\nKey observations:")
print("  - Arithmetic mean is too lenient: one bad metric is masked by good ones")
print("  - Simple product over-penalizes: 0.9 * 0.9 * 0.9 = 0.729 (too harsh)")
print("  - Harmonic mean over-penalizes the weakest component")
print("  - Geometric mean provides balanced penalty without being extreme")


print("\n" + "=" * 90)
print("TEST 5: MAPQ Phred normalization — Diminishing returns at high MAPQ")
print("=" * 90)
print(f"\n{'MAPQ':>6} {'P_error':>12} {'P_correct':>12} {'Interpretation':>35}")

for mq in [0, 1, 2, 3, 5, 10, 15, 20, 30, 40, 50, 60]:
    if mq == 0:
        pe = 1.0
    else:
        pe = math.pow(10, -mq/10.0)
    pc = phred_to_prob(mq)

    if mq == 0:
        interp = "No mapping confidence"
    elif mq < 5:
        interp = f"1 in {int(1/pe)} chance of error"
    elif mq < 20:
        interp = f"1 in {int(1/pe)} chance of error"
    else:
        interp = f"1 in {int(1/pe):,} chance of error"

    print(f"{mq:6d} {pe:12.6f} {pc:12.6f} {interp:>35}")

print("\nThis shows why Phred normalization is correct:")
print("  - MAPQ 10 vs 20 is a bigger jump in P_correct (0.90 vs 0.99) than")
print("    MAPQ 40 vs 60 (0.9999 vs 0.999999) — diminishing returns are natural")
print("  - A linear MAPQ/60 normalization would be WRONG: it would say MAPQ 30")
print("    is 'half as good' as MAPQ 60, when in reality both are near-perfect")


print("\n" + "=" * 90)
print("TEST 6: Edge cases")
print("=" * 90)
edge_cases = [
    ("All perfect",        1.00, 1.00, 60, "Should be ~1.0"),
    ("MAPQ=0",             0.95, 0.95,  0, "Should be 0.0 (unmappable)"),
    ("Zero map_ratio",     0.00, 0.95, 60, "Should be 0.0 (unmapped)"),
    ("Zero identity",      0.95, 0.00, 60, "Should be 0.0 (no match)"),
    ("All minimal",        0.01, 0.01,  1, "Near zero"),
    ("MAPQ=255 (unavail)", 0.95, 0.95, 255,"MAPQ unavailable in SAM spec"),
]

print(f"\n{'Case':>25} {'map_r':>7} {'ident':>7} {'MAPQ':>6} {'VeriScore':>10} {'Note':>35}")
for name, mr, ai, mq, note in edge_cases:
    vs = veriscore(mr, ai, mq)
    print(f"{name:>25} {mr:7.2f} {ai:7.2f} {mq:6d} {vs:10.4f} {note:>35}")

print("\nNote: MAPQ=255 means 'unavailable' in SAM spec. In PAF format from")
print("minimap2, MAPQ=255 is rare but should be handled. We cap MAPQ at 60")
print("in the implementation to avoid P_correct approaching 1.0 too closely.")


print("\n" + "=" * 90)
print("TEST 7: Ranking consistency — Does VeriScore rank SVtigs as a human expert would?")
print("=" * 90)

ranking_cases = [
    ("SVtig_A", 0.99, 0.99, 60, "Excellent: complete, accurate, unique"),
    ("SVtig_B", 0.95, 0.95, 40, "Very good: slightly less complete"),
    ("SVtig_C", 0.92, 0.90, 30, "Good: moderate quality"),
    ("SVtig_D", 0.88, 0.87, 20, "Borderline: just above thresholds"),
    ("SVtig_E", 0.95, 0.75, 60, "Problematic: good coverage but low identity"),
    ("SVtig_F", 0.70, 0.98, 60, "Problematic: high identity but partial mapping"),
    ("SVtig_G", 0.95, 0.95,  3, "Suspicious: good metrics but ambiguous placement"),
    ("SVtig_H", 0.60, 0.60, 10, "Poor: low across all metrics"),
]

print(f"\n{'Name':>10} {'map_r':>7} {'ident':>7} {'MAPQ':>6} {'VeriScore':>10}")
ranked = sorted(ranking_cases, key=lambda x: veriscore(x[1], x[2], x[3]), reverse=True)
for name, mr, ai, mq, desc in ranked:
    vs = veriscore(mr, ai, mq)
    print(f"{name:>10} {mr:7.2f} {ai:7.2f} {mq:6d} {vs:10.4f}   {desc}")

print("\nRanking matches expert intuition:")
print("  1. SVtig_A/B/C rank highest — good all-around quality")
print("  2. SVtig_D ranks above E/F — balanced metrics beat one strong + one weak")
print("  3. SVtig_E and SVtig_F rank similarly — different deficiencies, similar overall quality")
print("  4. SVtig_G penalized despite good map_ratio/identity — MAPQ matters for confidence")
print("  5. SVtig_H ranks lowest — poor across the board")


print("\n" + "=" * 90)
print("TEST 8: MAPQ Floor Analysis — Impact of floor=3 vs floor=0 (current)")
print("=" * 90)

def veriscore_floor(map_ratio, aln_identity, mapq, floor=3):
    """VeriScore with MAPQ floor."""
    effective_mapq = max(mapq, floor)
    p_correct = 1.0 - math.pow(10.0, -effective_mapq / 10.0)
    if map_ratio <= 0 or aln_identity <= 0:
        return 0.0
    return math.pow(map_ratio * aln_identity * p_correct, 1.0/3.0)

print("\n--- 8a: Same SVtig, varying MAPQ — floor=0 vs floor=3 ---")
print(f"{'MAPQ':>6} {'P_corr(raw)':>12} {'P_corr(fl3)':>12} {'VS(raw)':>10} {'VS(floor3)':>10} {'Delta':>8}")
for mq in [0, 1, 2, 3, 5, 10, 20, 30, 60]:
    vs_raw = veriscore(0.95, 0.95, mq)
    vs_fl = veriscore_floor(0.95, 0.95, mq, floor=3)
    pc_raw = phred_to_prob(mq) if mq > 0 else 0.0
    pc_fl = phred_to_prob(max(mq, 3))
    delta = vs_fl - vs_raw
    print(f"{mq:6d} {pc_raw:12.6f} {pc_fl:12.6f} {vs_raw:10.4f} {vs_fl:10.4f} {delta:+8.4f}")

print("\nKey insight: Floor only affects MAPQ 0-2. At MAPQ≥3 scores are identical.")
print("The floor converts 'no information' to 'uninformative prior' (P≈0.5).")

print("\n--- 8b: Ranking preservation — Does floor=3 change relative order? ---")
ranking_cases = [
    ("A: Perfect",           0.99, 0.99, 60),
    ("B: Good, MAPQ=0",     0.95, 0.95,  0),
    ("C: Good, MAPQ=3",     0.95, 0.95,  3),
    ("D: Good, MAPQ=10",    0.95, 0.95, 10),
    ("E: Mediocre, MAPQ=60",0.80, 0.80, 60),
    ("F: Mediocre, MAPQ=0", 0.80, 0.80,  0),
    ("G: Bad, MAPQ=60",     0.60, 0.60, 60),
    ("H: Bad, MAPQ=0",      0.60, 0.60,  0),
]

print(f"\n{'Case':>25} {'VS(raw)':>10} {'Rank':>5} {'VS(floor3)':>12} {'Rank':>5} {'Order preserved':>16}")
raw_scores = [(name, veriscore(mr, ai, mq)) for name, mr, ai, mq in ranking_cases]
fl_scores = [(name, veriscore_floor(mr, ai, mq, 3)) for name, mr, ai, mq in ranking_cases]
raw_ranked = sorted(raw_scores, key=lambda x: -x[1])
fl_ranked = sorted(fl_scores, key=lambda x: -x[1])

# Build rank maps
raw_rank = {name: i+1 for i, (name, _) in enumerate(raw_ranked)}
fl_rank = {name: i+1 for i, (name, _) in enumerate(fl_ranked)}

for name, mr, ai, mq in ranking_cases:
    vs_r = veriscore(mr, ai, mq)
    vs_f = veriscore_floor(mr, ai, mq, 3)
    rr = raw_rank[name]
    fr = fl_rank[name]
    preserved = "YES" if rr == fr else f"NO ({rr}->{fr})"
    print(f"{name:>25} {vs_r:10.4f} {rr:5d} {vs_f:12.4f} {fr:5d} {preserved:>16}")

print("\n--- 8c: Critical question — Does floor=3 create rank inversions? ---")
print("A rank inversion means: SVtig X was better than Y without floor, but worse with floor.")
print("This would be a mathematical problem that invalidates the metric.\n")

# Check all pairs
inversions = []
for i, (n1, mr1, ai1, mq1) in enumerate(ranking_cases):
    for j, (n2, mr2, ai2, mq2) in enumerate(ranking_cases):
        if i >= j:
            continue
        vs_r1 = veriscore(mr1, ai1, mq1)
        vs_r2 = veriscore(mr2, ai2, mq2)
        vs_f1 = veriscore_floor(mr1, ai1, mq1, 3)
        vs_f2 = veriscore_floor(mr2, ai2, mq2, 3)
        if (vs_r1 > vs_r2 and vs_f1 < vs_f2) or (vs_r1 < vs_r2 and vs_f1 > vs_f2):
            # But skip ties (both were 0.0)
            if vs_r1 == vs_r2 == 0.0:
                continue
            inversions.append((n1, n2, vs_r1, vs_r2, vs_f1, vs_f2))

if inversions:
    print(f"Found {len(inversions)} rank inversion(s):")
    for n1, n2, vr1, vr2, vf1, vf2 in inversions:
        print(f"  {n1} vs {n2}: raw={vr1:.4f}>{vr2:.4f}, floor={vf1:.4f} vs {vf2:.4f}")
else:
    print("NO rank inversions found among non-MAPQ=0 pairs.")
    print("(MAPQ=0 SVtigs that were tied at 0.0 now have distinct scores — this is")
    print(" expected and desirable: they were previously indistinguishable.)")

print("\n--- 8d: Sensitivity analysis — Different floor values ---")
print(f"\n{'Floor':>6} {'P_correct':>10} {'VS(0.95,0.95)':>14} {'VS(0.80,0.80)':>14} {'Interpretation':>30}")
for floor_val in [1, 2, 3, 5, 10]:
    pc = phred_to_prob(floor_val)
    vs_good = veriscore_floor(0.95, 0.95, 0, floor_val)
    vs_med = veriscore_floor(0.80, 0.80, 0, floor_val)
    if floor_val == 1:
        interp = "Too harsh (P=0.21)"
    elif floor_val == 2:
        interp = "Harsh (P=0.37)"
    elif floor_val == 3:
        interp = "Uninformative prior (P≈0.50)"
    elif floor_val == 5:
        interp = "Generous (P=0.68)"
    elif floor_val == 10:
        interp = "Too generous (P=0.90)"
    print(f"{floor_val:6d} {pc:10.4f} {vs_good:14.4f} {vs_med:14.4f} {interp:>30}")

print("\nFloor=3 is the natural choice:")
print("  - P_correct = 1 - 10^(-3/10) = 0.4988 ≈ 0.50")
print("  - This is exactly the maximum entropy (uninformative) prior")
print("  - Equivalent to saying 'this mapping has coin-flip confidence'")
print("  - Any higher floor would be over-confident given zero information")
print("  - Any lower floor would still over-penalize (P<0.37)")

print("\n--- 8e: Proof that floor=3 preserves monotonicity ---")
print("For VeriScore to be monotonic, increasing any input must increase output.")
print("The floor only changes behavior at MAPQ ∈ {0, 1, 2} → all become MAPQ=3.")
print("Within this range, all three MAPQ values produce the SAME score (no inversion).")
print("At MAPQ=3, the original formula takes over. Since the original is monotonic")
print("for MAPQ≥1, and the floor creates a flat region at [0,2], the combined")
print("function is monotonically non-decreasing everywhere. QED.\n")

# Verify monotonicity numerically
print("Numerical verification (map_ratio=0.95, identity=0.95):")
print(f"{'MAPQ':>6} {'VeriScore(floor=3)':>20} {'Monotonic':>10}")
prev_vs = 0.0
for mq in range(0, 61):
    vs = veriscore_floor(0.95, 0.95, mq, 3)
    mono = "YES" if vs >= prev_vs - 1e-10 else "NO !!!"
    if mq <= 5 or mq % 10 == 0:
        print(f"{mq:6d} {vs:20.6f} {mono:>10}")
    prev_vs = vs


print("\n" + "=" * 90)
print("MATHEMATICAL PROPERTIES SUMMARY")
print("=" * 90)
print("""
VeriScore = (R * I * C)^(1/3)

where:
  R = map_ratio       ∈ [0, 1]   (alignment completeness)
  I = aln_identity    ∈ [0, 1]   (alignment accuracy)
  C = 1-10^(-Q/10)   ∈ [0, 1)   (mapping confidence, Phred-transformed)

Properties:
  1. Range: [0, 1] — interpretable as overall quality fraction
  2. Monotonic: increasing in all three inputs
  3. Symmetric: no implicit weighting between components
  4. Zero-dominant: if ANY component is 0, VeriScore = 0
  5. Concave: diminishing returns — going from 0.8→0.9 helps more than 0.9→1.0
  6. Scale-invariant: geometric mean is the correct average for ratios
  7. Multiplicative penalty: a 10% drop in any metric reduces score by ~3.4%
     (since (0.9)^(1/3) ≈ 0.966)

Relationship to concordance:
  - Concordance is binary: concordant iff R >= T_r AND I >= T_i
  - VeriScore is continuous: provides ranking within concordant AND discordant groups
  - At default thresholds (0.85/0.85), concordant SVtigs with MAPQ≥20 have VeriScore≥0.87
  - VeriScore enables soft filtering: users can set a VeriScore threshold instead of
    or in addition to the binary concordance thresholds

Implementation note:
  - MAPQ is capped at 60 to avoid numerical issues (P_correct → 1.0)
  - MAPQ = 0 yields VeriScore = 0 (mapping has no confidence)
  - MAPQ = 255 (unavailable in SAM spec) is treated as 0 in PAF context
""")
