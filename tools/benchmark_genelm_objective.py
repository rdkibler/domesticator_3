#!/usr/bin/env python3
"""
Benchmark GeneLM TIS objective in DNAChisel optimization.

Compares optimization with and without GeneLM objective:
1. How effective it is at reducing off-target TIS
2. How much it slows down optimization
"""

import sys
import time
import random
import numpy as np
from pathlib import Path
from Bio import Seq

from dnachisel import (
    DnaOptimizationProblem,
    Location,
    reverse_translate,
    EnforceTranslation,
    AvoidHairpins,
    AvoidPattern,
    EnforceGCContent,
    MaximizeCAI,
)
from specifications.MinimizeNumKmers import MinimizeNumKmers
from specifications.MinimizeOffTargetTISConfidence import MinimizeOffTargetTISConfidence


# UTR sequences
UTR_5PRIME = "taatacgactcactataggggaattgtgagcggataacaattcccctctagaaataattttgtttaactttaagaaggagatatacat"
UTR_3PRIME = "tgagatccggctgctaacaaagcccgaaaggaagctgagttggctgctgccaccgctgagcaataactagc"


def generate_synthetic_protein(length, met_content=0.1, ag_content=0.6):
    """Generate a synthetic protein starting with MET and ending with LEHHHHHH."""
    # Amino acids with A/G-rich codons (high TIS potential)
    ag_rich_aa = ["M", "G", "A", "R", "N", "K"]

    # Normal amino acids
    all_aa = list("ACDEFGHIKLMNPQRSTVWY")

    # Build protein: MET + random middle + LEHHHHHH
    middle_len = length - 2 - 8  # MET + LEHHHHHH = 2 + 8 = 10
    middle = []
    for i in range(middle_len):
        if random.random() < met_content:
            middle.append("M")
        elif random.random() < ag_content:
            middle.append(random.choice(ag_rich_aa))
        else:
            middle.append(random.choice(all_aa))

    protein = ["M"] + middle + ["L", "E", "H", "H", "H", "H", "H", "H"]
    return "".join(protein)


def create_optimization_problem(dna_sequence, gene_start, gene_end, use_genelm=False):
    """Create a DnaOptimizationProblem with or without GeneLM objective."""
    location = Location(gene_start, gene_end)

    constraints = []
    objectives = []

    # Constraints
    constraints.append(EnforceTranslation(location=location))
    constraints.append(AvoidPattern("AAAAAA", location=location))
    constraints.append(AvoidPattern("TTTTTT", location=location))
    constraints.append(AvoidPattern("CCCCCC", location=location))
    constraints.append(AvoidPattern("GGGGGG", location=location))
    constraints.append(EnforceGCContent(mini=0.25, maxi=0.8, window=50, location=location))
    constraints.append(EnforceGCContent(mini=0.4, maxi=0.65, location=location))

    # Objectives
    objectives.append(MinimizeNumKmers(k=8, boost=10, location=location))
    objectives.append(AvoidHairpins(boost=1.0, location=location))
    objectives.append(MaximizeCAI(species="e_coli", boost=1.0, location=location))

    # GeneLM objective if requested
    if use_genelm:
        objectives.append(
            MinimizeOffTargetTISConfidence(
                intended_tis_sites=[(gene_start, gene_start + 3)],
                confidence_threshold=0.85,
                window_size=30,
                location=location,
                boost=5.0
            )
        )

    problem = DnaOptimizationProblem(
        sequence=dna_sequence,
        constraints=constraints,
        objectives=objectives,
        logger=None
    )

    return problem


def create_evaluation_problem(sequence, location):
    """Create a minimal DnaOptimizationProblem for evaluation."""
    return DnaOptimizationProblem(
        sequence=sequence,
        constraints=[],
        objectives=[],
        logger=None
    )


def main():
    print("=" * 80)
    print("GeneLM TIS Objective Benchmark")
    print("=" * 80)

    # Generate test protein sequences (fewer to speed up)
    random.seed(42)
    num_tests = 5  # Reduced from 10 for faster benchmark
    test_proteins = [
        generate_synthetic_protein(length=random.randint(50, 200))
        for _ in range(num_tests)
    ]

    print(f"\nGenerated {num_tests} test protein sequences:")
    for i, protein in enumerate(test_proteins):
        print(f"  Test {i+1}: Length={len(protein)}, num_met: {protein.count('M')}, num_ag_rich: {sum(protein.count(aa) for aa in ['M', 'G', 'A', 'R', 'N', 'K'])}")

    # Results storage
    results = []

    for i, protein_seq in enumerate(test_proteins):
        print(f"\n--- Test {i+1}/{num_tests} (length: {len(protein_seq)}) ---")

        # Create full DNA sequence with UTRs
        dna_no_utr = reverse_translate(protein_seq)
        full_sequence = UTR_5PRIME + dna_no_utr + UTR_3PRIME

        # Gene location (after 5' UTR)
        gene_start = len(UTR_5PRIME)
        gene_end = gene_start + len(dna_no_utr)
        location = Location(gene_start, gene_end)

        # Run without GeneLM
        print("  Running optimization without GeneLM...")
        problem_no_genelm = create_optimization_problem(
            full_sequence, gene_start, gene_end, use_genelm=False
        )

        start_time = time.time()
        problem_no_genelm.resolve_constraints_by_random_mutations()
        problem_no_genelm.optimize()
        no_genelm_time = time.time() - start_time

        no_genelm_seq = str(problem_no_genelm.sequence)

        # Run with GeneLM
        print("  Running optimization with GeneLM...")
        problem_genelm = create_optimization_problem(
            full_sequence, gene_start, gene_end, use_genelm=True
        )

        start_time = time.time()
        problem_genelm.resolve_constraints_by_random_mutations()
        problem_genelm.optimize()
        genelm_time = time.time() - start_time

        genelm_seq = str(problem_genelm.sequence)

        # Evaluate both with GeneLM spec to compare off-target scores
        print("  Evaluating off-target TIS scores...")

        # Create evaluation problems
        eval_problem_no_genelm = create_evaluation_problem(no_genelm_seq, location)
        eval_problem_genelm = create_evaluation_problem(genelm_seq, location)

        # Initialize spec on both problems
        spec = MinimizeOffTargetTISConfidence(
            intended_tis_sites=[(gene_start, gene_start + 3)],
            confidence_threshold=0.85,
            window_size=30,
            location=location,
            boost=1.0
        )
        spec.initialize_on_problem(eval_problem_no_genelm)
        eval_no_genelm = spec.evaluate(eval_problem_no_genelm)
        no_genelm_score = eval_no_genelm.score

        spec.initialize_on_problem(eval_problem_genelm)
        eval_genelm = spec.evaluate(eval_problem_genelm)
        genelm_score = eval_genelm.score

        # Count high-confidence off-target sites from the evaluation message
        import re
        no_genelm_match = re.search(r'(\d+) sites >=', eval_no_genelm.message)
        genelm_match = re.search(r'(\d+) sites >=', eval_genelm.message)
        no_genelm_high_conf = int(no_genelm_match.group(1)) if no_genelm_match else 0
        genelm_high_conf = int(genelm_match.group(1)) if genelm_match else 0

        # Store results
        results.append({
            "test": i + 1,
            "protein_length": len(protein_seq),
            "dna_length": len(full_sequence),
            "no_genelm_time": no_genelm_time,
            "no_genelm_high_conf": no_genelm_high_conf,
            "no_genelm_score": no_genelm_score,
            "genelm_time": genelm_time,
            "genelm_high_conf": genelm_high_conf,
            "genelm_score": genelm_score,
            "time_diff": genelm_time - no_genelm_time,
            "time_pct_diff": (genelm_time - no_genelm_time) / no_genelm_time * 100 if no_genelm_time > 0 else 0,
            "conf_reduction": no_genelm_high_conf - genelm_high_conf,
            "score_improvement": genelm_score - no_genelm_score
        })

        print(f"  No GeneLM: {no_genelm_time:.1f}s, {no_genelm_high_conf} high-conf TIS, score={no_genelm_score:.3f}")
        print(f"  GeneLM:    {genelm_time:.1f}s, {genelm_high_conf} high-conf TIS, score={genelm_score:.3f}")

    # Summary
    print("\n" + "=" * 80)
    print("Summary")
    print("=" * 80)

    avg_time_no = sum(r["no_genelm_time"] for r in results) / len(results)
    avg_time_genelm = sum(r["genelm_time"] for r in results) / len(results)
    avg_conf_no = sum(r["no_genelm_high_conf"] for r in results) / len(results)
    avg_conf_genelm = sum(r["genelm_high_conf"] for r in results) / len(results)
    avg_time_pct = sum(r["time_pct_diff"] for r in results) / len(results)
    avg_conf_reduction = sum(r["conf_reduction"] for r in results) / len(results)
    avg_score_no = sum(r["no_genelm_score"] for r in results) / len(results)
    avg_score_genelm = sum(r["genelm_score"] for r in results) / len(results)
    avg_score_improvement = sum(r["score_improvement"] for r in results) / len(results)

    print(f"\n{'Metric':<35} {'No GeneLM':>12} {'With GeneLM':>12} {'Diff':>12}")
    print("-" * 73)
    print(f"{'Optimization time (s)':<35} {avg_time_no:>12.1f} {avg_time_genelm:>12.1f} {avg_time_genelm - avg_time_no:>12.1f}")
    print(f"{'Time increase':<35} {'-':>12} {'-':>12} {avg_time_pct:>12.1f}%")
    print(f"{'Avg high-conf off-target TIS':<35} {avg_conf_no:>12.1f} {avg_conf_genelm:>12.1f} {avg_conf_reduction:>12.1f}")
    print(f"{'Avg off-target score':<35} {avg_score_no:>12.3f} {avg_score_genelm:>12.3f} {avg_score_genelm - avg_score_no:>12.3f}")

    print(f"\nGeneLM objective effectiveness:")
    print(f"  - Slows optimization by: {avg_time_pct:.1f}%")
    print(f"  - Reduces high-conf off-target TIS by: {avg_conf_reduction:.1f} sites")
    print(f"  - Improves off-target score (less negative = better): {avg_score_improvement:.3f}")

    # Save results to CSV
    import csv
    output_path = Path(__file__).parent / "genelm_objective_results.csv"
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    print(f"\nDetailed results saved to {output_path}")


if __name__ == "__main__":
    main()
