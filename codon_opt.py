#!/home/rdkibler/.conda/envs/michelangeneo/bin/python
# -*- coding: utf-8 -*-

###### finish imports ######
# Standard library imports
import copy
from mimetypes import init

# Third party imports
import dnachisel
from dnachisel import DnaOptimizationProblem, NoSolutionError
from dnachisel import DEFAULT_SPECIFICATIONS_DICT
from dnachisel import Location
from dnachisel import reverse_translate


from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
import numpy as np

# Local application imports
from tools import input_parsing, idt
from specifications.MinimizeNumKmers import MinimizeNumKmers

DEFAULT_SPECIFICATIONS_DICT["MinimizeNumKmers"] = MinimizeNumKmers


def optimize_single(
    amino_acid_sequence,
    kmers_weight=20,
    cai_weight=1.0,
    hairpins_weight=1.0,
    max_tries=10,
    species="e_coli",
):
    idt.user_info_file, idt.token_file = idt.use_dir("~/idt_credentials")
    idt_user_info = idt.get_user_info(idt.user_info_file)

    # doing this mostly
    # base_vector_record = input_parsing.load_vector_record("/home/rdkibler/projects/dom_dev/domesticator3/vectors/fragment.gb")
    # record = input_parsing.make_naive_vector_record_by_seq(base_vector_record,amino_acid_sequence)
    # initial_problem = DnaOptimizationProblem.from_record(record)
    naive_dna_sequence = reverse_translate(amino_acid_sequence)

    location = Location.from_biopython_location(
        FeatureLocation(0, len(amino_acid_sequence) * 3)
    )

    constraints = []
    objectives = []

    constraints.append(
        dnachisel.builtin_specifications.EnforceTranslation(location=location)
    )

    # A series of things that Genscript gets rid of
    constraints.append(
        dnachisel.builtin_specifications.AvoidPattern("CATATG", location=location)
    )
    constraints.append(
        dnachisel.builtin_specifications.AvoidPattern("CTCGAG", location=location)
    )
    constraints.append(
        dnachisel.builtin_specifications.AvoidPattern("GGAGG", location=location)
    )  # ecoli S-D
    constraints.append(
        dnachisel.builtin_specifications.AvoidPattern("TAAGGAG", location=location)
    )
    constraints.append(
        dnachisel.builtin_specifications.AvoidPattern("GCTGGTGG", location=location)
    )
    constraints.append(
        dnachisel.builtin_specifications.AvoidPattern("ATCTGTT", location=location)
    )
    constraints.append(
        dnachisel.builtin_specifications.AvoidPattern("GGRGGT", location=location)
    )

    constraints.append(
        dnachisel.builtin_specifications.AvoidPattern("AAAAAA", location=location)
    )
    constraints.append(
        dnachisel.builtin_specifications.AvoidPattern("TTTTTT", location=location)
    )
    constraints.append(
        dnachisel.builtin_specifications.AvoidPattern("CCCCCC", location=location)
    )
    constraints.append(
        dnachisel.builtin_specifications.AvoidPattern("GGGGGG", location=location)
    )

    # global
    constraints.append(
        dnachisel.builtin_specifications.EnforceGCContent(
            mini=0.4, maxi=0.65, location=location
        )
    )
    # local
    constraints.append(
        dnachisel.builtin_specifications.EnforceGCContent(
            mini=0.25, maxi=0.8, window=50, location=location
        )
    )
    # terminal
    # oops not implemented

    objectives.append(MinimizeNumKmers(k=8, boost=kmers_weight, location=location))
    objectives.append(
        dnachisel.builtin_specifications.AvoidHairpins(
            boost=hairpins_weight, location=location
        )
    )
    objectives.append(
        dnachisel.builtin_specifications.MaximizeCAI(
            species=species, boost=cai_weight, location=location
        )
    )

    initial_problem = DnaOptimizationProblem(
        naive_dna_sequence, constraints=constraints, objectives=objectives, logger=None
    )

    # initial_problem.constraints.extend(constraints)
    # initial_problem.objectives.extend(objectives)

    idt_threshold = 10.0
    num_iterations = 10
    solutions = []
    idt_scores = []
    solution_found = False
    for _ in range(num_iterations):
        if solution_found:
            break

        try:
            problem = copy.deepcopy(initial_problem)
            problem.resolve_constraints_by_random_mutations()
            problem.optimize()
            problem.resolve_constraints(final_check=True)
            solutions.append(problem)
        except NoSolutionError:

            initial_problem.max_random_iters += 1000
            continue

        i = len(solutions)
        solution = solutions[-1]
        seq = solution.sequence

        idt_token = idt.get_token(idt.token_file, idt_user_info)
        response = idt.query_complexity(seq, idt_token["access_token"])

        score_sum = 0
        for issue in response[0]:
            score_sum += issue["Score"]

        idt_scores.append((score_sum, response[0]))
        if score_sum < idt_threshold:
            solution_found = True
            break

    if len(solutions) == 0:
        raise NoSolutionError(
            f"no solution found for {amino_acid_sequence}", initial_problem
        )

    scores = [solution.objectives_evaluations().scores_sum() for solution in solutions]
    best_idx = np.argmin(scores)

    best_solution = solutions[best_idx]

    return best_solution.sequence, idt_scores[best_idx]


def main():
    import argparse
    from Bio import SeqIO

    parser = argparse.ArgumentParser(
        description="optimize sequences easily while checking against IDT"
    )
    parser.add_argument(
        "fasta_files", type=str, nargs="+", help="sequences to optimize in fasta format"
    )
    parser.add_argument(
        "--outpath",
        type=str,
        default="optimized_sequences.fasta",
        help="path to output optimized sequences",
    )
    args = parser.parse_args()

    for fasta_file in args.fasta_files:
        with open(fasta_file) as f:
            for record in SeqIO.parse(f, "fasta"):
                seq = str(record.seq)
                optimized_seq, idt_scores = optimize_single(seq)
                idt_sum, idt_breakdown = idt_scores
                print(idt_breakdown)

                with open(args.outpath, "a") as out:
                    out.write(f">{record.id}\n{optimized_seq}\n")


if __name__ == "__main__":
    main()
