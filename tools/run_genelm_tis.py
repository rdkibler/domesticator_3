#!/usr/bin/env python3
"""
Run GeneLM TIS predictor on a GenBank or FASTA file and print a summary of results.
"""

import argparse
import os
import sys
from pathlib import Path
from Bio import SeqIO

# Add GeneLM to path
SCRIPT_DIR = Path(__file__).resolve().parent
GENELM_RUN_SCRIPTS = SCRIPT_DIR / "GeneLM" / "run-as-script"
sys.path.insert(0, str(GENELM_RUN_SCRIPTS))

from genelm.core import AnnotatorPipeline
import torch
import logging
import uuid


def convert_gb_to_fasta(gb_path, output_path):
    """Convert GenBank file to FASTA format."""
    with open(gb_path, "r") as gb_file:
        records = list(SeqIO.parse(gb_file, "genbank"))

    with open(output_path, "w") as fasta_file:
        for record in records:
            SeqIO.write(record, fasta_file, "fasta")

    return records


def summarize_results(gff_path, input_file_path):
    """Summarize TIS prediction results."""
    print("=" * 60)
    print("GeneLM TIS Prediction Summary")
    print("=" * 60)

    # Read the input file to get annotations if available
    input_format = "genbank" if str(input_file_path).lower().endswith(".gb") else "fasta"
    input_records = {}
    input_records_by_desc = {}  # Also index by description for matching
    if input_format == "genbank":
        with open(input_file_path, "r") as f:
            for record in SeqIO.parse(f, "genbank"):
                # Get all features including CDS and gene
                all_features = [f for f in record.features if f.type in ("CDS", "gene")]
                input_records[record.id] = {
                    "features": all_features,
                    "seq_length": len(record.seq)
                }
                # Also index by description (like "pET29b")
                desc = record.description.split()[0] if record.description else record.id
                input_records_by_desc[desc] = input_records[record.id]
    else:
        with open(input_file_path, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                input_records[record.id] = {
                    "features": [],
                    "seq_length": len(record.seq)
                }
                # Index by short name
                short_name = record.id.split()[0]
                input_records_by_desc[short_name] = input_records[record.id]

    # Read GFF results
    predictions = []
    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("##") or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 9:
                predictions.append({
                    "seq_id": parts[0],
                    "source": parts[1],
                    "type": parts[2],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "score": parts[5],
                    "strand": parts[6],
                    "phase": parts[7],
                    "attributes": parts[8]
                })

    # Summary by sequence
    seq_stats = {}
    for pred in predictions:
        seq_id = pred["seq_id"]
        if seq_id not in seq_stats:
            seq_stats[seq_id] = {"tis_count": 0, "inside_annotation": 0, "outside_annotation": 0}

        seq_stats[seq_id]["tis_count"] += 1

        # Check if inside existing annotation
        if seq_id in input_records:
            features = input_records[seq_id].get("features", [])
            pred_start = pred["start"]
            pred_end = pred["end"]

            for feature in features:
                feat_start = feature.location.start.real
                feat_end = feature.location.end.real

                # Check for overlap
                if (pred_start <= feat_end and pred_end >= feat_start):
                    seq_stats[seq_id]["inside_annotation"] += 1
                    break
            else:
                seq_stats[seq_id]["outside_annotation"] += 1

    # Print results
    print(f"\nInput file: {input_file_path}")
    print(f"Input format: {input_format.upper()}")
    print(f"Total predictions: {len(predictions)}")

    for seq_id, stats in seq_stats.items():
        print(f"\nSequence: {seq_id}")
        if seq_id in input_records:
            print(f"  Sequence length: {input_records[seq_id]['seq_length']} bp")
            print(f"  Existing CDS features: {len(input_records[seq_id]['features'])}")
        print(f"  TIS sites predicted: {stats['tis_count']}")
        print(f"  Inside existing annotations: {stats['inside_annotation']}")
        print(f"  Outside existing annotations: {stats['outside_annotation']}")

    # Detailed predictions
    if predictions:
        print("\n" + "-" * 60)
        print("Detailed Predictions (TIS Refinement Results):")
        print("-" * 60)
        print(f"{'Seq ID':<20} {'Start':<10} {'End':<10} {'Strand':<8} {'Prob':<10} {'Status':<15}")
        print("-" * 60)

        for pred in predictions:
            # Parse attributes for probability
            attrs = pred["attributes"]
            prob = "N/A"
            if "Prob_cls1=" in attrs:
                prob = attrs.split("Prob_cls1=")[1].split(";")[0]

            status = "High-confidence"
            if pred["type"] == "CDS":
                status = "CDS"

            print(f"{pred['seq_id']:<20} {pred['start']:<10} {pred['end']:<10} "
                  f"{pred['strand']:<8} {prob:<10} {status:<15}")

    print("\n" + "=" * 60)
    print("Summary complete.")
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="Run GeneLM TIS predictor and summarize results"
    )
    parser.add_argument("input_file", help="Input GenBank (.gb) or FASTA file")
    parser.add_argument("--device", default="cpu", choices=["cpu", "gpu"],
                        help="Device to run on (default: cpu)")
    parser.add_argument("--output_dir", default=None,
                        help="Output directory for GFF file (default: GeneLM temp)")
    parser.add_argument("--keep_temp", action="store_true",
                        help="Keep temporary files")

    args = parser.parse_args()

    input_path = Path(args.input_file)
    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}")
        sys.exit(1)

    # Convert GenBank to FASTA if needed
    if str(input_path).lower().endswith(".gb"):
        print(f"Converting {input_path.name} to FASTA...")
        temp_fasta = input_path.with_suffix(".fasta")
        records = convert_gb_to_fasta(str(input_path), str(temp_fasta))
        input_path = temp_fasta
    else:
        records = list(SeqIO.parse(open(input_path), "fasta"))

    print(f"Processing {len(records)} sequence(s)...")
    print(f"Device: {args.device}")

    # Set device
    if args.device == "cpu":
        os.environ["CUDA_VISIBLE_DEVICES"] = ""

    # Initialize pipeline
    print("Loading GeneLM models (first run may download)...")
    annot = AnnotatorPipeline()

    # Run prediction
    print("Running annotation...")
    task_uuid = str(uuid.uuid4())
    tasks = {task_uuid: {"status": "Queued", "progress": 0, "result": None, "exec_state": {}}}

    result_path = annot.pipeline(
        input_path,
        "GFF",
        tasks,
        task_uuid,
        logging
    )

    if result_path is None:
        print("Error: Annotation failed.")
        sys.exit(1)

    print(f"\nResults saved to: {result_path}")

    # Summarize results
    summarize_results(str(result_path), args.input_file)

    # Cleanup temp FASTA if we converted
    if str(args.input_file).lower().endswith(".gb") and input_path.exists():
        input_path.unlink()


if __name__ == "__main__":
    main()
