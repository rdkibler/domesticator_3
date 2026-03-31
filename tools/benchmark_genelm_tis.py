#!/usr/bin/env python3
"""
Benchmark GeneLM TIS predictor speed on sequences of various lengths.
"""

import time
import random
import argparse
from pathlib import Path
import sys
import torch
from transformers import AutoModelForSequenceClassification, AutoTokenizer
import orfipy_core as oc
from Bio import Seq, SeqIO

# Add GeneLM to path
SCRIPT_DIR = Path(__file__).resolve().parent
GENELM_RUN_SCRIPTS = SCRIPT_DIR / "GeneLM" / "run-as-script"
sys.path.insert(0, str(GENELM_RUN_SCRIPTS))


def generate_random_dna(length):
    """Generate random DNA sequence."""
    return ''.join(random.choice('ACGT') for _ in range(length))


def get_potential_tis_sites(sequence):
    """Extract all potential TIS sites from sequence using orfipy."""
    seq_rc = str(Seq.Seq(sequence).reverse_complement())

    orfs = oc.start_search(
        sequence,                          # seq
        seq_rc,                            # seq_rc
        "temp",                            # seqname
        10,                                # minlen
        10000000,                          # maxlen
        'f',                               # strand (forward+reverse)
        ['TTG', 'CTG', 'ATG', 'GTG'],      # start codons
        ['TAA', 'TAG', 'TGA'],             # stop codons
        '1',                               # table
        True,                              # include_stop
        False,                             # partial3
        False,                             # partial5
        False,                             # between_stops
        False,                             # nested
        [False, False, True, False, False] # [bed12, bed, dna, rna, pep]
    )

    # Parse ORFs to get TIS sites
    tis_sites = []
    orf_data = orfs[2].strip().split('>')[1:] if orfs[2] else []

    for entry in orf_data:
        lines = entry.splitlines()
        header = lines[0]
        import re
        match = re.search(r'(\d+)-(\d+)\]\((\+|-)\)', header)
        if match:
            start = int(match.group(1))
            end = int(match.group(2))
            orf_strand = match.group(3)
            tis_sites.append((start, end, orf_strand))

    return tis_sites


def generate_kmer(sequence, k=6, overlap=1):
    """Generate k-mer representation of sequence."""
    return " ".join([sequence[j:j+k] for j in range(0, len(sequence) - k + 1, overlap)])


def extract_dna(sequence, start_pos=None, end_pos=None):
    """Extract DNA segment from sequence."""
    if start_pos is not None and end_pos is not None:
        return sequence[start_pos:end_pos]
    return sequence


def run_tis_prediction(model, tokenizer, sequence, tis_sites, window_size=30, device='cpu'):
    """Run TIS prediction on all sites in batch."""
    # Batch prepare all sequences
    batch_sequences = []
    for start, end, strand in tis_sites:
        window_start = max(0, start - window_size)
        window_end = min(len(sequence), start + window_size)
        window_seq = extract_dna(sequence, window_start, window_end)
        kmer_seq = generate_kmer(window_seq, k=6, overlap=1)
        batch_sequences.append(kmer_seq)

    if not batch_sequences:
        return []

    # Tokenize all at once
    inputs = tokenizer(
        batch_sequences,
        return_tensors="pt",
        padding=True,
        truncation=True
    )
    inputs = {key: val.to(device) for key, val in inputs.items()}

    # Batch prediction
    with torch.no_grad():
        outputs = model(**inputs)
        logits = outputs.logits
        probs = torch.softmax(logits, dim=-1).cpu().numpy()

    return probs


def benchmark_sequence_length(length, window_size=30, num_runs=3, device='cpu'):
    """Benchmark TIS prediction for a single sequence length."""
    # Generate random sequence
    sequence = generate_random_dna(length)

    # Get TIS sites
    tis_sites = get_potential_tis_sites(sequence)

    # Initialize model
    model_tis_checkpoint = "Genereux-akotenou/BacteriaTIS-DNABERT-K6-89M"
    model = AutoModelForSequenceClassification.from_pretrained(model_tis_checkpoint).to(device)
    tokenizer = AutoTokenizer.from_pretrained(model_tis_checkpoint)

    # Run multiple times and average
    times = []
    for _ in range(num_runs):
        start_time = time.time()
        probs = run_tis_prediction(model, tokenizer, sequence, tis_sites, window_size, device)
        end_time = time.time()
        times.append(end_time - start_time)

    avg_time = sum(times) / len(times)
    num_sites = len(tis_sites)

    return {
        'length': length,
        'num_tis_sites': num_sites,
        'avg_time': avg_time,
        'sites_per_second': num_sites / avg_time if num_sites > 0 else 0,
        'ms_per_site': (avg_time / num_sites * 1000) if num_sites > 0 else 0
    }


def main():
    parser = argparse.ArgumentParser(description='Benchmark GeneLM TIS predictor speed')
    parser.add_argument('--lengths', type=str, default='1000,5000,10000,20000,50000,100000',
                        help='Comma-separated sequence lengths to benchmark')
    parser.add_argument('--window-size', type=int, default=30, help='Window size around TIS')
    parser.add_argument('--num-runs', type=int, default=3, help='Number of runs per length')
    parser.add_argument('--device', type=str, default='cpu', choices=['cpu', 'gpu'],
                        help='Device to use (cpu or gpu)')
    parser.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility')

    args = parser.parse_args()

    random.seed(args.seed)
    device = torch.device('cuda' if args.device == 'gpu' else 'cpu')

    lengths = [int(x) for x in args.lengths.split(',')]

    print("=" * 80)
    print("GeneLM TIS Predictor Benchmark")
    print(f"Device: {device}")
    print(f"Window size: {args.window_size}")
    print(f"Runs per length: {args.num_runs}")
    print("=" * 80)

    results = []
    for length in lengths:
        result = benchmark_sequence_length(
            length,
            window_size=args.window_size,
            num_runs=args.num_runs,
            device=device
        )
        results.append(result)

        print(f"\nLength: {length:,} bp")
        print(f"  TIS sites found: {result['num_tis_sites']}")
        print(f"  Avg time: {result['avg_time']:.3f} s")
        print(f"  Sites/sec: {result['sites_per_second']:.1f}")
        print(f"  ms/site: {result['ms_per_site']:.2f}")

    # Summary table
    print("\n" + "=" * 80)
    print("Summary")
    print("=" * 80)
    print(f"{'Length':>12} {'TIS sites':>12} {'Time (s)':>12} {'Sites/sec':>12} {'ms/site':>12}")
    print("-" * 80)
    for r in results:
        print(f"{r['length']:>12,} {r['num_tis_sites']:>12} {r['avg_time']:>12.3f} "
              f"{r['sites_per_second']:>12.1f} {r['ms_per_site']:>12.2f}")


if __name__ == "__main__":
    main()
