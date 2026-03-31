import torch
from transformers import AutoModelForSequenceClassification, AutoTokenizer
from dnachisel import Specification, SpecEvaluation
from dnachisel.Location import Location
import numpy as np
from pathlib import Path
import sys
import re
from Bio import Seq, SeqIO

# Add GeneLM to path
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(__file__).resolve().parent.parent
GENELM_RUN_SCRIPTS = PROJECT_ROOT / "tools" / "GeneLM" / "run-as-script"
sys.path.insert(0, str(GENELM_RUN_SCRIPTS))

from genelm.core import AnnotatorPipeline
import orfipy_core as oc


class MinimizeOffTargetTISConfidence(Specification):
    """Minimize confidence in TIS predictions outside of intended TIS sites.

    This specification uses the GeneLM TIS predictor to score all potential
    TIS sites in the sequence and penalizes high-confidence predictions that
    fall outside of the intended TIS locations.

    The score is negative, with 0 being perfect (no off-target predictions
    with confidence above threshold).
    """

    best_possible_score = 0

    def __init__(self, intended_tis_sites=None, confidence_threshold=0.9,
                 window_size=30, location=None, boost=1.0):
        """
        Args:
            intended_tis_sites: List of (start, end) positions for intended TIS sites.
                If None, will use the location of the specification on the problem.
            confidence_threshold: Minimum confidence (probability) to count as a
                significant off-target prediction.
            window_size: Window size around each potential TIS site to evaluate.
            location: DNAChisel location specifying where this spec applies.
            boost: Multiplicative factor for the score.
        """
        self.intended_tis_sites = intended_tis_sites or []
        self.confidence_threshold = confidence_threshold
        self.window_size = window_size
        self.location = location
        self.boost = boost

        # Model is lazy-loaded when needed
        self._model = None
        self._tokenizer = None

    def _get_model(self):
        """Lazy load GeneLM TIS model."""
        if self._model is None:
            model_tis_checkpoint = "Genereux-akotenou/BacteriaTIS-DNABERT-K6-89M"
            device = torch.device("cpu")  # Use CPU for optimization
            self._model = AutoModelForSequenceClassification.from_pretrained(
                model_tis_checkpoint
            ).to(device)
            self._tokenizer = AutoTokenizer.from_pretrained(model_tis_checkpoint)
        return self._model, self._tokenizer

    def _generate_kmer(self, sequence, k=6, overlap=1):
        """Generate k-mer representation of sequence."""
        return " ".join([sequence[j:j+k] for j in range(0, len(sequence) - k + 1, overlap)])

    def _extract_dna(self, sequence, start_pos=None, end_pos=None):
        """Extract DNA segment from sequence."""
        if start_pos is not None and end_pos is not None:
            return sequence[start_pos:end_pos]
        return sequence

    def _get_all_potential_tis_sites(self, sequence):
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
            match = re.search(r'(\d+)-(\d+)\]\((\+|-)\)', header)
            if match:
                start = int(match.group(1))
                end = int(match.group(2))
                orf_strand = match.group(3)
                tis_sites.append((start, end, orf_strand))

        return tis_sites

    def _evaluate_tis_confidence_batch(self, sequence, tis_sites, intended_tis):
        """Evaluate TIS confidence using GeneLM model in batch."""
        model, tokenizer = self._get_model()
        device = torch.device("cpu")

        # Batch prepare all sequences
        batch_sequences = []
        for start, end, strand in tis_sites:
            window_start = max(0, start - self.window_size)
            window_end = min(len(sequence), start + self.window_size)
            window_seq = self._extract_dna(sequence, window_start, window_end)
            kmer_seq = self._generate_kmer(window_seq, k=6, overlap=1)
            batch_sequences.append(kmer_seq)

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

        # Build scores list
        scores = []
        for i, (start, end, strand) in enumerate(tis_sites):
            tis_confidence = float(probs[i][1])
            is_intended = any(
                abs(start - int_start) < self.window_size
                for int_start, _ in intended_tis
            )
            scores.append({
                'start': start,
                'end': end,
                'strand': strand,
                'confidence': tis_confidence,
                'is_intended': is_intended
            })

        return scores

    def initialize_on_problem(self, problem, role=None):
        return self._copy_with_full_span_if_no_location(problem)

    def evaluate(self, problem):
        """Evaluate the specification on the problem's sequence.

        Note: This evaluates on the CURRENT sequence (which may differ from
        the original during optimization). The TIS sites are re-extracted
        for each evaluation since the sequence may have mutated.
        """
        sequence = str(problem.sequence)

        # Set intended TIS if not specified
        if not self.intended_tis_sites and self.location:
            center = self.location.start + (self.location.length // 2)
            self.intended_tis_sites = [(center, center + 3)]

        # Re-extract TIS sites for current sequence
        current_tis = self._get_all_potential_tis_sites(sequence)

        # Evaluate confidence for each site
        scores = self._evaluate_tis_confidence_batch(
            sequence, current_tis, self.intended_tis_sites
        )

        # Calculate penalty for off-target high-confidence predictions
        total_penalty = 0.0
        off_target_count = 0

        for site in scores:
            if not site['is_intended'] and site['confidence'] >= self.confidence_threshold:
                excess_confidence = site['confidence'] - self.confidence_threshold
                total_penalty += excess_confidence
                off_target_count += 1

        # Score is negative (lower is worse), normalized by sequence length
        if len(sequence) > 0:
            score = -total_penalty * 100 / len(sequence)
        else:
            score = 0.0

        message = f"Off-target TIS: {off_target_count} sites >= {self.confidence_threshold:.0%}"
        message += f" (penalty: {total_penalty:.3f})"

        return SpecEvaluation(
            self,
            problem,
            score=score,
            locations=[self.location],
            message=message
        )

    def label_parameters(self):
        return [
            ("threshold", str(self.confidence_threshold)),
            ("window", str(self.window_size)),
        ]

    def short_label(self):
        return f"MinOffTIS {self.confidence_threshold}"

    def __str__(self):
        """String representation."""
        return f"MinimizeOffTargetTISConfidence(threshold={self.confidence_threshold})"
