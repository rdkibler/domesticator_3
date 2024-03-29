import sys

from Bio import SeqRecord, SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd


def iterate_circular_sequence(start_index: int, sequence, step: int = 3) -> str:
    """Returns a string

    A generator function which returns groups of characters in a string that
    wrap around at the end

    rdkibler 210320
    """

    # technically speaking there are two ways to do this. You could simply concat
    # the sequences and iterate the number of times equal to the length of the
    # starting sequence. That's memory intensive, but could be fast. Or you could
    # just keep track of indices and return slices. This is probably slower, but
    # less memory hungery. I'll use option 2 for now. This is better, too, if you
    # don't really expect to step through the entire sequence such as what we'll
    # use it for

    seq_len = len(sequence)
    for index in range(start_index, 2 * seq_len, step):
        yield sequence[index % seq_len : (index + step) % seq_len]


def find_polypeptides(vector_record) -> list:
    """Returns a list of Bio.SeqRecord.SeqRecords containing protein sequences

    Scans through the record looking for features of type "domesticator" and
    name "start". Scans though the sequence until it encounters a stop codon
    (one of TAA, TAG, or TGA), translates the strech between start and stop, and
    remembers it. Returns all products of this process

    rdkibler 210320
    """
    found_polypeptides = []

    for feature in vector_record.features:

        if feature.type == "domesticator" and feature.qualifiers["label"] == ["start"]:
            translation_start = feature.location.start
            cds_seq = ""
            for codon in iterate_circular_sequence(
                feature.location.start, vector_record.seq
            ):
                if codon in ["TAA", "TAG", "TGA"]:
                    break
                else:
                    # print(codon)
                    cds_seq += codon
            pp = SeqRecord.SeqRecord(cds_seq).translate()

            # look for all protein features between start and stop and add them to the name
            translation_length = len(cds_seq)
            translation_end = translation_start + translation_length
            name = ""
            features_for_naming = []

            for feature2 in vector_record.features:
                if (
                    feature2.type == "protein"
                    and feature2.location.start >= translation_start
                    and feature2.location.start < translation_end
                    and feature2.location.end > translation_start
                    and feature2.location.end <= translation_end + 3
                ):  # the +3 helps it detect stop codons
                    features_for_naming.append(feature2)

            features_for_naming.sort(key=lambda x: x.location.start)

            name = "-".join(
                [feat.qualifiers["label"][0] for feat in features_for_naming]
            )
            pp.name = name

            found_polypeptides.append(pp)
    return found_polypeptides


def get_params(records) -> pd.DataFrame:
    """Returns a DataFrame

    Retrieves data about each protein record in the list, such as molecular weight
    and extinction coefficient, and returns them as a pandas DataFrame

    rdkibler 210320
    """
    param_dict = {}

    for record in records:
        seq = str(record.seq)
        if "*" in seq:
            seq = seq.split("*")[0]
            print(
                f"Warning: The sequence for {record.id} contains a '*'. Everything after the '*' will be ignored",
                file=sys.stderr,
            )
        params = ProteinAnalysis(seq)
        mono_params = ProteinAnalysis(seq, True)
        reduced_ec, oxidized_ec = params.molar_extinction_coefficient()
        reduced_mec = reduced_ec / params.molecular_weight()
        oxidized_mec = oxidized_ec / params.molecular_weight()
        param_dict[record.name] = [
            seq,
            params.molecular_weight(),
            mono_params.molecular_weight(),
            params.isoelectric_point(),
            reduced_ec,
            reduced_mec,
            oxidized_ec,
            oxidized_mec,
        ]
    param_df = pd.DataFrame().from_dict(
        param_dict,
        orient="index",
        columns=[
            "seq",
            "average mass",
            "monoisotopic mass",
            "pI",
            "molar extinction coefficient (reduced)",
            "mass extinction coefficient (reduced)",
            "molar extinction coefficient (oxidized)",
            "mass extinction coefficient (oxidized)",
        ],
    )
    return param_df
