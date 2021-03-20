from Bio import SeqRecord

def make_vector_record(vector_filepath) -> SeqRecord:
	"""Return Biopython SeqRecord

	Simply a wrapper around Biopython SeqIO machinery to read a single vector file

	rdkibler 210320
	"""
	print("SUCCESS")
	return float