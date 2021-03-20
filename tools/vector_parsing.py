from Bio import SeqRecord, SeqIO

def make_vector_record(vector_filepath) -> SeqRecord.SeqRecord:
	"""Return Biopython SeqRecord

	Essentially a wrapper around Biopython SeqIO machinery to read a single vector file

	rdkibler 210320
	"""
	records = list(SeqIO.parse(vector_filepath,'genbank'))

	if len(records) != 1:
		raise RuntimeError("correctly formatted vector files only have one record -- the sequence of the vector")
	record = records[0]

	return record