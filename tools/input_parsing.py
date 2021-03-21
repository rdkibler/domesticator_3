import os
from typing import Generator

from Bio import SeqRecord, SeqIO
from Bio.Seq import Seq

from dnachisel import reverse_translate

def load_vector_record(vector_filepath) -> SeqRecord.SeqRecord:
	"""Return Biopython SeqRecord.SeqRecords

	Essentially a wrapper around Biopython SeqIO machinery to read a single vector file

	rdkibler 210320
	"""
	records = list(SeqIO.parse(vector_filepath,'genbank'))

	if len(records) != 1:
		raise RuntimeError("correctly formatted vector files only have one record -- the sequence of the vector")
	record = records[0]

	return record



def load_inserts(filenames) -> Generator[SeqRecord.SeqRecord,None,None]:
	""" yields Biopython SeqRecord.SeqRecords

	A generator function which always returns a list of SeqRecords holding DNA 
	encoding the protein sequences taken from the input files, but how it returns 
	them differes by input type.

	If the filename is a fasta file, then each line is returned as a separate record 
	(list of length 1)

	If the filename is a pdb file, then each chain is extracted and and they are 
	returned together

	"""

	for filename in filenames: 
		assert os.path.isfile(filename)
		basename = os.path.basename(filename)
		name,ext = os.path.splitext(basename)
		assert ext in ['.pdb','.fasta']

		if ext == '.fasta':
			for record in SeqIO.parse(filename, 'fasta'):
				orig_aa_seq = record.seq
				new_dna_seq = Seq(reverse_translate(record.seq))
				assert(orig_aa_seq == new_dna_seq.translate())
				yield [SeqRecord.SeqRecord(seq=new_dna_seq,id=record.id,name=record.name,description=record.description,annotations={"molecule_type": "DNA"})]
		else:
			records = []
			for record in SeqIO.parse(filename, "pdb-atom"):
				orig_aa_seq = record.seq
				new_dna_seq = Seq(reverse_translate(record.seq))
				assert(orig_aa_seq == new_dna_seq.translate())
				new_id = name + "_" + record.annotations['chain']
				new_name = name + "_" + record.annotations['chain']
				records.append(SeqRecord.SeqRecord(seq=new_dna_seq,id=new_id,name=new_name,description=record.description,annotations={"molecule_type": "DNA", "chain":record.annotations['chain']}))
			yield records



def get_insert_locations(record) -> dict:
	""" returns a dict of chain ID to dnachisel(?) locations

	This function looks through the record for features of type "domesticator"
	and name matching the form "!insert(?)" where "?" is a chain letter. Returns a
	dict of chain letter/ID to the location where it's supposed to go.
	"""
	raise NotImplementedError()


def put_insert_into_vector(insert_record,vector_record,location_ID:str) -> SeqRecord.SeqRecord:
	""" returns a SeqRecord

	This function will insert one insert_record into a copy of the vector_record at 
	the specified location (a single letter), making sure to adjust other annotations
	as needed.
	"""
	raise NotImplementedError()


def make_naive_vector_records(base_vector_record, protein_filepaths) -> Generator[SeqRecord.SeqRecord,None,None]:
	""" yields Biopython SeqRecord.SeqRecords

	A generator which loads the proteins, randomly reverse translates each one, inserts them into the vectors.
	and yields 

	"""
	print('hi')
	insert_locations = get_insert_locations(base_vector_record)
	print('how')
	for inserts in load_inserts(protein_filepaths):
		if len(insert_locations) == 1:
			for insert in inserts:
				yield put_insert_into_vector(insert,base_vector_record,"A")
		else:
			intermediate_vector_record = base_vector_record.copy()
			for insert in inserts:
				intermediate_vector_record = put_insert_into_vector(insert, intermediate_vector_record, insert.annotations['chain'])
			yield intermediate_vector_record


