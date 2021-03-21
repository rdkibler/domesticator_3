import copy
import os
from typing import Generator
import warnings

from Bio import SeqRecord, SeqIO, BiopythonParserWarning
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature

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

	rdkibler 210320
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
			with warnings.catch_warnings(): 
				warnings.filterwarnings("ignore", category=PDBConstructionWarning)
				warnings.filterwarnings("ignore", category=BiopythonParserWarning)
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

	This function looks through the record for features of type "misc_feature"
	and name matching the form "!insert(?)" where "?" is a chain letter. Returns a
	dict of chain letter/ID to the location where it's supposed to go.

	rdkibler 210320
	"""

	location_dict = {}

	for feature in record.features:
		if feature.type == "misc_feature" and feature.qualifiers['label'][0].startswith("!insert("):
			chain_letter = feature.qualifiers['label'][0].split("(")[1].split(")")[0]
			location = feature.location
			location_dict[chain_letter] = feature.location
	return location_dict





def replace_sequence_in_record(record, location, insert) -> SeqRecord.SeqRecord:
	""" Returns a modified seqrecord

	This function was borrowed from Domesticator1. It replaces the sequence in the record 
	and (importantly) adjusts all the locations of the annotations if possible. It returns
	the modified SeqRecord
	"""

	# I don't know if this is right. What does the strand number mean? --rdkibler 210320
	if location.strand >= 0:
		adjusted_seq = record.seq[:location.start] + insert.seq + record.seq[location.end:]
	else:
		adjusted_seq = record.seq[:location.start] + insert.reverse_complement().seq + record.seq[location.end:]

	record.seq = adjusted_seq

	seq_diff = len(insert) - len(location)
	orig_start = location.start
	orig_end = location.end

	processed_features = []

	#add a feature for the insert
	processed_features.append(SeqFeature(location=FeatureLocation(location.start,location.end + seq_diff, strand=location.strand), type="protein", qualifiers={'label':insert.id}))

	for feat in record.features:

		f_loc = feat.location

		loc_list = []

		for subloc in f_loc.parts:

			assert(subloc.start <= subloc.end) 
			
			#type 1: where the start and end are contained within the original location
			#-> do not add it to the processed_features list because I have no idea where it should go
			if subloc.start > location.start and subloc.start < location.end and subloc.end > location.start and subloc.end < location.end:
				continue

			#type 1b: where the start and end are the same which will happen a lot for storing constraints and objectives
			elif subloc.start == location.start and subloc.end == location.end:
				new_loc = FeatureLocation(location.start, location.end + seq_diff, strand=subloc.strand)

			#type 2: where they start or end inside the location
			#I assume that the total length of the annotation is important, so adjust the annotation to have the correct
			# length anchored outside of the insert, unless it'd extend through the insert
			elif subloc.start >= location.start and subloc.start <= location.end:
				#we already caught the case where it's fully within, so I know the end of the subloc extends after the end of the insert 
				new_loc = FeatureLocation(max(subloc.end + seq_diff - len(subloc),location.start), subloc.end + seq_diff, strand=subloc.strand)
				assert len(new_loc) == len(subloc)

			elif subloc.end >= location.start and subloc.end <= location.end:
				new_loc = FeatureLocation(subloc.start, min(subloc.end, location.end + seq_diff), strand=subloc.strand)
				assert len(new_loc) == len(subloc)
				
			#type 3: where they span the location 
			#-> keep the leftmost point same and add diff to rightmost. do not split
			elif location.start >= subloc.start and location.start <= subloc.end and location.end >= subloc.start and location.end <= subloc.end:
				new_loc = FeatureLocation(subloc.start, subloc.end + seq_diff, strand=subloc.strand)

			#type 4: where they start and end before location
			#-> add it to list unchanged
			elif subloc.start <= location.start and subloc.end <= location.start:
				new_loc = subloc

			#type 5: where they start and end after location
			#-> add diff to whole location
			elif subloc.start >= location.end and subloc.end >= location.end:
				new_loc = subloc + seq_diff

			loc_list.append(new_loc)
		
		#if the list is empty, it means that all the sublocs were contained within the insert
		if len(loc_list) > 0:
			feat.location = sum(loc_list)
			processed_features.append(feat)

	record.features = processed_features

	return record


def make_naive_vector_records(base_vector_record, protein_filepaths) -> Generator[SeqRecord.SeqRecord,None,None]:
	""" yields Biopython SeqRecord.SeqRecords

	A generator which loads the proteins, randomly reverse translates each one, inserts them into the vectors.
	and yields 

	rdkibler 210320
	"""

	insert_locations = get_insert_locations(base_vector_record)
	#print(insert_locations)
	for inserts in load_inserts(protein_filepaths):

		if len(insert_locations) == 1:
			for insert in inserts:
				intermediate_vector_record = copy.deepcopy(base_vector_record)
				intermediate_vector_record = replace_sequence_in_record(intermediate_vector_record, insert_locations[insert.annotations['chain']], insert)
				vec_name = intermediate_vector_record.name
				insert_name = insert.name
				intermediate_vector_record.name = f"{insert_name}__{vec_name}"
				yield 
		else:
			intermediate_vector_record = copy.deepcopy(base_vector_record)
			for insert in inserts:
				intermediate_insert_locations = get_insert_locations(intermediate_vector_record)
				intermediate_vector_record = replace_sequence_in_record(intermediate_vector_record, intermediate_insert_locations[insert.annotations['chain']], insert)

			vec_name = intermediate_vector_record.name
			insert_name = insert.name[:-2] #cuts off _A or whatever chain ID it is
			intermediate_vector_record.name = f"{insert_name}__{vec_name}"
			yield intermediate_vector_record


