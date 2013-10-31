import json
import sys
try:
	from Bio.Seq import Seq
	from Bio.Alphabet import IUPAC
except ImportError:
	print ("Need Biopython to use the IgBlast output parser class")

class cdr_analyzer():
	def __init__(self,initial_json_dictionary,query_seq):
		self.unmodified = initial_json_dictionary
		#framework1
		self.frames_and_cdrs = {}
		self.frames_and_cdrs_aa = {}

		#framework1
		self.framework1_frame_info = self.unmodified['alignment_summaries']['fr1_align']
		self.framework1_seq = Seq(query_seq[self.framework1_frame_info['from']-1:self.framework1_frame_info['to']-1], IUPAC.ambiguous_dna)
		self.framework1_seq_aa = self.framework1_seq.translate()

		
		#framework2
		self.framework2_frame_info = self.unmodified['alignment_summaries']['fr2_align']
		self.framework2_seq = Seq(query_seq[self.framework2_frame_info['from']-1:self.framework2_frame_info['to']-1], IUPAC.ambiguous_dna)
		self.framework2_seq_aa = self.framework2_seq.translate()
		
		#framework3
		self.framework3_frame_info = self.unmodified['alignment_summaries']['fr3_align']
		self.framework3_seq = Seq(query_seq[self.framework3_frame_info['from']-1:self.framework3_frame_info['to']-1], IUPAC.ambiguous_dna)
		self.framework3_seq_aa = self.framework3_seq.translate()

		#cdr1
		self.cdr1_frame_info = self.unmodified['alignment_summaries']['cdr1_align']
		self.cdr1_seq = Seq(query_seq[self.cdr1_frame_info['from']-1:self.cdr1_frame_info['to']-1], IUPAC.ambiguous_dna)
		self.cdr1_seq_aa = self.cdr1_seq.translate()

		#cdr2
		self.cdr2_frame_info = self.unmodified['alignment_summaries']['cdr2_align']
		self.cdr2_seq = Seq(query_seq[self.cdr2_frame_info['from']-1:self.cdr2_frame_info['to']-1], IUPAC.ambiguous_dna)
		self.cdr2_seq_aa = self.cdr2_seq.translate()

				
		
		#cdr3
		self.v_part_of_cdr3_info = self.unmodified['alignment_summaries']['cdr3_align']
		self.v_part_of_cdr3 = Seq(query_seq[self.v_part_of_cdr3_info['from']-1:self.v_part_of_cdr3_info['to']-1])
		self.v_d_junction = Seq(self.unmodified['v-d_junction'])
		self.d_region = Seq(self.unmodified['d_region'])
		self.d_j_junction = Seq(self.unmodified['d-j_junction'])
		self.j_region_info = self.unmodified['j_hits'][0]['rank_1']
		_begin = int(self.j_region_info['q._start']) - 1
		_end = int(self.j_region_info['q._end']) -1
		self.j_region = Seq(query_seq[_begin:_end])
		print self.j_region
		sys.exit()

		self.frames_and_cdrs = {'fw1':self.framework1_seq,'cdr1':self.cdr1_seq,'cdr3':self.v_part_of_cdr3,'fw2':self.framework2_seq,'cdr2':self.cdr2_seq,'fw3':self.framework3_seq}
		self.frames_and_cdrs_aa = {'fw1_aa':self.framework1_seq_aa,'cdr1_aa':self.cdr1_seq_aa,'fw2_aa':self.framework2_seq_aa,'cdr2_aa':self.cdr2_seq_aa,'fw3_aa':self.framework3_seq_aa}

		for i,j in zip(sorted(self.frames_and_cdrs),sorted(self.frames_and_cdrs_aa)):
			print i,self.frames_and_cdrs[i]
			print j,self.frames_and_cdrs_aa[j]
		sys.exit()

	def return_json_dict_with_cdr_analysis(self):
		print json.dumps(self.unmodified,sort_keys=1,indent=4)