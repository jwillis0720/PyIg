from arg_parse import blastargument_parser
import subprocess as sp
import multiprocessing as mp
import glob
import os
import output_parser
import sys
import gzip
from time import time
from split_fasta import split_fasta
try:
    import Bio
except ImportError("Trouble Installing BioPython:"):
    print ("Can't find BioPython Module in this path. PyIgBlast is dependent on Biopython")
    sys.exit()

class execute():
	'''A class encapsulating all the file handling'''
	def __init__(self):
		
		#set up instance
		self.arg_parser_instance = blastargument_parser()
		self.arg_dict = self.arg_parser_instance.return_parsed_args()
		
		#instance variables needed from argument parsed
		self.processors  = self.arg_parser_instance.return_number_procs() - 1
		_pool = mp.Pool(processes=self.processors)
		self.file_name = self.arg_dict['-query']
		self.file_path = os.path.dirname(self.file_name)
		self.file_base_name = os.path.basename(self.file_name).split('.')[0]
		
		#split fasta file up
		split_fasta(self.processors,self.file_name,suffix="tmp.fasta")
		self.split_up_starting_files = glob.glob(self.file_name+"*"+"tmp.fasta")

		#translation?
		self.translation_bool = self.arg_parser_instance.return_translation_bool()

		#output options
		self.zip_bool = self.arg_parser_instance.return_zip_bool()
		self.concat_bool = self.arg_parser_instance.return_concat_bool()
		self.json_bool = self.arg_parser_instance.return_json_bool()
		self.json_prefix = self.arg_parser_instance.return_json_prefix()
		self.output_prefix = self.arg_parser_instance.return_output_prefix()
		
		#run_protocol
		_pool.map(self._run_mp_and_delete,self.split_up_starting_files)

		print _pool, self.split_up_starting_files
		#concat
		if self.concat_bool:
			self._concat()

	def _concat(self):

		_json_file = self.json_prefix
		_out_file = self.output_prefix

		if self.zip_bool and self.json_bool:
			_zipped_and_json = glob.glob(self.file_name+"*.json.gz")
			with gzip.open(_json_file+".json.gz",'wb') as gf:
				for file in _zipped_and_json:
					f_in = gzip.open(file,'rb')
					gf.writelines(f_in)
					f_in.close()
					os.remove(file)

		elif self.json_bool and not self.zip_bool:
			_just_json = glob.glob(self.file_name+"*.json")
			with gzip.open(_json_file+".json",'w') as gf:
				for file in _just_json:
					f_in = open(file,'r')
					gf.writelines(f_in)
					f_in.close()
					os.remove(file)

		elif self.zip_bool and not self.json_bool:
			_zipped_only = glob.glob(self.file_name+"*.blast_out.gz")
			with gzip.open(_out_file+".blast_out.gz",'wb') as gf:
				for file in _zipped_only:
					f_in = gzip.open(file,'rb')
					gf.writelines(f_in)
					f_in.close()
					os.remove(file)

		elif not self.zip_bool and not self.json_bool:
			_blast_out = glob.glob(self.file_name+"*.blast_out")
			with open(_out_file+".blast_out",'w') as gf:
				for file in _blast_out:
					f_in = open(file,'r')
					gf.writelines(f_in)
					f_in.close()
					os.remove(file)




	def _run_mp_and_delete(self,_file):
		#file name outputs
		_blast_out = _file + ".blast_out"
		_json_out = _file + ".json"
		_zip_out = _file + ".json.gz"

		#set the filename in the instance:
		self.arg_dict['-query'] = _file
		self.arg_dict['-out'] = _blast_out

		#set up the command line
		_cline = self.arg_parser_instance.return_command_line_from_dict(self.arg_dict)

		if self.translation_bool:
			_cline.append("-show_translation")
		
		#and call igblast
		sp.call(_cline)

		#and now to zip up that instances file based on options
		#Case 1: We want json and we want it zipped
		if self.json_bool and self.zip_bool:
			output_parser.igblast_output(_blast_out,_json_out,gz=self.zip_bool)
			#Remove old garbage Case
			if self.concat_bool:
				os.remove(_blast_out)
				os.remove(_file)
		
		#Case 2: We want json, but don't want it zipped, and we want anciallary files:
		elif self.json_bool and not self.concat_bool:
			output_parser.igblast_output(_blast_out,_json_out,gz=self.zip_bool)

		#Case 3: We want everything zipped but not turned into json
		elif self.zip_bool and not self.json_bool:
			f_in = open(_blast_out,'rb')
			f_out = gzip.open(_zip_out,'wb')
			f_out.writelines(f_in)
			f_out.close()
			f_in.close()
			#get rid of garbage
			if concat_bool:
				os.remove(_file)
				os.remove(_blast_out)
		
		#Case 4: Only want blast output, but remove split fasta
		elif not self.zip_bool and not self.json_bool and self.concat_bool:
			os.remove(_file)




if __name__ == '__main__':
	ex = execute()