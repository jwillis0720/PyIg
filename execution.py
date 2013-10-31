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

arg_parser_instance = blastargument_parser()
arg_dict = arg_parser_instance.return_parsed_args()

def run_mp_and_delete(file_dict):

	#bools
	_json_bool = file_dict['json_bool']
	_concat_bool = file_dict['concat_bool']
	_zip_bool = file_dict['zip_bool']

	all_fasta = file_dict['all_fasta']

	#file name outputs
	_blast_out = file_dict['split_file'] + ".blast_out"
	_json_out = file_dict['split_file'] + ".json"
	_zip_out = file_dict['split_file'] + ".json.gz"

	_file = file_dict['split_file']

	#set the filename in the instance:
	arg_dict['-query'] = _file
	arg_dict['-out'] = _blast_out

	#set up the command line
	_cline = arg_parser_instance.return_command_line_from_dict(arg_dict)
	
	if file_dict['translation_bool']:
		_cline.append("-show_translation")

	#and call igblast
	sp.call(_cline)


	#and now to zip up that instances file based on options
	#Case 1: We want json and we want it zipped
	if _json_bool and _zip_bool:
		output_parser.igblast_output(_blast_out,_json_out,all_fasta,gz=_zip_bool)
		#Remove old garbage Case
		if _concat_bool:
			os.remove(_blast_out)
			os.remove(_file)

	elif _json_bool and _concat_bool:
		output_parser.igblast_output(_blast_out,_json_out,all_fasta,gz=_zip_bool)
		os.remove(_file)
		os.remove(_blast_out)
		

	#Case 3: We want json, but don't want it zipped, and we want anciallary files:
	elif _json_bool and not _concat_bool:
		output_parser.igblast_output(_blast_out,_json_out,all_fasta,gz=_zip_bool)

	#Case 4: We want everything zipped but not turned into json
	elif _zip_bool and not _json_bool:
		f_in = open(_blast_out,'rb')
		f_out = gzip.open(_zip_out,'wb')
		f_out.writelines(f_in)
		f_out.close()
		f_in.close()
		#get rid of garbage
		if concat_bool:
			os.remove(_file)
			os.remove(_blast_out)

	#Case 5: Only want blast output, but remove split fasta
	elif not _zip_bool and not _json_bool and _concat_bool:
		os.remove(_file)

def concat(_manager_dict):
		
		json_file = _manager_dict['json_prefix']
		out_file = _manager_dict['output_prefix']
		zip_bool = _manager_dict['zip_bool']
		json_bool = _manager_dict['json_bool']
		file_name = _manager_dict['file'].split('.fasta')[0]


		if zip_bool and json_bool:
			zipped_and_json = glob.glob(file_name+"*.json.gz")
			with gzip.open(json_file+".json.gz",'wb') as gf:
				for file in zipped_and_json:
					f_in = gzip.open(file,'rb')
					gf.writelines(f_in)
					f_in.close()
					os.remove(file)

		elif json_bool and not zip_bool:
			just_json = glob.glob(file_name+"*.json")
			with open(json_file+".json",'w') as gf:
				for file in just_json:
					f_in = open(file,'r')
					gf.writelines(f_in)
					f_in.close()
					os.remove(file)

		elif zip_bool and not json_bool:
			zipped_only = glob.glob(file_name+"*.blast_out.gz")
			with gzip.open(out_file+".blast_out.gz",'wb') as gf:
				for file in zipped_only:
					f_in = gzip.open(file,'rb')
					gf.writelines(f_in)
					f_in.close()
					os.remove(file)

		elif not zip_bool and not json_bool:
			blast_out = glob.glob(file_name+"*.blast_out")
			with open(out_file+".blast_out",'w') as gf:
				for file in blast_out:
					f_in = open(file,'r')
					gf.writelines(f_in)
					f_in.close()
					os.remove(file)


def execute():
	'''A function encapsulating all the file handling'''	
	#variables
	processors  = arg_parser_instance.return_number_procs() - 1
	pool = mp.Pool(processes=processors)
	file_name = arg_dict['-query']
	path = os.path.dirname(file_name)
	base_name = os.path.basename(file_name).split('.')[0]
	
	#split fasta file up
	all_fasta = split_fasta(processors,file_name,suffix=".tmp_fasta")
	split_up_starting_files=glob.glob(file_name.split('.fasta')[0]+"*"+".tmp_fasta")

	#translation?
	translation_bool = arg_parser_instance.return_translation_bool()
	
	#output options
	zip_bool = arg_parser_instance.return_zip_bool()
	concat_bool = arg_parser_instance.return_concat_bool()
	json_bool = arg_parser_instance.return_json_bool()
	json_prefix = arg_parser_instance.return_json_prefix()
	output_prefix = arg_parser_instance.return_output_prefix()
	
	#manager_dict
	_manager_list = []
	_manager_dict = {}
	for _file in split_up_starting_files:
		_manager_dict['file'] = file_name
		_manager_dict['split_file'] = _file
		_manager_dict['zip_bool'] = zip_bool
		_manager_dict['concat_bool'] = concat_bool
		_manager_dict['json_prefix'] = json_prefix
		_manager_dict['output_prefix'] = output_prefix
		_manager_dict['translation_bool'] = translation_bool
		_manager_dict['json_bool'] = json_bool
		_manager_dict['all_fasta'] = all_fasta
		_manager_list.append(_manager_dict)
		_manager_dict = {}

	#run_protocol
	#for i in _manager_list:
	#	run_mp_and_delete(i)
	pool.map(run_mp_and_delete,_manager_list)

	if concat_bool:
		concat(_manager_list[0])

if __name__ == '__main__':
	execute()