import subprocess as sp
import multiprocessing as mp
import glob
import os
import output_parser
import sys
import gzip
from time import time
import shutil
from split_fasta import split_fasta

try:
    import Bio
except ImportError("Trouble Installing BioPython:"):
    print("Can't find BioPython Module in this path. PyIgBlast is dependent on Biopython")
    sys.exit()


def run_mp_and_delete(manager):
    '''main method to run igblast through multiprocessing protocol,
    takes in a list of dictionaires each with a seperate set of arguments'''
    # blast options
    blast_options = manager['blast_options']

    # bools
    _zip_bool = manager['zip_bool']

    # file name outputs, these will all be temp files to be parsed later
    _file = manager['split_file']
    _blast_out = _file + ".blast_out"
    if not _zip_bool:
        _json_out = _file + ".json"
    else:
        _json_out = _file + ".json.gz"

    # set the filename in the instance:
    blast_options['-query'] = _file
    blast_options['-out'] = _blast_out

    # temporary path
    _temporary_path = manager['temporary_path']

    # output options
    _output_options = manager['output_options']

    # check on internal data
    _internal_data = manager['internal_data']
    _current_directory = os.getcwd()
    if not os.path.exists(os.path.join(_current_directory + os.path.basename(_internal_data))):
        shutil.copytree(_internal_data, os.getcwd())

    # set up the command line
    _cline = ["./igblastn"]  # we know its in this directory since we copied it here to make this executable
    for argument in blast_options:
        current_argument = [argument, blast_options[argument]]
        _cline += current_argument

   #_cline = " ".join(_cline)
    # print _cline
    # and call igblast
    sp.call(_cline)

    _output_type = manager['output_type']
    if _output_type == "blast_out":
        os.remove(_file)
    else:
        op = output_parser.igblast_output(_blast_out, _file, _temporary_path,
                                          _output_options, zip_bool=_zip_bool)
        op.parse_blast_file_to_type(_json_out, _output_type)
        os.remove(_file)
        os.remove(_blast_out)


def concat(_manager_dict):
    out_file = _manager_dict['output_file']
    file_type = _manager_dict['output_type']
    zip_bool = _manager_dict['zip_bool']
    file_names = os.path.dirname(_manager_dict['split_file']) \
        + "/" + os.path.basename(_manager_dict['non_split']).split('.')[0]
    if zip_bool and file_type == "json":
        zipped_and_json = glob.glob(file_names + "*.json.gz")
        with gzip.open(out_file + ".json.gz", 'wb') as gf:
            for file in zipped_and_json:
                f_in = gzip.open(file, 'rb')
                gf.writelines(f_in)
                f_in.close()
                os.remove(file)
    elif file_type == "json":
        just_json = glob.glob(file_names + "*.json")
        with open(out_file + ".json", 'w') as gf:
            for file in just_json:
                f_in = open(file, 'r')
                gf.writelines(f_in)
                f_in.close()
                os.remove(file)

    elif file_type == "csv":
        just_csv = glob.glob(file_names + "*.csv")
        with open(out_file + ".csv", 'w') as gf:
            for file in just_csv:
                for line in open(file):
                    gf.write(line)

                for files in just_csv[1:]:
                    f = open(files)
                    f.next()
                    for line in f:
                        gf.write(line)
                    f.close()
                # os.remove(file)


#     elif json_bool and not zip_bool:
#         just_json = glob.glob(file_name + "*.json")
#         with open(json_file + ".json", 'w') as gf:
#             for file in just_json:
#                 f_in = open(file, 'r')
#                 gf.writelines(f_in)
#                 f_in.close()
#                 os.remove(file)

#     elif zip_bool and not json_bool:
#         zipped_only = glob.glob(file_name + "*.blast_out.gz")
#         with gzip.open(out_file + ".blast_out.gz", 'wb') as gf:
#             for file in zipped_only:
#                 f_in = gzip.open(file, 'rb')
#                 gf.writelines(f_in)
#                 f_in.close()
#                 os.remove(file)

#     elif not zip_bool and not json_bool:
#         blast_out = glob.glob(file_name + "*.blast_out")
#         with open(out_file + ".blast_out", 'w') as gf:
#             for file in blast_out:
#                 f_in = open(file, 'r')
#                 gf.writelines(f_in)
#                 f_in.close()
#                 os.remove(file)


def execute(blast_options, outputoptions):
    '''A function that takes in and executes options from the gui widgets'''
    # variables
    processors = outputoptions['num_procs']
    pool = mp.Pool(processes=processors)
    file_name = outputoptions['pre_split_up_input']
    path = outputoptions['tmp_data_directory']

    # split fasta file up
    all_fasta = split_fasta(processors, path, file_name, suffix=".tmp_fasta")
    glob_path = path + os.path.basename(file_name).split('.fasta')[0] + "*.tmp_fasta"
    split_up_starting_files = glob.glob(glob_path)

    # output options
    zip_bool = outputoptions['zip_bool']
    output_file = outputoptions['final_outfile']

    # manager_dict
    _manager_list = []
    _manager_dict = {}
    for _file in split_up_starting_files:  # the full file name
        _manager_dict['non_split'] = file_name
        _manager_dict['split_file'] = _file
        _manager_dict['zip_bool'] = zip_bool
        _manager_dict['all_fasta'] = all_fasta
        _manager_dict['blast_options'] = blast_options  # all the blast options
        _manager_dict['internal_data'] = outputoptions['internal_data_directory']
        _manager_dict['output_type'] = outputoptions['output_type']
        _manager_dict['output_file'] = output_file
        _manager_dict['output_options'] = outputoptions['output_options']
        _manager_dict['temporary_path'] = path
        _manager_list.append(_manager_dict)
        _manager_dict = {}

    # run_protocol
    # for i in _manager_list:
    #     run_mp_and_delete(i)
    pool.map(run_mp_and_delete, _manager_list)
    concat(_manager_list[0])

if __name__ == '__main__':
    execute()
