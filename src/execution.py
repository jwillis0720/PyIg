#!/usr/bin/env python
import sys
import subprocess as sp
import multiprocessing as mp
import glob
import os
import gzip
from split_fasta import split_fasta
import time
import datetime
from distutils.dir_util import copy_tree as copytree
import output_parser
from arg_parse import argument_parser as ap


def run_mp_and_delete(manager):
    '''main method to run igblast through multiprocessing protocol,
    takes in a list of dictionaires each with a seperate set of arguments'''

    # bools
    _zip_bool = manager['zip_bool']
    _json_bool = manager['json_bool']
    _concat_bool = manager['concat_bool']
    _output_file = manager['output_prefix']
    #_output_options = manager['output_options']

    if _json_bool:
        _output_type = "json"
    else:
        _output_type = "csv"

    # file name outputs, these will all be temp files to be parsed later
    _file = manager['split_file']
    _blast_out = _file.split('.')[0] + ".blast_out"

    # temporary path
    _temporary_path = manager['tmp_path']

    # check on internal data
    _internal_data = manager['internal_data']
    _current_directory = os.getcwd()

    if not os.path.exists(os.path.join(_current_directory, os.path.basename(_internal_data))):
        print "Copying internal data to current directory"
        copytree(_internal_data, os.getcwd())

    # set up the command line
    _blast_options = manager['blast_options']

    # add executable
    _cline = [manager['executable']]
    # add all blast options
    for argument in _blast_options:
        arg = _blast_options[argument]
        current_argument = [argument, arg]
        _cline += current_argument

    # change them to strings
    _cline = [str(i) for i in _cline]
    # add query and output to command line
    _cline += ["-query", _file]
    _cline += ["-out", _blast_out]

    # run command line
    print "Running BLAST on processor {0} for split file {1}".format(manager['proc_number'], _file)
    sub = sp.Popen(_cline, stdout=sp.PIPE, stderr=sp.PIPE)

    # If we have stderr and stdout, lets print it
    stderr, stdout = sub.communicate()
    if stderr or stdout:
        print stderr, stdout
    sys.exit()

    # Now parse the output
    print "Parsing BLAST output to {0} on Processor {1}".format(_output_type, manager['proc_number'])
    op = output_parser.igblast_output(_blast_out, _file,
                                      _temporary_path, _output_options, zip_bool=_zip_bool)
    op.parse_blast_file_to_type(_output_file, _output_type)
    print "Done parsing {0} type".format(_output_type)
    if _concat_bool:
        print "Removing {0} and {1}".format(_file, _blast_out)
        os.remove(_file)
        os.remove(_blast_out)


def concat(_manager_dict):
    out_file = _manager_dict['output_file']
    file_type = _manager_dict['output_type']
    zip_bool = _manager_dict['zip_bool']

    file_names = os.path.dirname(_manager_dict['split_file']) \
        + "/" + os.path.basename(_manager_dict['non_split']).split('.')[0]

    marker = ""
    if zip_bool:
        marker = ".gz"
    print "Concatinating {0} files to {1}.{2}{3}".format(file_type, out_file, file_type, marker)

    if zip_bool and file_type == "json":
        zipped_and_json = glob.glob(file_names + "*.json.gz")
        with gzip.open(out_file + ".json.gz", 'wb') as gf:
            for file in zipped_and_json:
                f_in = gzip.open(file, 'rb')
                gf.writelines(f_in)
                f_in.close()
                os.remove(file)
                os.remove(file.split('.json.gz')[0] + '.db')

    elif file_type == "json" and not zip_bool:
        just_json = glob.glob(file_names + "*.json")
        with open(out_file + ".json", 'w') as gf:
            for file in just_json:
                f_in = open(file, 'r')
                gf.writelines(f_in)
                f_in.close()
                os.remove(file)
                os.remove(file.split('.json')[0] + '.db')

    elif zip_bool and file_type == "csv":
        csv_zip = glob.glob(file_names + "*.csv.gz")
        with gzip.open(out_file + ".csv.gz", 'wb') as gf:
            for file in csv_zip:
                for line in gzip.open(file, 'rb'):
                    gf.write(line)
                for files in csv_zip[1:]:
                    f = gzip.open(files, 'rb')
                    f.next()
                    for line in f:
                        gf.write(line)
                    f.close()
        for file in csv_zip:
            os.remove(file)
            os.remove(file.split('.csv.gz')[0] + '.db')

    elif file_type == "csv" and not zip_bool:
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
        for file in just_csv:
            os.remove(file)
            os.remove(file.split('.csv')[0] + '.db')

    elif file_type == "blast_out" and not zip_bool:
        blast_only = glob.glob(file_names + "*.blast_out")
        with open(out_file + ".blast_out", 'w') as gf:
            for file in blast_only:
                for line in open(file):
                    gf.write(line)
                os.remove(file)

    elif zip_bool and file_type == "blast_out":
        blast_only = glob.glob(file_names + "*.blast_out")
        with gzip.open(out_file + ".blast_out.gz", 'wb') as gf:
            for file in blast_only:
                for lines in open(file):
                    gf.write(lines)
                os.remove(file)


def execute(argument_class):
    '''A function that takes in and executes options from the gui widgets'''
    # variables
    ts = time.time()
    fomatted_time = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    print "Process Started {0}".format(fomatted_time)

    # query, path and internal database
    query_name = argument_class.get_query()
    tmp_path = argument_class.get_tmp_dir()
    processors = argument_class.get_procs()
    internal_data = argument_class.get_internal_directory()

    # output options
    zip_bool = argument_class.get_zip_bool()
    json_bool = argument_class.get_json_bool()
    concat_bool = argument_class.get_concat_bool()
    output_prefix = argument_class.get_output_prefix()

    # blast options - all the options specific to blast executable
    executable = argument_class.get_command()
    blast_options = argument_class.get_blast_options()

    # split fasta file up
    print "Splitting up file {0} into {1}".format(os.path.abspath(query_name), tmp_path)
    split_fasta(processors, tmp_path, query_name, suffix=".tmp_fasta")
    glob_path = os.path.join(tmp_path, os.path.basename(query_name).split('.')[0] + "*.tmp_fasta")

    # now grab all the temporary files in the temporary directory
    split_up_starting_files = glob.glob(glob_path)

    # manager_dict and list, holds all our values so we can pass them to varias processors
    _manager_list = []
    _manager_dict = {}
    for i, _file in enumerate(split_up_starting_files, start=1):  # the full file name
        _manager_dict['executable'] = executable
        _manager_dict['non_split'] = query_name
        _manager_dict['split_file'] = _file
        _manager_dict['zip_bool'] = zip_bool
        _manager_dict['json_bool'] = json_bool
        _manager_dict['concat_bool'] = concat_bool
        _manager_dict['output_prefix'] = output_prefix
        _manager_dict['tmp_path'] = tmp_path
        _manager_dict['internal_data'] = internal_data
        _manager_dict['blast_options'] = blast_options
        _manager_dict['proc_number'] = i
        _manager_list.append(_manager_dict)
        _manager_dict = {}

    # run_protocol
    for i in _manager_list:
        run_mp_and_delete(i)

    pool = mp.Pool(processes=processors)
    #pool.map(run_mp_and_delete, _manager_list)
    concat(_manager_list[0])
    print "Process is done"
    print "Took {0}".format(time.time() - ts)


if __name__ == '__main__':
    argument_class = ap()
    execute(argument_class)
