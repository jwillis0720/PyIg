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
    if not os.path.exists(os.path.join(_current_directory, os.path.basename(_internal_data))):
        print "Copying internal data to current directory"
        copytree(_internal_data, os.getcwd())

    # set up the command line
    _cline = [manager['executable']]  # we know its in this directory since we copied it here to make this executable
    for argument in blast_options:
        arg = blast_options[argument]
        if arg.startswith("C"):
            arg = '"' + arg + '"'
        current_argument = [argument, arg]
        _cline += current_argument

    print "Running BLAST on processor {0} for split file {1}".format(manager['proc_number'], _file)
    print " ".join(_cline)
    sub = sp.Popen(_cline, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = sub.communicate()

    # if we have output, lets print it
    if stdout or stderr:
        print stdout, stderr

    _output_type = manager['output_type']
    print "Parsing BLAST output to {0} on Processor {1}".format(_output_type, manager['proc_number'])

    if _output_type == "blast_out":
        os.remove(_file)
        print "Removing {0}".format(_file)
    else:
        print _blast_out, _file, _temporary_path
        op = output_parser.igblast_output(_blast_out, _file, _temporary_path,
                                          _output_options, gui=True, zip_bool=_zip_bool)
        op.parse_blast_file_to_type(_json_out, _output_type)
        print "Done parsing {0} type\nRemoving {1} and {2}".format(_output_type, _file, _blast_out)
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


def execute(blast_options, outputoptions):
    '''A function that takes in and executes options from the gui widgets'''
    # variables
    mp.freeze_support()
    ts = time.time()
    fomatted_time = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    print "Process Started {0}".format(fomatted_time)
    processors = outputoptions['num_procs']
    pool = mp.Pool(processes=processors)
    file_name = outputoptions['pre_split_up_input']
    path = outputoptions['tmp_data_directory']
    if not os.path.exists(path):
        msg = "{0} is not found, creating directory...".format(path)
        os.makedirs(os.path.abspath(path))
        print msg

    # split fasta file up
    all_fasta = split_fasta(processors, path, file_name, suffix=".tmp_fasta")
    glob_path = os.path.join(path, os.path.basename(file_name).split('.fasta')[0] + "*.tmp_fasta")

    print "Splitting up file {0} into {1}".format(file_name, path)
    split_up_starting_files = glob.glob(glob_path)

    # output options
    zip_bool = outputoptions['zip_bool']
    output_file = outputoptions['final_outfile']

    # manager_dict
    _manager_list = []
    _manager_dict = {}
    for i, _file in enumerate(split_up_starting_files, start=1):  # the full file name
        _manager_dict['executable'] = outputoptions['executable']
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
        _manager_dict['proc_number'] = i
        _manager_list.append(_manager_dict)
        _manager_dict = {}

    # run_protocol

    # for i in _manager_list:
     #   run_mp_and_delete(i)

    pool.map(run_mp_and_delete, _manager_list)
    concat(_manager_list[0])
    print "Process is done"
    print "Took {0}".format(time.time() - ts)


if __name__ == '__main__':
    execute()
