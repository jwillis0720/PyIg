#!/usr/bin/env python
import subprocess as sp
import multiprocessing as mp
import glob
import os
import gzip
import datetime
from shutil import copytree

#Non Standard Library
from pyig.backend import split_fasta
from pyig.backend import output_parser
from pyig.commandline import arg_parse

def run_mp_and_delete(manager):
    '''main method to run igblast through multiprocessing protocol,
    takes in a list of dictionaires each with a seperate set of arguments'''

    # bools
    _zip_bool = manager['zip_bool']
    _json_bool = manager['json_bool']
    _concat_bool = manager['concat_bool']
    _output_options = manager['output_options']

    if _json_bool:
        _output_type = "json"
    else:
        _output_type = "csv"

    # file name outputs, these will all be temp files to be parsed later
    _file = manager['split_file']
    _blast_out = _file.split('.')[0] + ".blast_out"
    _output_file = _file.split('.')[0] + "." + _output_type

    # temporary path
    _temporary_path = manager['tmp_path']

    # check on internal data
    _internal_data = manager['internal_data']
    _current_directory = os.getcwd()

    _species = manager['species']

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

    # Now parse the output
    print "Parsing BLAST output to {0} on Processor {1}".format(_output_type, manager['proc_number'])
    op = output_parser.igblast_output(_blast_out, _file,
                                      _temporary_path, _output_options, species=_species, gui=False, zip_bool=_zip_bool, germ_properties=manager['germ_properties'])
    op.parse_blast_file_to_type(_output_file, _output_type)
    print "Done parsing {0} type".format(_output_type)
    if _concat_bool:
        print "Removing {0} and {1}".format(_file, _blast_out)
        os.remove(_file)
        os.remove(_blast_out)


def concat(_manager_dict):
    out_file = _manager_dict['output_prefix']
    zip_bool = _manager_dict['zip_bool']
    json_bool = _manager_dict['json_bool']
    concat_bool = _manager_dict['concat_bool']

    # join the tmp file path with the query name to get all files that should be concatenated
    file_names = os.path.join(_manager_dict['tmp_path'],
                              os.path.basename(_manager_dict['non_split']).split('.')[0])

    if zip_bool and json_bool:
        zipped_and_json = glob.glob(file_names + "*.json.gz")
        with gzip.open(out_file + ".json.gz", 'wb') as gf:
            for file in zipped_and_json:
                f_in = gzip.open(file, 'rb')
                gf.writelines(f_in)
                f_in.close()
                if concat_bool:
                    os.remove(file)

    elif json_bool and not zip_bool:
        just_json = glob.glob(file_names + "*.json")
        with open(out_file + ".json", 'w') as gf:
            for file in just_json:
                f_in = open(file, 'r')
                gf.writelines(f_in)
                f_in.close()
                if concat_bool:
                    os.remove(file)

    elif zip_bool:
        csv_zip = glob.glob(file_names + "*.csv.gz")
        with gzip.open(out_file + ".csv.gz", 'wb') as gf:
            for line in gzip.open(csv_zip[0], 'rb'):
                    gf.write(line)
            for files in csv_zip[1:]:
                f = gzip.open(files, 'rb')
                f.next()
                for line in f:
                    gf.write(line)
                f.close()
        if concat_bool:
            for file in csv_zip:
                os.remove(file)

    else:
        just_csv = glob.glob(file_names + "*.csv")
        with open(out_file + ".csv", 'w') as gf:
            for line in open(just_csv[0]):
                    gf.write(line)
            for files in just_csv[1:]:
                f = open(files)
                f.next()
                for line in f:
                    gf.write(line)
                f.close()
        if concat_bool:
            for file in just_csv:
                os.remove(file)


def execute(argument_class):
    '''A function that takes in and executes options from the gui widgets'''
    # variables
    ts = datetime.time()
    #fomatted_time = datetime.datetime.fromtimestamp(float(ts)).strftime('%Y-%m-%d %H:%M:%S')
    print "Process Started {0}".format(ts)

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
    output_options = argument_class.get_output_options()

    # blast options - all the options specific to blast executable
    executable = argument_class.get_command()
    blast_options = argument_class.get_blast_options()

    # species
    species = argument_class.get_organism()


    #germ_properties_files
    germ_properties = argument_class.get_germ_file()

    # split fasta file up
    print "Splitting up file {0} into {1}".format(os.path.abspath(query_name), tmp_path)
    split_fasta.split_fasta(processors, tmp_path, query_name, suffix=".tmp_fasta")
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
        _manager_dict['output_options'] = output_options
        _manager_dict['species'] = species
        _manager_dict['proc_number'] = i
        _manager_dict['germ_properties'] = germ_properties
        _manager_list.append(_manager_dict)
        _manager_dict = {}

    # run_protocol
    for i in _manager_list:
        run_mp_and_delete(i)

    #pool = mp.Pool(processes=processors)
    #pool.map(run_mp_and_delete, _manager_list)
    concat(_manager_list[0])
    print "Process is done"
    print "Took {0}".format(datetime.time() - ts)
    os.removedirs(_manager_dict['tmp_path'])

def main():
    execute(arg_parse.argument_parser())

if __name__ == '__main__':
    main()
