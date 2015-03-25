#!/usr/bin/env python
import argparse
import os
from Bio import SeqIO
from multiprocessing import cpu_count


class PyIgArgumentParser():

    def __init__(self):
        self.arg_parse = argparse.ArgumentParser(
            prog="PyIg",
            description='''\
                    PyIg - Immunoglobulin and T-Cell receptor rearrangment software. \
                    It uses IgBLAST to call on V(D) and J genes. Then it recombines them to JSON format.\
                    PyIg is meant to be highly parralizable, so it uses multiple processors to call on \
                    multiple instances of BLAST making it ideal for high throughput sequencing. \
                    The JSON documents can be reconfigured to be uploaded to MySQL or NoSQL databases \
                    for efficient queries. In addition, you can use custom databases for V(D) and J gene lookup.
                    \nAuthor - Joran Willis
                     ''')
        necessary_arguments = self.arg_parse.add_argument_group(
            title="Necessary Arguments", description="Arguments that must be included")
        necessary_arguments.add_argument(
            'query', type=self._validate_fasta, metavar="query.fasta",
            help='The fasta file to be run through the protocol')

        type_arguments = self.arg_parse.add_argument_group(
            title="File paths and types", description="Database paths, search types")

        type_arguments.add_argument(
            '-d', '--database', default="/usr/local/pyig/data_dir/",
            type=str, help="If you used setup.py correctly,\
                this database should be /usr/local/pyig/database, but if you made custom changes to it, make sure to put in the path here")

        type_arguments.add_argument(
            '-r', '--receptor', default="Ig", choices=["Ig", "TCR"],
            help="The receptor you are analyzing, immunoglobulin or t cell receptor")

        type_arguments.add_argument(
            '-s', '--species', default='human', choices=['human', 'rabbit', 'mouse', 'rat', 'rhesus'],
            help='The Species you are analyzing')

        type_arguments.add_argument(
            '-c', '--chain', default='heavy', choices=['heavy', 'light', 'mixed'],
            help="The chain you want to analyze")

        blast_arguments = self.arg_parse.add_argument_group(
            title="BLAST Specific Arguments", description="Arguments Specific to IgBlast")

        blast_arguments.add_argument(
            "-nV", "--num_V_alignments", default="1", type=str,
            help="How many V genes do you want to match?")

        blast_arguments.add_argument(
            "-nD", "--num_D_alignments", default="1", type=str,
            help="How many D genes do you want to match?, does not apply for kappa and lambda")

        blast_arguments.add_argument(
            "-nJ", "--num_J_alignments", default="1", type=str,
            help="How many J genes do you want to match?")

        blast_arguments.add_argument(
            '-mD', "--minD", type=self._check_d_match_validity, default="5",
            help="The amount of nucleotide matches needed for a D gene match. >= 5 right now")

        blast_arguments.add_argument(
            "-x", '--executable',
            default="/usr/local/bin/igblastn",
            type=str,
            help="The location of IGBlastn, default is /usr/local/bin/igblastn if used setup.py")

        general_args = self.arg_parse.add_argument_group(
            title="General Arguments", description="Output and Miscellaneous Arguments")
        general_args.add_argument("-m", "--multi", default=cpu_count(), help="Multiprocess by the amount of CPUs you have. \
            Or you can enter a number or type 0 to turn it off")
        general_args.add_argument("-o", "--out", metavar="inputfile.json.gz",
                                  help="Output_file_name, defaults to inputfile.json.gz")
        general_args.add_argument("--debug", default=False, action="store_true",
                                  help="Debug mode, this will not delete the temporary blast files and will print some other useful things, like which regions did not parse")
        general_args.add_argument('--additional_field', type=self._additional_field_parse,
                                  help="A comma key,value pair for an additional field you want to add to the output json. Example \
                                  \n '--additional_field=donor,10` adds a donor field with value 10.")
        general_args.add_argument("-f", "--out-format", default="json", metavar="json",
                                  help="Output file format, defaults to json, but csv can be specified")

    def _check_d_match_validity(self, amount):
        if int(amount) >= 5:
            return str(amount)
        raise argparse.ArgumentTypeError(
            "The amount of D gene nucleotide matches must be >= 5, you have entered {0}".format(str(amount)))

    def _validate_fasta(self, text):
        try:
            SeqIO.parse(text, 'fasta').next()
            return text
        except StopIteration:
            raise argparse.ArgumentTypeError(
                "{0} is not fasta file".format(text))

    def _validate_path(self, path):
        if os.path.exists(os.path.abspath(path)):
            os.environ['IGDATA'] = os.path.abspath(path)
            return os.path.abspath(path)
        raise argparse.ArgumentTypeError(
            "{0} does not exist. Did you use setup.py correctly? Or do you have another location?".format(path))

    def _validate_executable(self, path):
        if os.path.exists(os.path.abspath(path)):
            return os.path.abspath(path)
        raise argparse.ArgumentTypeError(
            "{0} does not exists, please point to where igblastn is".format(path))

    def _additional_field_parse(self, keyvaluestring):
        '''
        I don't know what an exception will be, but have to have some kind of error checking
        '''
        try:
            keyvalue_split = keyvaluestring.split(',')
            return (keyvalue_split[0], keyvalue_split[1])
        except:
            raise argparse.ArgumentTypeError(
                "comma seperated error with {0}".format(keyvaluestring))

    def parse_arguments(self):
        arguments = self.arg_parse.parse_args()
        self._validate_path(arguments.database)
        self._validate_executable(arguments.executable)
        return arguments.__dict__


if __name__ == '__main__':
    print PyIgArgumentParser().parse_arguments()
