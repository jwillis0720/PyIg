from multiprocessing import cpu_count
import argparse
import sys
import textwrap
import os
import glob
from Bio import SeqIO


class PyIgArgumentParser():

    def __init__(self):
        self.arg_parse = argparse.ArgumentParser(
            prog="PyIg", formatter_class=argparse.RawTextHelpFormatter,
            description=textwrap.dedent('''\
                                                             PyIg
                     __________________________________________________________________________________________\n
                     PyIg is an easy igblast parser that fixes incorrect calls from  IGBLAST for nucleotides.
                     It gives a very nice comma seperated feedback, can connect to MongoDatabases, and use
                     multiproecessing to split up fasta file.
                     \nAuthor - Joran Willis
                     '''))
        necessary_arguments = self.arg_parse.add_argument_group(
            title="Necessary Arguments", description="Arguments that must be included")
        necessary_arguments.add_argument(
            'query', type=self._validate_fasta, metavar="query.fasta",
            help='The fasta file to be run through the protocol')

        type_arguments = self.arg_parse.add_argument_group(
            title="File paths and types", description="Database paths, search types")
        
        type_arguments.add_argument(
            '-d', '--database', default="/usr/local/pyig/data_dir/",
            type=self._validate_path, help="If you used setup.py correctly, this database should \
            be /usr/local/pyig/database, but if you made custom changes to it, make sure to put in the path here")

        type_arguments.add_argument(
            '-r','--receptor', default="Ig",choices=["Ig","TCR"],
            help="The receptor you are analyzing, immunoglobulin or t cell receptor")

        type_arguments.add_argument(
            '-s','--species', default='human',choices=['human','rabbit','mouse','rat','rhesus'],
            help='The Species you are analyzing')

        type_arguments.add_argument(
            '-c','--chain', default='heavy',choices=['heavy,kappa,lambda'],
            help="The chain you want to analyze")


        blast_arguments = self.arg_parse.add_argument_group(
            title="BLAST Specific Arguments", description="Arguments Specific to IgBlast")

        blast_arguments.add_argument(
            "-nV","--num_V_alignments",default="1",type=str,
            help="How many V genes do you want to match?")

        blast_arguments.add_argument(
            "-nD","--num_D_alignments",default="1",type=str,
            help="How many D genes do you want to match?, does not apply for kappa and lambda")

        blast_arguments.add_argument(
            "-nJ","--num_J_alignments",default="1",type=str,
            help="How many J genes do you want to match?")

        blast_arguments.add_argument(
            '-mD',"--minD",type=self._check_d_match_validity,default=5,
            help="The amount of nucleotide matches needed for a D gene match. >= 5 right now")

    
        blast_arguments.add_argument(
             "-x", '--executable',
             default="/usr/local/bin/igblastn",
             type=self._validate_executable,
             help="The location of IGBlastn, default is /usr/local/bin/igblastn if used setup.py")

        general_args = self.arg_parse.add_argument_group(
            title="General Arguments", description="Output and Miscellaneous Arguments")
        general_args.add_argument("-m","--multi",default=cpu_count(),help="Multiprocess by the amount of CPUs you have. \
            Or you can enter a number or type 0 to turn it off")

    def _check_d_match_validity(self, amount):
        if int(amount) >= 5:
            return str(amount)
        raise argparse.ArgumentTypeError("The amount of D gene nucleotide matches must be >= 5, you have entered {0}".format(str(amount)))


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
        raise argparse.ArgumentTypeError("{0} does not exist. Did you use setup.py correctly? Or do you have another location?".format(path))


    def _validate_executable(self, path):
        if os.path.exists(os.path.abspath(path)):
            return os.path.abspath(path)
        raise argparse.ArgumentTypeError("{0} does not exists, please point to where igblastn is".format(path))

    def parse_arguments(self):
        return self.arg_parse.parse_args().__dict__


if __name__ == '__main__':
   print PyIgArgumentParser().parse_arguments()
