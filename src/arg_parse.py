import argparse
import sys
import textwrap
from Bio import SeqIO
import os
import glob

from multiprocessing import cpu_count


class argument_parser():
  # call on all our file type parsers in the sequence_anlysis_method

    def __init__(self):
        self.db_directory = os.path.dirname(os.path.abspath(sys.argv[0]))

        """A customized argument parser that does a LOT of error checking"""
        self.parser = argparse.ArgumentParser(
            prog="igblast", formatter_class=argparse.RawTextHelpFormatter,
            description=textwrap.dedent('''\
                                                            PyIgBlast
                    __________________________________________________________________________________________\n
                    PyIgBlast calls upon BLAST for nucleotides. Uses multiproecessing to split up fasta file.
                    Parses the output to a csv/JSON and translates junctions using internal data structure.
                    \nAuthor - Joran Willis
                    '''))

        # Necessary Arguments
        neces = self.parser.add_argument_group(
            title='Necessary', description="These have to be included")
        # query
        neces.add_argument(
            "-q", "--query", metavar="query.fasta", required=True,
            type=self._check_if_fasta,
            help="The fasta file to be input into igBlast")

        # Database specific
        database_group = self.parser.add_argument_group(
            title="\nDatabase Paths")
        database_group.add_argument(
            "-d", "--db_path",
            type=self._check_if_db_exists,
            default=os.path.join(self.db_directory, "datafiles/database/"),
            help="The database path to the germline repertoire")

        # internal_data path
        database_group.add_argument(
            "-i", "--internal_data",
            default=os.path.join(self.db_directory, "datafiles/internal_data/"),
            type=self._check_if_db_exists,
            help="The database path to internal data repertoire")

        # auxilary group
        database_group.add_argument(
            "-a", "--aux_path",
            type=self._check_if_aux_path_exists,
            default=os.path.join(self.db_directory, "datafiles/optional_file/"),
            help="The auxilariay path that contains the frame origins of the germline genes for each repertoire.\n Helps produce translation and other metrics")

        # IGBlast Specific Options
        igspec = self.parser.add_argument_group(
            title="\nIgBlast Sprecific",
            description="IgBlast Specific Options with a Default")

        igspec.add_argument(
            '-y', '--type', default='Ig', choices=['Ig', 'TCR', 'custom'],
            help="Is this an IG or TCR recombination")

        igspec.add_argument(
            "-or", "--organism", default="human",
            choices=["human", "mouse"],
            help="The organism repeortire to blast against")

        igspec.add_argument(
            "-nV", "--num_v", default=1,
            type=int, help="How many V-genes to match?")

        igspec.add_argument(
            "-nD", "--num_d", default=1,
            type=int, help="How many D-genes to match?")

        igspec.add_argument(
            "-nJ", "--num_j", default=1,
            type=int, help="How many J-genes to match?")

        igspec.add_argument("-dgm", "--d_gene_matches",
                            default=5, type=int,
                            help="How many nuclodtieds in the D-gene must match to call it a hit")

        igspec.add_argument("-s", "--domain", default="imgt", choices=[
                            "imgt", "kabat"],
                            help="Which classification system do you want")

        # General Blast Settings
        general = self.parser.add_argument_group(
            title="\nGeneral Settings")

        general.add_argument(
            "-x", '--executable',
            default="/usr/bin/igblastn",
            type=self._check_if_executable_exists,
            help="The location of the executable, default is /usr/bin/igblastn")

        general.add_argument(
            "-o", "--out",
            help="output file prefix",
            default="igblast_out")

        general.add_argument(
            "-t", "--tmp",
            help="temporary directory to store files in.\nDefaults to ./tmp",
            type=self._check_if_tmp_exists,
            default="tmp")

        general.add_argument(
            "-e", "--e_value", type=str, default="1e-15",
            help="Real value for excpectation value threshold in blast.\nPut in scientific notation")

        general.add_argument(
            "-w", "--word_size", type=int, default=4,
            help="Word size for wordfinder algorithm")

        general.add_argument(
            "-pm", "--penalty_mismatch", type=int,
            default=-4, help="Penalty for nucleotide mismatch")

        general.add_argument(
            "-nP", "--num_procs",
            type=self._validate_cpu_count,
            default=cpu_count(),
            help="How many do you want to split the job across, default is the number of processors")

        formatter = self.parser.add_argument_group(
            title="Outputting Options")

        formatter.add_argument(
            "-op",
            "--output_options",
            type=self._validate_output_options,
            default="datafiles/output_options.txt",
            help="Open this file and comment out options you don't want in your final file.\nThe first column is the name of the option.\nThe second column is used by the parser and should not be changed.")

        formatter.add_argument(
            "-z",
            "--zip",
            default=False,
            action="store_true",
            help="Zip up all output files")

        formatter.add_argument(
            "-c", "--concatenate", default=True, action="store_false",
            help="Turn off automatic concatenation and deletion of temporary files.\nFiles are split up at the beginning to run across multiple processors")

        formatter.add_argument(
            "-j", "--json", action="store_true", default=False,
            help="Use the JSON output option that will format the text driven igblast output to a json document.\nDefaults to CSV")

        # return the arguments
        self.args = self.parser.parse_args()
        self._make_args_dict()

    # helper functions to validate arguments
    def _get_germline(self, vdj):
        if self.args.type.lower() == "ig":
            self.args.type = "Ig"
            path_to_data = os.path.join(self.args.db_path, self.args.type, self.args.organism)
        if self.args.type.lower() == "tcr":
            self.args.type = "TCR"
            path_to_data = os.path.join(self.args.db_path, self.args.type, self.args.organism)
        if self.args.type.lower() == "custom":
            path_to_data = os.path.join(self.args.db_path, self.args.type)
        files_in = glob.glob(os.path.join(path_to_data, "*"))
        common_prefix = os.path.basename(os.path.commonprefix(files_in))
        return os.path.join(self.args.db_path, self.args.type, self.args.organism, common_prefix + vdj)

    def _check_if_fasta(self, f_file):
        try:
            SeqIO.parse(f_file, "fasta").next()
            return f_file
        except StopIteration:
            msg = "{0} is not a fasta file\n".format(f_file)
            raise argparse.ArgumentTypeError(msg)

    def _validate_cpu_count(self, cpus):
        if int(cpus) > cpu_count():
            msg = "You have requested more processors than you have available\n\
            Currently have {0} available".format(cpu_count)
            raise argparse.ArgumentTypeError(msg)
        else:
            return int(cpus)

    def _check_if_executable_exists(self, x_path):
        if not os.path.exists(x_path):
            msg = "path to executable {0} does not exist, use -h for help\n".format(x_path)
            raise argparse.ArgumentTypeError(msg)
        if not os.access(x_path, os.R_OK):
            msg1 = "executable {0} does have not permission to run\n".format(x_path)
            raise argparse.ArgumentTypeError(msg1)
        else:
            return x_path

    def _check_if_db_exists(self, db_path):
        if os.path.exists(db_path):
            return db_path
        else:
            msg = "{0} path for does not exist for database\n".format(db_path)
            raise argparse.ArgumentTypeError(msg)

    def _check_if_tmp_exists(self, db_path):
        if os.path.exists(db_path):
            return db_path
        else:
            msg = "{0} path does not exist...making".format(db_path)
            print msg
            os.makedirs(db_path)
            return db_path

    def _check_if_aux_path_exists(self, aux_path):
        if os.path.exists(aux_path):
            return aux_path
        else:
            msg = "{0} path for aux files does not exist\n".format(aux_path)
            raise argparse.ArgumentTypeError(msg)

    def _validate_output_options(self, options_file):
        if os.path.exists(options_file):
            return os.path.abspath(options_file)
        else:
            msg = "{0} does not exists.\nThis file is necessary for the parser to run.".format(options_file)
            raise argparse.ArgumentTypeError(msg)

    def _make_args_dict(self):
        self.args_dict = {
            '-organism': self.args.organism,
            '-num_alignments_V': self.args.num_v,
            '-num_alignments_D': self.args.num_d,
            '-num_alignments_J': self.args.num_j,
            '-min_D_match': self.args.d_gene_matches,
            '-domain_system': self.args.domain,
            '-evalue': self.args.e_value,
            '-word_size': self.args.word_size,
            '-germline_db_V': self._get_germline("V"),
            '-germline_db_D': self._get_germline("D"),
            '-germline_db_J': self._get_germline("J"),
            '-ig_seqtype': self.args.type,
            '-D_penalty': self.args.penalty_mismatch,
            '-auxiliary_data': "{0}{1}_gl.aux".format(
                self.args.aux_path, self.args.organism)

        }

        # add special formatting instructions
        self.args.format_options = self._check_if_db_exists("datafiles/format_template.txt")
        formatting_titles = []
        for line in open(self.args.format_options).readlines():
            if line.startswith("#"):
                continue
            else:
                formatting_titles.append(line.split()[0])
        format = "7 " + " ".join(formatting_titles)
        self.args_dict['-outfmt'] = format

    def get_command(self):
        return self.args.executable

    def get_query(self):
        return self.args.query

    def get_procs(self):
        return self.args.num_procs

    def get_tmp_dir(self):
        return self.args.tmp

    def get_internal_directory(self):
        return self.args.internal_data

    def get_zip_bool(self):
        return self.args.zip

    def get_concat_bool(self):
        return self.args.concatenate

    def get_json_bool(self):
        return self.args.json

    def get_output_prefix(self):
        return self.args.out

    def get_blast_options(self):
        return self.args_dict

    def get_output_options(self):
        return self.args.output_options

    def get_organism(self):
        return self.args.organism
