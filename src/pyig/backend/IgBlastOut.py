import os
import tempfile
from pyig.backend.IgBlastOutSingle import IgBlastOutSingle


class IgBlastOut():

    '''
    This class handles the entire blast out handle that comes from each processor.
    That file is further divided up into individual blast hits that are then translated
    and joined

    Methods -
    IgBlastOut.set_seq_dictionary('dictionary with id value pairs')
    IgBlastOut.set_blast_output('outputfrom_blast.blast')
    IgBlastOut.set_species('human')
    IgBlastOutSingle.parse()

    IgBlastOutSingle.get_output_name() - returns the name of the json file after its been parsed
    '''

    def __init__(self, debug=False, out_format="json"):
        # output file to return
        self.parsed_output = tempfile.NamedTemporaryFile(suffix=".json", delete=False).name

        # set all these to empty so we can catch them if they are not set
        self.blast_output_handle = ""
        self.debug = debug
        self.seq_dictionary = {}
        self.species = ""
        self.additional_info = ""
        self.out_format = out_format

    def set_seq_dictionary(self, seq_dictionary):
        self.seq_dictionary = seq_dictionary

    def set_blast_output(self, blast_handle):
        self.blast_output_handle = blast_handle

    def set_species(self, species):
        self.species = species

    def set_additional_info(self, additional_info):
        self.additional_info = additional_info

    def get_output_name(self):
        # the only return method
        return self.parsed_output

    def parse(self):
        if not self.blast_output_handle or not self.species or not self.seq_dictionary:
            raise RuntimeError("blast output, species, or sequence dictionary is not set...see documentation")
        message = '''The sequecne ID from IgBlast and the fasta file are not finding each other. This is probably due to a 'strange' character in the fasta identifier. It is much easier to make every 'strange' character an underscore'''

        # focus lines is a set of lines that are each divided into each blast entry
        _focus_lines = []

        # get J gene HCDR3 stop postion in a properties file
        j_trans = {}
        j_lines = open(os.path.join(os.environ['IGDATA'], "germ_props", self.species, "properties.txt")).readlines()
        for line in j_lines:
            j_trans[line.split()[0]] = line.split()[1]

        # start parsing output
        with open(self.parsed_output, 'w') as out:
            for line in open(self.blast_output_handle):
                # IGBLASTN is the marker for each new entry
                if "IGBLASTN" in line:
                    # focus lines is intially entry so lets start filling it
                    if _focus_lines:
                        # Parse each individual entry
                        Single_Blast_Entry = IgBlastOutSingle(_focus_lines, j_trans, self.species, debug=self.debug)
                        Single_Blast_Entry.parse()
                        Single_Blast_Entry.get_id()
                        # set the actual nucleotide sequence with the lookup dictionary since its
                        # not in the blast output
                        try:
                            Single_Blast_Entry.set_seq(
                                self.seq_dictionary[Single_Blast_Entry.get_id()])
                        except KeyError:
                            raise KeyError(message)

                        # The method to join everything, could be called in single blast entries constructor
                        Single_Blast_Entry.join_and_translate()

                        if self.additional_info:
                            Single_Blast_Entry.set_additional_info(self.additional_info)

                        # get json or csv entry
                        if (self.out_format == "csv"):
                            out.write(Single_Blast_Entry.get_csv_entry())
                            out.write("\n")
                        else:
                            out.write(Single_Blast_Entry.get_json_entry())
                            out.write("\n")

                        # empty focus lines for next entry
                        _focus_lines = []
                    else:
                        # go to next line
                        continue
                else:
                    # fill all lines in between IGBLASTN
                    _focus_lines.append(line)

            # for last line - that can't be handled by our logic. do the same thing.
            if _focus_lines:
                Single_Blast_Entry = IgBlastOutSingle(_focus_lines, j_trans, self.species, debug=self.debug)
                Single_Blast_Entry.parse()
                Single_Blast_Entry.get_id()
                try:
                    Single_Blast_Entry.set_seq(self.seq_dictionary[Single_Blast_Entry.get_id()])
                except KeyError:
                    raise KeyError(message)

                if self.additional_info:
                    Single_Blast_Entry.set_additional_info(self.additional_info)

                Single_Blast_Entry.join_and_translate()
                if (self.out_format == "csv"):
                    out.write(Single_Blast_Entry.get_csv_entry())
                    out.write("\n")
                else:
                    out.write(Single_Blast_Entry.get_json_entry())
                    out.write("\n")
