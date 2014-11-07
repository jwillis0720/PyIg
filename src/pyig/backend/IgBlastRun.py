import os
import subprocess
import tempfile
from pyig.backend.IgBlastOut import IgBlastOut


class IgBlastRun():

    '''
    IgBlast single run is the class to call for a single IgBlast subprocess.
    This class is the most handy for multiprocessing but can be called alone
    for a single fasta file. It takes in two dictionaries. One dictionary is
    the arguments from the command line. The other is a dictionary of sequences
    with their respective ID.

    Examples:

    --Constructor:

    Ig_sr = IgBlast_SingleRun(argument_dictionary, sequence_dictionary)

    --Set the fasta file to be parsed. Not constant like the argument dictionary

    Ig_sr.set_query(fasta_file)

    --Run the process. Takes in single queue object. This is where the process will
    dump the output file

    Ig_sr.run_single_process(QueueObject)
    '''

    def __init__(self, arg_dict, sequence_dict, query):
        '''
        Constructor takes arguemnt dictionary and sequence dictionary
        '''

        # Set sequence dictionary to member
        self.seqs = sequence_dict

        # Debug mode - Bool
        self.debug = arg_dict['debug']

        # Fetch Argument Dictionary Values
        # IgBlast Specific
        self.executable = arg_dict['executable']
        self.minD = arg_dict['minD']
        self.numV = arg_dict['num_V_alignments']
        self.numD = arg_dict['num_D_alignments']
        self.numJ = arg_dict['num_J_alignments']
        self.species = arg_dict['species']
        self.receptor = arg_dict['receptor']
        self.temporary_output_file = tempfile.NamedTemporaryFile(
            suffix='.blast_out', delete=False).name

        # Database specific
        self.chain = arg_dict['chain']

        # First fetch the path to our data directory
        _path_to_data_base = os.path.join(
            os.environ['IGDATA'], self.receptor, self.chain, self.species)

        self.germline_v = os.path.join(_path_to_data_base, self.species + "_gl_V")
        self.germline_d = os.path.join(_path_to_data_base, self.species + "_gl_D")
        self.germline_j = os.path.join(_path_to_data_base, self.species + "_gl_J")
        self.auxilary_path = os.path.join(os.environ['IGDATA'], "aux", self.species + "_gl.aux")

        # Igblast specific but user can't change without hardcoding
        self.outfmt = "7"
        self.domain_system = "imgt"
        self.additional_info = arg_dict['additional_field']

        if self.debug:
            print "Setting query for {0}".format(query)

        self.query = query

    def _collect(self):
        '''Collect all blast arguments and put them in one list accessible by subproces.Popen'''
        arguments = [
            self.executable,
            '-min_D_match', self.minD,
            '-num_alignments_V', self.numV,
            '-num_alignments_D', self.numD,
            '-num_alignments_J', self.numJ,
            '-organism', self.species,
            '-ig_seqtype', self.receptor,
            '-germline_db_V', self.germline_v,
            '-germline_db_D', self.germline_d,
            '-germline_db_J', self.germline_j,
            '-auxiliary_data', self.auxilary_path,
            '-outfmt', self.outfmt,
            '-domain_system', self.domain_system,
            '-out', self.temporary_output_file,
            '-query', self.query]
        return arguments

    def run_single_process(self, queue):
        '''Call this method with a queue object to dump output file names too'''

        if self.debug:
            print "Running Process for {0}".format(self.query)
            for arg in self._collect():
                print arg,
            print "\nOutput File for igblastn is {0}".format(
                self.temporary_output_file)

        p = subprocess.Popen(self._collect(), stderr=subprocess.PIPE)
        stderr = p.communicate()
        if stderr[1]:
            raise RuntimeError("Error in calling Igblastn:\n\n {0}".format(stderr[1]))

        # until process is done until it moves on to the next line
        p.wait()

        # The IgBlast output class, set the blast output
        IgO = IgBlastOut(debug=self.debug)

        # Set the name of the blast output file
        # Put it as setters, but could put it in the constructor
        IgO.set_seq_dictionary(self.seqs)
        IgO.set_blast_output(self.temporary_output_file)
        IgO.set_species(self.species)
        if self.additional_info:
            IgO.set_additional_info(self.additional_info)

        print "Parsing IgBlast Output to human readable format.."
        # Where the magic happens. This class gets all the calls and sets the output dictionary we will use
        IgO.parse()

        # If we are debugging, we can keep temporary blast output and input
        if not self.debug:
            # remove blast
            os.remove(self.temporary_output_file)
            # remove input fasta
            os.remove(self.query)

        # The Queue Object can be passed around processors.
        # Here we put the name of the output json from Ig blast output parser
        queue.put(IgO.get_output_name())
