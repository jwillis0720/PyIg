from PyIgArgumentParser import PyIgArgumentParser as PyIg_Parse
import os
import subprocess 
from multiprocessing import Process
from split_fasta import split_fasta
import tempfile

arg_dict = PyIg_Parse().parse_arguments()
class IgBlast_SingleRun():
	'''IgBlast SubProcess for a Single Run -
	Use - IgBlast_SingRun(fasta_file).run()
	Just pass a single fasta file'''

	def __init__(self):
		'''Constructor takes fasta file'''
		self.ig_blast_executable = [arg_dict['executable']]
		self.minD = ["-min_D_match",str(arg_dict['minD'])]
		self.numV = ["-num_alignments_V",arg_dict['num_V_alignments']]
		self.numD = ["-num_alignments_D",arg_dict['num_D_alignments']]
		self.numJ = ["-num_alignments_J",arg_dict['num_J_alignments']]
		self.organism = ["-organism",arg_dict['species']]
		self.receptor = ["-ig_seqtype",arg_dict['receptor']]

		#Non Implemented
		self.chain = ["", arg_dict['chain']]

		#Not shown to user
		self.domain_system = ["-domain_system","imgt"]
		self.translate = ['-show_translation']

		#Add database and auxillary data
		_path_to_data_base = os.path.join(os.environ['IGDATA'],arg_dict['receptor'],arg_dict['species'])
		self.germline_v = ["-germline_db_V",os.path.join(_path_to_data_base,arg_dict['species']+"_gl_V")]
		self.germline_d = ["-germline_db_D",os.path.join(_path_to_data_base,arg_dict['species']+"_gl_D")]
		self.germline_j = ["-germline_db_J",os.path.join(_path_to_data_base,arg_dict['species']+"_gl_J")]
		self.auxilary_path = ['-auxiliary_data',os.path.join(os.environ['IGDATA'],"aux",arg_dict['species']+"_gl.aux")]


		#output_file for blast out
		self.output_file = ['-out', tempfile.NamedTemporaryFile(suffix='.blast_out')]
	
	def set_query(self,file):
		print "Setting query for {0}".format(file)
		self.query = ['-query',file]
		self.run_single_process()
		return self.output_file

	def _collect(self):
		'''Collect all blast arguments and put them in one list accessible by Popen'''
		return self.ig_blast_executable + \
			   self.minD + \
			   self.numV + \
			   self.numD + \
			   self.numJ + \
			   self.organism + \
			   self.receptor + \
			   self.domain_system + \
			   self.translate + \
			   self.germline_v + \
			   self.germline_d + \
			   self.germline_j + \
			   self.auxilary_path +\
			   self.query


	def run_single_process(self):
		print "Running Process for {0}".format(self.query[1])
		p = subprocess.Popen(self._collect())
		p.wait()
		
def run_all_process(fasta_file):
	num_procs = int(arg_dict['multi'])
	IgBlast = IgBlast_SingleRun()
	if num_procs > 1:
		jobs = []
		split_fasta_file_names = split_fasta(num_procs,fasta_file,delete=False)
		for name in split_fasta_file_names:
			p = Process(target=IgBlast.set_query,args=(name,))
			jobs.append(p)			
			p.start()

		for job,temp_file in zip(jobs,split_fasta_file_names):
			job.join()
			os.remove(temp_file)
	else:
		p = Process(target=IgBlast.set_query,args=(fasta_file,))
		p.start()

if __name__ == '__main__':
	run_all_process(arg_dict['query'])