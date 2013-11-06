import os, os.path, sys
import Tkinter
import ttk
import tkFileDialog as filedialog
import tkMessageBox
from Tkconstants import *
from Bio import SeqIO as so
from vertical_scroll import VerticalScrolledFrame as vsf

class pyigblast_gui():
	def __init__(self,root):
		#Initialization
		self.root = root 
		_program_name = sys.argv[0]
		_directory_name = os.path.dirname(os.path.abspath(_program_name))

		#argument dictionary we will pass to the arg parser eventually
		self.argument_dict = {
								'query':'',
								'database':	_directory_name+"/directory/",
								'in_data': _directory_name+"/internal_data/",
								'aux_data': _directory_name+"/optional_data/"}

		window_info = self.root.winfo_toplevel()
		window_info.wm_title('PyIgBLAST - GUI')
		window_info.geometry('1500x900+10+10')
		#creates main menu
		self.MainMenu()
		#creates main menu notebook inside
		self.TabNotebook()
	

	def MainMenu(self):
		main_menu = ttk.Frame(self.root)
		author_label = ttk.Label(main_menu,text="Jordan Willis")
		university_label = ttk.Label(main_menu,text="Vanderbilt University")
		exit_button = ttk.Button(main_menu,text="EXIT",command=lambda root=self.root:quit(root))
		author_label.pack(side=LEFT,padx=10,pady=3)
		exit_button.pack(side=LEFT,expand=YES,pady=3)
		university_label.pack(side=RIGHT,padx=10,pady=3)
		main_menu.pack(side=BOTTOM, fill=X)

	
	def TabNotebook(self):
		main_notebook_frame = ttk.Notebook(self.root, name='main_notebook')
		main_notebook_frame.enable_traversal()
		main_notebook_frame.pack(side=TOP,expand=1,fill=BOTH)
		self._create_files_and_directories(main_notebook_frame)
		self._create_readme(main_notebook_frame)


	def _create_files_and_directories(self,notebook_frame):
		#first tab
		f_and_d_frame = ttk.Frame(notebook_frame,name='f_and_d')
		fasta_input_frame = ttk.LabelFrame(f_and_d_frame)
		fasta_input_frame.pack(side=TOP,expand=0,fill=X,padx=10)
		self._make_fasta_entry(fasta_input_frame)
		
		#Set up Directory Frame and Labels
		directories_frame = ttk.LabelFrame(f_and_d_frame)
		directories_frame.pack(side=LEFT,expand=1,fill=BOTH)
		directory_label = ttk.Label(directories_frame,font=('Arial',20),text="Directories needed to run blast:")
		directory_label.pack(side=TOP,fill=X,padx=20,pady=10)
		self._set_up_directories(directories_frame)

		#place holder
		filler_frame = ttk.LabelFrame(f_and_d_frame,width=200)
		filler_frame.pack(side=LEFT,expand=1,fill=BOTH)

		#and add it to the big frame
		notebook_frame.add(f_and_d_frame,text="Files and Directories",underline=0,padding=2)

	def _set_up_directories(self,directories_frame):
		#blast directory
		blast_directories_frame = ttk.LabelFrame(directories_frame)
		blast_directories_frame.pack(side=TOP,fill=X,padx=10,expand=1)
		blast_directories_label = ttk.Label(blast_directories_frame,text="Compiled Blast Directory:",font=('Arial',16))
		blast_directories_label.pack(side=TOP,anchor=NW)
		blast_directory_entry = ttk.Entry(blast_directories_frame,width=10)
		blast_directory_entry.insert(END,self.argument_dict['database'])
		blast_directory_entry.pack(side=LEFT,expand=1,fill=X,pady=3)
		blast_directory_button = ttk.Button(blast_directories_frame,text="Browse...",
			command=lambda type='blast',entry=blast_directory_entry:self._enter_directory(type,entry))
		blast_directory_button.pack(side=LEFT,fill=X)
		
		#internal_blast directory
		internal_blast_directories_frame = ttk.LabelFrame(directories_frame)
		internal_blast_directories_frame.pack(side=TOP,fill=X,padx=10,expand=1)
		internal_blast_directories_label = ttk.Label(internal_blast_directories_frame,text="Internal Blast Directory:",font=('Arial',16))
		internal_blast_directories_label.pack(side=TOP,anchor=NW)
		internal_blast_directory_entry = ttk.Entry(internal_blast_directories_frame,width=10)
		internal_blast_directory_entry.insert(END,self.argument_dict['in_data'])
		internal_blast_directory_entry.pack(side=LEFT,expand=1,fill=X,pady=3)
		internal_blast_directory_button = ttk.Button(internal_blast_directories_frame,text="Browse...",
			command=lambda type='in_data',entry=internal_blast_directory_entry:self._enter_directory(type,entry))
		internal_blast_directory_button.pack(side=LEFT,fill=X)

		#auxilliary_blast directory
		aux_blast_directories_frame = ttk.LabelFrame(directories_frame)
		aux_blast_directories_frame.pack(side=TOP,fill=X,padx=10,expand=1)
		aux_blast_directories_label = ttk.Label(aux_blast_directories_frame,text="Auxillary Blast Directory:",font=('Arial',16))
		aux_blast_directories_label.pack(side=TOP,anchor=NW)
		aux_blast_directory_entry = ttk.Entry(aux_blast_directories_frame,width=10)
		aux_blast_directory_entry.insert(END,self.argument_dict['aux_data'])
		aux_blast_directory_entry.pack(side=LEFT,expand=1,fill=X,pady=3)
		aux_blast_directory_button = ttk.Button(aux_blast_directories_frame,text="Browse...",
			command=lambda type='aux_data',entry=aux_blast_directory_entry:self._enter_directory(type,entry))
		aux_blast_directory_button.pack(side=LEFT,fill=X)

	def _enter_directory(self,type,entry):
		direc = None
		tkMessageBox.showinfo("Consult:","Please Read the README tab before messing with this option")
		opts = {'title':"Select FASTA file to open...",
				'initialdir':entry.get()}
		direc = filedialog.askdirectory(**opts)
		if direc:
			entry.delete(0,END)
			entry.insert(END,fn)
			if type == 'blast':
				self.argument_dict['database'] = str(direc)
			elif type == 'in_data':
				self.argument_dict['in_data'] = str(direc)
			elif type == 'aux_data':
				self.argument_dict['aux_data'] = str(direct)	

	def _make_fasta_entry(self,fasta_input_frame):
		message = ttk.Label(fasta_input_frame,relief=FLAT, width=500, anchor=W,
								text='Enter the entry FASTA file here',font=('Arial',20))
		fasta_entry = ttk.Entry(fasta_input_frame,width=10)
		fasta_entry_button = ttk.Button(fasta_input_frame,text="Browse...",
			command=lambda entry=fasta_entry:self._enter_fasta(entry))
		message.pack(side=TOP,expand=1,fill=BOTH,padx=3,pady=3)
		fasta_entry.pack(side=LEFT,padx=3,expand=1,fill=X,pady=3)
		fasta_entry_button.pack(side=LEFT,fill=X)

	def _enter_fasta(self,entry):
		fn = None
		_not_fasta = True
		opts = {'title':"Select FASTA file to open...",
				'initialfile':entry.get()}
		while _not_fasta:
			fn = filedialog.askopenfilename(**opts)
			try:
				so.parse(str(fn),'fasta').next()
				_not_fasta = False
				if fn:
					entry.delete(0,END)
					entry.insert(END,fn)
					self.argument_dict['query'] = str(fn)
			except StopIteration:
				tkMessageBox.showwarning(
										"Open file",
										"Cannot open {0}, it is not a FASTA file\n".format(fn))

	def _create_readme(self,notebook_frame):
		readme_frame = ttk.Frame(notebook_frame,name="r_frame")
		#scroled_widget = ttk.ScrolledWindow(readme_frame,scrollbar='auto')
		#window = scroled.window
		vertical_scroll_frame = vsf(readme_frame)
		vertical_scroll_frame.pack(side=TOP,expand=1,fill=BOTH,anchor=NW)
		vsf_label = ttk.Label(vertical_scroll_frame.interior,text=open('README_gui.txt').readlines(),anchor=N)
		#scroled_widget.pack(side=TOP, fill=BOTH, expand=1)
		vsf_label.pack(side=TOP,fill=BOTH,expand=1,anchor=NW)
		notebook_frame.add(readme_frame,text="Readme",underline=0,padding=2)



if __name__ == '__main__':
	root = Tkinter.Tk()
	pyigblast_class = pyigblast_gui(root)
	root.mainloop()