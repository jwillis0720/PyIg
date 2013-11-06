import os, os.path, sys
import Tkinter
import ttk
import tkFileDialog as filedialog
from Tkconstants import *

class pyigblast_gui():
	def __init__(self,root):
		#Initialization
		self.root = root 
		self.exit = -1
		self.dir = None

		

		#local
		_program_name = sys.argv[0]
		_directory_name = os.path.dirname(os.path.abspath(_program_name))

		self.argument_dict = {
								'query':'',
								'database':	_directory_name+"/directory/"}

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
		fasta_input_frame.pack(side=TOP,expand=0,fill=BOTH,padx=10)
		self._make_fasta_entry(fasta_input_frame)
		notebook_frame.add(f_and_d_frame,text="Files and Directories",underline=0,padding=2)
		#self.sub_notebook = ttk.Notebook(f_and_d_page, ipadx=5, ipady=5,width=250)
		#self.sub_notebook.add('database',label="Germline Database", underline=0,
		#	createcmd=lambda self=self,nb=self.sub_notebook,name='database': self.sub_notebook_database(nb,name))
		#self.sub_notebook.add('i_database',label="Internal Database",underline=0,
		#	createcmd=lambda self=self,nb=self.sub_notebook,name='i_database': self.sub_notebook_i_database(nb,name))
		#self.sub_notebook.add('aux_database',label="Auxillary Database",underline=0,
	#		createcmd=lambda self=self,nb=self.sub_notebook,name='aux_database': self.sub_notebook_aux_database(nb,name))
	#	self.sub_notebook.pack(side=LEFT,expand=1,fill=BOTH)
	#	self.filler = ttk.Frame(f_and_d_page,padx=5,pady=5,width=500).pack(side=LEFT,fill=BOTH,expand=1)
		#self.input_frame.grid(in_=f_and_d_page,row=0,column=0,columnspan=2)
	
	# def sub_notebook_database(self,nb, name):
	#	notebook_frame = nb.page(name)
		# msg = ttk.Message(notebook_frame,relief=ttk.FLAT, width=500,
		# 	text='This directory contains all the germline files compiled to blast database formats. If you haven\'t created one or downloaded a new one. Leave it at the default')
		# self.blast_dirlist = ttk.DirList(notebook_frame,width=200,height=300)
		# self.blast_dirlist['command']
		# msg.pack(side=ttk.BOTTOM,padx=10,expand=Y,pady=10)
		# self.blast_dirlist.pack(side=ttk.BOTTOM, padx=10, pady=10,expand=Y)
		# current_db = ttk.LabelFrame(notebook_frame,label="Current Selected Database:",width=200,padx=10)
		# current_db.pack(side=BOTTOM,expand=Y)
		# current_db_label = ttk.Message(current_db,relief=SUNKEN,bg='red',text=self.argument_dict['database'],width=500)
		# current_db_label.pack(side=LEFT,pady=30, padx=20)
	
	# def get_blast_cwd(self,dirlist):
	# 	print "hello", dirlist

	# def sub_notebook_i_database(self,nb, name):
	# 	notebook_frame = nb.page(name)
	# 	msg = ttk.Message(notebook_frame,relief=ttk.FLAT, width=500, anchor=ttk.NW,
	# 		text='This directory contains all the germline files compiled to blast database formats. If you haven\'t created one or downloaded a new one. Leave it at the default')
	# 	self.internal_dirlist = ttk.DirList(notebook_frame)
	# 	msg.pack(side=ttk.TOP, expand=1, fill=ttk.BOTH, padx=3, pady=3)
	# 	self.internal_dirlist.pack(side=ttk.TOP,padx=20, pady=20)
	
	# def sub_notebook_aux_database(self,nb, name):
	# 	notebook_frame = nb.page(name)	
	# 	msg = ttk.Message(notebook_frame,relief=ttk.FLAT, width=500, anchor=ttk.NW,
	# 		text='This directory contains all the germline files compiled to blast database formats. If you haven\'t created one or downloaded a new one. Leave it at the default')
	# 	self.aux_dirlist = ttk.DirList(notebook_frame)
	# 	msg.pack(side=ttk.TOP, expand=1, fill=ttk.BOTH, padx=3, pady=3)
	# 	self.aux_dirlist.pack(side=ttk.TOP, padx=3, pady=3)


	def _make_fasta_entry(self,fasta_input_frame):
		message = ttk.Label(fasta_input_frame,relief=FLAT, width=500, anchor=W,
								text='Enter the entry FASTA file here',font=('Arial',16))
		fasta_entry = ttk.Entry(fasta_input_frame,width=10)
		fasta_entry_button = ttk.Button(fasta_input_frame,text="Browse...",
			command=lambda entry=fasta_entry:self._enter_fasta(entry))
		message.pack(side=TOP,expand=1,fill=BOTH,padx=3,pady=3)
		fasta_entry.pack(side=LEFT,padx=3,expand=1,fill=X,pady=3)
		fasta_entry_button.pack(side=LEFT,fill=X)

	def _enter_fasta(self,entry):
		fn = None
		opts = {'title':"Select FASTA file to open...",
				'initialfile':entry.get()}
		fn = filedialog.askopenfilename(**opts)
		if fn:
			entry.delete(0,END)
			entry.insert(END,fn)

	def _create_readme(self,notebook_frame):
		readme_frame = ttk.Frame(notebook_frame,name="r_frame")
		#scroled_widget = ttk.ScrolledWindow(readme_frame,scrollbar='auto')
		#window = scroled.window
		msg = ttk.Label(readme_frame,text=open('README.md').readlines(),anchor=N)
		#scroled_widget.pack(side=TOP, fill=BOTH, expand=1)
		msg.pack(side=TOP,fill=BOTH,expand=1)
		notebook_frame.add(readme_frame,text="Readme",underline=0,padding=2)

	def quitcmd (self):
		'''Quit our mainloop'''
		self.root.destory()




if __name__ == '__main__':
	root = Tkinter.Tk()
	pyigblast_class = pyigblast_gui(root)
	root.mainloop()