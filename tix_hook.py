import os, os.path, sys, Tix
from Tkconstants import *
import tkFileDialog
import traceback, tkMessageBox
from Tkinter import *

class pyigblast_gui():
	def __init__(self,root):
		#Initialization
		self.root = root 
		self.exit = -1
		self.dir = None

		self.argument_dict = {'query':''}

		#local
		_program_name = sys.argv[0]
		_directory_name = os.path.dirname(_program_name)

	def MainMenu(self):
		main_menu = Tix.Frame(self.root,bd=2,relief=RAISED)
		author_label = Tix.Label(main_menu,text="Jordan Willis")
		university_label = Tix.Label(main_menu,text="Vanderbilt University")
		exit_button = Tix.Button(main_menu,text="EXIT",command= lambda self=self:self.quitcmd())
		author_label.pack(side=LEFT)
		exit_button.pack(side=LEFT,expand=YES)
		university_label.pack(side=RIGHT)
		return main_menu

	
	def TabNotebook(self):
		notebook_frame = self.root
		notebook = Tix.NoteBook(notebook_frame, ipadx=5, ipady=5, bg='black')
		notebook.add('f_and_d', label="Files and Databases", underline=0,
			createcmd=lambda self=self,nb=notebook,name='f_and_d': self.files_and_directories(nb,name))
		notebook.add('readme', label="Usage", underline=0,
			createcmd=lambda self=self,nb=notebook,name='readme': self.readme(nb,name) )
		return notebook

	def files_and_directories(self,nb,name):
		f_and_d_page = nb.page(name)
		options = "label.padX4"
		self.input_frame = Tix.LabelFrame(f_and_d_page,options=options)
		self.input_frame.pack(side=TOP,expand=0,fill=BOTH)
		self.make_fasta_entry()
		
		#self.input_frame.grid(in_=f_and_d_page,row=0,column=0,columnspan=2)
	
	def make_fasta_entry(self):
		message = Tix.Message(self.input_frame,relief=Tix.FLAT, width=500, anchor=W,
								text='Enter the entry FASTA file here',font=('Arial',16))
		self.fasta_entry = Tix.FileEntry(self.input_frame, label="Select a FASTA file:",selectmode="normal")
		message.pack(side=TOP,expand=1,fill=BOTH,padx=3,pady=3)
		self.fasta_entry.pack(side=TOP,fill=X,padx=3,pady=3)
	
	def readme(self,nb, name):
		notebook_frame = nb.page(name)
		text = self.readmeText(notebook_frame)
		text.pack(side=TOP, fill=BOTH, expand=1)

	def readmeText(self,nf):
		notebook_frame = nf
		sw = Tix.ScrolledWindow(notebook_frame, scrollbar='auto')
		win = sw.window #what to put inside of the window
		msg = Tix.Message(win,
					  bd=0, width=400, anchor=N,
					  text=open('README.md').readlines())
		msg.pack(side=LEFT,fill=BOTH, padx=10, pady=10)
		return sw



	def build(self):
		window_info = self.root.winfo_toplevel()
		window_info.wm_title('PyIgBLAST - GUI')
		#if window_info <= 800:
		window_info.geometry('1500x900+10+10')
		frame1 = self.MainMenu()
		frame1.pack(side=BOTTOM, fill=X)
		frame2 = self.TabNotebook()
		frame2.pack(side=TOP,expand=1,fill=BOTH,padx=5,pady=5)
		window_info.wm_protocol("WM_DELETE_WINDOW", lambda self=self:self.quitcmd())

	def quitcmd (self):
		'''Quit our mainloop'''
		print self.fasta_entry['value']
		self.exit = 0

	def loop(self):
		while self.exit < 0:
			# There are 2 whiles here. The outer one lets you continue
			# after a ^C interrupt.
			try:
				# This is the replacement for _tkinter mainloop()
				# It blocks waiting for the next Tcl event using select.
				while self.exit < 0:
					self.root.tk.dooneevent(0)
			except SystemExit:
				# Tkinter uses SystemExit to exit
				self.exit = 1
				return
			except KeyboardInterrupt:
				if tkMessageBox.askquestion ('Interrupt', 'Really Quit?') == 'yes':
					# self.tk.eval('exit')
					self.exit = 1
					return
				continue
			except:
				# Otherwise it's some other error - be nice and say why
				t, v, tb = sys.exc_info()
				text = ""
				for line in traceback.format_exception(t,v,tb):
					text += line + '\n'
				try: tkMessageBox.showerror ('Error', text)
				except: pass
				self.exit = 1
				raise SystemExit, 1
	def destroy(self):
		self.root.destroy()

def RunMain(root):
	'''Method for calling on Class'''
	pyigblast_class = pyigblast_gui(root)
	pyigblast_class.build()
	pyigblast_class.loop()
	pyigblast_class.destroy()


if __name__ == '__main__':
	root = Tix.Tk()
	RunMain(root) #replaces loop with custom function