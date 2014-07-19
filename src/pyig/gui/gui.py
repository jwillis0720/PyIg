#!/usr/bin/env python
import os
import os.path
import sys
import collections
import time
import pickle

# TK standard library
import Tkinter
import ttk
import tkFileDialog as filedialog
import tkMessageBox
from Tkconstants import *

# Non Standard Library
from multiprocessing import cpu_count, freeze_support
import gui_execute
from output_tabs_checkboxes import all_checkboxes_dict
from Bio import SeqIO as so

class Checkbar():

    def __init__(self, parent=None, picks=[], startrow=0):
        self.vars = collections.OrderedDict()
        row_count = startrow
        col_count = 0
        for pick in picks:
            json_key = pick[2]
            var = Tkinter.IntVar()
            var.set(int(pick[1]))
            chk = ttk.Checkbutton(parent, onvalue=1, offvalue=0, text=pick[0], variable=var, command=lambda: self.state())
            if col_count % 5 == 0:
                row_count += 1
                col_count = 0
            chk.grid(row=row_count, column=col_count, sticky=NW, padx=2, pady=2)
            col_count += 1
            self.vars[pick[0]] = {"state": var, "json_key": json_key}
            parent.rowconfigure(row_count, weight=3)
            parent.columnconfigure(col_count, weight=3)

    def state(self):
        states_dict = collections.OrderedDict()
        for var in self.vars:
            states_dict[var] = {"state": self.vars[var]['state'].get(), "json_key": self.vars[var]['json_key']}
        return states_dict


class IORedirector(object):

    '''A general class for redirecting I/O to this Text widget.'''

    def __init__(self, text_area):
        self.text_area = text_area

    def flush(self):
        '''
        Frame sys.stdout's flush method.
        '''
        pass


class StdoutRedirector(IORedirector):

    '''A class for redirecting stdout to this Text widget.'''

    def write(self, str):
        self.text_area.insert(END, str)
        self.text_area.see(END)


class PyIg_gui():

    '''GUI Class, calls gui_execute function when it grabs all the input data'''

    def __init__(self, root):    
        # Initialization
        self.root = root

        # get directories
        data_name = os.path.join(os.path.dirname(gui_execute.__file__),'library_dir.txt')
        self._directory_name = pickle.load(open(data_name))
        #self._directory_name = os.path.dirname(os.path.abspath(sys.argv[0]))
        self._user_directory = os.path.expanduser("~")

        # this will later become a widget. I initialize it here because it needs to update
        self.output_entry = ""

        # Each output field is a term in the all_checkboxes dict found in another file for code tidiness
        self.all_checkboxes_dict = all_checkboxes_dict

        if os.name == "posix":
            self.os = "mac"
        else:
            self.os = "windows"


        # argument dictionary we will pass to the arg parser eventually
        self.argument_dict = {
            'query': '',
            'database': os.path.join(self._directory_name, "database"),
            'in_data': os.path.join(self._directory_name, "internal_data"),
            'aux_data': os.path.join(self._directory_name, "optional_file"),
            'output_file': os.path.join(self._user_directory, "pyigblast_output"),
            'tmp_data': os.path.join(self._user_directory, "pyigblast_temporary")}


        #hmmm how do i get back to /usr/local/bin dynamically
        if self.os == "mac":
            self.argument_dict['executable'] = os.path.join(self._directory_name,
                                                            "../../bin/igblastn")
        elif self.os == "windows":
            self.argument_dict['executable'] = os.path.join(self._directory_name,
                                                            "../../bin/igblast.exe")
        # Basic info about the window
        window_info = self.root.winfo_toplevel()
        window_info.wm_title('PyIg - GUI')

        self.MainMenu()

        '''Creates main menu - this will use pack geometry manager but some of the sub-frames
        will use the grid geometry manager...that may be confusing but it is needed since
        some of the options just take up too much space and need to be used by the grid manager'''
        # now created tabbed notebook that will house input and output inside the main menu
        self.TabNotebook()

    # Mainmenu Holder
    def MainMenu(self):
        '''The main menu has 3 buttons and 2 labels that. The remaining calls on the tab notebookis a tabbed notebook'''

        # Main Menu will have a frame
        main_menu = ttk.Frame(self.root)
        author_label = ttk.Label(main_menu, text="Jordan Willis")
        university_label = ttk.Label(main_menu, text="Vanderbilt University")

        # 3 buttons in the bottom of the rame
        exit_button = ttk.Button(
            main_menu, text="Exit", command=lambda root=self.root: root.destroy())
        refresh_button = ttk.Button(
            main_menu, text="Refresh", command=lambda self=self: self._update())
        run_button = ttk.Button(
            main_menu, text="Run", command=lambda self=self: self.execute_dummy())


        # Pack widgets
        author_label.pack(side=LEFT, fill=X, padx=10, pady=10)
        run_button.pack(side=LEFT, expand=1, fill=X, padx=10, pady=10)
        refresh_button.pack(side=LEFT, expand=1, fill=X, padx=10, pady=10)
        exit_button.pack(side=LEFT, expand=1, fill=X, padx=10, pady=10)
        university_label.pack(side=RIGHT, fill=X, padx=10, pady=10)
        main_menu.pack(side=BOTTOM, fill=X, pady=10)


    # Notebook holds different tabs of the output gui
    def TabNotebook(self):
        # top level method, creates notbook inside of root. Then we will fill each tab
        self.main_notebook_frame = ttk.Notebook(self.root, name='main_notebook')
        self.main_notebook_frame.enable_traversal()
        self.main_notebook_frame.pack(side=TOP, expand=1, fill=BOTH)
        self._create_files_and_directories(self.main_notebook_frame)
        self._create_format_output(self.main_notebook_frame)
        self._create_output_stream(self.main_notebook_frame)
        self._create_readme(self.main_notebook_frame)

    # Execution functions - Run Button
    def execute_dummy(self):
        '''A dummy function that switches to the output frame then executes blast'''
        self.root.focus()
        self.main_notebook_frame.select(2)
        self.root.after(100, self.execute)

    def execute(self):
        # We can print anything and it will be captured by standard out
        print "Starting Execution"

        if not self.argument_dict['query']:
            tkMessageBox.showwarning(
                "Input Missing",
                "At the very least we need a fasta\n")

        print self.j_gene_numb.get()
	# Grabs all information from the GUI output tabs that are checked that need will be parsed from the BLAST
        output_options = [self.general_class.state(),
                          self.nucleotide_class.state(),
                          self.amino_class.state(),
                          self.total_alignments_class.state(),
                          self.fw1_alignments_class.state(),
                          self.fw2_alignments_class.state(),
                          self.fw3_alignments_class.state(),
                          self.cdr1_alignments_class.state(),
                          self.cdr2_alignments_class.state(),
                          self.cdr3_alignments_class.state(),
                          self.v_hits_class.state(),
                          self.d_hits_class.state(),
                          self.j_hits_class.state()]

        # These go into blast
        blast_args_dict = {
            '-query': "",  # Will be put in the executable file
            '-organism': self.species_var.get(),
            '-num_alignments_V': self.v_gene_numb.get(),
            '-num_alignments_D': self.d_gene_numb.get(),
            '-num_alignments_J': self.j_gene_numb.get(),
            '-min_D_match': str(self.min_d_match.get()),
            '-D_penalty': str(int(self.penalty_mismatch.get())),
            '-domain_system': self.scheme_var.get(),
            '-out': "",
            '-evalue': str(self.evalue.get()),
            '-word_size': str(self.word_size.get()),
            '-max_target_seqs': str(500),
            '-germline_db_V': "{0}_gl_V".format(
                os.path.join(self.argument_dict['database'], "Ig", self.species_var.get(), self.species_var.get())),
            '-germline_db_D': "{0}_gl_D".format(
                os.path.join(self.argument_dict['database'], "Ig", self.species_var.get(), self.species_var.get())),
            '-germline_db_J': "{0}_gl_J".format(
                os.path.join(self.argument_dict['database'], "Ig", self.species_var.get(), self.species_var.get())),
            '-auxiliary_data': "{0}_gl.aux".format(
                os.path.join(self.argument_dict['aux_data'], self.species_var.get())),
            '-domain_system': self.scheme_var.get(),
            '-outfmt': str(7)  # it has to be this tab seperated format
        }

        # These go into the parser
        output_options_dict = {
            'executable': self.argument_dict['executable'],
            'final_outfile': self.output_file_entry.get().split('.')[0],
            'num_procs': self.proc_count.get(),
            'pre_split_up_input': self.argument_dict['query'],
            'zip_bool': self.zip_var.get(),
            'tmp_data_directory': self.argument_dict['tmp_data'],
            'internal_data_directory': self.argument_dict['in_data'],
            'output_type': self.output_type_var.get(),
            'output_options': output_options
        }

        # both dictionaries are put into execute
        gui_execute.g_execute(blast_args_dict, output_options_dict)

    # Internal GUI Functions
    def _create_files_and_directories(self, notebook_frame):
        '''This is the first frame that will house 3 subframes with fasta directories and options'''

        # First frame that contains input file
        f_and_d_frame = ttk.Frame(
            notebook_frame, name='f_and_d')
        message = ttk.Label(
            anchor=W, text="Enter FASTA File", font=("Arial", 20))
        fasta_input_frame = ttk.LabelFrame(
            f_and_d_frame, labelwidget=message)
        fasta_input_frame.pack(
            side=TOP, expand=0, fill=X, padx=10)

        # Put in fasta_entry frame to take in input
        self._make_fasta_entry(fasta_input_frame)

        # Second Frame - Set up directory frame within the tab that takes in all the directories needed to run blast
        directory_label = ttk.Label(
            font=('Arial', 20), text="Directories for BLAST:")
        directories_frame = ttk.LabelFrame(
            f_and_d_frame, labelwidget=directory_label)
        directories_frame.pack(
            side=LEFT, expand=1, fill=BOTH, pady=10, padx=10)

        # Now run functions to place directories within that frame
        self._set_up_directories(directories_frame)

        # On the other side of this tab we will put in the basic options
        option_label = ttk.Label(
            font=('Arial', 20), text="Options:")
        basic_options_frame = ttk.LabelFrame(
            f_and_d_frame, labelwidget=option_label)
        basic_options_frame.pack(
            side=LEFT, expand=1, fill=BOTH, pady=10)

        # Set up defaults for basic options
        self._set_up_basic_options(basic_options_frame)

        # Now add all frames into the notebook frame which is in turn inside main frame :)
        notebook_frame.add(
            f_and_d_frame, text="Input Options", underline=0, padding=2)

    def _make_fasta_entry(self, fasta_input_frame):
        fasta_entry = ttk.Entry(fasta_input_frame)
        fasta_entry_button = ttk.Button(fasta_input_frame, width=25, text="Browse...",
                                        command=lambda entry=fasta_entry: self._enter_fasta(entry))
        fasta_entry.pack(side=LEFT, padx=3, expand=1, fill=X, pady=3)
        fasta_entry_button.pack(side=LEFT, fill=X)

    def _enter_fasta(self, entry):
        '''inputs and validates fasta entry'''
        fn = None
        _not_fasta = True
        opts = {'title': "Select FASTA file to open...",
                'initialfile': entry.get(),
                'initialdir': self._user_directory}
        while _not_fasta:
            fn = filedialog.askopenfilename(**opts)
            try:
                so.parse(str(fn), 'fasta').next()
                _not_fasta = False
                if fn:
                    entry.delete(0, END)
                    entry.insert(END, fn)
                    self.argument_dict['query'] = str(fn)
            except StopIteration:
                tkMessageBox.showwarning(
                    "Open file",
                    "Cannot open {0}, it is not a FASTA file\n".format(fn))

    def _set_up_directories(self, directories_frame):
        # blast directory
        blast_directories_frame = ttk.LabelFrame(directories_frame)
        blast_directories_frame.pack(side=TOP, fill=X, padx=10, expand=1)
        blast_directories_label = ttk.Label(
            blast_directories_frame, text="Compiled Blast Directory:", font=('Arial', 16))
        blast_directories_label.pack(side=TOP, anchor=NW)
        blast_directory_entry = ttk.Entry(blast_directories_frame)
        blast_directory_entry.insert(END, self.argument_dict['database'])
        blast_directory_entry.pack(side=LEFT, expand=1, fill=X, pady=3)
        blast_directory_button = ttk.Button(
            blast_directories_frame, text="Browse...",
            command=lambda type='blast', entry=blast_directory_entry: self._enter_directory(type, entry))
        blast_directory_button.pack(side=LEFT, fill=X)

        # internal_blast directory
        internal_blast_directories_frame = ttk.LabelFrame(directories_frame)
        internal_blast_directories_frame.pack(
            side=TOP, fill=X, padx=10, expand=1)
        internal_blast_directories_label = ttk.Label(
            internal_blast_directories_frame, text="Internal Blast Directory:", font=('Arial', 16))
        internal_blast_directories_label.pack(side=TOP, anchor=NW)
        internal_blast_directory_entry = ttk.Entry(
            internal_blast_directories_frame)
        internal_blast_directory_entry.insert(
            END, self.argument_dict['in_data'])
        internal_blast_directory_entry.pack(
            side=LEFT, expand=1, fill=X, pady=3)
        internal_blast_directory_button = ttk.Button(
            internal_blast_directories_frame, text="Browse...",
            command=lambda type='in_data', entry=internal_blast_directory_entry: self._enter_directory(type, entry))
        internal_blast_directory_button.pack(side=LEFT, fill=X)

        # auxilliary_blast directory
        aux_blast_directories_frame = ttk.LabelFrame(directories_frame)
        aux_blast_directories_frame.pack(side=TOP, fill=X, padx=10, expand=1)
        aux_blast_directories_label = ttk.Label(
            aux_blast_directories_frame, text="Auxillary Blast Directory:", font=('Arial', 16))
        aux_blast_directories_label.pack(side=TOP, anchor=NW)
        aux_blast_directory_entry = ttk.Entry(
            aux_blast_directories_frame)
        aux_blast_directory_entry.insert(END, self.argument_dict['aux_data'])
        aux_blast_directory_entry.pack(side=LEFT, expand=1, fill=X, pady=3)
        aux_blast_directory_button = ttk.Button(
            aux_blast_directories_frame, text="Browse...",
            command=lambda type='aux_data', entry=aux_blast_directory_entry: self._enter_directory(type, entry))
        aux_blast_directory_button.pack(side=LEFT, fill=X)

        # auxilliary_blast directory
        tmp_directories_frame = ttk.LabelFrame(directories_frame)
        tmp_directories_frame.pack(side=TOP, fill=X, padx=10, expand=1)
        tmp_directories_label = ttk.Label(
            tmp_directories_frame, text="Directory to store temporary files in:",
            font=('Arial', 16))
        tmp_directories_label.pack(side=TOP, anchor=NW)
        tmp_directory_entry = ttk.Entry(tmp_directories_frame)
        tmp_directory_entry.insert(END, self.argument_dict['tmp_data'])
        tmp_directory_entry.pack(side=LEFT, expand=1, fill=X, pady=3)
        tmp_directory_button = ttk.Button(
            tmp_directories_frame, text="Browse...",
            command=lambda type='tmp_data', entry=tmp_directory_entry: self._enter_directory(type, entry))
        tmp_directory_button.pack(side=LEFT, fill=X)

    def _enter_directory(self, type, entry):
        direc = None
        tkMessageBox.showinfo(
            "Consult:", "Please Read the README tab before messing with this option")
        opts = {'title': "Select FASTA file to open...",
                'initialdir': entry.get()}
        direc = filedialog.askdirectory(**opts)
        if direc:
            entry.delete(0, END)
            entry.insert(END, direc)
            if type == 'blast':
                self.argument_dict['database'] = str(direc)
            elif type == 'in_data':
                self.argument_dict['in_data'] = str(direc)
            elif type == 'aux_data':
                self.argument_dict['aux_data'] = str(direc)
            elif type == 'tmp_data':
                self.argument_dict['tmp_data'] = str(direc)

    def _set_up_basic_options(self, basic_options_frame):
        # basic options
        #imgt or kabat
        scheme_label = ttk.Label(basic_options_frame, text="Scheme Output:", font=('Arial', 16))
        scheme_label.grid(row=0, column=0, sticky=W, padx=5, pady=5)
        self.scheme_var = Tkinter.StringVar()
        self.scheme_var.set("imgt")
        radio_button_imgt = ttk.Radiobutton(
            basic_options_frame, text="IMGT", variable=self.scheme_var, value="imgt")
        radio_button_kabat = ttk.Radiobutton(
            basic_options_frame, text="KABAT", variable=self.scheme_var, value="kabat")
        radio_button_imgt.grid(row=1, column=0, sticky=NW, padx=10, rowspan=3)
        radio_button_kabat.grid(row=1, column=1, sticky=NW, padx=10, rowspan=3)
        basic_options_frame.rowconfigure(1, weight=3)

        # species frame
        self.species_var = Tkinter.StringVar()
        # set initial value
        self.species_var.set("human")

        # label
        species_label = ttk.Label(
            basic_options_frame, text="Species Type:", font=('Arial', 16))
        species_label.grid(row=5, column=0, sticky=NW, padx=5, pady=5)

        # 4 buttons
        species_button_human = ttk.Radiobutton(
            basic_options_frame, text="HUMAN",
            variable=self.species_var, value="human")
        species_button_human.grid(row=6, column=0, sticky=NW, padx=10)
        species_button_mouse = ttk.Radiobutton(
            basic_options_frame, text="MOUSE",
            variable=self.species_var, value="mouse")
        species_button_mouse.grid(row=6, column=1, sticky=NW, padx=10)
        species_button_rabbit = ttk.Radiobutton(
            basic_options_frame, text="RABBIT",
            variable=self.species_var, value="rabbit", state="disabled")
        species_button_rabbit.grid(row=6, column=2, sticky=NW, padx=10)
        species_button_rat = ttk.Radiobutton(basic_options_frame, text="RAT",
                                             variable=self.species_var, value="rat", state="disabled")
        species_button_rat.grid(row=6, column=3, sticky=NW, padx=10)
        basic_options_frame.rowconfigure(6, weight=3)

        # output_type
        self.output_type_var = Tkinter.StringVar()
        self.output_type_var.set("csv")

        # label
        scheme_label = ttk.Label(
            basic_options_frame, text="Output Format:", font=('Arial', 16))
        scheme_label.grid(row=9, column=0, sticky=W, padx=5, pady=5)

        # button
        radio_button_json_output = ttk.Radiobutton(
            basic_options_frame, text="JSON", variable=self.output_type_var,
            command=lambda suffix="json", self=self: self._update_output(suffix), value="json")
        radio_button_json_output.grid(row=10, column=0, sticky=NW, padx=10)
        radio_button_csv_output = ttk.Radiobutton(
            basic_options_frame, text="CSV", variable=self.output_type_var,
            command=lambda suffix="csv", self=self: self._update_output(suffix), value="csv")
        radio_button_csv_output.grid(row=10, column=1, sticky=NW, padx=10)
        radio_button_raw_blast_output = ttk.Radiobutton(
            basic_options_frame, text="BLAST", variable=self.output_type_var,
            command=lambda suffix="blast_out", self=self: self._update_output(suffix), value="blast_out")
        radio_button_raw_blast_output.grid(row=10, column=2, sticky=NW, padx=10)
        basic_options_frame.rowconfigure(10, weight=3)

        # Number of V D and J genes
        nvdj_type_label = ttk.Label(
            basic_options_frame, text="VDJ:", font=('Arial', 16))
        nvdj_type_label.grid(row=13, column=0, sticky=NW, padx=5, pady=5)

        # Vgene
        self.v_gene_numb = Tkinter.StringVar()
        numbs = [i for i in xrange(1, 4)]
        v_gene_label = ttk.Label(basic_options_frame, text="V-Gene Matches")
        v_gene_combo = ttk.Combobox(
            basic_options_frame, values=numbs, textvariable=self.v_gene_numb)
        v_gene_combo.current(0)
        v_gene_label.grid(row=14, column=0, sticky=NW, padx=10)
        v_gene_combo.grid(row=15, column=0, sticky=NW, padx=10)

        self.d_gene_numb = Tkinter.StringVar()
        numbs = [i for i in xrange(1, 4)]
        d_gene_label = ttk.Label(basic_options_frame, text="D-Gene Matches")
        d_gene_combo = ttk.Combobox(
            basic_options_frame, values=numbs, textvariable=self.d_gene_numb)
        d_gene_combo.current(0)
        d_gene_label.grid(row=14, column=1, sticky=NW, padx=10)
        d_gene_combo.grid(row=15, column=1, sticky=NW, padx=10)

        self.j_gene_numb = Tkinter.StringVar()
        numbs = [i for i in xrange(1, 4)]
        j_gene_label = ttk.Label(basic_options_frame, text="J-Gene Matches")
        j_gene_combo = ttk.Combobox(
            basic_options_frame, values=numbs, textvariable=self.j_gene_numb)
        j_gene_combo.current(0)
        j_gene_label.grid(row=14, column=2, sticky=NW, padx=10)
        j_gene_combo.grid(row=15, column=2, sticky=NW, padx=10)
        basic_options_frame.rowconfigure(15, weight=3)

        # Blast options
        blast_options_label = ttk.Label(
            basic_options_frame, text="Blast Options:", font=('Arial', 16))
        blast_options_label.grid(row=17, column=0, sticky=NW, padx=5, pady=5)

        # initial values
        self.evalue = Tkinter.DoubleVar()
        self.evalue.set(1)
        self.word_size = Tkinter.IntVar()
        self.word_size.set(4)
        self.penalty_mismatch = Tkinter.DoubleVar()
        self.penalty_mismatch.set(-4)
        self.min_d_match = Tkinter.IntVar()
        self.min_d_match.set(5)
        self.proc_count = Tkinter.IntVar()
        self.proc_count.set(cpu_count())

        # evalue
        e_value_label = ttk.Label(
            basic_options_frame, text="e-Value Threshold")
        e_value_entry = ttk.Entry(basic_options_frame)
        e_value_entry.insert(0, self.evalue.get())
        e_value_entry.bind('<Return>', self._validate_e_value)
        e_value_entry.bind('<FocusOut>', self._validate_e_value)
        e_value_label.grid(row=18, column=0, sticky=NW, padx=10)
        e_value_entry.grid(row=19, column=0, sticky=NW, padx=10)

        # word size
        word_size_label = ttk.Label(
            basic_options_frame, text="Word Size")
        word_size_entry = ttk.Entry(basic_options_frame)
        word_size_entry.insert(0, self.word_size.get())
        word_size_entry.bind('<Return>', self._validate_word_value)
        word_size_entry.bind('<FocusOut>', self._validate_word_value)
        word_size_label.grid(row=18, column=1, sticky=NW, padx=10)
        word_size_entry.grid(row=19, column=1, sticky=NW, padx=10)

        penalty_mismatch_label = ttk.Label(
            basic_options_frame, text="Penalty Mismatch")
        penalty_mismatch_entry = ttk.Entry(basic_options_frame)
        penalty_mismatch_entry.insert(0, self.penalty_mismatch.get())
        penalty_mismatch_entry.bind(
            '<Return>', self._validate_penalty_mismatch_value)
        penalty_mismatch_entry.bind(
            '<FocusOut>', self._validate_penalty_mismatch_value)
        penalty_mismatch_label.grid(row=18, column=2, sticky=NW, padx=10)
        penalty_mismatch_entry.grid(row=19, column=2, sticky=NW, padx=10)

        # Min D Nucleotides
        min_d_match_label = ttk.Label(
            basic_options_frame, text="Minimal Number of D Nucleotides")
        min_d_match_entry = ttk.Entry(basic_options_frame)
        min_d_match_entry.insert(0, self.min_d_match.get())
        min_d_match_entry.bind(
            '<Return>', self._validate_min_match_value)
        min_d_match_entry.bind(
            '<FocusOut>', self._validate_min_match_value)
        min_d_match_label.grid(row=20, column=0, sticky=NW, padx=10)
        min_d_match_entry.grid(row=21, column=0, sticky=NW, padx=10)

        # how many cpus to use
        proc_count_label = ttk.Label(
            basic_options_frame, text="Processors")
        proc_count_entry = ttk.Entry(basic_options_frame)
        proc_count_entry.insert(0, self.proc_count.get())
        proc_count_entry.bind(
            '<Return>', self._validate_proc_count_value)
        proc_count_entry.bind(
            '<FocusOut>', self._validate_proc_count_value)
        proc_count_label.grid(row=20, column=1, sticky=NW, padx=10)
        proc_count_entry.grid(row=21, column=1, sticky=NW, padx=10)

        # zipped
        self.zip_var = Tkinter.IntVar()
        self.zip_var.set(0)
        zip_chk = ttk.Checkbutton(
            basic_options_frame, onvalue=1, offvalue=0, text="Zip output file",
            variable=self.zip_var, command=lambda: self.zip_var.get())
        zip_chk.grid(row=21, column=2, sticky=NW, padx=10)
        basic_options_frame.rowconfigure(21, weight=3)

    def _update_output(self, suffix):
        direct = os.path.dirname(self.output_file_entry.get())
        current = os.path.basename(self.output_file_entry.get()).split('.')[0]
        self.output_file_entry.delete(0, END)
        self.output_file_entry.insert(END, direct + "/" + current + "." + suffix)

    def _validate_e_value(self, event):
        entry_widget = event.widget
        content = entry_widget.get()
        try:
            if content.strip() == "":
                tkMessageBox.showwarning(
                    "Empty Field",
                    "Enter a value\n")
                entry_widget.insert(0, self.evalue.get())
                entry_widget.focus_set()
            elif float(content) < 0:
                tkMessageBox.showwarning(
                    "To Low",
                    "Value must be positive\n")
                entry_widget.delete(0, END)
                entry_widget.insert(END, self.evalue.get())
                entry_widget.focus_set()
            elif float(content) > 1000:
                tkMessageBox.showwarning(
                    "To High",
                    "Value must be less than 1000\n")
                entry_widget.delete(0, END)
                entry_widget.insert(END, self.evalue.get())
                entry_widget.focus_set()
            else:
                self.evalue.set(float(content))
        except ValueError:
            tkMessageBox.showwarning(
                "Not Float",
                "Value must be a Number\n")
            entry_widget.delete(0, END)
            entry_widget.insert(END, self.evalue.get())
            entry_widget.focus_set()

    def _validate_word_value(self, event):
        entry_widget = event.widget
        content = entry_widget.get()
        try:
            if content.strip() == "":
                tkMessageBox.showwarning(
                    "Empty Field",
                    "Enter a value\n")
                entry_widget.insert(0, self.word_size.get())
                entry_widget.focus_set()
            elif int(content) < 1:
                tkMessageBox.showwarning(
                    "To Low",
                    "Value must be greater than 1\n")
                entry_widget.delete(0, END)
                entry_widget.insert(END, self.word_size.get())
                entry_widget.focus_set()

            elif int(content) > 10:
                tkMessageBox.showwarning(
                    "To High",
                    "Value must be less than 10\n")
                entry_widget.delete(0, END)
                entry_widget.insert(END, self.word_size.get())
                entry_widget.focus_set()
            else:
                self.evalue.set(int(content))
        except ValueError:
            tkMessageBox.showwarning(
                "Not Int",
                "Value must be an integer\n")
            entry_widget.delete(0, END)
            entry_widget.insert(END, self.word_size.get())
            entry_widget.focus_set()

    def _validate_penalty_mismatch_value(self, event):
        entry_widget = event.widget
        content = entry_widget.get()
        try:
            if content.strip() == "":
                tkMessageBox.showwarning(
                    "Empty Field",
                    "Enter a value\n")
                entry_widget.insert(0, self.penalty_mismatch.get())
                entry_widget.focus_set()
            elif float(content) < -6.0:
                tkMessageBox.showwarning(
                    "To Low",
                    "Value must be greater than -6\n")
                entry_widget.delete(0, END)
                entry_widget.insert(END, self.penalty_mismatch.get())
                entry_widget.focus_set()

            elif float(content) > 0:
                tkMessageBox.showwarning(
                    "To High",
                    "Value must be less than 0\n")
                entry_widget.delete(0, END)
                entry_widget.insert(END, self.penalty_mismatch.get())
                entry_widget.focus_set()
            else:
                self.evalue.set(float(content))
        except ValueError:
            tkMessageBox.showwarning(
                "Not Float",
                "Value must be a Number\n")
            entry_widget.delete(0, END)
            entry_widget.insert(END, self.penalty_mismatch.get())
            entry_widget.focus_set()

    def _validate_min_match_value(self, event):
        entry_widget = event.widget
        content = entry_widget.get()
        try:
            if content.strip() == "":
                tkMessageBox.showwarning(
                    "Empty Field",
                    "Enter a value\n")
                entry_widget.insert(0, self.min_d_match.get())
                entry_widget.focus_set()
            elif int(content) < 5:
                tkMessageBox.showwarning(
                    "To Low",
                    "Value must be greater than 5\n")
                entry_widget.delete(0, END)
                entry_widget.insert(END, self.min_d_match.get())
                entry_widget.focus_set()

            elif int(content) > 15:
                tkMessageBox.showwarning(
                    "To High",
                    "Value must be less than 15\n")
                entry_widget.delete(0, END)
                entry_widget.insert(END, self.min_d_match.get())
                entry_widget.focus_set()
            else:
                self.evalue.set(int(content))
        except ValueError:
            tkMessageBox.showwarning(
                "Not Int",
                "Value must be an integer\n")
            entry_widget.delete(0, END)
            entry_widget.insert(END, self.min_d_match.get())
            entry_widget.focus_set()

    def _validate_proc_count_value(self, event):
        entry_widget = event.widget
        content = entry_widget.get()
        try:
            if content.strip() == "":
                tkMessageBox.showwarning(
                    "Empty Field",
                    "Enter a value\n")
                entry_widget.insert(0, self.proc_count.get())
                entry_widget.focus_set()
            elif int(content) < 0:
                tkMessageBox.showwarning(
                    "To Low",
                    "Value must be positive\n")
                entry_widget.delete(0, END)
                entry_widget.insert(END, self.proc_count.get())
                entry_widget.focus_set()

            elif int(content) > self.proc_count.get():
                tkMessageBox.showwarning(
                    "To High",
                    "Value is higher than the amount of processors you have\n")
                entry_widget.delete(0, END)
                entry_widget.insert(END, self.proc_count.get())
                entry_widget.focus_set()
            else:
                self.evalue.set(int(content))
        except ValueError:
            tkMessageBox.showwarning(
                "Not Int",
                "Value must be an integer\n")
            entry_widget.delete(0, END)
            entry_widget.insert(END, self.proc_count.get())
            entry_widget.focus_set()

    def _create_format_output(self, main_notebook_frame):
        '''second tab holds which output fields we want'''
        # second tab
        output_frame = ttk.Frame(main_notebook_frame)
        output_frame.pack(side=TOP, expand=1, fill=BOTH)
        self._make_output_file(output_frame)
        self._fill_output_format_tab(output_frame)
        main_notebook_frame.add(
            output_frame, text="Output Options", underline=0, padding=2)

    def _make_output_file(self, output_frame):
        output_file_frame = ttk.Frame(output_frame)
        output_file_frame.pack(side=TOP, fill=X, padx=2, pady=2)
        output_file_label = ttk.Label(
            output_file_frame, text='Enter Output File Name',
            font=('Arial', 20))
        output_file_label.pack(side=TOP, anchor=NW, padx=2, pady=2)

        self.output_file_entry = ttk.Entry(output_file_frame)
        self.output_file_entry.delete(0, END)
        self.output_file_entry.insert(END, self.argument_dict['output_file'] +
                                      "." + str(self.output_type_var.get()))
        output_file_entry_button = ttk.Button(
            output_file_frame, text="Browse...", width=25,
            command=lambda entry=self.output_file_entry: self._enter_output(entry))
        self.output_file_entry.pack(side=LEFT, padx=3, expand=1, fill=X)
        output_file_entry_button.pack(side=LEFT, fill=X)

    def _fill_output_format_tab(self, output_frame):
        # fill output_format_tab from another dictionary
        output_fields_label = ttk.Label(
            output_frame, text="Select Output Fields", font=('Arial', 20))
        output_fields_label.pack(side=TOP, anchor=NW, padx=5, pady=5)

        output_labels_nb_frame = ttk.Notebook(output_frame, name='output_labels')
        output_labels_nb_frame.enable_traversal()
        output_labels_nb_frame.pack(side=TOP, expand=1, fill=BOTH)

        # general frame output
        general_frame = ttk.Frame(output_frame)
        _general_options = self.all_checkboxes_dict['general']
        self.general_class = Checkbar(parent=general_frame, picks=[
                                      (i['formal'], i['default'], i['json_key']) for i in _general_options])
        output_labels_nb_frame.add(general_frame, text="General", underline=0, padding=2)

        # nucleotide
        nucleotide_frame = ttk.Frame(output_frame)
        _nucleotide_options = self.all_checkboxes_dict['nucleotide']
        self.nucleotide_class = Checkbar(parent=nucleotide_frame, picks=[
            (i['formal'], i['default'], i['json_key']) for i in _nucleotide_options])
        output_labels_nb_frame.add(nucleotide_frame, text="Nucleotide", underline=0, padding=2)

        # Translation Specific
        amino_frame = ttk.Frame(output_frame)
        _amino_options = self.all_checkboxes_dict['amino']
        self.amino_class = Checkbar(parent=amino_frame,
                                    picks=[(i['formal'], i['default'], i['json_key']) for i in _amino_options])
        output_labels_nb_frame.add(amino_frame, text="Amino Acid", underline=0, padding=2)

        # Alignment_frame
        alignment_frame = ttk.Frame(output_frame)
        _total_options = self.all_checkboxes_dict['total_alignments']
        self.total_alignments_class = Checkbar(parent=alignment_frame,
                                               picks=[(i['formal'], i['default'], i['json_key']) for i in _total_options], startrow=0)
        _fw1_options = self.all_checkboxes_dict['fw1_alignments']
        self.fw1_alignments_class = Checkbar(parent=alignment_frame,
                                             picks=[(i['formal'], i['default'], i['json_key']) for i in _fw1_options], startrow=1)

        _fw2_options = self.all_checkboxes_dict['fw2_alignments']
        self.fw2_alignments_class = Checkbar(parent=alignment_frame,
                                             picks=[(i['formal'], i['default'], i['json_key']) for i in _fw2_options], startrow=2)

        _fw3_options = self.all_checkboxes_dict['fw3_alignments']
        self.fw3_alignments_class = Checkbar(parent=alignment_frame,
                                             picks=[(i['formal'], i['default'], i['json_key']) for i in _fw3_options], startrow=3)

        _cdr1_options = self.all_checkboxes_dict['cdr1_alignments']
        self.cdr1_alignments_class = Checkbar(parent=alignment_frame,
                                              picks=[(i['formal'], i['default'], i['json_key']) for i in _cdr1_options], startrow=4)
        _cdr2_options = self.all_checkboxes_dict['cdr2_alignments']
        self.cdr2_alignments_class = Checkbar(parent=alignment_frame,
                                              picks=[(i['formal'], i['default'], i['json_key']) for i in _cdr2_options], startrow=5)

        _cdr3_options = self.all_checkboxes_dict['cdr3_alignments']
        self.cdr3_alignments_class = Checkbar(parent=alignment_frame,
                                              picks=[(i['formal'], i['default'], i['json_key']) for i in _cdr3_options], startrow=6)
        output_labels_nb_frame.add(alignment_frame, text="Alignment Frame", underline=0, padding=2)

        hits_frame = ttk.Frame(output_frame)
        _v_hits_options = self.all_checkboxes_dict['v_hits']
        self.v_hits_class = Checkbar(parent=hits_frame,
                                     picks=[(i['formal'], i['default'], i['json_key']) for i in _v_hits_options])
        _d_hits_options = self.all_checkboxes_dict['d_hits']
        self.d_hits_class = Checkbar(parent=hits_frame,
                                     picks=[(i['formal'], i['default'], i['json_key']) for i in _d_hits_options], startrow=2)
        _j_hits_options = self.all_checkboxes_dict['j_hits']
        self.j_hits_class = Checkbar(parent=hits_frame,
                                     picks=[(i['formal'], i['default'], i['json_key']) for i in _j_hits_options], startrow=4)
        output_labels_nb_frame.add(hits_frame, text="VDJ Hits", underline=0, padding=2)

        output_options_list = [self.general_class, self.nucleotide_class, self.amino_class, self.total_alignments_class,
                               self.fw1_alignments_class, self.fw2_alignments_class, self.fw3_alignments_class, self.cdr1_alignments_class, self.cdr2_alignments_class,
                               self.cdr3_alignments_class, self.v_hits_class, self.d_hits_class, self.j_hits_class]
        select_all_button = ttk.Button(output_frame, text="Select All", command=lambda o_options=output_options_list, self=self: self._select_all_outputs(o_options, 1))
        select_all_button.pack(side=LEFT, padx=2, pady=2)

        select_all_button = ttk.Button(output_frame, text="Deselect All", command=lambda o_options=output_options_list, self=self: self._select_all_outputs(o_options, 0))
        select_all_button.pack(side=LEFT, padx=2, pady=2)

    def _enter_output(self, entry):
        fo = None
        opts = {'title': "Select output file to open...",
                'initialdir': self._user_directory}
        fo = filedialog.asksaveasfilename(**opts)
        if fo:
            basename = os.path.splitext(str(fo))[0]
            out_file = basename + "." + self.output_type_var.get()
            entry.delete(0, END)
            entry.insert(END, out_file)
            self.argument_dict['output_file'] = out_file

    def _select_all_outputs(self, o_options, on):
        for instance in o_options:
            for options in instance.vars:
                instance.vars[options]['state'].set(on)

    def _create_output_stream(self, notebook_frame):
        self.create_output_frame = ttk.Frame(notebook_frame, name="o_frame")
        self.output_stream_text = Tkinter.Text(self.create_output_frame)
        self.output_stream_text.pack(side=LEFT, expand=1, fill=BOTH, anchor=NW)
        #comment out next two lines for regular output
        sys.stdout = StdoutRedirector(self.output_stream_text)
        sys.stderr = StdoutRedirector(self.output_stream_text)
        scroll = ttk.Scrollbar(self.create_output_frame)
        scroll.pack(side=RIGHT, fill=Y)

        scroll.config(command=self.output_stream_text.yview)
        self.output_stream_text.config(yscrollcommand=scroll.set)
        self.create_output_frame.pack(side=TOP, expand=1, fill=BOTH)
        self.output_stream_text.bind('<Key>', lambda e: 'break')
        # self.create_output_frame.update_idletasks()
        notebook_frame.add(self.create_output_frame, text="Output Stream", underline=0, padding=2)

    def _create_readme(self, notebook_frame):
        readme_frame = ttk.Frame(notebook_frame, name="r_frame")
        readme_text = Tkinter.Text(readme_frame)
        readme_text.pack(side=LEFT, expand=1, fill=BOTH, anchor=NW)
        scroll_bar = ttk.Scrollbar(readme_frame)
        scroll_bar.pack(side=RIGHT, fill=Y)
        scroll_bar.config(command=readme_text.yview)
        readme_text.config(yscrollcommand=scroll_bar.set)
        readme_frame.pack(side=TOP, expand=1, fill=BOTH)
        notebook_frame.update_idletasks()
        for line in open(os.path.join(self._directory_name, 'README_gui.txt')).readlines():
            readme_text.insert(END, line)
            readme_text.see(END)
            readme_text.update_idletasks()
        notebook_frame.add(readme_frame, text="Readme", underline=0, padding=2)

    def _update(self):
        import gui
        gui.main_refresh(self.root, gui)


def center(win):
    win.update_idletasks()
    frm_width = win.winfo_rootx() - win.winfo_x()
    win_width = win.winfo_width() + (frm_width * 2)
    titlebar_height = win.winfo_rooty() - win.winfo_y()
    win_height = win.winfo_height() + (titlebar_height + frm_width)
    x = (win.winfo_screenwidth() / 2) - (win_width / 2)
    y = (win.winfo_screenheight() / 2) - (win_height / 2)
    geom = (win.winfo_width(), win.winfo_height(), x, y)  # see note
    win.geometry('{0}x{1}+{2}+{3}'.format(*geom))


def main_refresh(root, gui):
    reload(gui)
    root.destroy()
    gui.main()


def main():
    root = Tkinter.Tk()
    PyIg_gui(root)
    center(root)
    root.mainloop()
    root.update_idletasks()


if __name__ == '__main__':
    if os.name == "nt":
        freeze_support()
        main()
    else:
        main()
