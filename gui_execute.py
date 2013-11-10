import os
import os.path
import sys
import Tkinter
import ttk
import tkFileDialog as filedialog
import tkMessageBox
from Tkconstants import *
from Bio import SeqIO as so
from vertical_scroll import VerticalScrolledFrame as vsf
from output_tabs_checkboxes import all_checkboxes
from multiprocessing import cpu_count


class pyigblast_gui():

    '''gui class, will do everything including calling on executables'''

    def __init__(self, root):
        # Initialization
        self.root = root
        _program_name = sys.argv[0]
        self._directory_name = os.path.dirname(os.path.abspath(_program_name))
        self._user_directory = os.path.expanduser("~")
        # argument dictionary we will pass to the arg parser eventually
        self.argument_dict = {
            'query': '',
            'database': self._directory_name + "/directory/",
            'in_data': self._directory_name + "/internal_data/",
            'aux_data': self._directory_name + "/optional_data/",
            'output_file': self._user_directory + "/pyigblast_output",
            'tmp_data': self._user_directory + "/pyigblast_temporary/"}
        window_info = self.root.winfo_toplevel()
        window_info.wm_title('PyIgBLAST - GUI')
        window_info.geometry('1500x900+10+10')
        # creates main menu
        self.MainMenu()
        # creates main menu notebook inside
        self.TabNotebook()

    def MainMenu(self):
        main_menu = ttk.Frame(self.root)
        author_label = ttk.Label(main_menu, text="Jordan Willis")
        university_label = ttk.Label(main_menu, text="Vanderbilt University")
        exit_button = ttk.Button(
            main_menu, text="Exit",
            command=lambda root=self.root: root.destroy())
        refresh_button = ttk.Button(
            main_menu, text="Refresh", command=lambda self=self: self._update())
        author_label.pack(side=LEFT, padx=10, pady=3)
        refresh_button.pack(side=LEFT, padx=10, expand=YES, pady=3)
        exit_button.pack(side=LEFT, expand=YES, pady=3)
        university_label.pack(side=RIGHT, padx=10, pady=3)
        main_menu.pack(side=BOTTOM, fill=X)

    def TabNotebook(self):
        main_notebook_frame = ttk.Notebook(self.root, name='main_notebook')
        main_notebook_frame.enable_traversal()
        main_notebook_frame.pack(side=TOP, expand=1, fill=BOTH)
        self._create_files_and_directories(main_notebook_frame)
        self._create_format_output(main_notebook_frame)
        self._create_readme(main_notebook_frame)

    def _create_files_and_directories(self, notebook_frame):
        # This is the first tab that houses files, directories and Options
        # First frame in the notebook
        f_and_d_frame = ttk.Frame(notebook_frame, name='f_and_d')
        fasta_input_frame = ttk.LabelFrame(f_and_d_frame)
        fasta_input_frame.pack(side=TOP, expand=0, fill=X, padx=10)
        # put in fasta_entry frame to take in input
        self._make_fasta_entry(fasta_input_frame)

        # Set up directory frame within the tab that takes in all the
        # directories needed to run run blast
        directories_frame = ttk.LabelFrame(f_and_d_frame)
        directories_frame.pack(side=LEFT, expand=1, fill=BOTH)
        directory_label = ttk.Label(directories_frame, font=('Arial', 20),
                                    text="Directories needed to run blast:")
        directory_label.pack(side=TOP, fill=X, padx=20, pady=10)

        # set now run functions to place directories within that frame
        self._set_up_directories(directories_frame)

        # On the other side of this tab we will put in the basic options
        # including output type and blast
        basic_options_frame = ttk.LabelFrame(f_and_d_frame)
        basic_options_frame.pack(side=LEFT, fill=BOTH, padx=5)
        self._set_up_basic_options(basic_options_frame)

        # and add it to the big frame
        notebook_frame.add(
            f_and_d_frame, text="Input Options", underline=0, padding=2)

    def _make_fasta_entry(self, fasta_input_frame):
        message = ttk.Label(
            fasta_input_frame, relief=FLAT, width=500, anchor=W,
            text='Enter the entry FASTA file here', font=('Arial', 20))
        fasta_entry = ttk.Entry(fasta_input_frame, width=10)
        fasta_entry_button = ttk.Button(fasta_input_frame, text="Browse...",
                                        command=lambda entry=fasta_entry: self._enter_fasta(entry))
        message.pack(side=TOP, expand=1, fill=BOTH, padx=3, pady=3)
        fasta_entry.pack(side=LEFT, padx=3, expand=1, fill=X, pady=3)
        fasta_entry_button.pack(side=LEFT, fill=X)

    def _enter_fasta(self, entry):
        fn = None
        _not_fasta = True
        opts = {'title': "Select FASTA file to open...",
                'initialfile': entry.get()}
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
            entry.insert(END, fn)
            if type == 'blast':
                self.argument_dict['database'] = str(direc)
            elif type == 'in_data':
                self.argument_dict['in_data'] = str(direc)
            elif type == 'aux_data':
                self.argument_dict['aux_data'] = str(direct)
            elif type == 'tmp_data':
                self.argument_dict['tmp_data'] = str(direct)

    def _set_up_basic_options(self, basic_options_frame):
        # basic options

        #imgt or kabat
        scheme_frame = ttk.LabelFrame(basic_options_frame)
        scheme_frame.pack(side=TOP, fill=X, expand=1, padx=5, pady=5)
        self._set_up_scheme_frame(scheme_frame)

        # heavy or light chain
        chain_type_frame = ttk.LabelFrame(basic_options_frame)
        chain_type_frame.pack(side=TOP, fill=X, padx=5, pady=5)
        self._set_up_chain_type_frame(chain_type_frame)

        # heavy or light chain
        species_type_frame = ttk.LabelFrame(basic_options_frame)
        species_type_frame.pack(side=TOP, fill=X, expand=1, padx=5, pady=5)
        self._set_up_species_type_frame(species_type_frame)

        # output_type
        output_type_frame = ttk.LabelFrame(basic_options_frame)
        output_type_frame.pack(side=TOP, fill=X, expand=1, padx=5, pady=5)
        self._set_up_output_type_frame(output_type_frame)

        # Number of V D and J genes
        nvdj_type_frame = ttk.LabelFrame(
            basic_options_frame, text="VDJ Matches")
        nvdj_type_frame.pack(side=TOP, fill=X, expand=1, padx=5, pady=5)
        self._set_up_nvdj_type_frame(nvdj_type_frame)

        # Blast options
        blast_options_frame = ttk.LabelFrame(
            basic_options_frame, text="Blast Options")
        blast_options_frame.pack(side=TOP, fill=X, expand=1, padx=5, pady=5)
        self._set_up_blast_options(blast_options_frame)

        # zipped
        zip_bool_type_frame = ttk.LabelFrame(
            basic_options_frame, text="Zip output file")
        zip_bool_type_frame.pack(side=BOTTOM, anchor=SW, padx=5, pady=5)
        self.zip_var = Tkinter.IntVar()
        self.zip_var.set(1)
        zip_chk = ttk.Checkbutton(
            zip_bool_type_frame, onvalue=1, offvalue=0, text="Zip output file",
            variable=self.zip_var, command=lambda: self.zip_var.get())
        zip_chk.pack(side=TOP, anchor=SW)

    def _set_up_scheme_frame(self, scheme_frame):
        self.scheme_var = Tkinter.StringVar()
        self.scheme_var.set("imgt")
        scheme_label = ttk.Label(
            scheme_frame, text="Scheme output:", font=('Arial', 16))
        scheme_label.pack(side=TOP, anchor=NW)
        radio_button_imgt = ttk.Radiobutton(
            scheme_frame, text="IMGT", variable=self.scheme_var, value="imgt")
        radio_button_imgt.pack(side=LEFT, fill=X, expand=1)
        radio_button_kabat = ttk.Radiobutton(
            scheme_frame, text="KABAT", variable=self.scheme_var, value="kabat")
        radio_button_kabat.pack(side=LEFT, fill=X, expand=1)

    def _set_up_chain_type_frame(self, chain_type_frame):
        self.chain_var = Tkinter.StringVar()
        self.chain_var.set("heavy")
        chain_label = ttk.Label(
            chain_type_frame, text="Chain Type:", font=('Arial', 16))
        chain_label.pack(side=TOP, anchor=NW)
        chain_button_heavy = ttk.Radiobutton(
            chain_type_frame, text="HEAVY", variable=self.chain_var, value="heavy")
        chain_button_heavy.pack(side=LEFT, fill=X, expand=1)
        chain_button_light = ttk.Radiobutton(
            chain_type_frame, text="LIGHT", variable=self.chain_var, value="light")
        chain_button_light.pack(side=LEFT, fill=X, expand=1)

    def _set_up_species_type_frame(self, species_type_frame):
        self.species_var = Tkinter.StringVar()
        self.species_var.set("human")
        species_label = ttk.Label(
            species_type_frame, text="Species Type:", font=('Arial', 16))
        species_label.pack(side=TOP, anchor=NW)
        species_button_human = ttk.Radiobutton(
            species_type_frame, text="HUMAN",
            variable=self.species_var, value="human")
        species_button_human.pack(side=LEFT, fill=X, expand=1)
        species_button_mouse = ttk.Radiobutton(
            species_type_frame, text="MOUSE",
            variable=self.species_var, value="mouse")
        species_button_mouse.pack(side=LEFT, fill=X, expand=1)
        species_button_rabbit = ttk.Radiobutton(
            species_type_frame, text="RABBIT",
            variable=self.species_var, value="rabbit", state="disabled")
        species_button_rabbit.pack(side=LEFT, fill=X, expand=1)
        species_button_rat = ttk.Radiobutton(species_type_frame, text="RAT",
                                             variable=self.species_var, value="rat", state="disabled")
        species_button_rat.pack(side=LEFT, fill=X, expand=1)

    def _set_up_output_type_frame(self, output_type_frame):
        self.output_type_var = Tkinter.StringVar()
        self.output_type_var.set("json")
        scheme_label = ttk.Label(
            output_type_frame, text="Output format:", font=('Arial', 16))
        scheme_label.pack(side=TOP, anchor=NW)
        radio_button_json_output = ttk.Radiobutton(
            output_type_frame, text="JSON", variable=self.output_type_var, value="json")
        radio_button_json_output.pack(side=LEFT, fill=X, expand=1)
        radio_button_csv_output = ttk.Radiobutton(
            output_type_frame, text="CSV", variable=self.output_type_var, value="csv")
        radio_button_csv_output.pack(side=LEFT, fill=X, expand=1)
        radio_button_raw_blast_output = ttk.Radiobutton(
            output_type_frame, text="BLAST", variable=self.output_type_var, value="blast_out")
        radio_button_raw_blast_output.pack(side=LEFT, fill=X, expand=1)

    def _set_up_nvdj_type_frame(self, nvdj_type_frame):
        self.v_gene_numb = ""
        numbs = [i for i in xrange(1, 4)]
        v_gene_label = ttk.LabelFrame(nvdj_type_frame, text="V-Gene Matches")
        v_gene_combo = ttk.Combobox(
            v_gene_label, values=numbs, textvariable=self.v_gene_numb)
        v_gene_combo.current(0)
        v_gene_label.pack(side=LEFT, expand=1, pady=5, padx=20, fill=X)
        v_gene_combo.pack(side=TOP, expand=1, pady=5, padx=10, fill=X)

        self.d_gene_numb = ""
        numbs = [i for i in xrange(1, 4)]
        d_gene_label = ttk.LabelFrame(nvdj_type_frame, text="D-Gene Matches")
        d_gene_combo = ttk.Combobox(
            d_gene_label, values=numbs, textvariable=self.d_gene_numb)
        d_gene_combo.current(0)
        d_gene_label.pack(side=LEFT, expand=1, pady=5, padx=20, fill=X)
        d_gene_combo.pack(side=TOP, expand=1, pady=5)

        self.j_gene_numb = ""
        numbs = [i for i in xrange(1, 4)]
        j_gene_label = ttk.LabelFrame(nvdj_type_frame, text="J-Gene Matches")
        j_gene_combo = ttk.Combobox(
            j_gene_label, values=numbs, textvariable=self.j_gene_numb)
        j_gene_combo.current(0)
        j_gene_label.pack(side=LEFT, expand=1, pady=5, padx=20, fill=X)
        j_gene_combo.pack(side=TOP, expand=1, pady=5, padx=10, fill=X)

    def _set_up_blast_options(self, blast_options_frame):
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
        e_value_label = ttk.LabelFrame(
            blast_options_frame, text="e-Value Threshold")
        e_value_entry = ttk.Entry(e_value_label)
        e_value_entry.insert(0, self.evalue.get())
        e_value_entry.bind('<Return>', self._validate_e_value)
        e_value_entry.bind('<FocusOut>', self._validate_e_value)
        e_value_label.pack(side=LEFT, expand=1, pady=5, padx=5, fill=X)
        e_value_entry.pack(side=TOP, expand=1, pady=5, padx=5, fill=X)

        # word size
        word_size_label = ttk.LabelFrame(
            blast_options_frame, text="Word Size")
        word_size_entry = ttk.Entry(word_size_label)
        word_size_entry.insert(0, self.word_size.get())
        word_size_entry.bind('<Return>', self._validate_word_value)
        word_size_entry.bind('<FocusOut>', self._validate_word_value)
        word_size_label.pack(side=LEFT, expand=1, pady=5, padx=5, fill=X)
        word_size_entry.pack(side=TOP, expand=1, pady=5, padx=5, fill=X)

        penalty_mismatch_label = ttk.LabelFrame(
            blast_options_frame, text="Penalty Mismatch")
        penalty_mismatch_entry = ttk.Entry(penalty_mismatch_label)
        penalty_mismatch_entry.insert(0, self.penalty_mismatch.get())
        penalty_mismatch_entry.bind(
            '<Return>', self._validate_penalty_mismatch_value)
        penalty_mismatch_entry.bind(
            '<FocusOut>', self._validate_penalty_mismatch_value)
        penalty_mismatch_label.pack(side=LEFT, expand=1, pady=5,
                                    padx=5, fill=X)
        penalty_mismatch_entry.pack(side=TOP, expand=1, pady=5,
                                    padx=5, fill=X)
        # Min D Nucleotides
        min_d_match_label = ttk.LabelFrame(
            blast_options_frame, text="Minimal Number of D Nucleotides")
        min_d_match_entry = ttk.Entry(min_d_match_label)
        min_d_match_entry.insert(0, self.min_d_match.get())
        min_d_match_entry.bind(
            '<Return>', self._validate_min_match_value)
        min_d_match_entry.bind(
            '<FocusOut>', self._validate_min_match_value)
        min_d_match_label.pack(side=LEFT, expand=1, pady=5, padx=5, fill=X)
        min_d_match_entry.pack(side=TOP, expand=1, pady=5, padx=5, fill=X)

        # how many cpus to use
        proc_count_label = ttk.LabelFrame(
            blast_options_frame, text="Processors")
        proc_count_entry = ttk.Entry(proc_count_label)
        proc_count_entry.insert(0, self.proc_count.get())
        proc_count_entry.bind(
            '<Return>', self._validate_proc_count_value)
        proc_count_entry.bind(
            '<FocusOut>', self._validate_proc_count_value)
        proc_count_label.pack(side=LEFT, expand=1, pady=5, padx=5, fill=X)
        proc_count_entry.pack(side=TOP, expand=1, pady=5, padx=5, fill=X)

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
        # secon tab
        output_frame = ttk.LabelFrame(main_notebook_frame)
        output_frame.pack(side=TOP, expand=1, fill=BOTH, padx=10, pady=10)
        self._make_output_file(output_frame)
        self._fill_output_format_tab(output_frame)
        main_notebook_frame.add(
            output_frame, text="Output Options", underline=0, padding=2)

    def _make_output_file(self, output_frame):
        output_file_frame = ttk.LabelFrame(output_frame)
        output_file_frame.pack(side=TOP, fill=X, padx=10, pady=10)
        output_file_label = ttk.Label(
            output_file_frame, text='Enter Output File Name',
            font=('Arial', 26))
        output_file_label.pack(side=TOP, anchor=NW, padx=5, pady=5)

        output_file_entry = ttk.Entry(output_file_frame, width=10)
        output_file_entry.delete(0, END)
        output_file_entry.insert(END, self.argument_dict['output_file'] +
                                 "." + str(self.output_type_var.get()))
        output_file_entry_button = ttk.Button(
            output_file_frame, text="Browse...",
            command=lambda entry=output_file_entry: self._enter_output(entry))
        output_file_entry.pack(side=LEFT, padx=3, expand=1, fill=X, pady=3)
        output_file_entry_button.pack(side=LEFT, fill=X)

    def _fill_output_format_tab(self, output_frame):
        # fill output_format_tab from another dictionary
        output_fields_frame = ttk.LabelFrame(output_frame)
        output_fields_top_label = ttk.Label(
            output_fields_frame, text="Select Output Fields", font=('Arial', 20))
        output_fields_top_label.pack(side=TOP, padx=3, pady=3, anchor=NW)
        output_fields_frame.pack(
            side=TOP, fill=BOTH, expand=1, padx=10, pady=10)

        # general options
        general_frame = ttk.LabelFrame(output_fields_frame)
        general_label = ttk.Label(
            general_frame, text="General fields:", font=('Arial', 20))
        general_label.pack(side=TOP, anchor=NW, padx=5, pady=5)
        _general_options = all_checkboxes['general']
        general_class = Checkbar(parent=general_frame, picks=[
                                 (i['formal'], i['default']) for i in _general_options])
        general_frame.pack(side=TOP, expand=1, fill=X, padx=10)

        # nucleotide
        nucleotide_frame = ttk.LabelFrame(output_fields_frame)
        nucleotide_label = ttk.Label(nucleotide_frame, font=('Arial', 20),
                                     text="Nucleotide Specific:")
        nucleotide_label.pack(side=TOP, anchor=NW, padx=5, pady=5)
        _nucleotide_options = all_checkboxes['nucleotide']
        nucleotide_class = Checkbar(parent=nucleotide_frame, picks=[
                                    (i['formal'], i['default']) for i in _nucleotide_options])
        nucleotide_frame.pack(side=TOP, fill=X, expand=1, padx=10)

        # Translation Specific
        amino_frame = ttk.LabelFrame(output_fields_frame)
        amino_label = ttk.Label(amino_frame, font=('Arial', 20),
                                text="Translation Specific:")
        amino_label.pack(side=TOP, anchor=NW, padx=5, pady=5)
        _amino_options = all_checkboxes['amino']
        amino_class = Checkbar(parent=amino_frame, picks=[
                               (i['formal'], i['default']) for i in _amino_options])
        amino_frame.pack(side=TOP, expand=1, fill=X, padx=10)

    def _enter_output(self, entry):
        fo = None
        opts = {'title': "Select FASTA file to open...",
                'initialfile': entry.get(),
                'initialdir': self._user_directory}
        fo = filedialog.asksaveasfilename(**opts)
        if fo:
            entry.delete(0, END)
            entry.insert(END, fn)
            self.argument_dict['output_file'] = str(fn)

    def _create_readme(self, notebook_frame):
        readme_frame = ttk.Frame(notebook_frame, name="r_frame")
        # scroled_widget = ttk.ScrolledWindow(readme_frame,scrollbar='auto')
        # window = scroled.window
        vertical_scroll_frame = vsf(readme_frame)
        vertical_scroll_frame.pack(side=TOP, expand=1, fill=BOTH, anchor=NW)
        vsf_label = ttk.Label(
            vertical_scroll_frame.interior, text=open(self._directory_name + '/README_gui.txt').readlines(), anchor=N)
        # scroled_widget.pack(side=TOP, fill=BOTH, expand=1)
        vsf_label.pack(side=TOP, fill=BOTH, expand=1, anchor=NW)
        notebook_frame.add(readme_frame, text="Readme", underline=0, padding=2)

    def _update(self):
        import gui_execute
        gui_execute.main_refresh(self.root, gui_execute)


def main_refresh(root, gui_execute):
    reload(gui_execute)
    root.destroy()
    gui_execute.main()


def main():
    root = Tkinter.Tk()
    pyigblast_class = pyigblast_gui(root)
    root.mainloop()


class Checkbar():

    def __init__(self, parent=None, picks=[], side=LEFT, anchor=W):
        self.vars = {}
        for pick in picks:
            var = Tkinter.IntVar()
            var.set(int(pick[1]))
            chk = ttk.Checkbutton(parent, onvalue=1, offvalue=0, text=pick[
                                  0], variable=var, command=lambda: self.state())
            chk.pack(side=side, anchor=anchor, expand=YES)
            self.vars[pick[0]] = var

    def state(self):
        states_dict = {}
        for var in self.vars:
            states_dict[var] = self.vars[var].get()
        print states_dict
        return states_dict

if __name__ == '__main__':
    main()
