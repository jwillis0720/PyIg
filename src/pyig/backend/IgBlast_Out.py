import sys
import json
import gzip
#from cdr_analyzer import cdr_analyzer
import shelve
from Bio import SeqIO
import os.path
from collections import OrderedDict
import csv
import tempfile
from TranslateAndJoin import TranslateAndJoin
from pyig.backend.DefaultOrderedDict import DefaultOrderedDict


class SingleOutput_Entry():

    def __init__(self, entry, j_trans, species, debug=False):
        self.entry = entry
        self.j_trans = j_trans
        self.species = species
        self.debug = debug
        self.output = OrderedDict((('Sequence Id',  ""),
                                 ('Query Sequence', ""),
            ('Chain type', ""),
            ('Format Type', ""),
            ('Species', self.species),
            ('Top V Hit', ""),
            ('Top D Hit', ""),
            ('Top J Hit', ""),
            ('Productive', "False"),
            ('Productive CDR3', 'False'),
            ('Strand', ""),
            ('Framework 1 Nucleotides', ""),
            ('Framework 2 Nucleotides', ""),
            ('Framework 3 Nucleotides', ""),
            ('Framework 4 Nucleotides', ""),
            ('CDR1 Nucleotides', ""),
            ('CDR2 Nucleotides', ""),
            ('CDR3 Nucleotides', ""),
            ('Framework 1 AA', ""),
            ('Framework 2 AA', ""),
            ('Framework 3 AA', ""),
            ('Framework 4 AA', ""),
            ('Framework 1 AA Length', ""),
            ('Framework 2 AA Length', ""),
            ('Framework 3 AA Length', ""),
            ('Framework 4 AA Length', ""),
            ('CDR1 AA', ""),
            ('CDR2 AA', ""),
            ('CDR3 AA', ""),
            ('CDR1 AA Length', ""),
            ('CDR2 AA Length', ""),
            ('CDR3 AA Length', ""),
            ('Total V Alignment Matches', ""),
            ('Total V Alignment Mismatches', ""),
            ('Total V Alignment Length', ""),
            ('Total V Alignment Gaps', ""),
            ('Total V Alignment Identity', ""),
            ("FW1 Alignment From", ""),
            ("FW1 Alignment To", "",),
            ('FW1 Alignment Matches', ""),
            ('FW1 Alignment Mismatches', ""),
            ('FW1 Alignment Length', ""),
            ('FW1 Alignment Gaps', ""),
            ('FW1 Alignment Identity', ""),
            ("FW2 Alignment From", ""),
            ("FW2 Alignment To", ""),
            ('FW2 Alignment Matches', ""),
            ('FW2 Alignment Mismatches', ""),
            ('FW2 Alignment Length', ""),
            ('FW2 Alignment Gaps', ""),
            ('FW2 Alignment Identity', ""),
            ("FW3 Alignment From", ""),
            ("FW3 Alignment To", ""),
            ('FW3 Alignment Matches', ""),
            ('FW3 Alignment Mismatches', ""),
            ('FW3 Alignment Length', ""),
            ('FW3 Alignment Gaps', ""),
            ('FW3 Alignment Identity', ""),
            ('CDR1 Alignment From', ""),
            ('CDR1 Alignment To', ""),
            ('CDR1 Alignment Matches', ""),
            ('CDR1 Alignment Mismatches', ""),
            ('CDR1 Alignment Length', ""),
            ('CDR1 Alignment Gaps', ""),
            ('CDR1 Alignment Identity', ""),
            ('CDR2 Alignment From', ""),
            ('CDR2 Alignment To', ""),
            ('CDR2 Alignment Matches', ""),
            ('CDR2 Alignment Mismatches', ""),
            ('CDR2 Alignment Length', ""),
            ('CDR2 Alignment Gaps', ""),
            ('CDR2 Alignment Identity', ""),
            ("CDR3 Alignment From", ""),
            ("CDR3 Alignment To", ""),
            ('CDR3 Alignment Matches', ""),
            ('CDR3 Alignment Mismatches', ""),
            ('CDR3 Alignment Length', ""),
            ('CDR3 Alignment Gaps', ""),
            ('CDR3 Alignment Identity', ""),
            ('Junction V-End', ""),
            ('V-D Junction', ""),
            ('Junction D-Gene', ""),
            ('D-J Junction', ""),
            ('Junction J-Start', ""),
            ('D or J Junction', ""),
            ('Junction Merged', "")))

        # Title fields, these are the four sections it divides up to.
        # The fields are essentially the header to each section.
        # We could hard code this, but considering that it will be unique to each entry
        # We should parse it for each entry
        self.rearrangment_summary_titles = ""
        self.junction_detail_titles = ""
        self.alignment_summary_titles = ""
        self.hit_fields = ""

        # Here is the meat - all the good stuff goes into here
        self.rearrangment_summary = ""
        self.junction_detail = ""
        self.cdr1_alignment_summary = ""
        self.cdr2_alignment_summary = ""
        self.cdr3_alignment_summary = ""
        self.fr1_alignment_summary = ""
        self.fr2_alignment_summary = ""
        self.fr3_alignment_summary = ""
        self.total_alignment_summary = ""

        self.framework1_align = {}
        self.framework2_align = {}
        self.framework3_align = {}
        self.framework4_align = {}
        self.cdr1_align = {}
        self.cdr2_align = {}
        self.cdr3_align = {}
        self.total_v_align = {}

        # VD and J hits
        self.hits_v = []  # to be parsed in another function
        self.hits_d = []  # to be parsed in another function
        self.hits_j = []  # to be parsed in another function

    def get_json_entry(self):
        self.json = json.dumps(self.output, indent=4)
        return self.json

    def get_id(self):
        return self.output['Sequence Id']

    def set_seq(self, sequence):
        self.output['Query Sequence'] = str(sequence).upper()

    def parse(self):

        _rearrangment_breaker = False
        _junction_breaker = False
        _fields_breaker = False

        # Query Name
        for line in self.entry:
            if "Query" in line:
                self.output['Sequence Id'] = line.split(":")[1].strip()

            if "Domain classification requested:" in line:
                # imgt or kabat
                self.output['Format Type'] = line.split(":")[1].strip()

            if "rearrangement summary" in line:
                '''Okay we found rearrangement summary fields, that means
                the next line will be the actually rearrangment summary.
                So we turn the "breaker on" The rearrangment is just basic info
                returned like the top matches and if was productive etc'''
                self.rearrangment_summary_titles = line.strip().split("(")[2].split(")")[0].split(",")
                _rearrangment_breaker = True
                continue

            if _rearrangment_breaker:
                '''the meat of the rearranment, right after the field titles'''
                self.rearrangment_summary = line.strip().split("\t")
                _rearrangment_breaker = False

            if "junction details" in line:
                '''We found the junction detail line, get the fields,
                and the next line will be the actual junction details
                The junction details are the nucleotides of the Vend, VD,D,Jstart'''
                self.junction_detail_titles = line.strip().split(
                    "(")[2].split(")")[0].split(",")
                _junction_breaker = True
                continue

            if _junction_breaker:
                '''The meat of the junction details, right after the fields'''
                self.junction_detail = line.strip().split("\t")
                _junction_breaker = False

            if "Alignment summary" in line:
                '''Alignment summaries include the regions for CDR and framework,
                probably the most useful as it contains the areas in our sequence
                which match each regions CDR1,2,3 FW1,2,3 and the total (which is just the v portion)'''
                self.alignment_summary_titles = line.strip().split(
                    "(")[1].split(")")[0].split(",")

            '''Ok, in some blast versions, the line starts with just the field we are looking for,
            In other versions it has a '-', so I will just check for both'''
            # check for hypen
            try:
                if line.split()[0].split('-')[0] == "FR1":
                    self.fr1_alignment_summary = line.strip().split()[1:]
                if line.split()[0].split('-')[0] == "CDR1":
                    self.cdr1_alignment_summary = line.strip().split()[1:]
                if line.split()[0].split('-')[0] == "FR2":
                    self.fr2_alignment_summary = line.strip().split()[1:]
                if line.split()[0].split('-')[0] == "CDR2":
                    self.cdr2_alignment_summary = line.strip().split()[1:]
                if line.split()[0].split('-')[0] == "FR3":
                    self.fr3_alignment_summary = line.strip().split()[1:]
                if line.split()[0].split('-')[0] == "CDR3":
                    self.cdr3_alignment_summary = line.strip('\t').split()[1:]
                if line.split()[0].split('-')[0] == "Total":
                    self.total_alignment_summary = line.strip().split()[1:]
            except IndexError:
                pass

            # else the line will just start with the region of interest
            if line.startswith("FWR1"):
                self.fr1_alignment_summary = line.strip().split()[1:]
            if line.startswith("CDR1"):
                self.cdr1_alignment_summary = line.strip().split()[1:]
            if line.startswith("FWR2"):
                self.fr2_alignment_summary = line.strip().split()[1:]
            if line.startswith("CDR2"):
                self.cdr2_alignment_summary = line.strip().split()[1:]
            if line.startswith("FWR3"):
                self.fr3_alignment_summary = line.strip().split()[1:]
            if line.startswith("CDR3"):
                self.cdr3_alignment_summary = line.strip().split('\t')[1:]
            if line.startswith("Total"):
                self.total_alignment_summary = line.strip().split()[1:]

            # Finally parse VDJ hits
            if "# Fields:" in line:
                self.hit_titles = line.strip().split(":")[1].split(",")
                _fields_breaker = True
            if _fields_breaker:
                '''vdj hits have to be a list, since there can be more than one depending
                what the user asked for'''
                if line.startswith("V"):
                    self.hits_v.append(line)
                elif line.startswith("D"):
                    self.hits_d.append(line)
                elif line.startswith("J"):
                    self.hits_j.append(line)

        self.parse_rearranment()
        self.parse_junction()
        self.parse_fw1_align()
        self.parse_fw2_align()
        self.parse_fw3_align()
        self.parse_cdr1_align()
        self.parse_cdr2_align()
        self.parse_cdr3_align()
        self.parse_total_v_align()
        self.parse_v_hits()
        self.parse_d_hits()
        self.parse_j_hits()

    def parse_rearranment(self):
        '''parse the rearrangement summary (just the basic statistics)portion of the blast hit'''
        for title, value in zip(self.rearrangment_summary_titles, self.rearrangment_summary):

            # Retitle the IgBlast titles to what we like :)
            if title.strip() == "stop codon":
                title = "Stop Codon"
            if title.strip() == "Top V gene match":
                title = "Top V Hit"
            if title.strip() == "Top D gene match":
                title = "Top D Hit"
            if title.strip() == "Top J gene match":
                title = "Top J Hit"
            if title.strip() == "Top ":
                title = "Top J Hit"

            if len(value.split(',')) > 1:
                # sometimes there are multiple values for the rearranment values
                self.output[title.strip()] = tuple(value.split(','))
            else:
                self.output[title.strip()] = value

    def parse_junction(self):
        '''Return a dictionary of the junction that is the
        nuceotide sequence of teh CDR3 junctions'''

        self.junction_merged = ""
        for title, value in zip(self.junction_detail_titles, self.junction_detail):
            if title.strip() == "V end":
                title = "Junction V-End"
            if title.strip() == "V-D junction":
                title = "V-D Junction"
            if title.strip() == "D region":
                title = "Junction D-Gene"
            if title.strip() == "D-J junction":
                title = "D-J Junction"
            if title.strip() == "J start":
                title = "Junction J-Start"
            if "(" in value:
                self.output["D or J Junction"] = value.split(
                    "(")[1].split(")")[0]
                self.junction_merged += self.output["D or J Junction"]
            else:
                self.output[title.strip()] = value
                if value != "N/A":
                    self.junction_merged += value

            self.output['Junction Merged'] = self.junction_merged
            # right now we have junction merged - which is an incomplete HCDR3 because
            # they only give you half the J region...We will fix that shortly

    # Now lets start parsing the alignment summaries into each sections FW1,2,3,4,CDR1,2,3
    def parse_fw1_align(self):
        for title, value in zip(self.alignment_summary_titles, self.fr1_alignment_summary):
            if title.strip() == 'percent identity':
                title = 'Identity'
            try:
                self.output["FW1 Alignment " + title.strip().capitalize()] = float(value)
            except ValueError:
                self.output["FW1 Alignment " + title.strip().capitalize()] = value

    def parse_fw2_align(self):
        for title, value in zip(self.alignment_summary_titles, self.fr2_alignment_summary):
            if title.strip() == 'percent identity':
                title = 'Identity'
            try:
                self.output["FW2 Alignment " + title.strip().capitalize()] = float(value)
            except ValueError:
                self.output["FW2 Alignment " + title.strip().capitalize()] = value

    def parse_fw3_align(self):
        for title, value in zip(self.alignment_summary_titles, self.fr3_alignment_summary):
            if title.strip() == 'percent identity':
                title = 'Identity'
            try:
                self.output["FW3 Alignment " + title.strip().capitalize()] = float(value)
            except ValueError:
                self.output["FW3 Alignment " + title.strip().capitalize()] = value

    def parse_cdr1_align(self):
        for title, value in zip(self.alignment_summary_titles, self.cdr1_alignment_summary):
            if title.strip() == 'percent identity':
                title = 'Identity'
            try:
                self.output['CDR1 Alignment ' + title.strip().capitalize()] = float(value)
            except ValueError:
                self.output['CDR1 Alignment ' + title.strip().capitalize()] = value

    def parse_cdr2_align(self):
        for title, value in zip(self.alignment_summary_titles, self.cdr2_alignment_summary):
            if title.strip() == 'percent identity':
                title = 'Identity'
            try:
                self.output['CDR2 Alignment ' + title.strip().capitalize()] = float(value)
            except ValueError:
                self.output['CDR2 Alignment ' + title.strip().capitalize()] = value

    def parse_cdr3_align(self):
        for title, value in zip(self.alignment_summary_titles, self.cdr3_alignment_summary):
            if title.strip() == 'percent identity':
                title = 'Identity'
            try:
                self.output['CDR3 Alignment ' + title.strip().capitalize()] = float(value)
            except ValueError:
                self.output['CDR3 Alignment ' + title.strip().capitalize()] = value

    def parse_total_v_align(self):
        '''The total alignment is just the Whole V region region which is pretty nice'''
        for title, value in zip(self.alignment_summary_titles, self.total_alignment_summary):
            if title.strip() == 'percent identity':
                title = 'Identity'
            if title.strip() == 'from' or title.strip() == 'to':
                continue
            try:
                self.output["Total V Alignment " + title.strip().capitalize()] = float(value)
            except ValueError:
                self.output["Total V Alignment " + title.strip().capitalize()] = value

    def parse_v_hits(self):
        for rank, entry in enumerate(self.hits_v, start=1):
            _entry_dict = OrderedDict()
            for value, title in zip(entry.split()[1:], self.hit_titles):
                # sometimes there is nothing there so it can't cast to a float
                try:
                    _entry_dict['V-Gene Rank_' + str(rank) + " " + title.strip(
                    ).capitalize().replace("%", "Percent")] = float(value)
                except ValueError:
                    _entry_dict['V-Gene Rank_' + str(
                        rank) + " " + title.strip().capitalize().replace("%", "Percent")] = value
            self.output['V-Gene Rank_' + str(rank)] = _entry_dict

    def parse_d_hits(self):
        for rank, entry in enumerate(self.hits_d, start=1):
            _entry_dict = OrderedDict()
            for value, title in zip(entry.split()[1:], self.hit_titles):
                # sometimes there is nothing there so it can't cast to a float
                try:
                    _entry_dict['D-Gene Rank_' + str(rank) + " " + title.strip(
                    ).capitalize().replace("%", "Percent")] = float(value)
                except ValueError:
                    _entry_dict['D-Gene Rank_' + str(
                        rank) + " " + title.strip().capitalize().replace("%", "Percent")] = value
            self.output['D-Gene Rank_' + str(rank)] = _entry_dict

    def parse_j_hits(self):
        for rank, entry in enumerate(self.hits_j, start=1):
            _entry_dict = OrderedDict()
            for value, title in zip(entry.split()[1:], self.hit_titles):
                # sometimes there is nothing there so it can't cast to a float
                try:
                    _entry_dict['J-Gene Rank_' + str(rank) + " " + title.strip(
                    ).capitalize().replace("%", "Percent")] = float(value)
                except ValueError:
                    _entry_dict['J-Gene Rank_' + str(
                        rank) + " " + title.strip().capitalize().replace("%", "Percent")] = value
            self.output['J-Gene Rank_' + str(rank)] = _entry_dict

    def join_and_translate(self):
        TranslateAndJoin(self)
        # Find out if the join is productive:
        if "*" not in self.output['CDR3 AA']:
            self.output['Productive CDR3'] = "True"


class IgBlast_Out():

    def __init__(self, debug=False):
        self.blast_output_handle = ""
        self.parsed_output = tempfile.NamedTemporaryFile(suffix=".json", delete=False).name
        self.debug = debug

    def set_seq_dictionary(self, seq_dictionary):
        self.seq_dictionary = seq_dictionary

    def set_blast_output(self, blast_handle):
        self.blast_output_handle = blast_handle

    def set_input_query(self, query):
        self.input_query = {}
        for entry in SeqIO.parse(query, 'fasta'):
            self.input_query[entry.id] = entry.seq

    def set_species(self, species):
        self.species = species

    def get_output_name(self):
        return self.parsed_output

    def parse(self):
        _focus_lines = []
        j_trans = {}
        j_lines = open(os.path.join(os.environ['IGDATA'], "germ_props", self.species, "properties.txt")).readlines()
        for line in j_lines:
            j_trans[line.split()[0]] = line.split()[1]
        with open(self.parsed_output, 'w') as out:
            for line in open(self.blast_output_handle):
                if "IGBLASTN" in line:
                    if _focus_lines:
                        Single_Blast_Entry = SingleOutput_Entry(_focus_lines, j_trans, self.species, debug=self.debug)
                        Single_Blast_Entry.parse()
                        Single_Blast_Entry.get_id()
                        Single_Blast_Entry.set_seq(
                            self.seq_dictionary[Single_Blast_Entry.get_id()])
                        Single_Blast_Entry.join_and_translate()
                        out.write(Single_Blast_Entry.get_json_entry())
                        out.write("\n")
                        _focus_lines = []
                    else:
                        continue
                else:
                    _focus_lines.append(line)
            # for last line -
            if _focus_lines:
                Single_Blast_Entry = SingleOutput_Entry(_focus_lines, j_trans, self.species, debug=self.debug)
                Single_Blast_Entry.parse()
                Single_Blast_Entry.get_id()
                Single_Blast_Entry.set_seq(self.seq_dictionary[Single_Blast_Entry.get_id()])
                Single_Blast_Entry.join_and_translate()
                out.write(Single_Blast_Entry.get_json_entry())
                out.write("\n")