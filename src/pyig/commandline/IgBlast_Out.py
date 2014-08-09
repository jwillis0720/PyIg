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
from pyig.backend.DefaultOrderedDict import DefaultOrderedDict


class SingleOutput_Entry():

    def __init__(self,entry,species):
        self.entry = entry 
        self.species = species
        self.output =  OrderedDict((('Sequence Id' ,  ""),
                        ('Query Sequence' , ""),
                        ('Chain type' , ""),
                        ('Format Type', ""),
                        ('Top V Hit', ""),
                        ('Top D Hit', ""),
                        ('Top J Hit', ""),
                        ('Productive', "False"),
                        ('Productive CDR3', 'False'),
                        ('Strand', ""),
                        ('Framework 1 Nucleotides',""),
                        ('Framework 2 Nucleotides',""),
                        ('Framework 3 Nucleotides',""),
                        ('Framework 4 Nucleotides',""),
                        ('CDR1 Nucleotides', ""),
                        ('CDR2 Nucleotides', ""),
                        ('CDR3 Nucleotides', ""),
                        ('Framework 1 AA' , ""),
                        ('Framework 2 AA' , ""),
                        ('Framework 3 AA' , ""),
                        ('Framework 4 AA' , ""),
                        ('Framework 1 AA Length', ""),
                        ('Framework 2 AA Length', ""),
                        ('Framework 3 AA Length', ""),
                        ('Framework 4 AA Length', ""),
                        ('CDR 1 AA' , ""),
                        ('CDR 2 AA' , ""),
                        ('CDR 3 AA' , ""),
                        ('CDR 1 AA Length' , ""),
                        ('CDR 2 AA Length' , ""),
                        ('CDR 3 AA Length' , ""),
                        ('Total V Alignment Matches' , ""),
                        ('Total V Alignment Mismatches' , ""),
                        ('Total V Alignment Length' , ""),
                        ('Total V Alignment Gaps' , ""),
                        ('Total V Alignment Identity' , ""),
                        ('FW1 Alignment Matches' , ""),
                        ('FW1 Alignment Mismatches' , ""),
                        ('FW1 Alignment Length' , ""),
                        ('FW1 Alignment Gaps' , ""),
                        ('FW1 Alignment Identity' , ""),
                        ('FW2 Alignment Matches' , ""),
                        ('FW2 Alignment Mismatches' , ""),
                        ('FW2 Alignment Length' , ""),
                        ('FW2 Alignment Gaps' , ""),
                        ('FW2 Alignment Identity' , ""),
                        ('FW3 Alignment Matches' , ""),
                        ('FW3 Alignment Mismatches' , ""),
                        ('FW3 Alignment Length' , ""),
                        ('FW3 Alignment Gaps' , ""),
                        ('FW3 Alignment Identity' , ""),
                        ('CDR1 Alignment Matches' , ""),
                        ('CDR1 Alignment Mismatches' , ""),
                        ('CDR1 Alignment Length' , ""),
                        ('CDR1 Alignment Gaps' , ""),
                        ('CDR1 Alignment Identity' , ""),
                        ('CDR2 Alignment Matches' , ""),
                        ('CDR2 Alignment Mismatches' , ""),
                        ('CDR2 Alignment Length' , ""),
                        ('CDR2 Alignment Gaps' , ""),
                        ('CDR2 Alignment Idenity' , ""),
                        ('CDR3 Alignment Matches' , ""),
                        ('CDR3 Alignment Mismatches' , ""),
                        ('CDR3 Alignment Length' , ""),
                        ('CDR3 Alignment Gaps' , ""),
                        ('CDR3 Alignment Identity' , ""),
                        # ('Top V gene match' , ""),
                        # ('V-Gene Mismatches' , ""),
                        # ('V-Gene Percent Identity' , ""),
                        # ('V-Gene Gaps' , ""),
                        # ('V-Gene e-Value' , ""),
                        # ('V-Gene Bit Score' , ""),
                        # ('V-Gene Alignment Length' , ""),
                        # ('Top D gene match' , ""),
                        # ('D-Gene Mismatches' , ""),
                        # ('D-Gene Percent Identity' , ""),
                        # ('D-Gene Gaps' , ""),
                        # ('D-Gene e-Value' , ""),
                        # ('D-Gene Bit Score' , ""),
                        # ('D-Gene Alignment Length' , ""),
                        # ('Top J gene match' , ""),
                        # ('J-Gene Mismatches' , ""),
                        # ('J-Gene Percent Identity' , ""),
                        # ('J-Gene Gaps' , ""),
                        # ('J-Gene e-Value' , ""),
                        # ('J-Gene Bit Score' , ""),
                        # ('J-Gene Alignment Length' , ""),
                        ('Junction V-End',""),
                        ('V-D Junction',""),    
                        ('Junction D-Gene',""),
                        ('D-J Junction',""),
                        ('Junction J-Start',""),
                        ('D or J Junction',"")))

         # Where to end translation for CDR3 loops specified in an output file
        self.end_cdr3_dicts = os.path.join(os.environ['IGDATA'],"germ_props",self.species,"proerties.txt")

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
        self.json = json.dumps(self.output,indent=4)
        return self.json
    
    def get_id(self):
        return self.output['Sequence Id']

    def set_seq(self, sequence):
        self.output['Query Sequence'] = str(sequence).upper()

    def parse(self):

        _rearrangment_breaker = False
        _junction_breaker = False
        _fields_breaker = False 
       
       #Query Name 
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

        
            #Finally parse VDJ hits
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


    def parse_rearranment(self):
        '''parse the rearrangement summary (just the basic statistics)portion of the blast hit'''
        for title, value in zip(self.rearrangment_summary_titles, self.rearrangment_summary):
            
            #Retitle the IgBlast titles to what we like :)
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
                #sometimes there are multiple values for the rearranment values
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

            #right now we have junction merged - which is an incomplete HCDR3 because 
            #they only give you half the J region...We will fix that shortly


    #Now lets start parsing the alignment summaries into each sections FW1,2,3,4,CDR1,2,3
    def parse_fw1_align(self):
        for title, value in zip(self.alignment_summary_titles, self.fr1_alignment_summary):
            if title.strip() == 'percent identity':
                title = 'Identity'
            self.output["FW1 Alignment "+title.strip().capitalize()] = value
            
    def parse_fw2_align(self):
        for title, value in zip(self.alignment_summary_titles, self.fr2_alignment_summary):
            if title.strip() == 'percent identity':
                title = 'Identity'
            self.output["FW2 Alignment "+title.strip().capitalize()] = value   

    def parse_fw3_align(self):
        for title, value in zip(self.alignment_summary_titles, self.fr3_alignment_summary):
            if title.strip() == 'percent identity':
                title = 'Identity'
            self.output["FW3 Alignment "+title.strip().capitalize()] = value   

    def parse_cdr1_align(self):
        for title, value in zip(self.alignment_summary_titles, self.cdr1_alignment_summary):
            if title.strip() == 'percent identity':
                title = 'Identity'
            self.output['CDR1 Alignment '+ title.strip().capitalize()] = value

    def parse_cdr2_align(self):
        for title, value in zip(self.alignment_summary_titles, self.cdr2_alignment_summary):
            if title.strip() == 'percent identity':
                title = 'Identity'
            self.output['CDR2 Alignment '+ title.strip().capitalize()] = value

    def parse_cdr3_align(self):
        for title, value in zip(self.alignment_summary_titles, self.cdr3_alignment_summary):
            if title.strip() == 'percent identity':
                title = 'Identity'
            self.output['CDR3 Alignment '+ title.strip().capitalize()] = value

    def parse_total_v_align(self):
        '''The total alignment is just the Whole V region region which is pretty nice'''
        for title, value in zip(self.alignment_summary_titles, self.total_alignment_summary):
            if title.strip() == 'percent identity':
                title = 'Identity'
            if title.strip() == 'from' or title.strip() == 'to':
                continue
            self.output["Total V Alignment "+ title.strip().capitalize()] = value


    def parse_v_hits(self):       
        for rank,entry in enumerate(self.hits_v,start=1):
            _entry_dict = OrderedDict()
            for value, title in zip(entry.split()[1:], self.hit_titles):
                # sometimes there is nothing there so it can't cast to a float
                try:
                    _entry_dict['V-Gene Rank_'+ str(rank) +" "+title.strip().capitalize().replace("%", "Percent")] = float(value)
                except ValueError:
                    _entry_dict['V-Gene Rank_'+ str(rank) +" "+title.strip().capitalize().replace("%", "Percent")] = value 
            self.output['V-Gene Rank_' + str(rank)] = _entry_dict


#     def generate_blast_dict(self):
#         '''THE blast dict is a dictionary of dictionaries containing all of
#         our information about the hit, this will be filled by several methods,
#         each of which return their own dictionaries'''
#         blast_dict = {}
#         blast_dict[self.query] = {
#             "domain_classification": self.domain_classification,
#             "raw_seq": self.get_raw_seq(),
#             "rearrangement": self.parse_rearranment(),
#             "junction": self.parse_junction(),
#             "fr1_align": self.parse_fr1_align(),
#             "fr2_align": self.parse_fr2_align(),
#             "fr3_align": self.parse_fr3_align(),
#             "cdr1_align": self.parse_cdr1_align(),
#             "cdr2_align": self.parse_cdr2_align(),
#             "cdr3_align": self.parse_cdr3_align(),
#             "total_align": self.parse_total_align(),
#             "v_hits": self.parse_v_hits(),
#             "d_hits": self.parse_d_hits(),
#             "j_hits": self.parse_j_hits()
#         }

#         # add the analysis
#         blast_dict = cdr_analyzer(blast_dict, self.end_cdr3_dicts).return_modified_dict()
#         return blast_dict

#     def get_raw_seq(self):
#         return self.raw_seqs_handle[self.query].upper()

#     def parse_junction(self):
#         '''Return a dictionary of the junction that is the
#         nuceotide sequence of teh CDR3 junctions'''
#         _return_dict = {}
#         self.junction_together = ""
#         for title, value in zip(self.junction_detail_titles, self.junction_detail):
#             if "(" in value:
#                 _return_dict["d-or-j_junction"] = value.split(
#                     "(")[1].split(")")[0]
#                 self.junction_together += value.split(
#                     "(")[1].split(")")[0]
#             else:
#                 _return_dict[title.strip().lower().replace(" ", "_")] = value
#                 if value != "N/A":
#                     self.junction_together += value
#         return _return_dict

#     # Now lets start parsing the alignment summaries into each sections FW1,2,3,4,CDR1,2,3
#     def parse_fr1_align(self):
#         _return_dict = {}
#         for title, value in zip(self.alignment_summary_titles, self.fr1_alignment_summary):
#             try:
#                 _return_dict[title] = int(value)
#             except ValueError:
#                 try:
#                     _return_dict[title] = float(value)
#                 except ValueError:
#                     _return_dict[title] = value
#         return _return_dict

#     def parse_fr2_align(self):
#         _return_dict = {}
#         for title, value in zip(self.alignment_summary_titles, self.fr2_alignment_summary):
#             try:
#                 _return_dict[title] = int(value)
#             except ValueError:
#                 try:
#                     _return_dict[title] = float(value)
#                 except ValueError:
#                     _return_dict[title] = value
#         return _return_dict

#     def parse_fr3_align(self):
#         _return_dict = {}
#         for title, value in zip(self.alignment_summary_titles, self.fr3_alignment_summary):
#             try:
#                 _return_dict[title] = int(value)
#             except ValueError:
#                 try:
#                     _return_dict[title] = float(value)
#                 except ValueError:
#                     _return_dict[title] = value
#         return _return_dict

#     def parse_cdr1_align(self):
#         _return_dict = {}
#         for title, value in zip(self.alignment_summary_titles, self.cdr1_alignment_summary):
#             try:
#                 _return_dict[title] = int(value)
#             except ValueError:
#                 try:
#                     _return_dict[title] = float(value)
#                 except ValueError:
#                     _return_dict[title] = value
#         return _return_dict

#     def parse_cdr2_align(self):
#         _return_dict = {}
#         for title, value in zip(self.alignment_summary_titles, self.cdr2_alignment_summary):
#             try:
#                 _return_dict[title] = int(value)
#             except ValueError:
#                 try:
#                     _return_dict[title] = float(value)
#                 except ValueError:
#                     _return_dict[title] = value
#         return _return_dict

#     def parse_cdr3_align(self):
#         _return_dict = {}
#         for title, value in zip(self.alignment_summary_titles, self.cdr3_alignment_summary):
#             try:
#                 _return_dict[title] = int(value)
#             except ValueError:
#                 try:
#                     _return_dict[title] = float(value)
#                 except ValueError:
#                     _return_dict[title] = value
#         return _return_dict

#     def parse_total_align(self):
#         '''The total alignment is just the Whole region which is pretty nice'''
#         _return_dict = {}
#         for title, value in zip(self.alignment_summary_titles, self.total_alignment_summary):
#             try:
#                 _return_dict[title] = int(value)
#             except ValueError:
#                 try:
#                     _return_dict[title] = float(value)
#                 except ValueError:
#                     _return_dict[title] = value
#         return _return_dict

#     def parse_v_hits(self):
#         _return_dict = {}
#         rank = 1
#         for entry in self.hits_v:
#             _entry_dict = {}
#             for value, title in zip(entry.split()[1:], self.hit_fields):
#                 # sometimes there is nothing there so it can't cast to a float
#                 try:
#                     _entry_dict[title.strip().replace(' ', '_').replace("%", "percent")] = float(value)
#                 except:
#                     _entry_dict[title.strip().replace(' ', '_').replace("%", "percent")] = value
#             _return_dict['rank_' + str(rank)] = _entry_dict
#             rank += 1
#         # return a v_hit dictionary that is unsorted, but we kept the rank
#         # by how it appeared in the balst out
#         return _return_dict

#     def parse_d_hits(self):
#         _return_dict = {}
#         rank = 1
#         for entry in self.hits_d:
#             _entry_dict = {}
#             #_entry_dict["rank"] = int(rank)
#             for value, title in zip(entry.split()[1:], self.hit_fields):
#                 try:
#                     _entry_dict[title.strip().replace(' ', '_').replace("%", "percent")] = float(value)
#                 except:
#                     _entry_dict[title.strip().replace(' ', '_').replace("%", "percent")] = value
#             _return_dict['rank_' + str(rank)] = _entry_dict
#             rank += 1
#         return _return_dict

#     def parse_j_hits(self):
#         _return_dict = {}
#         rank = 1
#         for entry in self.hits_j:
#             _entry_dict = {}
#             for value, title in zip(entry.split()[1:], self.hit_fields):
#                 try:
#                     _entry_dict[title.strip().replace(' ', '_').replace("%", "percent")] = float(value)
#                 except:
#                     _entry_dict[title.strip().replace(' ', '_').replace("%", "percent")] = value
#             _return_dict['rank_' + str(rank)] = _entry_dict
#             rank += 1
#         return _return_dict



class IgBlast_Out():

    '''Pass the whole output to this class and it will deal with it
    Args:
    blast file - from output format 7 of blastign
    general_options - The options you want to output for general options passed as a dictionary
    nuc_options - The nucleotide options you want passed as a dictionary
    aa_options - The
    output - string for the filename to the output
    '''

    def __init__(self):
        self.blast_output_handle = ""
        self.parsed_output = tempfile.NamedTemporaryFile(suffix=".json",delete=False).name

    def set_seq_dictionary(self, seq_dictionary):
        self.seq_dictionary = seq_dictionary

    def set_blast_output(self,blast_handle):
        self.blast_output_handle = blast_handle

    def set_input_query(self,query):
        self.input_query = {}
        for entry in SeqIO.parse(query, 'fasta'):
            self.input_query[entry.id] = entry.seq

    def get_output_name(self):
        return self.parsed_output

    def parse(self):
        _focus_lines = []
        with open(self.parsed_output, 'w') as out:
            for line in open(self.blast_output_handle):
                if "IGBLASTN" in line:
                    if _focus_lines:
                            Single_Blast_Entry = SingleOutput_Entry(_focus_lines,'human')
                            Single_Blast_Entry.parse()
                            Single_Blast_Entry.get_id()
                            Single_Blast_Entry.set_seq(
                                self.seq_dictionary[Single_Blast_Entry.get_id()])
                            out.write(Single_Blast_Entry.get_json_entry())
                            out.write("\n")
                            _focus_lines = []                    
                    else:
                         continue
                else:
                     _focus_lines.append(line)
            #for last line -
            if _focus_lines:
                Single_Blast_Entry = SingleOutput_Entry(_focus_lines,'human')
                Single_Blast_Entry.parse()
                Single_Blast_Entry.get_id()
                Single_Blast_Entry.set_seq(self.seq_dictionary[Single_Blast_Entry.get_id()])
                out.write(Single_Blast_Entry.get_json_entry())
                out.write("\n")
        

#         # The blast file that was output, format it to a handle
#         self.blast_file_handle = ""
        
#         # the file to use that will tell us where the CDR3 ends, pass this to the CDR analyzer
#         #self.end_dict = {}
#         # these are dictionaries with the options we want to output. See output_tabs checkbozes

#         # bool to tell us if we want the file zipped
#         #self.zip = kwargs['zip_bool']
#         #self.option_list = options
#         # asks if we are parsing to the gui or not
#         #self.gui = kwargs['gui']
#         #self.species = kwargs['species']
#         # where you want to store the output before it gets concat
#         #self.temporary_directory = temporary_directory

#         # the input file given to this blast instance, used to get the raw seq
#         #self.input_file = input_file

#         # get back the db handle to the raw seqs
#         #self.raw_seqs_db_handle = self.get_raw_seqs_db(os.path.basename(blast_file).split('.')[0])

#         #self.germ_file = kwargs['germ_properties']

#         #try:
#          #   for line in open(self.germ_file).readlines():
#          #       line_split = line.split()
#          #       self.end_dict[line_split[0]] = int(line_split[1])
#         #except IOError:
#          #   print "Cant open germ_properties.txt.\nNeed this file to process CDR3 regions"
#          #   sys.exit()

#     def parse_blast_file_to_type(self, output_file, o_type="json"):
#         self.header_written = False
#         focus_lines = []
#         output_file = output_file.split('.')[0] + "." + o_type
#         if self.zip:
#             json_output_file_handle = gzip.open(output_file + ".gz", 'wb')
#         else:
#             json_output_file_handle = open(output_file, 'w')

#         with json_output_file_handle as openfile:
#             for line in self.blast_file_handle:
#                 if "IGBLASTN" in line:
#                     if focus_lines:
#                         sbe = single_blast_entry(focus_lines, self.raw_seqs_db_handle, self.end_dict)
#                         blast_dictionary = sbe.generate_blast_dict()  # return single blast entry
#                         if o_type == "json":
#                             json_document_trimmed = trim_json(blast_dictionary, self.option_list, gui=self.gui)  # trim the json according to input
#                             openfile.write(json.dumps(json_document_trimmed))
#                             openfile.write("\n")
#                               # write it out
#                             focus_lines = []
#                         if o_type == "csv":
#                             csv_document_trimmed = trim_json(blast_dictionary, self.option_list, gui=self.gui)  # trim the json according to input
#                             dw = csv.DictWriter(openfile, fieldnames=csv_document_trimmed.keys(), delimiter=",", quotechar="|", quoting=csv.QUOTE_MINIMAL)
#                             if not self.header_written:
#                                 dw.writeheader()
#                                 self.header_written = True
#                             dw.writerow(csv_document_trimmed)
#                             focus_lines = []
#                     else:
#                         continue
#                 else:
#                     focus_lines.append(line)
#             # and do it for the end too

#             sbe = single_blast_entry(focus_lines, self.raw_seqs_db_handle, self.end_dict)
#             blast_dictionary = sbe.generate_blast_dict()  # return single blast entry
#             if o_type == "json":
#                 json_document_trimmed = trim_json(blast_dictionary, self.option_list, gui=self.gui)  # trim the json according to input
#                 openfile.write(json.dumps(json_document_trimmed))
#                 openfile.write("\n")  # write it out
#             elif o_type == "csv":
#                 csv_document_trimmed = trim_json(blast_dictionary, self.option_list, gui=self.gui)
#                 dw = csv.DictWriter(openfile, fieldnames=csv_document_trimmed.keys(), delimiter=",", quotechar="|", quoting=csv.QUOTE_MINIMAL)
#                 if not self.header_written:
#                     dw.writeheader()
#                 dw.writerow(csv_document_trimmed)

#         self.raw_seqs_db_handle.close()
#         self.blast_file_handle.close()

#     def get_raw_seqs_db(self, name):
#         self.db_name = os.path.join(self.temporary_directory, name)
#         shelf = shelve.open(self.db_name)
#         for sequence_object in SeqIO.parse(self.input_file, 'fasta'):
#             shelf[str(sequence_object.id)] = str(sequence_object.seq)
#         return shelf


# class single_blast_entry():

#     '''The helper class to parse an individual blast result'''

#     def __init__(self, single_entry, raw_seqs_handle, end_cdr3_dicts):
#         # some breaker options that will tell us when to go on to the next section
#         _rearrangment_breaker = False
#         _junction_breaker = False
#         _fields_breaker = False

#         # the database handle
#         self.raw_seqs_handle = raw_seqs_handle

#         # Where to end translation for CDR3 loops specified in an output file
#         self.end_cdr3_dicts = end_cdr3_dicts

#         # initalize the fields, some can be empty on crappy reads
#         # basic
#         self.query = ""
#         self.full_query_seq = ""
#         self.domain_classification = ""

#         # title fields, these are the four sections it divides up to.
#         # The fields are essentially the header to each section.
#         # We could hard code this, but considering that it will be unique to each entry
#         # We should parse it for each entry
#         self.rearrangment_summary_titles = ""
#         self.alignment_summary_titles = ""
#         self.junction_detail_titles = ""
#         self.hit_fields = ""

#         # Here is the meat.
#         self.rearrangment_summary = ""
#         self.junction_detail = ""
#         self.cdr1_alignment_summary = ""
#         self.cdr2_alignment_summary = ""
#         self.cdr3_alignment_summary = ""
#         self.fr1_alignment_summary = ""
#         self.fr2_alignment_summary = ""
#         self.fr3_alignment_summary = ""
#         self.total_alignment_summary = ""

#         # hits
#         self.hits_v = []  # to be parsed in another function
#         self.hits_d = []  # to be parsed in another function
#         self.hits_j = []  # to be parsed in another function

#         for line in single_entry:
#             if "Query" in line:
#                 # this will be the id field which blast calls the quey
#                 self.query = line.split(":")[1].strip()

#             if "Domain classification requested:" in line:
#                 # imgt or kabat
#                 self.domain_classification = line.split(":")[1].strip()

#             if "rearrangement summary" in line:
#                 '''Okay we found rearrangement summary fields, that means
#                 the next line will be the actually rearrangment summary.
#                 So we turn the "breaker on" The rearrangment is just basic info
#                 returned like the top matches and if was productive etc'''
#                 self.rearrangment_summary_titles = line.strip().split(
#                     "(")[2].split(")")[0].split(",")
#                 self.rearrangment_summary_titles = [x.strip().lower().replace(" ", "_") for x in self.rearrangment_summary_titles]
#                 _rearrangment_breaker = True
#                 continue

#             if _rearrangment_breaker:
#                 # the meat of the rearranment, right after the fields
#                 self.rearrangment_summary = line.strip().split("\t")
#                 _rearrangment_breaker = False

#             if "junction details" in line:
#                 '''We found the junction detail line, get the fields,
#                 and the next line will be the actual junction details
#                 The junction details are the nucleotides of the Vend, VD,D,Jstart'''
#                 self.junction_detail_titles = line.strip().split(
#                     "(")[2].split(")")[0].split(",")
#                 _junction_breaker = True
#                 continue

#             if _junction_breaker:
#                 '''The meat of the junction details, right after the fields'''
#                 self.junction_detail = line.strip().split("\t")
#                 _junction_breaker = False

#             if "Alignment summary" in line:
#                 '''Alignment summaries include the regions for CDR and framework,
#                 probably the most useful as it contains the areas in our sequence
#                 which match each regions CDR1,2,3 FW1,2,3 and the total (which is just the v portion)'''
#                 self.alignment_summary_titles = line.strip().split(
#                     "(")[1].split(")")[0].split(",")
#                 # strip off white marks in the fields
#                 self.alignment_summary_titles = [x.strip().replace(" ", "_") for x in self.alignment_summary_titles]
#             '''Ok, in some blast versions, the line starts with just the field we are looking for,
#             In other versions it has a '-', so I will just check for both'''

#             # check for hypen
#             try:
#                 if line.split()[0].split('-')[0] == "FR1":
#                     self.fr1_alignment_summary = line.strip().split()[1:]
#                 if line.split()[0].split('-')[0] == "CDR1":
#                     self.cdr1_alignment_summary = line.strip().split()[1:]
#                 if line.split()[0].split('-')[0] == "FR2":
#                     self.fr2_alignment_summary = line.strip().split()[1:]
#                 if line.split()[0].split('-')[0] == "CDR2":
#                     self.cdr2_alignment_summary = line.strip().split()[1:]
#                 if line.split()[0].split('-')[0] == "FR3":
#                     self.fr3_alignment_summary = line.strip().split()[1:]
#                 if line.split()[0].split('-')[0] == "CDR3":
#                     self.cdr3_alignment_summary = line.strip('\t').split()[1:]
#                 if line.split()[0].split('-')[0] == "Total":
#                     self.total_alignment_summary = line.strip().split()[1:]
#             except IndexError:
#                 pass

#             # else the line will just start with the region of interest
#             if line.startswith("FWR1"):
#                 self.fr1_alignment_summary = line.strip().split()[1:]
#             if line.startswith("CDR1"):
#                 self.cdr1_alignment_summary = line.strip().split()[1:]
#             if line.startswith("FWR2"):
#                 self.fr2_alignment_summary = line.strip().split()[1:]
#             if line.startswith("CDR2"):
#                 self.cdr2_alignment_summary = line.strip().split()[1:]
#             if line.startswith("FWR3"):
#                 self.fr3_alignment_summary = line.strip().split()[1:]
#             if line.startswith("CDR3"):
#                 self.cdr3_alignment_summary = line.strip().split('\t')[1:]
#             if line.startswith("Total"):
#                 self.total_alignment_summary = line.strip().split()[1:]

#             '''Last but not least, we can pass all the VDJ hits along with all the
#             info about them that blast gives, the VDJ hits can give a ton of information
#             about the junction, the identity to each gene, the frame, gaps, evalue, etc,
#             these are ranked and put into a list'''
#             if "# Fields:" in line:
#                 self.hit_fields = line.strip().split(":")[1].split(",")
#                 _fields_breaker = True
#             if _fields_breaker:
#                 '''vdj hits have to be a list, since there can be more than one depending
#                 what the user asked for'''
#                 if line.startswith("V"):
#                     self.hits_v.append(line)
#                 elif line.startswith("D"):
#                     self.hits_d.append(line)
#                 elif line.startswith("J"):
#                     self.hits_j.append(line)
#         '''now that I got all the blast output parsed,
#         lets put it into the requested format, before we do that however lets make
#         one dictionary'''

#     def generate_blast_dict(self):
#         '''THE blast dict is a dictionary of dictionaries containing all of
#         our information about the hit, this will be filled by several methods,
#         each of which return their own dictionaries'''
#         blast_dict = {}
#         blast_dict[self.query] = {
#             "domain_classification": self.domain_classification,
#             "raw_seq": self.get_raw_seq(),
#             "rearrangement": self.parse_rearranment(),
#             "junction": self.parse_junction(),
#             "fr1_align": self.parse_fr1_align(),
#             "fr2_align": self.parse_fr2_align(),
#             "fr3_align": self.parse_fr3_align(),
#             "cdr1_align": self.parse_cdr1_align(),
#             "cdr2_align": self.parse_cdr2_align(),
#             "cdr3_align": self.parse_cdr3_align(),
#             "total_align": self.parse_total_align(),
#             "v_hits": self.parse_v_hits(),
#             "d_hits": self.parse_d_hits(),
#             "j_hits": self.parse_j_hits()
#         }

#         # add the analysis
#         blast_dict = cdr_analyzer(blast_dict, self.end_cdr3_dicts).return_modified_dict()
#         return blast_dict

#     def get_raw_seq(self):
#         return self.raw_seqs_handle[self.query].upper()

#     def parse_rearranment(self):
#         '''parse the rearrangement summary (just the basic statistics)
#         portion of the blast hit'''
#         _return_dict = {}
#         for title, value in zip(self.rearrangment_summary_titles, self.rearrangment_summary):
#             if len(value.split(',')) > 1:
#                 '''sometimes there are multiple values for the rearranment values'''
#                 _return_dict[title.strip()] = tuple(value.split(','))
#             else:
#                 _return_dict[title.strip()] = value
#         return _return_dict

#     def parse_junction(self):
#         '''Return a dictionary of the junction that is the
#         nuceotide sequence of teh CDR3 junctions'''
#         _return_dict = {}
#         self.junction_together = ""
#         for title, value in zip(self.junction_detail_titles, self.junction_detail):
#             if "(" in value:
#                 _return_dict["d-or-j_junction"] = value.split(
#                     "(")[1].split(")")[0]
#                 self.junction_together += value.split(
#                     "(")[1].split(")")[0]
#             else:
#                 _return_dict[title.strip().lower().replace(" ", "_")] = value
#                 if value != "N/A":
#                     self.junction_together += value
#         return _return_dict

#     # Now lets start parsing the alignment summaries into each sections FW1,2,3,4,CDR1,2,3
#     def parse_fr1_align(self):
#         _return_dict = {}
#         for title, value in zip(self.alignment_summary_titles, self.fr1_alignment_summary):
#             try:
#                 _return_dict[title] = int(value)
#             except ValueError:
#                 try:
#                     _return_dict[title] = float(value)
#                 except ValueError:
#                     _return_dict[title] = value
#         return _return_dict

#     def parse_fr2_align(self):
#         _return_dict = {}
#         for title, value in zip(self.alignment_summary_titles, self.fr2_alignment_summary):
#             try:
#                 _return_dict[title] = int(value)
#             except ValueError:
#                 try:
#                     _return_dict[title] = float(value)
#                 except ValueError:
#                     _return_dict[title] = value
#         return _return_dict

#     def parse_fr3_align(self):
#         _return_dict = {}
#         for title, value in zip(self.alignment_summary_titles, self.fr3_alignment_summary):
#             try:
#                 _return_dict[title] = int(value)
#             except ValueError:
#                 try:
#                     _return_dict[title] = float(value)
#                 except ValueError:
#                     _return_dict[title] = value
#         return _return_dict

#     def parse_cdr1_align(self):
#         _return_dict = {}
#         for title, value in zip(self.alignment_summary_titles, self.cdr1_alignment_summary):
#             try:
#                 _return_dict[title] = int(value)
#             except ValueError:
#                 try:
#                     _return_dict[title] = float(value)
#                 except ValueError:
#                     _return_dict[title] = value
#         return _return_dict

#     def parse_cdr2_align(self):
#         _return_dict = {}
#         for title, value in zip(self.alignment_summary_titles, self.cdr2_alignment_summary):
#             try:
#                 _return_dict[title] = int(value)
#             except ValueError:
#                 try:
#                     _return_dict[title] = float(value)
#                 except ValueError:
#                     _return_dict[title] = value
#         return _return_dict

#     def parse_cdr3_align(self):
#         _return_dict = {}
#         for title, value in zip(self.alignment_summary_titles, self.cdr3_alignment_summary):
#             try:
#                 _return_dict[title] = int(value)
#             except ValueError:
#                 try:
#                     _return_dict[title] = float(value)
#                 except ValueError:
#                     _return_dict[title] = value
#         return _return_dict

#     def parse_total_align(self):
#         '''The total alignment is just the Whole region which is pretty nice'''
#         _return_dict = {}
#         for title, value in zip(self.alignment_summary_titles, self.total_alignment_summary):
#             try:
#                 _return_dict[title] = int(value)
#             except ValueError:
#                 try:
#                     _return_dict[title] = float(value)
#                 except ValueError:
#                     _return_dict[title] = value
#         return _return_dict

#     def parse_v_hits(self):
#         _return_dict = {}
#         rank = 1
#         for entry in self.hits_v:
#             _entry_dict = {}
#             for value, title in zip(entry.split()[1:], self.hit_fields):
#                 # sometimes there is nothing there so it can't cast to a float
#                 try:
#                     _entry_dict[title.strip().replace(' ', '_').replace("%", "percent")] = float(value)
#                 except:
#                     _entry_dict[title.strip().replace(' ', '_').replace("%", "percent")] = value
#             _return_dict['rank_' + str(rank)] = _entry_dict
#             rank += 1
#         # return a v_hit dictionary that is unsorted, but we kept the rank
#         # by how it appeared in the balst out
#         return _return_dict

#     def parse_d_hits(self):
#         _return_dict = {}
#         rank = 1
#         for entry in self.hits_d:
#             _entry_dict = {}
#             #_entry_dict["rank"] = int(rank)
#             for value, title in zip(entry.split()[1:], self.hit_fields):
#                 try:
#                     _entry_dict[title.strip().replace(' ', '_').replace("%", "percent")] = float(value)
#                 except:
#                     _entry_dict[title.strip().replace(' ', '_').replace("%", "percent")] = value
#             _return_dict['rank_' + str(rank)] = _entry_dict
#             rank += 1
#         return _return_dict

#     def parse_j_hits(self):
#         _return_dict = {}
#         rank = 1
#         for entry in self.hits_j:
#             _entry_dict = {}
#             for value, title in zip(entry.split()[1:], self.hit_fields):
#                 try:
#                     _entry_dict[title.strip().replace(' ', '_').replace("%", "percent")] = float(value)
#                 except:
#                     _entry_dict[title.strip().replace(' ', '_').replace("%", "percent")] = value
#             _return_dict['rank_' + str(rank)] = _entry_dict
#             rank += 1
#         return _return_dict


# def trim_json(blast_dictionary, all_options, gui=False):
#     '''Our Main Function that will return a json type document'''
#     # to be converted to a json document
#     json_dictionary = OrderedDict()
#     if gui:
#         for options in all_options:
#             json_dictionary.update(add_entries(blast_dictionary, json_dictionary, options, gui))
#     else:
#         json_dictionary.update(add_entries(blast_dictionary, json_dictionary, all_options, gui))

#     return json_dictionary


# def add_entries(blast_dictionary, dictionary, options, gui):
#     if gui:
#         for gentry in options:
#             formal = gentry
#             state = options[gentry]['state']
#             json_key = options[gentry]['json_key']
#             split_keys = json_key.split('.')
#             length_of_key = len(split_keys)
#             if state == 0:
#                 continue
#             try:
#                 if length_of_key == 1:
#                     dictionary[formal] = blast_dictionary[json_key]
#                 elif length_of_key == 2:
#                     dictionary[formal] = blast_dictionary[split_keys[0]][split_keys[1]]
#                 elif length_of_key == 3:
#                     dictionary[formal] = blast_dictionary[split_keys[0]][split_keys[1]][split_keys[2]]
#             except KeyError:
#                 dictionary[formal] = ""

#     else:
#         for entry in open(options).readlines():
#             formal = entry.split(',')[0]
#             json_key = entry.split(',')[1].strip()
#             split_keys = json_key.split('.')
#             length_of_key = len(split_keys)
#             if formal.startswith('#'):
#                 continue
#             try:
#                 if length_of_key == 1:
#                     dictionary[formal] = blast_dictionary[json_key]
#                 elif length_of_key == 2:
#                     dictionary[formal] = blast_dictionary[split_keys[0]][split_keys[1]]
#                 elif length_of_key == 3:
#                     dictionary[formal] = blast_dictionary[split_keys[0]][split_keys[1]][split_keys[2]]
#             except KeyError:
#                 dictionary[formal] = ""

#     return dictionary


# def trim_csv(blast_dictionary, general_options, nuc_options, aa_options):
#     pass

# if __name__ == '__main__':
#     to_convert = sys.argv[1]
#     raw_sequences = sys.argv[2]
#     temporary = "."
#     from output_tabs_checkboxes import all_checkboxes as ac
#     igo = igblast_output(to_convert, raw_sequences, temporary, ac, zip_bool=False)
#     igo.parse_blast_file_to_type("testing.csv", o_type="csv")
