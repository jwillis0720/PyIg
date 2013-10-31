import sys
import json
import gzip
from cdr_analyzer import cdr_analyzer

try:
	from Bio.Seq import Seq
	from Bio.Alphabet import IUPAC
except ImportError:
	print ("Need Biopython to use the IgBlast output parser class")


class igblast_output():

    '''Pass the whole output to this class and it will deal with it
    Args:
    file - string for the filename to parse
    output - string for the filename to the output
    '''

    def __init__(self, file, output_file, fasta_files,gz=False):
        breaker = True
        query_holder = []
        self.fasta_files = fasta_files
        _end_dict = {}
        try:
            for line in open('human_germ_properties.txt').readlines():
                line_split = line.split()
                _end_dict[line_split[0]] = int(line_split[1])
        except IOError:
            print "Cant open human_germ_properties.txt, need this file to process CDR3 regions"
            sys.exit()
        if gz:
            z = gzip.open(output_file+".gz",'wb')
        else:
            z = open(output_file,'w')
        with z as f:
            for line in open(file):
                if "IGBLASTN" in line:
                    breaker = False
                    if query_holder:
                        f.write(
                            single_blast_entry(query_holder,self.fasta_files,_end_dict).return_json_document())
                        f.write("\n")
                    query_holder = []
                    continue
                if not breaker:
                    query_holder.append(line)
            f.write(single_blast_entry(query_holder,self.fasta_files,_end_dict).return_json_document())
            f.write("\n")

class single_blast_entry():
    '''The helper class to parse an individual blast result'''

    def __init__(self, query,fasta_files,end_translation_dictionaries):
        _rearrangment_breaker = False
        _junction_breaker = False
        _fields_breaker = False

        self.fasta_files = fasta_files
        #Where to end translation for CDR3 loops specified in an output file
        self.end_translation_dictionaries = end_translation_dictionaries

        # initalize the fields, some can be empty on crappy reads
        # basic
        self.query = ""
        self.full_query_seq = ""
        self.domain_classification = ""

        # title fields
        self.rearrangment_summary_titles = ""
        self.alignment_summary_titles = ""
        self.junction_detail_titles = ""
        self.hit_fields = ""

        # summaries
        self.rearrangment_summary = ""
        self.junction_detail = ""
        self.cdr1_alignment_summary = ""
        self.cdr2_alignment_summary = ""
        self.cdr3_alignment_summary = ""
        self.fr1_alignment_summary = ""
        self.fr2_alignment_summary = ""
        self.fr3_alignment_summary = ""
        self.total_alignment_summary = ""

        # hits
        self.hits_v = []  # to be parsed in another function
        self.hits_d = []  # to be parsed in another function
        self.hits_j = []  # to be parsed in another function

        # will hold everything
        self.blast_dict = {}
        for line in query:
            if "Query" in line:
                self.query = line.split(":")[1].strip()
            if "Domain classification requested:" in line:
                self.domain_classification = line.split(":")[1].strip()
            if "rearrangement summary" in line:
                self.rearrangment_summary_titles = line.strip().split(
                    "(")[2].split(")")[0].split(",")
                _rearrangment_breaker = True
                continue
            if _rearrangment_breaker:
                self.rearrangment_summary = line.strip().split("\t")
                _rearrangment_breaker = False
            if "junction details" in line:
                self.junction_detail_titles = line.strip().split(
                    "(")[2].split(")")[0].split(",")
                _junction_breaker = True
                continue
            if _junction_breaker:
                self.junction_detail = line.strip().split("\t")
                _junction_breaker = False
            if "Alignment summary" in line:
                self.alignment_summary_titles = line.strip().split(
                    "(")[1].split(")")[0].split(",")
                self.alignment_summary_titles = self.alignment_summary_titles
                self.alignment_summary_titles = [x.strip() for x in self.alignment_summary_titles]
            try:
                if line.split()[0].split('-')[0] == "FWR1":
                    self.fr1_alignment_summary = line.strip().split()[1:]
                if line.split()[0].split('-')[0] == "CDR1":
                    self.cdr1_alignment_summary = line.strip().split()[1:]
                if line.split()[0].split('-')[0] == "FWR2":
                    self.fr2_alignment_summary = line.strip().split()[1:]
                if line.split()[0].split('-')[0] == "CDR2":
                    self.cdr2_alignment_summary = line.strip().split()[1:]
                if line.split()[0].split('-')[0] == "FWR3":
                    self.fr3_alignment_summary = line.strip().split()[1:]
                if line.split()[0].split('-')[0] == "CDR3":
                    self.cdr3_alignment_summary = line.strip('\t').split()[1:]
                if line.split()[0].split('-')[0] == "Total":
                    self.total_alignment_summary = line.strip().split()[1:]
            except IndexError:
                pass
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
            if "# Fields:" in line:
                self.hit_fields = line.strip().split(":")[1].split(",")
                _fields_breaker = True
            if _fields_breaker:
                if line.startswith("V"):
                    self.hits_v.append(line)
                elif line.startswith("D"):
                    self.hits_d.append(line)
                elif line.startswith("J"):
                    self.hits_j.append(line)
        self.full_query_seq = self.fasta_files[self.query]
        self.process()

    def process(self):
        self.blast_dict[
            self.query] = {"domain_classification": self.domain_classification,
                           "rearrangement": self.parse_rearranment(),
                           "junction": self.parse_junction(),
                           "fr1_align": self.parse_fr1_align(),
                           "fr2_align": self.parse_fr2_align(),
                           "fr3_align": self.parse_fr3_align(),
                           "cdr1_align": self.parse_cdr1_align(),
                           "cdr2_align": self.parse_cdr2_align(),
                           "cdr3_align": self.parse_cdr3_align(),
                           "total_align": self.parse_total_align(),
                           "v_hits": self.parse_v_hits(),
                           "d_hits": self.parse_d_hits(),
                           "j_hits": self.parse_j_hits()
                           }
        self.translate_junction()

    def translate_junction(self):
    	self.cdr3_partial = ""
    	if self.junction_together:
    		coding_region = Seq(self.junction_together,IUPAC.ambiguous_dna)
    		self.cdr3_partial = str(coding_region.translate())

    def parse_rearranment(self):
        _return_dict = {}
        for title, value in zip(self.rearrangment_summary_titles, self.rearrangment_summary):
            if len(value.split(',')) > 1:
                #cast multiple entries for tuple, makes them easier for json
                _return_dict[title.strip()] = tuple(value.split(','))
            else:
                _return_dict[title.strip()] = value
        return _return_dict

    def parse_junction(self):
        _return_dict = {}
        self.junction_together = ""
        for title, value in zip(self.junction_detail_titles, self.junction_detail):
            if "(" in value:
                _return_dict["d-or-j_junction"] = value.split(
                    "(")[1].split(")")[0]
                self.junction_together += value.split(
                    "(")[1].split(")")[0]
            else:
                _return_dict[title.strip().lower().replace(" ", "_")] = value
               	if value != "N/A":
               		self.junction_together += value
        return _return_dict

    def parse_fr1_align(self):
            _return_dict = {}
            for title, value in zip(self.alignment_summary_titles, self.fr1_alignment_summary):
                try:
                    _return_dict[title] = int(value)
                except ValueError:
                    try:
                        _return_dict[title] = float(value)
                    except ValueError:
                        _return_dict[title] = value
            return _return_dict

    def parse_fr2_align(self):
            _return_dict = {}
            for title, value in zip(self.alignment_summary_titles, self.fr2_alignment_summary):
                try:
                    _return_dict[title] = int(value)
                except ValueError:
                    try:
                        _return_dict[title] = float(value)
                    except ValueError:
                        _return_dict[title] = value
            return _return_dict

    def parse_fr3_align(self):
        _return_dict = {}
        for title, value in zip(self.alignment_summary_titles, self.fr3_alignment_summary):
            try:
                _return_dict[title] = int(value)
            except ValueError:
                try:
                    _return_dict[title] = float(value)
                except ValueError:
                    _return_dict[title] = value
        return _return_dict

    def parse_cdr1_align(self):
        _return_dict = {}
        for title, value in zip(self.alignment_summary_titles, self.cdr1_alignment_summary):
            try:
                _return_dict[title] = int(value)
            except ValueError:
                try:
                    _return_dict[title] = float(value)
                except ValueError:
                    _return_dict[title] = value
        return _return_dict
    
    def parse_cdr2_align(self):
        _return_dict = {}
        for title, value in zip(self.alignment_summary_titles, self.cdr2_alignment_summary):
            try:
                _return_dict[title] = int(value)
            except ValueError:
                try:
                    _return_dict[title] = float(value)
                except ValueError:
                    _return_dict[title] = value
        return _return_dict
 
    def parse_cdr3_align(self):
        _return_dict = {}
        for title, value in zip(self.alignment_summary_titles, self.cdr3_alignment_summary):
            try:
                _return_dict[title] = int(value)
            except ValueError:
                try:
                    _return_dict[title] = float(value)
                except ValueError:
                    _return_dict[title] = value
        return _return_dict

    def parse_total_align(self):
        _return_dict = {}
        for title, value in zip(self.alignment_summary_titles, self.total_alignment_summary):
            try:
                _return_dict[title] = int(value)
            except ValueError:
                try:
                    _return_dict[title] = float(value)
                except ValueError:
                    _return_dict[title] = value
        return _return_dict

    def parse_v_hits(self):
        _return_dict = {}
        rank = 1
        for entry in self.hits_v:
            _entry_dict = {}
            #_entry_dict['rank'] = int(rank)
            for value, title in zip(entry.split()[1:], self.hit_fields):
                try:
                    _entry_dict[title.strip().replace(' ', '_')] = float(value)
                except:
                    _entry_dict[title.strip().replace(' ', '_')] = value
            _return_dict['rank_' + str(rank)] = _entry_dict
            rank += 1
        return _return_dict

    def parse_d_hits(self):
        _return_dict = {}
        rank = 1
        for entry in self.hits_d:
            _entry_dict = {}
            #_entry_dict["rank"] = int(rank)
            for value, title in zip(entry.split()[1:], self.hit_fields):
                try:
                    _entry_dict[title.strip().replace(' ', '_')] = float(value)
                except:
                    _entry_dict[title.strip().replace(' ', '_')] = value
            _return_dict['rank_' + str(rank)] = _entry_dict
            rank += 1
        return _return_dict

    def parse_j_hits(self):
        _return_dict = {}
        rank = 1
        for entry in self.hits_j:
            _entry_dict = {}
            #_entry_dict['rank'] = int(rank)
            for value, title in zip(entry.split()[1:], self.hit_fields):
                try:
                    _entry_dict[title.strip().replace(' ', '_')] = float(value)
                except:
                    _entry_dict[title.strip().replace(' ', '_')] = value
            _return_dict['rank_' + str(rank)] = _entry_dict
            rank += 1
        return _return_dict

    def return_json_document(self):
        '''Our Main Function that will return a json type document'''
        # to be converted to a json document
        self.json_dictionary = {}
        
        #hits arrays if we have more than one hit we kept in the blast query
        self.v_hits_array = []
        self.j_hits_array = []
        self.d_hits_array = []

        self.json_dictionary = {
            "_id": self.query,
            "format": self.blast_dict[self.query]['domain_classification']
        }

        # Most important should be considered individually
        try:
            self.json_dictionary["top_v"] = self.blast_dict[
                self.query]['rearrangement']['Top V gene match']
        except KeyError:
            self.json_dictionary["top_v"] = "N/A"
        try:
            self.json_dictionary["top_d"] = self.blast_dict[
                self.query]['rearrangement']['Top D gene match']
        except KeyError:
            self.json_dictionary["top_d"] = "N/A"
        try:
            self.json_dictionary["top_j"] = self.blast_dict[
                self.query]['rearrangement']['Top J gene match']
        except KeyError:
            self.json_dictionary["top_j"] = "N/A",
        try:
            self.json_dictionary["strand"] = self.blast_dict[
                self.query]['rearrangement']['Strand']
        except KeyError:
            self.json_dictionary["strand"] = "N/A"
        try:
            self.json_dictionary["chain_type"] = self.blast_dict[
                self.query]['rearrangement']['Chain type']
        except KeyError:
            self.json_dictionary["chain_type"] = "N/A"
        try:
            self.json_dictionary["stop_codon"] = self.blast_dict[
                self.query]['rearrangement']['stop codon']
        except KeyError:
            self.json_dictionary["stop_codon"] = "N/A"
        try:
            self.json_dictionary["productive"] = self.blast_dict[
                self.query]['rearrangement']['Productive']
        except KeyError:
            self.json_dictionary["productive"] = "N/A"
        try:
            self.json_dictionary["in_frame"] = self.blast_dict[
                self.query]['rearrangement']['V-J frame']
        except KeyError:
            self.json_dictionary["in_frame"] = "N/A"

        # add junctions. won't modify key names this time
        try:
            for junction_entry in self.blast_dict[self.query]['junction']:
                self.json_dictionary[junction_entry] = self.blast_dict[
                    self.query]['junction'][junction_entry]
        except KeyError:
            for junction_title in self.junction_detail_titles:
                self.json_dictionary[junction_title] = "N/A"

        # alignment_summary will be empty if it is empty, no need for try and
        self.alignment_summaries = {}
        if self.blast_dict[self.query]['fr1_align']:
            self.alignment_summaries['fr1_align'] = self.blast_dict[self.query]['fr1_align']

        else:
            self.alignment_summaries['fr1_align'] = 'N/A'

        if self.blast_dict[self.query]['cdr1_align']:
            self.alignment_summaries['cdr1_align'] = self.blast_dict[self.query]['cdr1_align']

        else:
            self.alignment_summaries['cdr1_align'] = 'N/A'
        
        if self.blast_dict[self.query]['fr2_align']:
            self.alignment_summaries['fr2_align'] = self.blast_dict[self.query]['fr2_align']

        else:
            self.alignment_summaries['fr2_align'] = 'N/A'

        if self.blast_dict[self.query]['cdr2_align']:
            self.alignment_summaries['cdr2_align'] = self.blast_dict[self.query]['cdr2_align']

        else:
            self.alignment_summaries['cdr2_align'] = 'N/A'
        
        if self.blast_dict[self.query]['fr3_align']:
            self.alignment_summaries['fr3_align'] = self.blast_dict[self.query]['fr3_align']

        else:
            self.alignment_summaries['fr3_align'] = 'N/A'
        
        if self.blast_dict[self.query]['cdr3_align']:
            self.alignment_summaries['cdr3_align'] = self.blast_dict[self.query]['cdr3_align']

        else:
            self.alignment_summaries['cdr3_align'] = 'N/A'
        
        if self.blast_dict[self.query]['total_align']:
            self.alignment_summaries['total_align'] = self.blast_dict[self.query]['total_align']
         
        else:
            self.alignment_summaries['total_align'] = 'N/A'

        self.json_dictionary['alignment_summaries'] = self.alignment_summaries

        # vhits
        try:
            for rank in sorted(self.blast_dict[self.query]['v_hits']):
                self.v_hits_array.append(
                    {rank: self.blast_dict[self.query]['v_hits'][rank]})
        except ValueError:
            self.v_hits_array = "N/A"

        # dhits
        try:
            for rank in sorted(self.blast_dict[self.query]['d_hits']):
                self.d_hits_array.append(
                    {rank: self.blast_dict[self.query]['d_hits'][rank]})
        except ValueError:
            self.d_hits_array = "N/A"

        # jhits
        try:
            for rank in sorted(self.blast_dict[self.query]['j_hits']):
                self.j_hits_array.append(
                    {rank: self.blast_dict[self.query]['j_hits'][rank]})
        except KeyError:
            self.j_hits_array = "N/A"

        self.json_dictionary["v_hits"] = self.v_hits_array
        self.json_dictionary["d_hits"] = self.d_hits_array
        self.json_dictionary['j_hits'] = self.j_hits_array

        if self.json_dictionary["productive"].lower()  == "yes":
        	self.json_dictionary["partial_cdr3_aa"] = self.cdr3_partial

        self.json_dictionary = cdr_analyzer(self.json_dictionary,self.full_query_seq,self.end_translation_dictionaries).return_json_dict_with_cdr_analysis() 


        # convert dictionary to json object
        self.json = json.dumps(self.json_dictionary, sort_keys=1)

        # and finally return the object
        return self.json

if __name__ == '__main__':
	to_convert = sys.argv[1]
	igblast_output(to_convert,"test_out.json")
