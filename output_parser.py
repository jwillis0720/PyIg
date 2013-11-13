import sys
import json
import gzip
from cdr_analyzer import cdr_analyzer

# try:
# except ImportError:
#    print("Need Biopython to use the IgBlast output parser class")


class igblast_output():

    '''Pass the whole output to this class and it will deal with it
    Args:
    blast file - from output format 7 of blastign
    general_options - The options you want to output for general options passed as a dictionary
    nuc_options - The nucleotide options you want passed as a dictionary
    aa_options - The
    output - string for the filename to the output
    '''

    def __init__(self, blast_file, general_options,
                 nuc_options, aa_options, zip_bool=False):

        # The blast file that was output, format it to a handle
        self.blast_file_handle = open(blast_file)
        # the file to use that will tell us where the CDR3 ends, pass this to the CDR analyzer
        self.end_dict = {}
        # these are dictionaries with the options we want to output. See output_tabs checkbozes
        self.general_options = general_options
        self.nuc_options = nuc_options
        self.aa_options = aa_options
        # bool to tell us if we want the file zipped
        self.zip = zip_bool
        try:
            for line in open('junctional_data/human_germ_properties.txt').readlines():
                line_split = line.split()
                self.end_dict[line_split[0]] = int(line_split[1])
        except IOError:
            print "Cant open human_germ_properties.txt, \
            need this file to process CDR3 regions"
            sys.exit()

    def parse_blast_file_to_type(self, output_file, o_type="json"):
        focus_lines = []
        if self.zip:
            json_output_file_handle = gzip.open(output_file + ".gz", 'wb')
        else:
            json_output_file_handle = open(output_file, 'w')
        with json_output_file_handle as openfile:
            for line in self.blast_file_handle:
                if "IGBLASTN" in line:
                    if focus_lines:
                        sbe = single_blast_entry(focus_lines, self.end_dict)
                        blast_dictionary = sbe.generate_blast_dict()  # return single blast entry
                        if o_type == "json":
                            json_document_trimmed = trim_json(blast_dictionary,
                                                              self.general_options, self.nuc_options, self.aa_options)  # trim the json according to input
                            openfile.write(json_document_trimmed + "\n")  # write it out
                            focus_lines = []
                        if o_type == "csv":
                            csv_document_trimmed = trim_csv(blast_dictionary,
                                                            self.general_options, self.nuc_options, self.aa_options)  # trim the json according to input
                            openfile.write(csv_document_trimmed + "\n")  # write it out
                            focus_lines = []
                    else:
                        continue
                else:
                    focus_lines.append(line)
            # and do it for the end too
            sbe = single_blast_entry(focus_lines, self.end_dict)
            blast_dictionary = sbe.generate_blast_dict()  # return single blast entry
            if o_type == "json":
                json_document_trimmed = trim_json(blast_dictionary,
                                                  self.general_options, self.nuc_options, self.aa_option)  # trim the json according to input
                openfile.write(json_document_trimmed + "\n")  # write it out
            elif o_type == "csv":
                csv_document_trimmed = trim_csv(blast_dictionary,
                                                self.general_options, self.nuc_options, self.aa_option)  # trim the json according to input
                openfile.write(csv_document_trimmed + "\n")  # write it out


class single_blast_entry():

    '''The helper class to parse an individual blast result'''

    def __init__(self, single_entry, end_cdr3_dicts):
        # some breaker options that will tell us when to go on to the next section
        _rearrangment_breaker = False
        _junction_breaker = False
        _fields_breaker = False

        # Where to end translation for CDR3 loops specified in an output file
        self.end_cdr3_dicts = end_cdr3_dicts

        # initalize the fields, some can be empty on crappy reads
        # basic
        self.query = ""
        self.full_query_seq = ""
        self.domain_classification = ""

        # title fields, these are the four sections it divides up to.
        # The fields are essentially the header to each section.
        # We could hard code this, but considering that it will be unique to each entry
        # We should parse it for each entry
        self.rearrangment_summary_titles = ""
        self.alignment_summary_titles = ""
        self.junction_detail_titles = ""
        self.hit_fields = ""

        # Here is the meat.
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

        for line in single_entry:
            if "Query" in line:
                # this will be the id field which blast calls the quey
                self.query = line.split(":")[1].strip()

            if "Domain classification requested:" in line:
                # imgt or kabat
                self.domain_classification = line.split(":")[1].strip()

            if "rearrangement summary" in line:
                '''Okay we found rearrangement summary fields, that means
                the next line will be the actually rearrangment summary.
                So we turn the "breaker on" The rearrangment is just basic info
                returned like the top matches and if was productive etc'''
                self.rearrangment_summary_titles = line.strip().split(
                    "(")[2].split(")")[0].split(",")
                self.rearrangment_summary_titles = [x.strip().lower().replace(" ", "_") for x in self.rearrangment_summary_titles]
                _rearrangment_breaker = True
                continue

            if _rearrangment_breaker:
                # the meat of the rearranment, right after the fields
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
                # strip off white marks in the fields
                self.alignment_summary_titles = [x.strip() for x in self.alignment_summary_titles]
            '''Ok, in some blast versions, the line starts with just the field we are looking for,
            In other versions it has a '-', so I will just check for both'''

            # check for hypen
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

            '''Last but not least, we can pass all the VDJ hits along with all the
            info about them that blast gives, the VDJ hits can give a ton of information
            about the junction, the identity to each gene, the frame, gaps, evalue, etc,
            these are ranked and put into a list'''
            if "# Fields:" in line:
                self.hit_fields = line.strip().split(":")[1].split(",")
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
        '''now that I got all the blast output parsed,
        lets put it into the requested format, before we do that however lets make
        one dictionary'''

    def generate_blast_dict(self):
        '''THE blast dict is a dictionary of dictionaries containing all of
        our information about the hit, this will be filled by several methods,
        each of which return their own dictionaries'''
        blast_dict = {}
        blast_dict[self.query] = {
            "domain_classification": self.domain_classification,
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

        #blast_dict = cdr_analyzer()

        return blast_dict

    def parse_rearranment(self):
        '''parse the rearrangement summary (just the basic statistics)
        portion of the blast hit'''
        _return_dict = {}
        for title, value in zip(self.rearrangment_summary_titles, self.rearrangment_summary):
            if len(value.split(',')) > 1:
                '''sometimes there are multiple values for the rearranment values'''
                _return_dict[title.strip()] = tuple(value.split(','))
            else:
                _return_dict[title.strip()] = value
        return _return_dict

    def parse_junction(self):
        '''Return a dictionary of the junction that is the
        nuceotide sequence of teh CDR3 junctions'''
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

    # Now lets start parsing the alignment summaries into each sections FW1,2,3,4,CDR1,2,3
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
        '''The total alignment is just the Whole region which is pretty nice'''
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
            for value, title in zip(entry.split()[1:], self.hit_fields):
                # sometimes there is nothing there so it can't cast to a float
                try:
                    _entry_dict[title.strip().replace(' ', '_')] = float(value)
                except:
                    _entry_dict[title.strip().replace(' ', '_')] = value
            _return_dict['rank_' + str(rank)] = _entry_dict
            rank += 1
        # return a v_hit dictionary that is unsorted, but we kept the rank
        # by how it appeared in the balst out
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


def trim_json(blast_dictionary, general_options, nuc_options, aa_options):
        '''Our Main Function that will return a json type document'''
        # to be converted to a json document
        json_dictionary = {}

        # hits arrays if we have more than one hit we kept in the blast query
        query = blast_dictionary.keys()[0]

        blast_dictionary = blast_dictionary[query]
        for general_entry in general_options:
            default = general_entry['default']
            key = general_entry['json_key']
            formal = general_entry['formal']
            if default == 0:
                continue
            if key == "_id":
                json_dictionary["_id"] = query
                continue
            split_keys = key.split('.')
            length_of_key = len(split_keys)
            if length_of_key == 1:
                json_dictionary[formal] = blast_dictionary[key]
            elif length_of_key == 2:
                json_dictionary[formal] = blast_dictionary[split_keys[0]][split_keys[1]]
            elif length_of_key == 3:
                json_dictionary[formal] = blast_dictionary[split_keys[0]][split_keys[1]][split_keys[2]]

        # Most important should be considered individually
        # try:
        #     self.json_dictionary["top_v"] = self.blast_dict[
        #         self.query]['rearrangement']['Top V gene match']
        # except KeyError:
        #     self.json_dictionary["top_v"] = "N/A"
        # try:
        #     self.json_dictionary["top_d"] = self.blast_dict[
        #         self.query]['rearrangement']['Top D gene match']
        # except KeyError:
        #     self.json_dictionary["top_d"] = "N/A"
        # try:
        #     self.json_dictionary["top_j"] = self.blast_dict[
        #         self.query]['rearrangement']['Top J gene match']
        # except KeyError:
        #     self.json_dictionary["top_j"] = "N/A",
        # try:
        #     self.json_dictionary["strand"] = self.blast_dict[
        #         self.query]['rearrangement']['Strand']
        # except KeyError:
        #     self.json_dictionary["strand"] = "N/A"
        # try:
        #     self.json_dictionary["chain_type"] = self.blast_dict[
        #         self.query]['rearrangement']['Chain type']
        # except KeyError:
        #     self.json_dictionary["chain_type"] = "N/A"
        # try:
        #     self.json_dictionary["stop_codon"] = self.blast_dict[
        #         self.query]['rearrangement']['stop codon']
        # except KeyError:
        #     self.json_dictionary["stop_codon"] = "N/A"
        # try:
        #     self.json_dictionary["productive"] = self.blast_dict[
        #         self.query]['rearrangement']['Productive']
        # except KeyError:
        #     self.json_dictionary["productive"] = "N/A"
        # try:
        #     self.json_dictionary["in_frame"] = self.blast_dict[
        #         self.query]['rearrangement']['V-J frame']
        # except KeyError:
        #     self.json_dictionary["in_frame"] = "N/A"

        # add junctions. won't modify key names this time
        # try:
        #     for junction_entry in self.blast_dict[self.query]['junction']:
        #         self.json_dictionary[junction_entry] = self.blast_dict[
        #             self.query]['junction'][junction_entry]
        # except KeyError:
        #     for junction_title in self.junction_detail_titles:
        #         self.json_dictionary[junction_title] = "N/A"

        # alignment_summary will be empty if it is empty, no need for try and
        # self.alignment_summaries = {}
        # if self.blast_dict[self.query]['fr1_align']:
        #     self.alignment_summaries[
        #         'fr1_align'] = self.blast_dict[self.query]['fr1_align']

        # else:
        #     self.alignment_summaries['fr1_align'] = 'N/A'

        # if self.blast_dict[self.query]['cdr1_align']:
        #     self.alignment_summaries[
        #         'cdr1_align'] = self.blast_dict[self.query]['cdr1_align']

        # else:
        #     self.alignment_summaries['cdr1_align'] = 'N/A'

        # if self.blast_dict[self.query]['fr2_align']:
        #     self.alignment_summaries[
        #         'fr2_align'] = self.blast_dict[self.query]['fr2_align']

        # else:
        #     self.alignment_summaries['fr2_align'] = 'N/A'

        # if self.blast_dict[self.query]['cdr2_align']:
        #     self.alignment_summaries[
        #         'cdr2_align'] = self.blast_dict[self.query]['cdr2_align']

        # else:
        #     self.alignment_summaries['cdr2_align'] = 'N/A'

        # if self.blast_dict[self.query]['fr3_align']:
        #     self.alignment_summaries[
        #         'fr3_align'] = self.blast_dict[self.query]['fr3_align']

        # else:
        #     self.alignment_summaries['fr3_align'] = 'N/A'

        # if self.blast_dict[self.query]['cdr3_align']:
        #     self.alignment_summaries[
        #         'cdr3_align'] = self.blast_dict[self.query]['cdr3_align']

        # else:
        #     self.alignment_summaries['cdr3_align'] = 'N/A'

        # if self.blast_dict[self.query]['total_align']:
        #     self.alignment_summaries['total_align'] = self.blast_dict[
        #         self.query]['total_align']

        # else:
        #     self.alignment_summaries['total_align'] = 'N/A'

        # self.json_dictionary['alignment_summaries'] = self.alignment_summaries

        # vhits
        # try:
        #     for rank in sorted(self.blast_dict[self.query]['v_hits']):
        #         self.v_hits_array.append(
        #             {rank: self.blast_dict[self.query]['v_hits'][rank]})
        # except ValueError:
        #     self.v_hits_array = "N/A"

        # dhits
        # try:
        #     for rank in sorted(self.blast_dict[self.query]['d_hits']):
        #         self.d_hits_array.append(
        #             {rank: self.blast_dict[self.query]['d_hits'][rank]})
        # except ValueError:
        #     self.d_hits_array = "N/A"

        # jhits
        # try:
        #     for rank in sorted(self.blast_dict[self.query]['j_hits']):
        #         self.j_hits_array.append(
        #             {rank: self.blast_dict[self.query]['j_hits'][rank]})
        # except KeyError:
        #     self.j_hits_array = "N/A"

        # self.json_dictionary["v_hits"] = self.v_hits_array
        # self.json_dictionary["d_hits"] = self.d_hits_array
        # self.json_dictionary['j_hits'] = self.j_hits_array

        # if self.json_dictionary["productive"].lower() == "yes":
        #     self.json_dictionary["partial_cdr3_aa"] = self.cdr3_partial

        # self.json_dictionary = cdr_analyzer(
        #     self.json_dictionary, self.full_query_seq, self.end_translation_dictionaries).return_json_dict_with_cdr_analysis()

        # convert dictionary to json object
        # self.json = json.dumps(self.json_dictionary, sort_keys=1)

        # and finally return the object
        # return self.json


def trim_csv(blast_dictionary, general_options, nuc_options, aa_options):
    pass

if __name__ == '__main__':
    to_convert = sys.argv[1]
    from output_tabs_checkboxes import all_checkboxes as ac
    igo = igblast_output(to_convert, ac['general'], ac['nucleotide'], ac['amino'], zip_bool=True)
    igo.parse_blast_file_to_type("testing.json.gz", o_type="json")
