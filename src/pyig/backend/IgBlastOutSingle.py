import json
from collections import OrderedDict
from pyig.backend.TranslateAndJoin import TranslateAndJoin


class IgBlastOutSingle():

    '''
    Class for parsing one entry of IgBlastn. This should be subclassed from IgBlastOut. @Todo

    The only method you should use is the parse method after using the constructor.
    '''

    def __init__(self, entry, j_trans, species, debug=False):
        '''
        entry - an iterator containing the lines from blast output
        j_trans - a dictionary with the end of the cdr3 J gene positions
        species - the species we are parsing
        '''
        self.entry = entry
        self.j_trans = j_trans
        self.species = species
        self.debug = debug

        # The main output - this should be the __repr__ function when overloading @TODO
        # Set up with all the arguments you should have except the VDJ hits which are added later
        # This way the JSON can be pretty consistent
        self.output = OrderedDict((('Sequence Id', ""),
                                   ('Query Sequence', ""),
                                   ('Chain type', ""),
                                   ('Format Type', ""),
                                   ('Species', self.species),
                                   ('Top V Hit', ""),
                                   ('Top V Family', ""),
                                   ('Top V Gene', ""),
                                   ('Top D Hit', ""),
                                   ('Top D Family', ""),
                                   ('Top D Gene', ""),
                                   ('Top J Hit', ""),
                                   ('Top J Gene', ""),
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
                                   ('Junction Merged', ""),
                                   ('Isotype', "")))

        # Title fields, these are the four sections it divides up to.
        # The fields are essentially the header to each section.
        # We could hard code this, but considering that it will be unique to each entry
        # We should parse it for each entry
        # See the actual blast output file using --debug and it will make a lot more sense
        self.rearrangement_summary_titles = ""
        self.junction_detail_titles = ""
        self.alignment_summary_titles = ""
        self.hit_fields = ""

        # Here is the meat - all the good stuff goes into here, these are zipped with the
        # titles in the methods
        self.rearrangement_summary = ""
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

        # VD and J hits, can be empty if nothing is found
        # Can be varying size depending on what the user asks for
        # Always take top ranked one for the VDJ assignment
        self.hits_v = []  # to be parsed in another function
        self.hits_d = []  # to be parsed in another function
        self.hits_j = []  # to be parsed in another function

    def get_json_entry(self):
        '''Dumps JSON text from this output'''
        #@Todo - Can change around indent if you want, sometimes mongo complains
        temp_dict = OrderedDict()
        for k, v in self.output.iteritems():
            if v != "":
                temp_dict[k] = v
        self.output = temp_dict
        self.json = json.dumps(self.output, indent=4)
        return self.json

    def get_id(self):
        '''Get the sequence id returned by blast - should match the fasta id'''
        return self.output['Sequence Id']

    def set_seq(self, sequence):
        '''
        Set the sequence from the fasta since
        blast does not keep the input query sequence in the output
        '''
        self.output['Query Sequence'] = str(sequence).upper()

    def set_additional_info(self, additional_info):
        'Set tuple as additional info'
        for field in xrange(0,len(additional_info),2):
          self.output[additional_info[field]] = additional_info[field+1]

    def parse(self):
        '''The method that iterates through the output lines'''

        # tells the iterator when it has reached the end of a certain section
        _rearrangement_breaker = False
        _junction_breaker = False
        _fields_breaker = False

        # Query Name
        for line in self.entry:
                # Starting to iterate
            if "Query" in line:
                self.output['Sequence Id'] = line.split(":")[1].strip()

            if "Domain classification requested:" in line:
                # imgt or kabat
                self.output['Format Type'] = line.split(":")[1].strip()

            if "rearrangement summary" in line:
                '''Okay we found rearrangement summary fields, that means
                the next line will be the actually rearrangement summary.
                So we turn the "breaker on" The rearrangement is just basic info
                returned like the top matches and if was productive etc'''
                self.rearrangement_summary_titles = line.strip().split("(")[2].split(")")[0].split(",")
                _rearrangement_breaker = True
                continue

            if _rearrangement_breaker:
                '''the meat of the rearranment, right after the field titles'''
                self.rearrangement_summary = line.strip().split("\t")

                # shut off rearrangement breaker to break out of this section
                _rearrangement_breaker = False

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

            # Finally parse VDJ hits - it can be found with the # Fields tag
            if "# Fields:" in line:
                self.hit_titles = line.strip().split(":")[1].split(",")
                _fields_breaker = True
            if _fields_breaker:
                '''
                VDJ hits have to be a list
                since there can be more than
                one depending what the user asked for'''

                if line.startswith("V"):
                    self.hits_v.append(line)
                elif line.startswith("D"):
                    self.hits_d.append(line)
                elif line.startswith("J"):
                    self.hits_j.append(line)

        # Now put together titles and of the actual content into the output
        self._parse_rearranment()
        self._parse_junction()
        self._parse_fw1_align()
        self._parse_fw2_align()
        self._parse_fw3_align()
        self._parse_cdr1_align()
        self._parse_cdr2_align()
        self._parse_cdr3_align()
        self._parse_total_v_align()
        self._parse_v_hits()
        self._parse_d_hits()
        self._parse_j_hits()
        self._get_families_and_genes()

    def _get_families_and_genes(self):
        '''parse the top v hits into familes and genes'''
        self.output['Top V Family'] = self.output['Top V Hit'].split('-')[0]
        self.output['Top V Gene'] = self.output['Top V Hit'].split('*')[0]
        self.output['Top D Family'] = self.output['Top D Hit'].split('-')[0]
        self.output['Top D Gene'] = self.output['Top D Hit'].split('*')[0]
        self.output['Top J Gene'] = self.output['Top D Hit'].split('*')[0]

    def _parse_rearranment(self):
        '''parse the rearrangement summary (just the basic statistics) portion of the blast hit'''
        for title, value in zip(self.rearrangement_summary_titles, self.rearrangement_summary):

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

    def _parse_junction(self):
        '''parse the part that has the nucleotides of the junction'''
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
            # they only give you half the J region...We will fix that shortly with translate and join

    # Now lets start parsing the alignment summaries into each sections FW1,2,3,4,CDR1,2,3
    def _parse_fw1_align(self):
        for title, value in zip(self.alignment_summary_titles, self.fr1_alignment_summary):
            if title.strip() == 'percent identity':
                title = 'Identity'
            try:
                self.output["FW1 Alignment " + title.strip().capitalize()] = float(value)
            except ValueError:
                self.output["FW1 Alignment " + title.strip().capitalize()] = value

    def _parse_fw2_align(self):
        for title, value in zip(self.alignment_summary_titles, self.fr2_alignment_summary):
            if title.strip() == 'percent identity':
                title = 'Identity'
            try:
                self.output["FW2 Alignment " + title.strip().capitalize()] = float(value)
            except ValueError:
                self.output["FW2 Alignment " + title.strip().capitalize()] = value

    def _parse_fw3_align(self):
        for title, value in zip(self.alignment_summary_titles, self.fr3_alignment_summary):
            if title.strip() == 'percent identity':
                title = 'Identity'
            try:
                self.output["FW3 Alignment " + title.strip().capitalize()] = float(value)
            except ValueError:
                self.output["FW3 Alignment " + title.strip().capitalize()] = value

    def _parse_cdr1_align(self):
        for title, value in zip(self.alignment_summary_titles, self.cdr1_alignment_summary):
            if title.strip() == 'percent identity':
                title = 'Identity'
            try:
                self.output['CDR1 Alignment ' + title.strip().capitalize()] = float(value)
            except ValueError:
                self.output['CDR1 Alignment ' + title.strip().capitalize()] = value

    def _parse_cdr2_align(self):
        for title, value in zip(self.alignment_summary_titles, self.cdr2_alignment_summary):
            if title.strip() == 'percent identity':
                title = 'Identity'
            try:
                self.output['CDR2 Alignment ' + title.strip().capitalize()] = float(value)
            except ValueError:
                self.output['CDR2 Alignment ' + title.strip().capitalize()] = value

    def _parse_cdr3_align(self):
        for title, value in zip(self.alignment_summary_titles, self.cdr3_alignment_summary):
            if title.strip() == 'percent identity':
                title = 'Identity'
            try:
                self.output['CDR3 Alignment ' + title.strip().capitalize()] = float(value)
            except ValueError:
                self.output['CDR3 Alignment ' + title.strip().capitalize()] = value

    def _parse_total_v_align(self):
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

    # Here are the individual hits table which can be longer than one
    # Go through and rank them by what appears first
    # Append them to a list - it dumps to json as a nice nested array
    def _parse_v_hits(self):
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
            # Adding a whole dictionary that makes it nested
            self.output['V-Gene Rank_' + str(rank)] = _entry_dict

    def _parse_d_hits(self):
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
            # Adding a whole dictionary that makes it nested
            self.output['D-Gene Rank_' + str(rank)] = _entry_dict

    def _parse_j_hits(self):
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
             # Adding a whole dictionary that makes it nested
            self.output['J-Gene Rank_' + str(rank)] = _entry_dict

    def join_and_translate(self):
        '''
        Pass the self class to join and translate so
        join and translate have access to all the member
        variables
        '''
        TranslateAndJoin(self)
        # Finally
        # Find out if the join is productive by seeing if a stop codon (*) is in there
        if "*" not in self.output['CDR3 AA']:
            self.output['Productive CDR3'] = "True"
