import os
import sys
import warnings
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
warnings.filterwarnings('ignore')

class TranslateAndJoin():

    def __init__(self, IgO):
        self.species = 'human'
        self.IgO = IgO
        self.output = IgO.output
        self.sequence = self.output['Query Sequence']
        self.seq_id = self.output['Sequence Id']
        self.debug = IgO.debug

        try:
            self.framework1_set_and_translate()
        except (KeyError, ValueError):
            if self.debug:
                print "FW1 does not exist for {0}".format(self.seq_id)
            pass

        try:
            self.framework2_set_and_translate()
        except (KeyError, ValueError):
            if self.debug:
                print "FW2 does not exist for {0}".format(self.seq_id)
            pass

        try:
            self.framework3_set_and_translate()
        except (KeyError, ValueError):
            if self.debug:
                print "FW3 does not exist for {0}".format(self.seq_id)
            pass

        try:
            self.CDR1_set_and_translate()
        except (KeyError, ValueError):
            if self.debug:
                print "CDR1 does not exist for {0}".format(self.seq_id)
            pass

        try:
            self.CDR2_set_and_translate()
        except (KeyError, ValueError):
            if self.debug:
                print "CDR2 does not exist for {0}".format(self.seq_id)
            pass

        try:
            self.CDR3_set_and_translate()
        except (KeyError, ValueError):
            if self.debug:
                print "CDR3 does not exist for {0}".format(self.seq_id)
            pass
        try:
            self.FW4_set_and_translate()
        except (KeyError, ValueError):
            if self.debug:
                print "FW4 does not exist for {0}".format(self.seq_id)
            pass

    def framework1_set_and_translate(self):
        _from = int(self.output["FW1 Alignment From"])
        _to = int(self.output["FW1 Alignment To"])
        self.output['Framework 1 Nucleotides'] = self.sequence[_from - 1:_to]
        self.output['Framework 1 AA'] = str(
            Seq(self.output['Framework 1 Nucleotides'], IUPAC.ambiguous_dna).translate())
        self.output['Framework 1 AA Length'] = len(self.output['Framework 1 AA'])

    def framework2_set_and_translate(self):
        _from = int(self.output["FW2 Alignment From"])
        _to = int(self.output["FW2 Alignment To"])
        self.output['Framework 2 Nucleotides'] = self.sequence[_from - 1:_to]
        self.output['Framework 2 AA'] = str(
            Seq(self.output['Framework 2 Nucleotides'], IUPAC.ambiguous_dna).translate())
        self.output['Framework 2 AA Length'] = len(self.output['Framework 2 AA'])

    def framework3_set_and_translate(self):
        _from = int(self.output["FW3 Alignment From"])
        _to = int(self.output["FW3 Alignment To"])
        self.output['Framework 3 Nucleotides'] = self.sequence[_from - 1:_to]
        self.output['Framework 3 AA'] = str(
            Seq(self.output['Framework 3 Nucleotides'], IUPAC.ambiguous_dna).translate())
        self.output['Framework 3 AA Length'] = len(self.output['Framework 3 AA'])

    def CDR1_set_and_translate(self):
        _from = int(self.output["CDR1 Alignment From"])
        _to = int(self.output["CDR1 Alignment To"])
        self.output['CDR1 Nucleotides'] = self.sequence[_from - 1:_to]
        self.output['CDR1 AA'] = str(
            Seq(self.output['CDR1 Nucleotides'], IUPAC.ambiguous_dna).translate())
        self.output['CDR1 AA Length'] = len(self.output['CDR1 AA'])

    def CDR2_set_and_translate(self):
        _from = int(self.output["CDR2 Alignment From"])
        _to = int(self.output["CDR2 Alignment To"])
        self.output['CDR2 Nucleotides'] = self.sequence[_from - 1:_to]
        self.output['CDR2 AA'] = str(Seq(self.output['CDR2 Nucleotides'], IUPAC.ambiguous_dna).translate())
        self.output['CDR2 AA Length'] = len(self.output['CDR2 AA'])

    def CDR3_set_and_translate(self):
        j_gene = self.output['J-Gene Rank_1']['J-Gene Rank_1 Subject id']
        try:
            j_length = int(self.IgO.j_trans[j_gene])
        except KeyError:
            raise Exception(
                "Germline Gene {0} not found in CDR3 Lenth Description, this needs to be set to determine when the J reading frame ends for the CDR3, refer to documentation for setting this gene".format(j_gene))
            sys.exit(1)

        j_begin = int(self.output['J-Gene Rank_1']['J-Gene Rank_1 Q. start']) - 1
        j_end = int(self.output['J-Gene Rank_1']['J-Gene Rank_1 S. start'])
        j_end = (j_length - j_end) + j_begin

        if self.output['CDR3 Alignment From']:
            v_part_from = int(self.output['CDR3 Alignment From']) - 1
            v_part_to = int(self.output['CDR3 Alignment To'])
            v_part_nuc = self.sequence[v_part_from - 1:v_part_to]

        else:
            v_part_from = self.sequence.index(self.IgO.junction_merged[5:]) - 1

        self.output['CDR3 Nucleotides'] = self.sequence[v_part_from:j_end]
        self.output['CDR3 AA'] = str(Seq(self.output['CDR3 Nucleotides'], IUPAC.ambiguous_dna).translate())
        self.output['CDR3 AA Length'] = len(self.output['CDR3 AA'])

    def FW4_set_and_translate(self):
        _from = self.sequence.index(self.output['CDR3 Nucleotides']) + len(self.output['CDR3 Nucleotides'])
        _to = int(self.output['J-Gene Rank_1']['J-Gene Rank_1 Q. end'])
        self.output['Framework 4 Nucleotides'] = self.sequence[_from - 2:_to]
        self.output['Framework 4 AA'] = str(
            Seq(self.output['Framework 4 Nucleotides'], IUPAC.ambiguous_dna).translate())
        self.output['Framework 4 AA Length'] = len(self.output['Framework 4 AA'])
