import os
import sys
import warnings

warnings.filterwarnings('ignore')

try:
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC
except ImportError:
    print("Need Biopython to use the IgBlast output parser class")


class TranslateAndJoin():

    def __init__(self,IgO):
        self.species = 'human'
        self.IgO = IgO      
        self.output = IgO.output
        self.sequence = self.output['Query Sequence']
        self.framework1_set_and_translate()
        self.framework2_set_and_translate()
        self.framework3_set_and_translate()
        self.CDR1_set_and_translate()
        self.CDR2_set_and_translate()
        self.CDR3_set_and_translate()

    def framework1_set_and_translate(self):
        _from = int(self.output["FW1 Alignment From"])
        _to = int(self.output["FW1 Alignment To"])
        self.output['Framework 1 Nucleotides'] = self.sequence[_from-1:_to]
        self.output['Framework 1 AA'] = str(Seq(self.output['Framework 1 Nucleotides'], IUPAC.ambiguous_dna).translate())
        self.output['Framework 1 AA Length'] = len(self.output['Framework 1 AA'])


    def framework2_set_and_translate(self):
        _from = int(self.output["FW2 Alignment From"])
        _to = int(self.output["FW2 Alignment To"])
        self.output['Framework 2 Nucleotides'] = self.sequence[_from-1:_to]
        self.output['Framework 2 AA'] = str(Seq(self.output['Framework 2 Nucleotides'], IUPAC.ambiguous_dna).translate())
        self.output['Framework 2 AA Length'] = len(self.output['Framework 2 AA'])


    def framework3_set_and_translate(self):
        _from = int(self.output["FW3 Alignment From"])
        _to = int(self.output["FW3 Alignment To"])
        self.output['Framework 3 Nucleotides'] = self.sequence[_from-1:_to]
        self.output['Framework 3 AA'] = str(Seq(self.output['Framework 3 Nucleotides'], IUPAC.ambiguous_dna).translate())
        self.output['Framework 3 AA Length'] = len(self.output['Framework 3 AA'])


    def CDR1_set_and_translate(self):
        _from = int(self.output["CDR1 Alignment From"])
        _to = int(self.output["CDR1 Alignment To"])
        self.output['CDR1 Nucleotides'] = self.sequence[_from-1:_to]
        self.output['CDR1 AA'] = str(Seq(self.output['CDR1 Nucleotides'], IUPAC.ambiguous_dna).translate())
        self.output['CDR1 AA Length'] = len(self.output['CDR1 AA'])

    def CDR2_set_and_translate(self):
        _from = int(self.output["CDR2 Alignment From"])
        _to = int(self.output["CDR2 Alignment To"])
        self.output['CDR2 Nucleotides'] = self.sequence[_from-1:_to]
        self.output['CDR2 AA'] = str(Seq(self.output['CDR2 Nucleotides'], IUPAC.ambiguous_dna).translate())
        self.output['CDR2 AA Length'] = len(self.output['CDR2 AA'])


    def CDR3_set_and_translate(self):
        v_part_from = int(self.output['CDR3 Alignment From'])
        v_part_to = int(self.output['CDR3 Alignment To'])
        v_part_nuc = self.sequence[v_part_from-1:v_part_to]
        junction_merged = self.IgO.junction_merged
       
        try: 
            j_gene =  self.output['J-Gene Rank_1']['J-Gene Rank_1 Subject id']
            j_begin = int(self.output['J-Gene Rank_1']['J-Gene Rank_1 Q. start']) - 1
            j_length = int(self.IgO.j_trans[j_gene])
            j_end = int(self.output['J-Gene Rank_1']['J-Gene Rank_1 S. start'])
            j_end = (j_length - j_end) + j_begin
            self.output['CDR3 Nucleotides'] = self.sequence[v_part_from-1:j_end]
            self.output['CDR3 AA'] = str(Seq(self.output['CDR3 Nucleotides'], IUPAC.ambiguous_dna).translate())
            self.output['CDR3 AA Length'] = len(self.output['CDR3 AA'])
        except KeyError:
            print self.sequence



    #     '''Now the J region, which will be the beginning of the V region of the
    #     CDR3 to the end of the translated region of the J gene, the translated region
    #     is not implicitly coded by igblastn so we have to manually look at it from a file'''
    #     # j_region
    #     try:
    #         self.j_region_info = self.unmodified['j_hits']['rank_1']
    #         # where does the query start on the j gene
    #         _begin = int(self.j_region_info['q._start']) - 1
    #         # name of jgene
    #         self.j_gene = self.j_region_info['subject_id']

    #         # where does the CDR3 loop stop in the j_gene stop translating for
    #         # this particular gene
    #         self.j_gene_end = self.end_translation_dictionaries[self.j_gene]

    #         # and make it relative to the sequence read
    #         _end = (self.j_gene_end - int(
    #             self.j_region_info['s._start'])) + _begin

    #         # j_region
    #         self.j_region = Seq(query_seq[_begin:_end])
    #         self.full_cdr3 = Seq(
    #             query_seq[self.v_part_of_cdr3_info['from'] - 1:_end])
    #         self.full_cdr3_aa = self.full_cdr3.translate()

    #         # get the rest of the sequence until the end
    #         self.framework4_seq = Seq(query_seq[_end - 2:])
    #         self.framework4_seq_aa = self.framework4_seq.translate()
    #     except (TypeError, KeyError, IndexError):
    #         self.j_gene = ""
    #         self.j_region = ""
    #         self.full_cdr3 = ""
    #         self.full_cdr3_aa = ""
    #         self.framework4_seq = ""
    #         self.framework4_seq_aa = ""

    #     self.frames_and_cdrs = {
    #         'fw1': self.framework1_seq,
    #         'cdr1': self.cdr1_seq,
    #         'cdr3': self.full_cdr3,
    #         'fw2': self.framework2_seq,
    #         'cdr2': self.cdr2_seq,
    #         'fw3': self.framework3_seq,
    #         'fw4': self.framework4_seq}
    #     self.frames_and_cdrs_aa = {
    #         'fw1_aa': self.framework1_seq_aa,
    #         'cdr1_aa': self.cdr1_seq_aa,
    #         'fw2_aa': self.framework2_seq_aa,
    #         'cdr2_aa': self.cdr2_seq_aa,
    #         'fw3_aa': self.framework3_seq_aa,
    #         'cdr3_aa': self.full_cdr3_aa,
    #         'fw4_aa': self.framework4_seq_aa}

    #     self.func_seq = self.framework1_seq + self.cdr1_seq + \
    #         self.framework2_seq + self.cdr2_seq + \
    #         self.framework3_seq + self.full_cdr3
    #     self.func_seq_aa = self.framework1_seq_aa + self.cdr1_seq_aa + \
    #         self.framework2_seq_aa + self.cdr2_seq_aa + \
    #         self.framework3_seq_aa + self.full_cdr3_aa
    #     self.unmodified['functional_seq'] = str(
    #         self.func_seq).upper() + str(self.framework4_seq).upper()
    #     self.unmodified['functional_seq_aa'] = str(
    #         self.func_seq_aa) + str(self.framework4_seq_aa)

    #     if str(self.func_seq_aa).find('*') != -1 or not str(self.func_seq_aa):
    #         self.unmodified['productive'] = 'False'
    #     else:
    #         self.unmodified['productive'] = 'True'

    #     if str(self.full_cdr3_aa).find('*') != -1 or not str(self.full_cdr3_aa):
    #         self.unmodified['productive_cdr3'] = 'False'
    #     else:
    #         self.unmodified['productive_cdr3'] = 'True'

    #     for i, j in zip(self.frames_and_cdrs, self.frames_and_cdrs_aa):
    #         self.frames_and_cdrs[i] = str(self.frames_and_cdrs[i]).upper()
    #         self.frames_and_cdrs_aa[j] = str(
    #             self.frames_and_cdrs_aa[j]).upper()

    #     self.unmodified['cdr3_aa_length'] = len(str(self.full_cdr3_aa))
    #     self.unmodified['cdr1_aa_length'] = len(str(self.cdr1_seq_aa))
    #     self.unmodified['cdr2_aa_length'] = len(str(self.cdr2_seq_aa))
    #     self.unmodified['fw1_aa_length'] = len(str(self.framework1_seq_aa))
    #     self.unmodified['fw2_aa_length'] = len(str(self.framework2_seq_aa))
    #     self.unmodified['fw3_aa_length'] = len(str(self.framework3_seq_aa))
    #     self.unmodified['fw4_aa_length'] = len(str(self.framework4_seq_aa))

    #     self.unmodified['regions'] = self.frames_and_cdrs
    #     self.unmodified['regions_aa'] = self.frames_and_cdrs_aa
    #     self.unmodified['raw_seq'] = str(query_seq).upper()

    # def return_modified_dict(self):
    #     return self.unmodified
