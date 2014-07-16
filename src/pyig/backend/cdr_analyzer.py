import sys
import warnings

warnings.filterwarnings('ignore')

try:
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC
except ImportError:
    print("Need Biopython to use the IgBlast output parser class")


class cdr_analyzer():

    def __init__(self, initial_blast_dictionary, end_translation_dictionaries):
        query = initial_blast_dictionary.keys()[0]

        self.unmodified = initial_blast_dictionary[query]
        self.end_translation_dictionaries = end_translation_dictionaries

        # the raw sequence we want to query
        self.unmodified['_id'] = query
        query_seq = self.unmodified['raw_seq']

        # framework1
        self.frames_and_cdrs = {}
        self.frames_and_cdrs_aa = {}

        # framework1
        self.framework1_frame_info = self.unmodified['fr1_align']
        try:
            self.framework1_seq = Seq(
                query_seq[self.framework1_frame_info['from'] - 1:self.framework1_frame_info['to']], IUPAC.ambiguous_dna)
            self.framework1_seq_aa = self.framework1_seq.translate()
        except (TypeError, KeyError):
            self.framework1_seq = ""
            self.framework1_seq_aa = ""

        # framework2
        self.framework2_frame_info = self.unmodified['fr2_align']
        try:
            self.framework2_seq = Seq(
                query_seq[self.framework2_frame_info['from'] - 1:self.framework2_frame_info['to']], IUPAC.ambiguous_dna)
            self.framework2_seq_aa = self.framework2_seq.translate()
        except (TypeError, KeyError):
            self.framework2_seq = ""
            self.framework2_seq_aa = ""

        # framework3
        self.framework3_frame_info = self.unmodified['fr3_align']
        try:
            self.framework3_seq = Seq(
                query_seq[self.framework3_frame_info['from'] - 1:self.framework3_frame_info['to']], IUPAC.ambiguous_dna)
            self.framework3_seq_aa = self.framework3_seq.translate()
        except (TypeError, KeyError):
            self.framework3_seq = ""
            self.framework3_seq_aa = ""

        # cdr1
        self.cdr1_frame_info = self.unmodified['cdr1_align']
        try:
            self.cdr1_seq = Seq(
                query_seq[self.cdr1_frame_info['from'] - 1:self.cdr1_frame_info['to']], IUPAC.ambiguous_dna)
            self.cdr1_seq_aa = self.cdr1_seq.translate()
        except (TypeError, KeyError):
            self.cdr1_seq = ""
            self.cdr1_seq_aa = ""

        # cdr2
        self.cdr2_frame_info = self.unmodified['cdr2_align']
        try:
            self.cdr2_seq = Seq(
                query_seq[self.cdr2_frame_info['from'] - 1:self.cdr2_frame_info['to']], IUPAC.ambiguous_dna)
            self.cdr2_seq_aa = self.cdr2_seq.translate()
        except (TypeError, KeyError):
            self.cdr2_seq = ""
            self.cdr2_seq_aa = ""

        # cdr3 v_part
        self.v_part_of_cdr3_info = self.unmodified['cdr3_align']
        try:
            self.v_part_of_cdr3 = Seq(
                query_seq[self.v_part_of_cdr3_info['from'] - 1:self.v_part_of_cdr3_info['to']])
        except (KeyError, TypeError):
            self.v_part_of_cdr3 = ""

        # junctional info
        junctions = self.unmodified['junction']
        try:
            self.v_d_junction = Seq(junctions['v-d_junction'])
        except (TypeError, KeyError):
            self.v_d_junction = ""

        try:
            self.d_region = Seq(junctions['d_region'])
        except (TypeError, KeyError):
            self.d_region = ""

        try:
            self.d_j_junction = Seq(junctions['d-j_junction'])
        except (TypeError, KeyError):
            self.d_j_junction = ""

        try:
            self.d_or_j_junction = Seq(junctions['d-or-j_junction'])
        except (TypeError, KeyError):
            self.d_or_j_junction = ""

        '''Now the J region, which will be the beginning of the V region of the
        CDR3 to the end of the translated region of the J gene, the translated region
        is not implicitly coded by igblastn so we have to manually look at it from a file'''
        # j_region
        try:
            self.j_region_info = self.unmodified['j_hits']['rank_1']
            # where does the query start on the j gene
            _begin = int(self.j_region_info['q._start']) - 1
            # name of jgene
            self.j_gene = self.j_region_info['subject_id']

            # where does the CDR3 loop stop in the j_gene stop translating for
            # this particular gene
            self.j_gene_end = self.end_translation_dictionaries[self.j_gene]

            # and make it relative to the sequence read
            _end = (self.j_gene_end - int(
                self.j_region_info['s._start'])) + _begin

            # j_region
            self.j_region = Seq(query_seq[_begin:_end])
            self.full_cdr3 = Seq(
                query_seq[self.v_part_of_cdr3_info['from'] - 1:_end])
            self.full_cdr3_aa = self.full_cdr3.translate()

            # get the rest of the sequence until the end
            self.framework4_seq = Seq(query_seq[_end - 2:])
            self.framework4_seq_aa = self.framework4_seq.translate()
        except (TypeError, KeyError, IndexError):
            self.j_gene = ""
            self.j_region = ""
            self.full_cdr3 = ""
            self.full_cdr3_aa = ""
            self.framework4_seq = ""
            self.framework4_seq_aa = ""

        self.frames_and_cdrs = {
            'fw1': self.framework1_seq,
            'cdr1': self.cdr1_seq,
            'cdr3': self.full_cdr3,
            'fw2': self.framework2_seq,
            'cdr2': self.cdr2_seq,
            'fw3': self.framework3_seq,
            'fw4': self.framework4_seq}
        self.frames_and_cdrs_aa = {
            'fw1_aa': self.framework1_seq_aa,
            'cdr1_aa': self.cdr1_seq_aa,
            'fw2_aa': self.framework2_seq_aa,
            'cdr2_aa': self.cdr2_seq_aa,
            'fw3_aa': self.framework3_seq_aa,
            'cdr3_aa': self.full_cdr3_aa,
            'fw4_aa': self.framework4_seq_aa}

        self.func_seq = self.framework1_seq + self.cdr1_seq + \
            self.framework2_seq + self.cdr2_seq + \
            self.framework3_seq + self.full_cdr3
        self.func_seq_aa = self.framework1_seq_aa + self.cdr1_seq_aa + \
            self.framework2_seq_aa + self.cdr2_seq_aa + \
            self.framework3_seq_aa + self.full_cdr3_aa
        self.unmodified['functional_seq'] = str(
            self.func_seq).upper() + str(self.framework4_seq).upper()
        self.unmodified['functional_seq_aa'] = str(
            self.func_seq_aa) + str(self.framework4_seq_aa)

        if str(self.func_seq_aa).find('*') != -1 or not str(self.func_seq_aa):
            self.unmodified['productive'] = 'False'
        else:
            self.unmodified['productive'] = 'True'

        if str(self.full_cdr3_aa).find('*') != -1 or not str(self.full_cdr3_aa):
            self.unmodified['productive_cdr3'] = 'False'
        else:
            self.unmodified['productive_cdr3'] = 'True'

        for i, j in zip(self.frames_and_cdrs, self.frames_and_cdrs_aa):
            self.frames_and_cdrs[i] = str(self.frames_and_cdrs[i]).upper()
            self.frames_and_cdrs_aa[j] = str(
                self.frames_and_cdrs_aa[j]).upper()

        self.unmodified['cdr3_aa_length'] = len(str(self.full_cdr3_aa))
        self.unmodified['cdr1_aa_length'] = len(str(self.cdr1_seq_aa))
        self.unmodified['cdr2_aa_length'] = len(str(self.cdr2_seq_aa))
        self.unmodified['fw1_aa_length'] = len(str(self.framework1_seq_aa))
        self.unmodified['fw2_aa_length'] = len(str(self.framework2_seq_aa))
        self.unmodified['fw3_aa_length'] = len(str(self.framework3_seq_aa))
        self.unmodified['fw4_aa_length'] = len(str(self.framework4_seq_aa))

        self.unmodified['regions'] = self.frames_and_cdrs
        self.unmodified['regions_aa'] = self.frames_and_cdrs_aa
        self.unmodified['raw_seq'] = str(query_seq).upper()

    def return_modified_dict(self):
        return self.unmodified
