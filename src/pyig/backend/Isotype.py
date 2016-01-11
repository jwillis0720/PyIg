from Bio import SeqIO
from Bio import pairwise2


def get_isotype(sequence):
    seq = sequence
    alignments = {}
    alignments['IgG_Rev'] = pairwise2.align.localms(str(seq), "GCCTCCACCAAGGGCCCATCG", 2, -1, -1, -1)
    alignments['IgG_Rev1'] = pairwise2.align.localms(str(seq), "GCCTCCACCAAGGGCCCATCT", 2, -1, -1, -1)
    alignments['IgG_Rev2'] = pairwise2.align.localms(str(seq), "GCTTCCACCAAGGGCCCATCG", 2, -1, -1, -1)
    alignments['IgG_Rev3'] = pairwise2.align.localms(str(seq), "GCTTCCACCAAGGGCCCATCT", 2, -1, -1, -1)
    alignments['IgM_Rev'] = pairwise2.align.localms(str(seq), "GGAGTGCATCCGCCCCAACC", 2, -1, -1, -1)

    best = max(alignments.iterkeys(), key=(lambda key: alignments[key][0][2]))
    if alignments[best][0][2] < 33:
        return "Ambiguous"
    isotype = best.split('_')[0]
    return isotype
