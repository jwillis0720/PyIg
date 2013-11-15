# output checkboxes


general = [
    {'default': 1, 'formal': 'Sequence ID', 'json_key': '_id'},
    {'default': 1, 'formal': 'Input Sequence', 'json_key': 'raw_seq'},
    {'default': 1, 'formal': 'Chain Type', 'json_key': 'rearrangement.chain_type'},
    {'default': 1, 'formal': 'Format Type', 'json_key': 'domain_classification'},
    {'default': 1, 'formal': 'Query Sequence', 'json_key': 'functional_seq'},
    {'default': 1, 'formal': 'Query Translated', 'json_key': 'functional_seq_aa'},
    {'default': 1, 'formal': 'Top V Hit', 'json_key': 'rearrangement.top_v_gene_match'},
    {'default': 1, 'formal': 'Top D Hit', 'json_key': 'rearrangement.top_d_gene_match'},
    {'default': 1, 'formal': 'Top J Hit', 'json_key': 'rearrangement.top_j_gene_match'},
    {'default': 1, 'formal': 'Productive', 'json_key': 'productive'},
    {'default': 1, 'formal': 'Productive CDR3', 'json_key': 'productive_cdr3'},
    {'default': 1, 'formal': 'Strand', 'json_key': 'rearrangement.strand'}
]

nucleotide = [
    {'default': 1, 'formal': 'Framework 1 Nuc.', 'json_key': 'regions.fw1'},
    {'default': 1, 'formal': 'Framework 2 Nuc.', 'json_key': 'regions.fw2'},
    {'default': 1, 'formal': 'Framework 3 Nuc.', 'json_key': 'regions.fw3'},
    {'default': 1, 'formal': 'Framework 4 Nuc.', 'json_key': 'regions.fw4'},
    {'default': 1, 'formal': 'CDR1 Nuc.', 'json_key': 'regions.cdr1'},
    {'default': 1, 'formal': 'CDR2 Nuc.', 'json_key': 'regions.cdr2'},
    {'default': 1, 'formal': 'CDR3 Nuc.', 'json_key': 'regions.cdr3'}]

amino = [
    {'default': 1, 'formal': 'Framework 1 AA', 'json_key': 'regions_aa.fw1_aa'},
    {'default': 1, 'formal': 'Framework 2 AA', 'json_key': 'regions_aa.fw2_aa'},
    {'default': 1, 'formal': 'Framework 3 AA', 'json_key': 'regions_aa.fw3_aa'},
    {'default': 1, 'formal': 'Framework 4 AA', 'json_key': 'regions_aa.fw4_aa'},
    {'default': 1, 'formal': 'Framework 1 AA Length', 'json_key': 'fw1_aa_length'},
    {'default': 1, 'formal': 'Framework 2 AA Length', 'json_key': 'fw2_aa_length'},
    {'default': 1, 'formal': 'Framework 3 AA Length', 'json_key': 'fw3_aa_length'},
    {'default': 1, 'formal': 'Framework 4 AA Length', 'json_key': 'fw4_aa_length'},
    {'default': 1, 'formal': 'CDR1 AA', 'json_key': 'regions_aa.cdr1_aa'},
    {'default': 1, 'formal': 'CDR2 AA', 'json_key': 'regions_aa.cdr2_aa'},
    {'default': 1, 'formal': 'CDR3 AA', 'json_key': 'regions_aa.cdr3_aa'},
    {'default': 1, 'formal': 'CDR1 AA Length', 'json_key': 'cdr1_aa_length'},
    {'default': 1, 'formal': 'CDR2 AA Length', 'json_key': 'cdr2_aa_length'},
    {'default': 1, 'formal': 'CDR3 AA Length', 'json_key': 'cdr3_aa_length'}
]

total_alignments = [
    {'default': 1, 'formal': 'Total Alignment Matches', 'json_key': 'total_align.matches'},
    {'default': 1, 'formal': 'Total Alignment Mismatches', 'json_key': 'total_align.mismatches'},
    {'default': 1, 'formal': 'Total Alignment Length', 'json_key': 'total_align.length'},
    {'default': 1, 'formal': 'Total Alignment Gaps', 'json_key': 'total_align.gaps'},
    {'default': 1, 'formal': 'Total Alignment Percent Identity', 'json_key': 'total_align.percent_identity'}
]

fw1_alignments = [
    {'default': 0, 'formal': 'FW1 Alignment Matches', 'json_key': 'fr1_align.matches'},
    {'default': 0, 'formal': 'FW1 Alignment Mismatches', 'json_key': 'fr1_align.mismatches'},
    {'default': 0, 'formal': 'FW1 Alignment Length', 'json_key': 'fr1_align.length'},
    {'default': 0, 'formal': 'FW1 Alignment Gaps', 'json_key': 'fr1_align.gaps'},
    {'default': 0, 'formal': 'FW1 Identity', 'json_key': 'fr1_align.percent_identity'}
]

fw2_alignments = [
    {'default': 0, 'formal': 'FW2 Alignment Matches', 'json_key': 'fr2_align.matches'},
    {'default': 0, 'formal': 'FW2 Alignment Mismatches', 'json_key': 'fr2_align.mismatches'},
    {'default': 0, 'formal': 'FW2 Alignment Length', 'json_key': 'fr2_align.length'},
    {'default': 0, 'formal': 'FW2 Alignment Gaps', 'json_key': 'fr2_align.gaps'},
    {'default': 0, 'formal': 'FW2 Identity', 'json_key': 'fr2_align.percent_identity'}
]

fw3_alignments = [
    {'default': 0, 'formal': 'FW3 Alignment Matches', 'json_key': 'fr3_align.matches'},
    {'default': 0, 'formal': 'FW3 Alignment Mismatches', 'json_key': 'fr3_align.mismatches'},
    {'default': 0, 'formal': 'FW3 Alignment Length', 'json_key': 'fr3_align.length'},
    {'default': 0, 'formal': 'FW3 Alignment Gaps', 'json_key': 'fr3_align.gaps'},
    {'default': 0, 'formal': 'FW3 Identity', 'json_key': 'fr3_align.percent_identity'}
]

cdr1_alignments = [
    {'default': 0, 'formal': 'CDR1 Alignment Matches', 'json_key': 'cdr1_align.matches'},
    {'default': 0, 'formal': 'CDR1 Alignment Mismatches', 'json_key': 'cdr1_align.mismatches'},
    {'default': 0, 'formal': 'CDR1 Alignment Length', 'json_key': 'cdr1_align.length'},
    {'default': 0, 'formal': 'CDR1 Alignment Gaps', 'json_key': 'cdr1_align.gaps'},
    {'default': 0, 'formal': 'CDR1 Alignment Identity', 'json_key': 'cdr1_align.percent_identity'}
]

cdr2_alignments = [
    {'default': 0, 'formal': 'CDR2 Alignment Matches', 'json_key': 'cdr2_align.matches'},
    {'default': 0, 'formal': 'CDR2 Alignment Mismatches', 'json_key': 'cdr2_align.mismatches'},
    {'default': 0, 'formal': 'CDR2 Alignment Length', 'json_key': 'cdr2_align.length'},
    {'default': 0, 'formal': 'CDR2 Alignment Gaps', 'json_key': 'cdr2_align.gaps'},
    {'default': 0, 'formal': 'CDR2 Alignment Idenity', 'json_key': 'cdr2_align.percent_identity'}
]

cdr3_alignments = [
    {'default': 0, 'formal': 'CDR3 Alignment Matches', 'json_key': 'cdr3_align.matches'},
    {'default': 0, 'formal': 'CDR3 Alignment Mismatches', 'json_key': 'cdr3_align.mismatches'},
    {'default': 0, 'formal': 'CDR3 Alignment Length', 'json_key': 'cdr3_align.length'},
    {'default': 0, 'formal': 'CDR3 Alignment Gaps', 'json_key': 'cdr3_align.gaps'},
    {'default': 0, 'formal': 'CDR3 Alignment Identity', 'json_key': 'cdr3_align.percent_identity'}
]

v_hits = [
    {'default': 0, 'formal': 'V-Gene', 'json_key': 'v_hits.rank_1.subject_id'},
    {'default': 0, 'formal': 'V-Gene Mismatches', 'json_key': 'v_hits.rank_1.mismatches'},
    {'default': 0, 'formal': 'V-Gene Percent Identity', 'json_key': 'v_hits.rank_1.percent_identity'},
    {'default': 0, 'formal': 'V-Gene Percent Gaps', 'json_key': 'v_hits.rank_1.gaps'},
    {'default': 0, 'formal': 'V-Gene Percent e-Value', 'json_key': 'v_hits.rank_1.evalue'},
    {'default': 0, 'formal': 'V-Gene Percent Bit Score', 'json_key': 'v_hits.rank_1.bit_score'},
    {'default': 0, 'formal': 'V-Gene Percent Alignment Length', 'json_key': 'v_hits.rank_1.alignment_length'}
]

d_hits = [
    {'default': 0, 'formal': 'D-Gene', 'json_key': 'd_hits.rank_1.subject_id'},
    {'default': 0, 'formal': 'D-Gene Mismatches', 'json_key': 'd_hits.rank_1.mismatches'},
    {'default': 0, 'formal': 'D-Gene Percent Identity', 'json_key': 'd_hits.rank_1.percent_identity'},
    {'default': 0, 'formal': 'D-Gene Percent Gaps', 'json_key': 'd_hits.rank_1.gaps'},
    {'default': 0, 'formal': 'D-Gene Percent e-Value', 'json_key': 'd_hits.rank_1.evalue'},
    {'default': 0, 'formal': 'D-Gene Percent Bit Score', 'json_key': 'd_hits.rank_1.bit_score'},
    {'default': 0, 'formal': 'D-Gene Percent Alignment Length', 'json_key': 'd_hits.rank_1.alignment_length'}
]

j_hits = [
    {'default': 0, 'formal': 'J-Gene', 'json_key': 'j_hits.rank_1.subject_id'},
    {'default': 0, 'formal': 'J-Gene Mismatches', 'json_key': 'j_hits.rank_1.mismatches'},
    {'default': 0, 'formal': 'J-Gene Percent Identity', 'json_key': 'j_hits.rank_1.percent_identity'},
    {'default': 0, 'formal': 'J-Gene Percent Gaps', 'json_key': 'j_hits.rank_1.gaps'},
    {'default': 0, 'formal': 'J-Gene Percent e-Value', 'json_key': 'j_hits.rank_1.evalue'},
    {'default': 0, 'formal': 'J-Gene Percent Bit Score', 'json_key': 'j_hits.rank_1.bit_score'},
    {'default': 0, 'formal': 'J-Gene Percent Alignment Length', 'json_key': 'j_hits.rank_1.alignment_length'}
]

all_checkboxes = [general + nucleotide + amino + total_alignments + fw1_alignments + fw2_alignments + fw3_alignments +
                  cdr1_alignments + cdr2_alignments + cdr3_alignments + v_hits + d_hits + j_hits]
all_checkboxes_dict = {
    'general': general,
    'nucleotide': nucleotide,
    'amino': amino,
    'total_alignments': total_alignments,
    'fw1_alignments': fw1_alignments,
    'fw2_alignments': fw2_alignments,
    'fw3_alignments': fw3_alignments,
    'cdr1_alignments': cdr1_alignments,
    'cdr2_alignments': cdr2_alignments,
    'cdr3_alignments': cdr3_alignments,
    "v_hits": v_hits,
    "d_hits": d_hits,
    "j_hits": j_hits
}
