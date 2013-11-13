# output checkboxes


general = [
    {'default': 1, 'formal': 'Sequence ID', 'json_key': '_id'},
    {'default': 0, 'formal': 'Input Sequence', 'json_key': 'raw_seq'},
    {'default': 1, 'formal': 'Chain Type', 'json_key': 'rearrangement.chain_type'},
    {'default': 1, 'formal': 'Format Type', 'json_key': 'domain_classification'},
    {'default': 0, 'formal': 'Query Sequence', 'json_key': 'full_seq'},
    {'default': 0, 'formal': 'Query Translated', 'json_key': 'full_seq_aa'},
    {'default': 1, 'formal': 'Top V Hit', 'json_key': 'rearrangement.top_v_gene_match'},
    {'default': 1, 'formal': 'Top D Hit', 'json_key': 'rearrangement.top_d_gene_match'},
    {'default': 1, 'formal': 'Top J Hit', 'json_key': 'rearrangement.top_j_gene_match'},
    {'default': 0, 'formal': 'Productive', 'json_key': 'productive'},
    {'default': 0, 'formal': 'Productive CDR3', 'json_key': 'productive_cdr3'}
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
    {'default': 1, 'formal': 'CDR1 AA', 'json_key': 'regions_aa.cdr1_aa'},
    {'default': 1, 'formal': 'CDR2 AA', 'json_key': 'regions_aa.cdr2_aa'},
    {'default': 1, 'formal': 'CDR3 AA', 'json_key': 'regions_aa.cdr3_aa'}
]

all_checkboxes = {
    'general': general, 'nucleotide': nucleotide, 'amino': amino,
    'v_genes': {}, 'd_genes': {}, 'j_genes': {}, 'junctional': {}
}
