PyIgBlast
=========

PyIgBlast - Open source parser to call IgBlast and parse results for high-throughput sequencing.
Uses Python multi-processing to get around bottlenecks of IgBlast multi-threading.
Parses blast output to deliminated files (csv,json) for uploading to databases.

#ONLY BENCHMARKED WITH HUMAN HEAVY SO FAR

Requires
=========

1.   >Python2.7
2.   [BioPython - Python tools for biological computations](http://biopython.org/wiki/Download)
3.   [IgBlastn - BLAST algorithm for analyzing immunoglobulin repertoires](ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/)

Usage
========

    pyigblast -h

Required Arguments
--------

    -q --query
      The fasta file to blast

    ex.
    pyigblast -q query.fasta -d database/ -i internal/


Database Paths
---------

These options are not required because they default to the databases under the included datafiles directory

    -d --db_path
        The database path to the germline repertoire
    -i --internal_data
        The database path to the internal data repertoire
    -a --aux_path
        The auxilary path to the optional database

    ex.
    pyigblast -q query.fasta -d datafiles/database/ -i datafiles/internal -a datafiles/optional_files/


IgBlast Specific Arguments
--------
Arguments for the IgBLAST algorithm, all of which have defaults

    -or --organism
        The organism immunoglobulin repertoire to blast against
    -nV --num_v
        How many V-genes to return
    -nD --num_d
        How many D-genes to return
    -nJ --num_j
        How many J-genes to return
    -dgm --d_gene_matches
        How many nucleotides in the D-gene must match to call it a hit
    -s --domain
        Classification system to use

     ex.
     pyigblast -q query.fasta -d datafiles/database/ -i datafiles/internal -a datafiles/optional_files/ -or human -nV 2 -nD 2 -nJ 1 -dgm 5 -s imgt


General Settings
----------
Arguments for the PyIg parser. All have defaults.

    -x --executable
        The location of the executable if it is not under /usr/bin/igblastn
    -o --out
        The outfile prefix
    -t --tmp
        The directory that we will store temporary files in like split fastas and blast output
    -e --evalue
        Real value for the expectation value threshold in BLAST.
    -w --word_size
        Word size for the BLAST algorithm
    -pm --penalty_mismatch
        The mismatch penalty in the Dgene
    -nP --num_procs
        How many processors to use.

    ex.
    PyIg -q query.fasta -d datafiles/database/ -i datafiles/internal -a datafiles/optional_files/ -or human -nV 2 -nD 2 -nJ 1 -dgm 5 -s imgt -x /usr/local/bin/igblastn -o Myprettyoutput -t /tmp/temp -e 1E15 -w 4 -pm -5 -np 3


Output Options
----------
The output options. Defaults to a csv file with common output options

    -op --output_options
        Open this file and comment out options you don't want in your final output. The output option default are noncommented out under datafiles/output_options.txt Don't change the second column
    -z --zip
        Zip up all output files
    -c --concatenate
        Turn off automatic concatenation of result files. Pyigblast splits up files across processors,
        if you want them to be put back together.
    -j, --json
        Use the JSON output option that will format the text driven igblast output to a json document

    ex.
    PyIg -q query.fasta -d datafiles/database/ -i datafiles/internal -a datafiles/optional_files/ -x /bin/igblastn -o Myspiffyoutput -e 1E-15
     -w 4 -pm 3 -rm 5 -or human -nV 2 -nD 2 -nJ 1 -dgm 5 -s imgt -op datafiles/output_options.txt -z -c -j

Output Options File
---------------
The content of the output options file, just comment out to remove from output, again don't change the text:

    #Name,internal parser name
    Sequence ID,_id
    Input Sequence,raw_seq
    Chain Type,rearrangement.chain_type
    Format Type,domain_classification
    Query Sequence,functional_seq
    Query Translated,functional_seq_aa
    Top V Hit,rearrangement.top_v_gene_match
    Top D Hit,rearrangement.top_d_gene_match
    Top J Hit,rearrangement.top_j_gene_match
	Productive,productive
	Productive CDR3,productive_cdr3
	Strand,rearrangement.strand
	Framework 1 Nuc.,regions.fw1
	Framework 2 Nuc.,regions.fw2
	Framework 3 Nuc.,regions.fw3
	Framework 4 Nuc.,regions.fw4
	CDR1 Nuc.,regions.cdr1
	CDR2 Nuc.,regions.cdr2
	CDR3 Nuc.,regions.cdr3
	Framework 1 AA,regions_aa.fw1_aa
	Framework 2 AA,regions_aa.fw2_aa
	Framework 3 AA,regions_aa.fw3_aa
	Framework 4 AA,regions_aa.fw4_aa
	Framework 1 AA Length,fw1_aa_length
	Framework 2 AA Length,fw2_aa_length
	Framework 3 AA Length,fw3_aa_length
	Framework 4 AA Length,fw4_aa_length
	CDR1 AA,regions_aa.cdr1_aa
	CDR2 AA,regions_aa.cdr2_aa
	CDR3 AA,regions_aa.cdr3_aa
	CDR1 AA Length,cdr1_aa_length
	CDR2 AA Length,cdr2_aa_length
	CDR3 AA Length,cdr3_aa_length
	Total Alignment Matches,total_align.matches
	Total Alignment Mismatches,total_align.mismatches
	Total Alignment Length,total_align.length
	Total Alignment Gaps,total_align.gaps
	Total Alignment Percent Identity,total_align.percent_identity
	#FW1 Alignment Matches,fr1_align.matches
	#FW1 Alignment Mismatches,fr1_align.mismatches
	#FW1 Alignment Length,fr1_align.length
	#FW1 Alignment Gaps,fr1_align.gaps
	#FW1 Identity,fr1_align.percent_identity
	#FW2 Alignment Matches,fr2_align.matches
	#FW2 Alignment Mismatches,fr2_align.mismatches
	#FW2 Alignment Length,fr2_align.length
	#FW2 Alignment Gaps,fr2_align.gaps
	#FW2 Identity,fr2_align.percent_identity
	#FW3 Alignment Matches,fr3_align.matches
	#FW3 Alignment Mismatches,fr3_align.mismatches
	#FW3 Alignment Length,fr3_align.length
	#FW3 Alignment Gaps,fr3_align.gaps
	#FW3 Identity,fr3_align.percent_identity
	#CDR1 Alignment Matches,cdr1_align.matches
	#CDR1 Alignment Mismatches,cdr1_align.mismatches
	#CDR1 Alignment Length,cdr1_align.length
	#CDR1 Alignment Gaps,cdr1_align.gaps
	#CDR1 Alignment Identity,cdr1_align.percent_identity
	#CDR2 Alignment Matches,cdr2_align.matches
	#CDR2 Alignment Mismatches,cdr2_align.mismatches
	#CDR2 Alignment Length,cdr2_align.length
	#CDR2 Alignment Gaps,cdr2_align.gaps
	#CDR2 Alignment Idenity,cdr2_align.percent_identity
	#CDR3 Alignment Matches,cdr3_align.matches
	#CDR3 Alignment Mismatches,cdr3_align.mismatches
	#CDR3 Alignment Length,cdr3_align.length
	#CDR3 Alignment Gaps,cdr3_align.gaps
	#CDR3 Alignment Identity,cdr3_align.percent_identity
	#V-Gene,v_hits.rank_1.subject_id
	#V-Gene Mismatches,v_hits.rank_1.mismatches
	#V-Gene Percent Identity,v_hits.rank_1.percent_identity
	#V-Gene Gaps,v_hits.rank_1.gaps
	#V-Gene e-Value,v_hits.rank_1.evalue
	#V-Gene Bit Score,v_hits.rank_1.bit_score
	#V-Gene Alignment Length,v_hits.rank_1.alignment_length
	#D-Gene,d_hits.rank_1.subject_id
	#D-Gene Mismatches,d_hits.rank_1.mismatches
	#D-Gene Percent Identity,d_hits.rank_1.percent_identity
	#D-Gene Gaps,d_hits.rank_1.gaps
	#D-Gene e-Value,d_hits.rank_1.evalue
	#D-Gene Bit Score,d_hits.rank_1.bit_score
	#D-Gene Alignment Length,d_hits.rank_1.alignment_length
	#J-Gene,j_hits.rank_1.subject_id
	#J-Gene Mismatches,j_hits.rank_1.mismatches
	#J-Gene Percent Identity,j_hits.rank_1.percent_identity
	#J-Gene Gaps,j_hits.rank_1.gaps
	#J-Gene e-Value,j_hits.rank_1.evalue
	#J-Gene Bit Score,j_hits.rank_1.bit_score
	#J-Gene Alignment Length,j_hits.rank_1.alignment_length
