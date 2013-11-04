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
    -d --db_path 
      The database path from igblast
    -i --internal_data 
      The internal database from igblast
      
    pyigblast -q query.fasta -d database/ -i internal/

Recommended Arguments
---------

    -a --aux_path 
        The auxiliray path that contains the frame origins of the germline gene
    -x --executable 
        If igblastn is not in /usr/bin , specifiy path to executable
        
    pyigblast -q query.fasta -d database/ -i internal -a optional/ -x /bin/igblastn
    
BLAST Specfic Arguments
---------
Have defaults

    -o --out 
        The output filename
    -e --e_value 
        Expectation value threshold for blast.
    -w --word_size 
        Word size for wordmatch algorithm
    -pm --penalty_mistmatch
        Penalty for nucleotide mismatch
    -rm --reward_match
        Reward for nucleotide match
    
    pyigblast -q query.fasta -d database/ -i internal -a optional/ -x /bin/igblastn -o blast_output -e 1E-15 
    -w 4 -pm 3 -rm 5

IgBlast Specific Arguments
--------
Have defaults

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
    -sT --show_translation
        Show translation of alignment
     
     pyigblast -q query.fasta -d database/ -i internal -a optional/ -x /bin/igblastn -o blast_output -e 1E-15 
     -w 4 -pm 3 -rm 5 -or human -nV 2 -nD 2 -nJ 1 -dgm 5 -s imgt -sT

Formatting
----------
Have defaults

    -f --format_options
        Default is a tab seperated format of :
            qseqid sseqid pident length mismatch gapopen qstart qend sstart send
        Or:
            Pass format_file.txt and uncomment out metrics to return
    -z --zip
        Zip up all output files
    -c --concatenate
        Turn off automatic concatenation of result files. Pyigblast splits up files across processors, 
        if you want them to be put back together.
    -j, --json
        Use the JSON output option that will format the text driven igblast output to a json document
    -jp --json_prefix
        The prefix for json_output files
    
    pyigblast -q query.fasta -d database/ -i internal -a optional/ -x /bin/igblastn -o blast_output -e 1E-15 
     -w 4 -pm 3 -rm 5 -or human -nV 2 -nD 2 -nJ 1 -dgm 5 -s imgt -sT -f format_templat.txt -z -c -j 
     -jp my_json_file.txt

Formatting File
---------------
The format file, just comment out to remove from output:

     qseqid # Query Seq-id
    #qgi # Query GI
    #qacc # Query accesion
    #qaccver # Query accesion.version
    qlen # Query sequence length
    sseqid # Subject Seq-id
    #sallseqid # All subject Seq-id(s), separated by a ';'
    #sgi # Subject GI
    #sallgi # All subject GIs
    sacc # Subject accession
    #saccver # Subject accession.version
    #sallacc # All subject accessions
    slen # Subject sequence length
    qstart # Start of alignment in query
    qend # End of alignment in query
    sstart # Start of alignment in subject
    send # End of alignment in subject
    qseq # Aligned part of query sequence
    sseq # Aligned part of subject sequence
    evalue # Expect value
    bitscore # Bit score
    score # Raw score
    length # Alignment length
    pident # Percentage of identical matches
    nident # Number of identical matches
    mismatch # Number of mismatches
    positive # Number of positive-scoring matches
    gapopen # Number of gap openings
    gaps # Total number of gaps
    ppos # Percentage of positive-scoring matche
    #frames # Query and subject frames separated by a '/'
    qframe # Query frame
    sframe # Subject frame
    #btop # Blast traceback operations (BTOP)
    #staxids # unique Subject Taxonomy ID(s), separated by a ';'(in numerical order)
    #sscinames # unique Subject Scientific Name(s), separated by a ';'
    #scomnames # unique Subject Common Name(s), separated by a ';'
    #sblastnames # unique Subject Blast Name(s), separated by a ';'(in alphabetical order)
    #sskingdoms # unique Subject Super Kingdom(s), separated by a ';'(in alphabetical order) 
    #stitle # Subject Title
    #salltitles # All Subject Title(s), separated by a '<>'
    #sstrand # Subject Strand
    #qcovs # Query Coverage Per Subject
    #qcovhsp # Query Coverage Per HSP
