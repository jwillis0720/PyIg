PyIgBlast
=========

PyIgBlast - Open source parser to call IgBlast and parse results for high-throughput sequencing. 
Uses Python multi-processing to get around bottlenecks of IgBlast multi-threading. 
Parses blast output to deliminated files (csv,json) for uploading to databases. 


Requires
=========

1.   >Python2.7
2.   [BioPython - Python tools for biological computations](http://biopython.org/wiki/Download)
3.   [IgBlastn - BLAST algorithm for analyzing immunoglobulin repertoires](ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/)

Usage
========
    
    pyigblast -h 

Required
--------

    -q --query 
      The fasta file to blast
    -d --db_path 
      The database path from igblast
    -i --internal_data 
      The internal database from igblast
      
    pyigblast -q query.fasta -d database/ -i internal/
















Major commits
=========

8/6/2013 Full functionality with multiprocessing and a formatting template

8/7/2013 Can run on blue cluster

8/9/2013 Parses output from igblast into python classes

8/16/2013 Parses to JSON, zips, and cleans up directories
