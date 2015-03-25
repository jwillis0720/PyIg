PyIg
=========

PyIg - Immunoglobulin and T-Cell receptor rearrangment software. It uses IgBLAST to call on V(D) and J genes. Then it recombines them to JSON format. PyIg is meant to be highly parralizable, so it uses multiple processors to call on multiple instances of BLAST making it ideal for high throughput sequencing. The JSON documents can be reconfigured to be uploaded to MySQL or NoSQL databases for efficient queries. In addition, you can use custom databases for V(D) and J gene lookup.

Requires
=========

1.   >Python2.7
2.   [BioPython - Python tools for biological computations](http://biopython.org/wiki/Download)


Installation
=========

    sudo python setup.py install

Usage
========

    /usr/local/bin/PyIg -h

Required Positional Arguments
--------


The fasta file to be run through PyIg.

    query.fasta

    ex.
    PyIg -q query.fasta


Database
--------

The setup script should of placed this in /usr/local/pyig/. If you have it located somewhere else (not recommended), use this option.

    -d DATABASE, --database DATABASE



Types
--------

Arguments for the database to query

    -r {Ig,TCR}, --receptor {Ig,TCR}
                    The receptor you are analyzing, immunoglobulin or T cell Receptor - defaults to Ig
    -s {human,rabbit,mouse,rat,rhesus}, --species {human,rabbit,mouse,rat,rhesus}
                    The Species you are analyzing, defaults to human
    -c {heavy,light, mixed}, --chain {heavy,light, mixed}
                    The chain you want to analyze, defaults to heavy

                    ex.

                    /usr/local/bin/PyIg -s mouse -c light mouse_light_chains.fasta

BLAST Specific Arguments:
--------

 Arguments Specific to IgBlast

    -nV NUM_V_ALIGNMENTS, --num_V_alignments NUM_V_ALIGNMENTS
                         How many V genes do you want to match?
    -nD NUM_D_ALIGNMENTS, --num_D_alignments NUM_D_ALIGNMENTS
                         How many D genes do you want to match?, does not apply for kappa and lambda
    -nJ NUM_J_ALIGNMENTS, --num_J_alignments NUM_J_ALIGNMENTS
                         How many J genes do you want to match?
    -mD MIND, --minD MIND
                        The amount of nucleotide matches needed for a D gene match. >= 5 right now
    -x EXECUTABLE, --executable EXECUTABLE
                         The location of IGBlastn, default is /usr/local/bin/igblastn if used setup.py



General Arguments
--------

Output and Miscellaneous Arguments

     -m MULTI, --multi MULTI
          Multiprocess by the amount of CPUs you have. Or you can enter a number or type 0 to turn it off
     -o OUT, --out OUT     Output_file_name
     --debug Debug mode, this will not delete the temporary blast files and will print some other useful things, like which regions did not parse


For details on using the database, creating your own database, and developing, see the [Wiki](https://github.com/jwillis0720/PyIg/wiki).
