* PyIg - Python Immunoglobulin Blast

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but without any warranty.

* Third-Party Software Used by PyIg

IgBLASTN - Nucleotide-Nucleotide BLAST for immunoglobulin sequences 2.2.28+

Released under public domain as "United Stated Government work" under terms of the United States
copyright act. It can not be copyrighted. See:

Ye, J, Ma N, Madden TL, Ostell JM, 2013 Nuc. Acid Research IgBLAST: an immunoglobulin variable domain sequence analysis tool

-------------------------------------------------------------------------------------------------------------

*Basic usage

All that is required is a FASTA file containing any number of sequences. The rest of the arguments are defaulted detailed as followed. It will output one file in either JSON or CSV formats. The blast output is raw igblastn output and is not parsed.

*Input FASTA file

Can be one or many fasta entries to blast across. Will get one row or entry in the output for every sequence submitted.

*Databases

The database directories are used by blast and the parsing scripts to scan junctions and assign germline genes.

-Compiled BLAST Directory

The BLAST formatted germline files to BLAST against. To use your own, you would need to use 'makeblastdb' from the BLAST software suite to format sequences for BLAST.

-Internal Data Directory

Used by BLAST to assign frameworks and cdrs for each gene.

-Auxillary

Contains coding frames for the JH gnes. Necessary for BLAST translation, however, the end parser puts together the junctions and retranslates. So it's not necessary.

*Options

-Scheme

Should the numbering and selection scheme use IMGT or KABAT.

-Species

Which species to BLAST against. No need to change the database, just give this option and it will take care of it for you.

-Output format

JSON (JavaScript Object Notation) - This output format is useful if you want to upload to databases like MONGO. It contains all the fields and values within the file. If you don't know that you need this, just use CSV.

CSV (Comma seperated values) - All the fields are on the first line, followed by the data. This can be opened in EXCEL or uploaded to a relational database like mySQL.

BLAST output - This is just the raw output from blast if you wanted it. Basically gives junctional regions but is not parsed or organized.

-VDJ

How many V,D,J genes to show in each output. A number > 1 is only relevant to JSON output since it is the only thing that can handle multiple entries of the same object (nested). If you choose 3, it will rank them by sequence identity, and put them in order in the JSON. If you use CSV, it will only give back rank 1.

-BLAST Options

e-value - The minimal expected value to use for BLAST to consider a gene a match. Ex. For 10, you would see 10 matches in the db for a given score just by chance.

word size - BLAST breaks up the search into "words", what word size should it be used.

penalty mismatch - The penalty of a score for having a mismatch

minimal D nucleotides - The minimal number of matching D nucleotides needed for BLAST to consider it a match.

processors - How many processors to use for blasting. PyIg will break up the input fasta and submit it to each processor. It automatically garbage collects

zip files - For large files, should they be zipped up at the end

**
http://stackoverflow.com/questions/3333334/stdout-to-tkinter-gui