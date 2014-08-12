import os
import os.path
import sys
from tempfile import NamedTemporaryFile
from Bio import SeqIO


def replace_non_ascii(string_characters):
    acceptable = '-_0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    for i in set(filter(lambda x: x not in acceptable, string_characters)):
        s = string_characters.replace(i, '_')
        string_characters = s
    return string_characters


def split_fasta(num_procs, fasta_file, suffix=".tmp_fasta", delete=False):
    '''Split the file name by the number of processors you specify
     arguments -
     num_procs - The amount of processors you are running
     suffix - gives the fasta files an identifiable suffix to let you know they are  split
     '''

    # return all fasta so we can get the raw sequence from it which blast does not provide.
    parent_file = []
    print "Counting entries in fasta files {0}".format(os.path.abspath(fasta_file))

    for i, j in enumerate(SeqIO.parse(fasta_file, 'fasta')):
        if i % 10000 == 0 and i != 0:
            print "coutned {0} entries".format(i)
        parent_file.append(j)
    print "{0} entry in fasta file".format(len(parent_file))

    files_per_temporary_file = float(len(parent_file)) / float(num_procs)
    list_of_temporary_files = []

    # carries all of our records to be written
    joiner = []
    num = 1

    # enumerate huge fasta
    for record in parent_file:
        # append records to our list holder
        joiner.append(">" + replace_non_ascii(record.id) + "\n" + str(record.seq))
        # if we have reached the maximum numbers to be in that file, write
        if num > files_per_temporary_file:
            joiner.append("")
            file_name = NamedTemporaryFile(suffix=suffix, delete=delete)
            with open(file_name.name, 'w') as f:
                f.write("\n".join(joiner))
            list_of_temporary_files.append(file_name.name)
            joiner = []
            num = 1
        else:
            num += 1

    # if joiner still has stuff in it
    if joiner:
        # for left over fasta entries, very important or else they will
        # just hang out in limbo
        joiner.append("")
        file_name = NamedTemporaryFile(suffix=suffix, delete=delete)
        with open(file_name.name, 'w') as f:
            f.write("\n".join(joiner))
        list_of_temporary_files.append(file_name.name)
    return list_of_temporary_files

if __name__ == '__main__':
    blash = split_fasta(sys.argv[1], os.path.dirname(os.path.abspath(sys.argv[2])), sys.argv[2])
