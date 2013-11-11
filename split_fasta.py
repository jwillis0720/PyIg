import Bio.SeqIO
import os
import os.path
import sys


def split_fasta(num_procs, path, file_name, suffix=".tmp_fasta"):
    '''Split the file name by the number of processors you specify
     arguments -
     num_procs - The amount of processors you are running
     path - the output path you want  the files dumped in
     suffix - gives the fasta files an identifiable suffix to let you know they are  split
     '''

    # return all fasta so we can get the raw sequence from it which blast does not provide.
    #@todo, find a way to get the memory down in these functions
    all_fasta = {}
    if not os.path.exists(path):
        os.makedirs(path)

    length_parent_file = 0
    parent_file = []
    file_prefix = os.path.basename(file_name).split('.fasta')[0]
    print "Counting entries in fasta files..."
    for i, j in enumerate in Bio.SeqIO.parse(file_name, 'fasta'):
        if i % 10000 == 0:
            print "coutned {0} entries".format(i)
        length_parent_file += i
        parent_file.append(j)
    print "{0} in fasta file".format(length_parent_file)

    files_per_tmp = float(length_parent_file) / float(num_procs)
    print "{0} processors, blasting {1} entries per processor".format(num_procs, files_per_tmp)

    # carries all of our records to be written
    joiner = []
    file_counter = 1
    num = 1

    # enumerate huge fasta
    for record in parent_file:
        # append records to our list holder
        joiner.append(">" + record.id + "\n" + str(record.seq))
        all_fasta[record.id] = str(record.seq)
        # if we have reached the maximum numbers to be in that file, write
        # to a file, and then clear
        if num > files_per_tmp:
            joiner.append("")
            with open(path + file_prefix + str(file_counter) + suffix, 'w') as f:
                f.write("\n".join(joiner))
            # change file name,clear record holder, and change the file
            # count
            joiner = []
            file_counter += 1
            num = 1
        else:
            num += 1
    if joiner:
        # for left over fasta entries, very important or else they will
        # just hang out in limbo
        joiner.append("")
        with open(path + file_prefix + str(file_counter) + suffix, 'w') as f:
            f.write("\n".join(joiner))
    return all_fasta

if __name__ == '__main__':
    split_fasta(sys.argv[1], os.path.dirname(os.path.abspath(sys.argv[2])), sys.argv[2])
