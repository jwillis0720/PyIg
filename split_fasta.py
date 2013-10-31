import Bio.SeqIO

def split_fasta(num_procs,file_name,suffix=".tmp_fasta"):
        '''Split the file name by the number of processors you specify

        arguments (num_procs - The amount of processors you are running, lets the function know how much to split up the file into with one file per processor)'''

        all_fasta = {}
        file_prefix = file_name.split('.fasta')[0]
        print "Counting entries in fasta files..."
        parent_file = list(Bio.SeqIO.parse(file_name, 'fasta'))
        length_parent_file = len(parent_file)
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
                with open(file_prefix + str(file_counter) + suffix, 'w') as f:
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
            with open(file_prefix + str(file_counter) + suffix, 'w') as f:
                f.write("\n".join(joiner))
        return all_fasta
