
from Bio import SeqIO as so
import sys
def write_file(parent_file,how_many):
        #get file_counter and base name of fasta_file
        parent_file_base_name = parent_file.split(".")[0]
        counter = 1

        #our first file name
        file = parent_file_base_name + "_" + str(counter) + ".fasta"

        #carries all of our records to be written
        joiner = []
        #enumerate huge fasta
        for num,record in enumerate(so.parse(parent_file, "fasta"),start=1):

            #append records to our list holder
            joiner.append(">" + record.id + "\n" + str(record.seq))

            #if we have reached the maximum numbers to be in that file, write to a file, and then clear
            #record holder
            if num % how_many == 0:
                joiner.append("")
                with open(file,'w') as f:
                    f.write("\n".join(joiner))

                #change file name,clear record holder, and change the file count
                counter += 1
                file = parent_file_base_name + "_" + str(counter) + ".fasta"
                joiner = []
        if joiner:
            joiner.append("")
            with open(file,'w') as f:
                f.write("\n".join(joiner))

if __name__ == "__main__":
    try:
        parent_file = sys.argv[1]
        how_many = int(sys.argv[2])
    except IndexError:
        print "fasta_splitter.py parent_file.fasta how_many_per_file"
        sys.exit()
    write_file(parent_file,how_many)
    print "fasta_splitter.py is finished."
