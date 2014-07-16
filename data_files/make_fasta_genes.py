from Bio import SeqIO as so
import sys
in_file = sys.argv[1]
out_file = sys.argv[2]

fi = [i for i in so.parse(in_file,'fasta')]

with open(out_file,'w') as f:
	for line in fi:
		id = line.id.split('|')[1]
		seq = str(line.seq).replace('.','').upper()
		str_to_write = ">{0}\n{1}\n".format(id,seq)
		f.write(str_to_write)
