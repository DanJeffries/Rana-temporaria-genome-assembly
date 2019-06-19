### cline version
from Bio import SeqIO
from Bio import Seq
import sys
import gzip

## define recplacement function

def harder(myrecord):

    new_seq = []

    for i in myrecord.seq:
        if i.isupper():
            new_seq.append(i)
        elif i.islower():
            new_seq.append("N")

    myrecord.seq = Seq.Seq("".join(new_seq))

    return myrecord


## Parse the command line options

if len(sys.argv) != 2:
	sys.exit("USAGE:    harder.py < full/path/to/genome.fa >")

genome_path = sys.argv[1]

if genome_path.endswith("gz"):
    genome = SeqIO.parse(gzip.open(genome_path, 'r'), "fasta")
else:
    genome = SeqIO.parse(open(genome_path, 'r'), "fasta")

## Output path is the same place as the input file

outpath = "%s/%s" % (genome_path.rpartition("/")[0], "hardmasked.fa")

hardmasked = open(outpath, 'w')

counter = 0
printcounter = 0
                     
for record in genome:
    
    counter += 1
    printcounter += 1
    
    if printcounter == 500:
        print "%s Sequences processed" % counter
        printcounter = 0
    
    new_record = harder(record)
    
    SeqIO.write(new_record, hardmasked, "fasta")
                     
hardmasked.close()
                     
print "Hard masked genome is here %s" % outpath
