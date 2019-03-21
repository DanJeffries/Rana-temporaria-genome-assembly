import gzip
import sys

def fasta_maka(whitey, cat, out = None):

    """
    This script will make a .fasta file from the RADtag sequences found in a
    Stacks catalog.tags.tsv file. It requires a whitelist of catalog tag IDs.

    USAGE:
	python  Tags_2_fasta_with_whitelist.py  /full/path/to/whitelist.txt  /full/path/to/catalog.tags.tsv  [full/path/to/out.fasta]

    OUT:
        A .fasta file contaning the nucleotide sequences for each tag in the whitelist. The Tag_ID will become the sequence
        header. 

    """

    import sys
    import gzip
    
    if whitey is str:
        loci = open(whitey, 'r').readlines()
    elif isinstance(whitey, (list, set)):
        loci = whitey
    else:
        sys.exit("Unknown whitelist format - expected a python list or a file path")
        
    if cat.endswith("gz"):
        tags = gzip.open(cat, 'r').readlines()
    else:
        tags = open(cat, 'r').readlines()

    ## Pull out the locus ID's from the whitelist

    Loc_IDs = []
    for locus in loci:
        if locus.startswith("compli"):
            Loc_id = locus.split("_")[1]
        else:
            Loc_id = locus.split("_")[0]
        Loc_IDs.append(Loc_id)
    
    print "Number of tags in whitelist:",len(Loc_IDs)

    ## Write the fasta
    
    if not out == None:
        fasta = open(out, 'w')
        outpath = out
    else:
        fasta = open("%s/%s" % (cat.rpartition('/')[0], 'Whitelist_tags.fa'), 'w')
        outpath = "%s/%s" % (cat.rpartition('/')[0], 'Whitelist_tags.fa')

    count = 0
    for line in tags:
        if 'consensus' in line:
            Tag_ID = line.split()[2]
            if Tag_ID in Loc_IDs:
                count+=1
                fasta.write('>'+ Tag_ID +'\n'+line.split()[8]+'\n')
                
    print count, "sequences written to", outpath

    fasta.close()



if len(sys.argv) < 3:
    sys.exit("\n## Not enough arguments ##\n\n%s" % fasta_maka.__doc__)

elif len(sys.argv) > 3:
    sys.exit("\n## Too many arguments ## \n\n%s" % fasta_maka.__doc__)

else:
    fasta_maka(sys.argv[1], sys.argv[2])

