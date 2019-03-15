import vcf
import sys

vcf_path = sys.argv[1]
popmap = open(sys.argv[2], 'r').readlines()


## vcf module and Stacks disagree on the format of the Alelle depth field, so here I just alter the header in the VCF, otherwise the module throws an error.
## altered vcf file is saved in same location as original 

alteredvcfpath = "%s%s" % (vcf_path, ".altered")

oldvcf = open(vcf_path, 'r').readlines()
alteredvcf = open(alteredvcfpath, 'w')        

for line in oldvcf:
    
    if "Allele Depth" not in line:            
        alteredvcf.write(line)
    elif "Allele Depth" in line:            
        line = '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allele Depth">\n'
        alteredvcf.write(line)
        
alteredvcf.close()    
altered_vcf = open(alteredvcfpath, 'r')    
myvcf = vcf.Reader(altered_vcf)


altered_vcf = open(alteredvcfpath, 'r')    
myvcf = vcf.Reader(altered_vcf)

per_sample_dict = {}

loc_order_list = []
sample_order_list = []


for record in myvcf:
    #print record
    loc_id = "%s_%s" % (record.ID, record.POS)
    
    loc_order_list.append(loc_id) ## as dictionaries don't preserve order, doing that here with a list.
    
    ## Change the genotpe format
    
    for sample in record:
        if sample.sample not in per_sample_dict:
            per_sample_dict[sample.sample] = {}
            sample_order_list.append(sample.sample)
        
        if sample['GT'] == None:
            new_gt = "0 0"
        elif sample['GT'] == "0/0":
            new_gt = "1 1"
        elif sample['GT'] == "0/1":
            new_gt = "1 2"
        elif sample['GT'] == "1/0":
            new_gt = "2 1"
        elif sample['GT'] == "1/1":
            new_gt = "2 2"
        
        per_sample_dict[sample.sample][loc_id] = new_gt
            
sample_loc_lists = {}     
sample_lines = {}

## put the line together for each sample

for sample in sample_order_list:
    sample_loc_lists[sample] = []
    
    
    for loc in loc_order_list:
        sample_loc_lists[sample].append(per_sample_dict[sample][loc])

    sample_lines[sample] = "\t".join(sample_loc_lists[sample])
    
    

## get samples to keep from the popmap file

samples = []

for line in popmap:
    samples.append(line.split()[0])

print "\nKeeping %s samples\n" % len(samples)
    
## lastly, output the kept samples in new format. 

outfile = open("%s_Lepmap_input.dat" % sys.argv[1].rpartition(".")[0], 'w')

for sample in samples:
    outfile.write("%s\t%s\n" % (sample, sample_lines[sample]))
                  
outfile.close()

print "All done, your shiny new LepMap input file is here: %s_Lepmap_input.dat" % sys.argv[1].rpartition(".")[0]
