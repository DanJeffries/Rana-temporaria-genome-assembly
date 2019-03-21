import sys

def LepMap_2_Stacks_IDs(Map_path, VCF_path):

    """
    This script will replace the marker IDs in the map file outputted by the LepMap3 "OrderMarkers2" module
    with those from a VCF. This script was written based on VCFs output by Stacks 1.48.
    
    USAGE:
    
        LepMap_2_Stacks_IDs.py  /full/path/to/map_file  /full/path/to/VCF
    
    OUPUTS:
    
        Will create a new map file with a similar name with LepMap IDs replaced by Stacks IDs
    
    """
    
    
    ## Make VCF ID dictionary

    VCF_ID_dict = {}

    tag_index = 1 ## keeps track of the order of loci in the vcf. This is essentially the LepMap ID. 

    for line in open(VCF_path, 'r').readlines():
        if not line.startswith("#"):        
            tag_ID = line.split()[2]

            VCF_ID_dict[str(tag_index)] = tag_ID

            tag_index += 1
            
            
    # Make new file to write the new IDs to - will be the same name as the input, just with "RealIDs" in it.

    new_map_file_path = "%s_%s.%s" % (Map_path.rpartition(".")[0], "RealIDs", Map_path.rpartition(".")[2])
    New_map_file = open(new_map_file_path, 'w')

    ## For each marker in the map file, rewrite the line to the new map file with the Stacks ID instead. 
    
    for line in open(Map_path, 'r').readlines():

        if line.startswith("#"):
            
            New_map_file.write(line)
            
        else:
            
            LepMap_ID = line.split()[0]
            Stacks_ID = VCF_ID_dict[LepMap_ID]  ## getting stacks ID here
            rest_of_line = "\t".join(line.split("\t")[1:])
            New_map_file.write("%s\t%s" % (Stacks_ID, rest_of_line))
            
    New_map_file.close()
    
    print "\nDone :) \n\nYour new map file is here: %s\n" % new_map_file_path


if len(sys.argv) < 3:
    sys.exit("\n## Not enough arguments ##\n\n%s" % LepMap_2_Stacks_IDs.__doc__)

elif len(sys.argv) > 3:
    sys.exit("\n## Too many arguments ## \n\n%s" %LepMap_2_Stacks_IDs.__doc__)

else:
    LepMap_2_Stacks_IDs(sys.argv[1], sys.argv[2])

