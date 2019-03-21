import sys

def Genetic_mapper_in_maka(LepMap_outs_path):
    
    """
    Genetic_mapper_in_maka   < full/path/to/Lepmap_file >
    
    Makes two input files for genetic_mapper.pl, one for the male and one for the female map
    
    """
    

    LEPmap_outs = open(LepMap_outs_path, 'r').readlines()
    
    male_out_path = "%s_MALE_genetic_mapper.dat" % LepMap_outs_path.rpartition(".")[0]
    female_out_path = "%s_FEMALE_genetic_mapper.dat" % LepMap_outs_path.rpartition(".")[0]

    male_genetic_mapper_input = open(male_out_path, 'w')
    male_genetic_mapper_input.write("ID\tLG\tPOS\tLOD\n")
    female_genetic_mapper_input = open(female_out_path, 'w')
    female_genetic_mapper_input.write("ID\tLG\tPOS\tLOD\n")

    for line in LEPmap_outs:
        if line.startswith("#***"):
            LG = line.split()[3]
            #print LG
        elif not line.startswith("#") and "duplicate" not in line:
            marker_ID = line.split()[0]
            male_pos = line.split()[1]
            female_pos = line.split()[2]

            male_genetic_mapper_input.write("%s\t%s\t%s\n" % (marker_ID,LG,male_pos))
            female_genetic_mapper_input.write("%s\t%s\t%s\n" % (marker_ID,LG,female_pos))

    male_genetic_mapper_input.close()
    female_genetic_mapper_input.close()



if len(sys.argv) != 2:
    sys.exit(Genetic_mapper_in_maka.__doc__)
    
else:
    Genetic_mapper_in_maka(sys.argv[1])
