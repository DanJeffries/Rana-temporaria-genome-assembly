library(ggplot2)
library(reshape2)
install.packages("cowplot")
library(cowplot)

##########################################################
###                                                    ###
###   Here I will plot the output of the IDB module    ###
###   of LepMap3. For the individual family files I    ###
###   split files manually in the terminal.            ###
###                                                    ###
##########################################################


### All families ###

ibd = read.delim("/home/djeffrie/Data/Genomes/Rtemp_hybrid/ALLMAPS_2019/SbfI/LepMap3_2019/INPUTS/ibd_only.txt")

p <- ggplot(ibd, aes(Indiv_1, Indiv_2)) + 
  geom_tile(aes(fill = IBD), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p

### ST01 only ###

ibd_ST01 = read.delim("/home/djeffrie/Data/Genomes/Rtemp_hybrid/ALLMAPS_2019/SbfI/LepMap3_2019/INPUTS/ibd_ST01_only.txt")

p_ST01 <- ggplot(ibd_ST01, aes(Indiv_1, Indiv_2)) +
  geom_tile(aes(fill = IBD), colour = "white") +
  scale_fill_gradient(low= "white", high = "steelblue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p_ST01

### ST16 only ###

ibd_ST16 = read.delim("/home/djeffrie/Data/Genomes/Rtemp_hybrid/ALLMAPS_2019/SbfI/LepMap3_2019/INPUTS/ibd_ST16_only.txt")

p_ST16 <- ggplot(ibd_ST16, aes(Indiv_1, Indiv_2)) +
  geom_tile(aes(fill = IBD), colour = "white") +
  scale_fill_gradient(low= "white", high = "steelblue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p_ST16

### ST18 only ###

ibd_ST18 = read.delim("/home/djeffrie/Data/Genomes/Rtemp_hybrid/ALLMAPS_2019/SbfI/LepMap3_2019/INPUTS/ibd_ST18_only.txt")

p_ST18 <- ggplot(ibd_ST18, aes(Indiv_1, Indiv_2)) +
  geom_tile(aes(fill = IBD), colour = "white") +
  scale_fill_gradient(low= "white", high = "steelblue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10), axis.text.y = element_text(size = 10))
p_ST18

### ST43 only ###

ibd_ST43 = read.delim("/home/djeffrie/Data/Genomes/Rtemp_hybrid/ALLMAPS_2019/SbfI/LepMap3_2019/INPUTS/ibd_ST43_only.txt")

p_ST43 <- ggplot(ibd_ST43, aes(Indiv_1, Indiv_2)) +
  geom_tile(aes(fill = IBD), colour = "white") +
  scale_fill_gradient(low= "white", high = "steelblue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10), axis.text.y = element_text(size = 10))
p_ST43

### ST45 only ###

ibd_ST45 = read.delim("/home/djeffrie/Data/Genomes/Rtemp_hybrid/ALLMAPS_2019/SbfI/LepMap3_2019/INPUTS/ibd_ST45_only.txt")

p_ST45 <- ggplot(ibd_ST45, aes(Indiv_1, Indiv_2)) +
  geom_tile(aes(fill = IBD), colour = "white") +
  scale_fill_gradient(low= "white", high = "steelblue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10), axis.text.y = element_text(size = 10))
p_ST45

### ST99 only ###

ibd_ST99 = read.delim("/home/djeffrie/Data/Genomes/Rtemp_hybrid/ALLMAPS_2019/SbfI/LepMap3_2019/INPUTS/ibd_ST99_only.txt")

p_ST99 <- ggplot(ibd_ST99, aes(Indiv_1, Indiv_2)) +
  geom_tile(aes(fill = IBD), colour = "white") +
  scale_fill_gradient(low= "white", high = "steelblue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10), axis.text.y = element_text(size = 10))
p_ST99

## Can plot to grid like this, but need to change font sizes.

plot_grid(p_ST01,p_ST16,p_ST18,p_ST43,p_ST45,p_ST99, cnol = 2)
