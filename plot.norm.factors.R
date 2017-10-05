library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(broom)
library(genefilter)

setwd("/Users/dimitras/Documents/dimitra/Workspace/RNA-Seq/plot.norm.factors/")

#NOTE: Remove the lines below "-------" that explain the headers in the input files.

####################################################
#####################GENE LEVEL#####################
####################################################
gnf <- read.table("gene_normalization_factors.txt", header = TRUE, check.names = FALSE)

#total number of reads
gnf.m = gnf %>%
  gather(Measure,total.number.reads,totalnumreads)

gnf.m %>%   
  ggplot(aes(x=sample, y=total.number.reads/1000000, group=Measure, color=Measure)) + 
  # scale_y_log10() +
  geom_line(size=2) +
  theme_bw(base_size = 50) +
  theme(legend.position = "none") +
  labs(y="Total number of reads in millions", x="Samples")

ggsave("total.number.of.reads.plot.pdf", height = 35, width = 70, units ="cm")

#all percentages measures
library(RColorBrewer)
mypalette = c("purple","turquoise2","violet", "coral1", "darkolivegreen3", "firebrick2", "yellow2", "royalblue1")

gnf.m = gnf %>%
  gather(Measure,count,`%U`:`%senseGeneU`)
 
gnf.m %>%   
  ggplot(aes(x=sample, y=count, group=Measure, color=Measure)) + 
  geom_line(size=2) +
  scale_color_manual(values = mypalette) +
  theme_bw(base_size = 50) +
  theme(legend.position="bottom", legend.direction = "horizontal", legend.title = element_blank(), 
        legend.text = element_text(size = 40), legend.key.height = unit(2,"cm")) +
  # theme(legend.position = c(0.92,0.60), legend.title = element_blank(),
  #       legend.text = element_text(size = 40), legend.key.height = unit(2,"cm")) +
  labs(y="%count", x="Samples")

ggsave("gene.norm.factors.comboplot.pdf", height = 30, width = 70, units ="cm")



####################################################
######################EIJ LEVEL#####################
####################################################
enf <- read.table("exon-intron-junction_normalization_factors.txt", header = TRUE, check.names = FALSE)

#all percentages measures
mypalette = c("purple","turquoise2","violet", "coral1", "darkolivegreen3", "firebrick2", "yellow2", "royalblue1")

enf.m = enf %>%
  gather(Measure,count,`%exonicU_sense`:`%exon_inconsistentU`)

enf.m %>%   
  ggplot(aes(x=sample, y=count, group=Measure, color=Measure)) + 
  geom_line(size=2) +
  scale_color_manual(values = mypalette) +
  theme_bw(base_size = 50) +
  theme(legend.position="bottom", legend.direction = "horizontal", legend.title = element_blank(), 
        legend.text = element_text(size = 40), legend.key.height = unit(2,"cm")) +
  # theme(legend.position = c(0.92,0.60), legend.title = element_blank(),
  #       legend.text = element_text(size = 40), legend.key.height = unit(2,"cm")) +
  labs(y="%count", x="Samples")

ggsave("eij.norm.factors.comboplot.pdf", height = 30, width = 70, units ="cm")
