# Title     : to plot genome organization of SARS-CoV-2
# Objective : TODO
# Created by: yiquan
# Created on: 12/5/20
library(ggplot2)
library(gggenes)
library(cowplot)
library(ggpubr)
library("gridExtra")
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape2)
library(stringr)
library(dplyr)
library(crayon)

CoVgenome <- data.frame(
 molecule = "",
 Genome = "",
 start = c(266,21563,26245,26523,28274),
 end = c(21555,25384,26472,27191,29533),
 gene = c("ORFab","Spike","E","M","N")
)

classifying_mut_type <- function(Con){
  if (grepl('\\+',Con)){
    return ("insertion")
    }
  else if (grepl('\\-',Con)){
    return ("deletion")
    }
  else{
    return ("mismatch")
    }
  }

format_varfreq <- function(VarFreq){
  return (as.numeric(str_replace(VarFreq, '%', ''))/100)
  }



snp_table <- read_tsv('results/all_variants.snp') %>%
               mutate(mut_type=mapply(classifying_mut_type, Cons)) %>%
               mutate(VarFreq=mapply(format_varfreq, VarFreq)) %>%
               mutate(Sample=factor(Sample, levels=unique(rev(Sample)))) %>%
               select(Sample, Position, Cons, mut_type, VarFreq)


colorscale  <- brewer.pal(8,"Accent")
textsize <- 13
p <- ggplot(snp_table,aes(y=Sample,x=Position,color=mut_type)) +
  geom_point(aes(fill=mut_type,size=VarFreq),color='black',shape=21,stroke=0.3) +
  scale_size_continuous(range = c(0,3.5)) +
  scale_fill_manual(values=alpha(colorscale, 0.5)) +
  theme_cowplot(20) +
  theme(axis.title=element_text(size=textsize,face="bold"),
        axis.text=element_text(size=textsize,face="bold"),
        axis.text.x=element_text(angle = 0, hjust = 0.5,size=textsize, vjust=0.5,face="bold"),
        legend.key.size=unit(0.05,'in'),
        legend.title=element_blank(),
        legend.text=element_text(size=textsize,face="bold"),
        legend.position='right') +
  labs(y=expression(""),x=expression(bold("Genome Position")))

genomeplot <- ggplot(CoVgenome,
                     aes(xmin= start, xmax= end, y= "",fill = gene, label=gene)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  geom_gene_label(align = "centre") +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3")+theme_genes()+
  theme(axis.text.x=element_text(size=13,face="bold"))


png("graph/SARS-CoV-2_mut_freq_genome.png",
    width = 26, height = 30, units = "cm", res = 500)
t <- ggarrange(genomeplot, p + font("xy.text", size = 13)+ font("x", size = 13), heights = c(0.1,0.9),
                    ncol = 1, nrow = 2, align = "v")
annotate_figure(t,
                top = text_grob("Visualizing SARS-CoV-2 Genome", color = "black", face = "bold", size = 14),
                bottom = text_grob(" anything could be here ", color = "blue",
                                   hjust = 1, x = 1, face = "italic", size = 10),
                left = text_grob("anything", color = "green", rot = 90),
                fig.lab = "Figure 1", fig.lab.face = "bold"
                )
dev.off()