#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape2)
library(stringr)
library(dplyr)
require(cowplot)

plot_mut_freq_on_genome <- function(snp_table, graphname){
  colorscale  <- brewer.pal(8,"Accent")
  textsize <- 6
  p <- ggplot(snp_table,aes(y='P1',x=Position,color=mut_type)) +
         geom_point(aes(fill=mut_type,size=VarFreq),color='black',shape=21,stroke=0.3) +
         scale_size_continuous(range = c(0,3.5)) +
         scale_fill_manual(values=alpha(colorscale, 0.8)) +
         theme_cowplot(12) +
         theme(axis.title=element_text(size=textsize,face="bold"),
               axis.text=element_text(size=textsize,face="bold"),
               axis.text.x=element_text(angle = 0, hjust = 0.5,size=textsize, vjust=0.5,face="bold"),
               legend.key.size=unit(0.05,'in'),
               legend.title=element_blank(),
               legend.text=element_text(size=textsize,face="bold"),
               legend.position='right') +
         labs(y=expression(""),x=expression(bold("genome position")))
  ggsave(graphname, p, height=1.2, width=3.5)
  }

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

snp_table <- read_tsv('results/variants.snp') %>%
               mutate(mut_type=mapply(classifying_mut_type, Cons)) %>%
               mutate(VarFreq=mapply(format_varfreq, VarFreq)) %>%
               select(Position, Cons, mut_type, VarFreq)
plot_mut_freq_on_genome(snp_table, 'graph/mut_freq_genome.png')
