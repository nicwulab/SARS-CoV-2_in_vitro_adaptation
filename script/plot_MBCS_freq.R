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

plot_freq <- function(freq_table, graphname){
  textsize <- 7
  p <- ggplot(freq_table,aes(x=MBCS,y=freq),fill='black') +
         geom_bar(stat="identity", position=position_dodge(), width=0.5, fill='black') +
         theme_cowplot(12) +
         theme(axis.title=element_text(size=textsize,face="bold"),
               axis.text=element_text(size=textsize,face="bold"),
			   axis.text.x=element_text(angle = 90, hjust = 1,size=textsize, vjust=0.5,face="bold"),
			   legend.key.size=unit(0.12,'in'),
			   legend.title=element_blank(),
			   legend.text=element_text(size=textsize,face="bold"),
			   legend.position='right') +
         labs(y=expression(bold("Frequency")),x=expression())
  ggsave(graphname, p, height=1.5, width=1.5)
  }

renaming <- function(v){
  if (v=='-----'){return ('deletion')}
  else {return (v)}
  }

wrapper <- function(ID){
  MBCS_freq_table <- read_tsv(paste('Cov2_fastq/P',ID,'_S',ID,'/MBCS_freq.tsv',sep='')) %>%
                       mutate(MBCS=mapply(renaming, MBCS)) %>%
                       mutate(MBCS=factor(MBCS, levels=rev(MBCS))) 
  plot_freq(MBCS_freq_table, paste('graph/MBCS_freq_P',ID,'_S',ID,'.png',sep=''))
  }

for (ID in seq(1, 38)){
  if (ID==6){next}
  wrapper(ID)
  }

