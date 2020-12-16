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
         scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0), labels=c("0.0","0.2","0.4","0.6","0.8","1.0"), limits=c(0,1)) +
         labs(y=expression(bold("Frequency")),x=expression())
  ggsave(graphname, p, height=1.5, width=1.5)
  }

MBCS_levels <- c('RRARS','RRARG','RRAHS','LRARS','deletion')
MBCS_freq_table <- read_tsv("results/all_MBCS_freq.tsv") %>%
					 mutate(MBCS=factor(MBCS, levels=MBCS_levels)) 
samples <- unique(MBCS_freq_table$sampleID)
for (sample in samples){
  sample_table <- MBCS_freq_table %>%
                    filter(sampleID==sample)
  plot_freq(sample_table, paste('graph/MBCS_freq_',sample,'.png',sep=''))
  }
