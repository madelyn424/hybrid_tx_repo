---
title: "Hybrid transcript validation notebook"
output: html_notebook
---

#install packages
install.packages("tidyverse")
install.packages("hrbrthemes")
install.packages("viridis")

#Load packages
library(tidyverse)
library(hrbrthemes)
library(viridis)

# Inhibition graph --------------------------------------------------------

#8/29/2023 to do to get standard error calculated:

#1. Get means among replicates in columns for each transcript type AND total counts
#2. Get st.error among replicates in columns for each transcript type AND total counts

library(matrixStats)
df_replication_2nd_try <- read.csv("df_replication_2nd_try.csv")

df_replication_2nd_try$mean_3LTR <- rowMeans(df_replication_2nd_try[,c('LTR3_human_sequences_Replicate_1', 'LTR3_human_sequences_Replicate_2', 'LTR3_human_sequences_Replicate_3', 'LTR3_human_sequences_Replicate_4')],na.rm=TRUE)

df_replication_2nd_try$stdev_3LTR <- rowSds(as.matrix(df_replication_2nd_try[c('LTR3_human_sequences_Replicate_1', 'LTR3_human_sequences_Replicate_2', 'LTR3_human_sequences_Replicate_3', 'LTR3_human_sequences_Replicate_4')]),na.rm=TRUE)

df_replication_2nd_try <- df_replication_2nd_try %>%
  mutate(SE_3LTR = rowSds(as.matrix(df_replication_2nd_try[c('LTR3_human_sequences_Replicate_1', 'LTR3_human_sequences_Replicate_2', 'LTR3_human_sequences_Replicate_3', 'LTR3_human_sequences_Replicate_4')]),na.rm=TRUE)/sqrt(df_replication_2nd_try[c('Number_of_replicates')]))
df_replication_2nd_try <- df_replication_2nd_try %>% as.data.frame()

df_replication_2nd_try$mean_5LTR <- rowMeans(df_replication_2nd_try[,c('LTR5_human_sequences_Replicate_1', 'LTR5_human_sequences_Replicate_2', 'LTR5_human_sequences_Replicate_3', 'LTR5_human_sequences_Replicate_4')])

df_replication_2nd_try$stdev_5LTR <- rowSds(as.matrix(df_replication_2nd_try[c('LTR5_human_sequences_Replicate_1', 'LTR5_human_sequences_Replicate_2', 'LTR5_human_sequences_Replicate_3', 'LTR5_human_sequences_Replicate_4')]),na.rm=TRUE)

df_replication_2nd_try <- df_replication_2nd_try %>%
  mutate(SE_5LTR = rowSds(as.matrix(df_replication_2nd_try[c('LTR5_human_sequences_Replicate_1', 'LTR5_human_sequences_Replicate_2', 'LTR5_human_sequences_Replicate_3', 'LTR5_human_sequences_Replicate_4')]),na.rm=TRUE)/sqrt(df_replication_2nd_try[c('Number_of_replicates')]))
df_replication_2nd_try <- df_replication_2nd_try %>% as.data.frame()

df_replication_2nd_try$mean_HIV_only <- rowMeans(df_replication_2nd_try[,c('HIV_only_sequences_Replicate_1', 'HIV_only_sequences_Replicate_2', 'HIV_only_sequences_Replicate_3', 'HIV_only_sequences_Replicate_4')])

df_replication_2nd_try$stdev_HIV_only <- rowSds(as.matrix(df_replication_2nd_try[c('HIV_only_sequences_Replicate_1', 'HIV_only_sequences_Replicate_2', 'HIV_only_sequences_Replicate_3', 'HIV_only_sequences_Replicate_4')]),na.rm=TRUE)

df_replication_2nd_try <- df_replication_2nd_try %>%
  mutate(SE_HIV_only = rowSds(as.matrix(df_replication_2nd_try[c('HIV_only_sequences_Replicate_1', 'HIV_only_sequences_Replicate_2', 'HIV_only_sequences_Replicate_3', 'HIV_only_sequences_Replicate_4')]),na.rm=TRUE)/sqrt(df_replication_2nd_try[c('Number_of_replicates')]))
df_replication_2nd_try <- df_replication_2nd_try %>% as.data.frame()

df_replication_2nd_try$Total_sequences_Replicate_1 <- rowSums(df_replication_2nd_try[c('LTR3_human_sequences_Replicate_1', 'LTR5_human_sequences_Replicate_1', 'HIV_only_sequences_Replicate_1')])

df_replication_2nd_try$Total_sequences_Replicate_2 <- rowSums(df_replication_2nd_try[c('LTR3_human_sequences_Replicate_2', 'LTR5_human_sequences_Replicate_2', 'HIV_only_sequences_Replicate_2')])

df_replication_2nd_try$Total_sequences_Replicate_3 <- rowSums(df_replication_2nd_try[c('LTR3_human_sequences_Replicate_3', 'LTR5_human_sequences_Replicate_3', 'HIV_only_sequences_Replicate_3')])

df_replication_2nd_try$Total_sequences_Replicate_4 <- rowSums(df_replication_2nd_try[c('LTR3_human_sequences_Replicate_4', 'LTR5_human_sequences_Replicate_4', 'HIV_only_sequences_Replicate_4')])

df_replication_2nd_try$mean_Total_sequences <- rowMeans(df_replication_2nd_try[,c('Total_sequences_Replicate_1', 'Total_sequences_Replicate_2', 'Total_sequences_Replicate_3', 'Total_sequences_Replicate_4')])

df_replication_2nd_try$Stdev_Total_sequences <- rowSds(as.matrix(df_replication_2nd_try[c('Total_sequences_Replicate_1', 'Total_sequences_Replicate_2', 'Total_sequences_Replicate_3', 'Total_sequences_Replicate_4')]),na.rm=TRUE)

df_replication_2nd_try <- df_replication_2nd_try %>%
  mutate(SE_Total_sequences = rowSds(as.matrix(df_replication_2nd_try[c('Total_sequences_Replicate_1', 'Total_sequences_Replicate_2', 'Total_sequences_Replicate_3', 'Total_sequences_Replicate_4')]),na.rm=TRUE)/sqrt(df_replication_2nd_try[c('Number_of_replicates')]))
df_replication_2nd_try <- df_replication_2nd_try %>% as.data.frame()

df_replication_2nd_try$Total_hybrid_transcript_sequences_Replicate_1 <- rowSums(df_replication_2nd_try[c('LTR3_human_sequences_Replicate_1', 'LTR5_human_sequences_Replicate_1')])

df_replication_2nd_try$Total_hybrid_transcript_sequences_Replicate_2 <- rowSums(df_replication_2nd_try[c('LTR3_human_sequences_Replicate_2', 'LTR5_human_sequences_Replicate_2')])

df_replication_2nd_try$Total_hybrid_transcript_sequences_Replicate_3 <- rowSums(df_replication_2nd_try[c('LTR3_human_sequences_Replicate_3', 'LTR5_human_sequences_Replicate_3')])

df_replication_2nd_try$Total_hybrid_transcript_sequences_Replicate_4 <- rowSums(df_replication_2nd_try[c('LTR3_human_sequences_Replicate_4', 'LTR5_human_sequences_Replicate_4')])

df_replication_2nd_try$mean_Total_hybrid_transcript_sequences <- rowMeans(df_replication_2nd_try[,c('Total_hybrid_transcript_sequences_Replicate_1', 'Total_hybrid_transcript_sequences_Replicate_2', 'Total_hybrid_transcript_sequences_Replicate_3', 'Total_hybrid_transcript_sequences_Replicate_4')])

df_replication_2nd_try$Stdev_Total_hybrid_transcript_sequences <- rowSds(as.matrix(df_replication_2nd_try[c('Total_hybrid_transcript_sequences_Replicate_1', 'Total_hybrid_transcript_sequences_Replicate_2', 'Total_hybrid_transcript_sequences_Replicate_3', 'Total_hybrid_transcript_sequences_Replicate_4')]),na.rm=TRUE)

df_replication_2nd_try <- df_replication_2nd_try %>%
  mutate(SE_Total_sequences = rowSds(as.matrix(df_replication_2nd_try[c('Total_hybrid_transcript_sequences_Replicate_1', 'Total_hybrid_transcript_sequences_Replicate_2', 'Total_hybrid_transcript_sequences_Replicate_3', 'Total_hybrid_transcript_sequences_Replicate_4')]),na.rm=TRUE)/sqrt(df_replication_2nd_try[c('Number_of_replicates')]))
df_replication_2nd_try <- df_replication_2nd_try %>% as.data.frame()


#3. Do this all in counts/ug

#covert 3LTR into counts/ug
df_replication_2nd_try <- df_replication_2nd_try %>%
  mutate(LTR3_human_sequences_per_ug_Replicate_1 = (df_replication_2nd_try[c('LTR3_human_sequences_Replicate_1')])/(df_replication_2nd_try[c('Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA')]))

df_replication_2nd_try <- df_replication_2nd_try %>%
  mutate(LTR3_human_sequences_per_ug_Replicate_2 = (df_replication_2nd_try[c('LTR3_human_sequences_Replicate_2')])/(df_replication_2nd_try[c('Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA')]))

df_replication_2nd_try <- df_replication_2nd_try %>%
  mutate(LTR3_human_sequences_per_ug_Replicate_3 = (df_replication_2nd_try[c('LTR3_human_sequences_Replicate_3')])/(df_replication_2nd_try[c('Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA')]))

df_replication_2nd_try <- df_replication_2nd_try %>%
  mutate(LTR3_human_sequences_per_ug_Replicate_4 = (df_replication_2nd_try[c('LTR3_human_sequences_Replicate_4')])/(df_replication_2nd_try[c('Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA')]))

#Convert 5LTRMSD  into transcripts/ug

df_replication_2nd_try <- df_replication_2nd_try %>%
  mutate(LTR5_human_sequences_per_ug_Replicate_1 = (df_replication_2nd_try[c('LTR5_human_sequences_Replicate_1')])/(df_replication_2nd_try[c('Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA')]))

df_replication_2nd_try <- df_replication_2nd_try %>%
  mutate(LTR5_human_sequences_per_ug_Replicate_2 = (df_replication_2nd_try[c('LTR5_human_sequences_Replicate_2')])/(df_replication_2nd_try[c('Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA')]))

df_replication_2nd_try <- df_replication_2nd_try %>%
  mutate(LTR5_human_sequences_per_ug_Replicate_3 = (df_replication_2nd_try[c('LTR5_human_sequences_Replicate_3')])/(df_replication_2nd_try[c('Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA')]))

df_replication_2nd_try <- df_replication_2nd_try %>%
  mutate(LTR5_human_sequences_per_ug_Replicate_4 = (df_replication_2nd_try[c('LTR5_human_sequences_Replicate_4')])/(df_replication_2nd_try[c('Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA')]))

#HIV-only

#df_replication_2nd_try <- df_replication_2nd_try %>%
  #mutate(HIV_only_sequences_per_ug_Replicate_1 = (df_replication_2nd_try[c('HIV_only_sequences_Replicate_1')])/(df_replication_2nd_try[c('Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA')]))

#df_replication_2nd_try <- df_replication_2nd_try %>%
  #mutate(HIV_only_sequences_per_ug_Replicate_2 = (df_replication_2nd_try[c('HIV_only_sequences_Replicate_2')])/(df_replication_2nd_try[c('Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA')]))

#df_replication_2nd_try <- df_replication_2nd_try %>%
  #mutate(HIV_only_sequences_per_ug_Replicate_3 = (df_replication_2nd_try[c('HIV_only_sequences_Replicate_3')])/(df_replication_2nd_try[c('Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA')]))

#df_replication_2nd_try <- df_replication_2nd_try %>%
  #mutate(HIV_only_sequences_per_ug_Replicate_4 = (df_replication_2nd_try[c('HIV_only_sequences_Replicate_4')])/(df_replication_2nd_try[c('Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA')]))

#Total sequences/ug
df_replication_2nd_try <- df_replication_2nd_try %>%
  mutate(Total_hybrid_transcript_sequences_per_ug_Replicate_1 = (df_replication_2nd_try[c('Total_hybrid_transcript_sequences_Replicate_1')])/(df_replication_2nd_try[c('Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA')]))

df_replication_2nd_try <- df_replication_2nd_try %>%
  mutate(Total_hybrid_transcript_sequences_per_ug_Replicate_2 = (df_replication_2nd_try[c('Total_hybrid_transcript_sequences_Replicate_2')])/(df_replication_2nd_try[c('Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA')]))

df_replication_2nd_try <- df_replication_2nd_try %>%
  mutate(Total_hybrid_transcript_sequences_per_ug_Replicate_3 = (df_replication_2nd_try[c('Total_hybrid_transcript_sequences_Replicate_3')])/(df_replication_2nd_try[c('Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA')]))

df_replication_2nd_try <- df_replication_2nd_try %>%
  mutate(Total_hybrid_transcript_sequences_per_ug_Replicate_4 = (df_replication_2nd_try[c('Total_hybrid_transcript_sequences_Replicate_4')])/(df_replication_2nd_try[c('Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA')]))

#Mean counts/ug

df_replication_2nd_try$mean_3LTR_sequences_per_ug <- rowMeans(df_replication_2nd_try[,c('LTR3_human_sequences_per_ug_Replicate_1', 'LTR3_human_sequences_per_ug_Replicate_2', 'LTR3_human_sequences_per_ug_Replicate_3', 'LTR3_human_sequences_per_ug_Replicate_4')],na.rm=TRUE)

df_replication_2nd_try$stdev_3LTR_sequences_per_ug <- rowSds(as.matrix(df_replication_2nd_try[c('LTR3_human_sequences_per_ug_Replicate_1', 'LTR3_human_sequences_per_ug_Replicate_2', 'LTR3_human_sequences_per_ug_Replicate_3', 'LTR3_human_sequences_per_ug_Replicate_4')]),na.rm=TRUE)

df_replication_2nd_try <- df_replication_2nd_try %>%
  mutate(SE_3LTR_sequences_per_ug = rowSds(as.matrix(df_replication_2nd_try[c('LTR3_human_sequences_per_ug_Replicate_1', 'LTR3_human_sequences_per_ug_Replicate_2', 'LTR3_human_sequences_per_ug_Replicate_3', 'LTR3_human_sequences_per_ug_Replicate_4')]),na.rm=TRUE)/sqrt(df_replication_2nd_try[c('Number_of_replicates')]))
df_replication_2nd_try <- df_replication_2nd_try %>% as.data.frame()

df_replication_2nd_try$mean_5LTR_sequences_per_ug <- rowMeans(df_replication_2nd_try[,c('LTR5_human_sequences_per_ug_Replicate_1', 'LTR5_human_sequences_per_ug_Replicate_2', 'LTR5_human_sequences_per_ug_Replicate_3', 'LTR5_human_sequences_per_ug_Replicate_4')])

df_replication_2nd_try$stdev_5LTR_sequences_per_ug <- rowSds(as.matrix(df_replication_2nd_try[c('LTR5_human_sequences_per_ug_Replicate_1', 'LTR5_human_sequences_per_ug_Replicate_2', 'LTR5_human_sequences_per_ug_Replicate_3', 'LTR5_human_sequences_per_ug_Replicate_4')]),na.rm=TRUE)

df_replication_2nd_try <- df_replication_2nd_try %>%
  mutate(SE_5LTR_sequences_per_ug = rowSds(as.matrix(df_replication_2nd_try[c('LTR5_human_sequences_per_ug_Replicate_1', 'LTR5_human_sequences_per_ug_Replicate_2', 'LTR5_human_sequences_per_ug_Replicate_3', 'LTR5_human_sequences_per_ug_Replicate_4')]),na.rm=TRUE)/sqrt(df_replication_2nd_try[c('Number_of_replicates')]))
df_replication_2nd_try <- df_replication_2nd_try %>% as.data.frame()

#df_replication_2nd_try$mean_HIV_only_sequences_per_ug <- rowMeans(df_replication_2nd_try[,c('HIV_only_sequences_per_ug_Replicate_1', 'HIV_only_sequences_per_ug_Replicate_2', 'HIV_only_sequences_per_ug_Replicate_3', 'HIV_only_sequences_per_ug_Replicate_4')])

#df_replication_2nd_try$stdev_HIV_only_sequences_per_ug <- rowSds(as.matrix(df_replication_2nd_try[c('HIV_only_sequences_per_ug_Replicate_1', 'HIV_only_sequences_per_ug_Replicate_2', 'HIV_only_sequences_per_ug_Replicate_3', 'HIV_only_sequences_per_ug_Replicate_4')]),na.rm=TRUE)

#df_replication_2nd_try <- df_replication_2nd_try %>%
  #mutate(SE_HIV_only_sequences_per_ug = rowSds(as.matrix(df_replication_2nd_try[c('HIV_only_sequences_per_ug_Replicate_1', 'HIV_only_sequences_per_ug_Replicate_2', 'HIV_only_sequences_per_ug_Replicate_3', 'HIV_only_sequences_per_ug_Replicate_4')]),na.rm=TRUE)/sqrt(df_replication_2nd_try[c('Number_of_replicates')]))
#df_replication_2nd_try <- df_replication_2nd_try %>% as.data.frame()

df_replication_2nd_try$Total_hybrid_transcript_sequences_per_ug_Replicate_1 <- rowSums(df_replication_2nd_try[c('LTR3_human_sequences_per_ug_Replicate_1', 'LTR5_human_sequences_per_ug_Replicate_1')])

df_replication_2nd_try$Total_hybrid_transcript_sequences_per_ug_Replicate_2 <- rowSums(df_replication_2nd_try[c('LTR3_human_sequences_per_ug_Replicate_2', 'LTR5_human_sequences_per_ug_Replicate_2')])

df_replication_2nd_try$Total_hybrid_transcript_sequences_per_ug_Replicate_3 <- rowSums(df_replication_2nd_try[c('LTR3_human_sequences_per_ug_Replicate_3', 'LTR5_human_sequences_per_ug_Replicate_3')])

df_replication_2nd_try$Total_hybrid_transcript_sequences_per_ug_Replicate_4 <- rowSums(df_replication_2nd_try[c('LTR3_human_sequences_per_ug_Replicate_4', 'LTR5_human_sequences_per_ug_Replicate_4')])

df_replication_2nd_try$mean_Total_hybrid_transcript_sequences_per_ug <- rowMeans(df_replication_2nd_try[,c('Total_hybrid_transcript_sequences_per_ug_Replicate_1', 'Total_hybrid_transcript_sequences_per_ug_Replicate_2', 'Total_hybrid_transcript_sequences_per_ug_Replicate_3', 'Total_hybrid_transcript_sequences_per_ug_Replicate_4')],na.rm=TRUE)

df_replication_2nd_try$Stdev_Total_hybrid_transcript_sequences_per_ug <- rowSds(as.matrix(df_replication_2nd_try[c('Total_hybrid_transcript_sequences_per_ug_Replicate_1', 'Total_hybrid_transcript_sequences_per_ug_Replicate_2', 'Total_hybrid_transcript_sequences_per_ug_Replicate_3', 'Total_hybrid_transcript_sequences_per_ug_Replicate_4')]),na.rm=TRUE)

df_replication_2nd_try <- df_replication_2nd_try %>%
  mutate(SE_Total_hybrid_transcript_sequences_per_ug = rowSds(as.matrix(df_replication_2nd_try[c('Total_hybrid_transcript_sequences_per_ug_Replicate_1', 'Total_hybrid_transcript_sequences_per_ug_Replicate_2', 'Total_hybrid_transcript_sequences_per_ug_Replicate_3', 'Total_hybrid_transcript_sequences_per_ug_Replicate_4')]),na.rm=TRUE)/sqrt(df_replication_2nd_try[c('Number_of_replicates')]))
df_replication_2nd_try <- df_replication_2nd_try %>% as.data.frame()

unlist(df_replication_2nd_try)

df_replication_2nd_try <- df_replication_2nd_try %>% as.data.frame()

#4. Redo graphs with standard error instead of standard deviation 

ggplot(df_replication_2nd_try, aes(x=PID, y=mean_Total_hybrid_transcript_sequences_per_ug, fill=as.character(Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA))) +
  geom_bar (data = subset(df_replication_2nd_try, PID %in% c("PIC1709-Exp1","PIC1709-Exp2","PIC1352","PIC1373","PIC1858")), stat="identity", position=position_dodge2(preserve="single"), width = NULL) +
  ylab("Mean Hybrid Transcripts/Ug")+
  geom_errorbar(data = subset(df_replication_2nd_try, PID %in% c("PIC1709-Exp1","PIC1709-Exp2","PIC1352","PIC1373","PIC1858")), aes(ymin=(mean_Total_hybrid_transcript_sequences_per_ug-SE_Total_hybrid_transcript_sequences_per_ug$Number_of_replicates), ymax=(mean_Total_hybrid_transcript_sequences_per_ug+SE_Total_hybrid_transcript_sequences_per_ug$Number_of_replicates)), position=position_dodge2(preserve="single"))+
  scale_fill_manual(name="Ug RNA Input", labels = c("0.13","0.25","0.5","1","1.9","2","2.9","3"),      values = c("#F8766D","#00BE67","#00BFC4","#00A9FF","#C77CFF","#C77CFF","#FF61CC","#FF61CC"))+
  theme(axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), axis.line.x = element_line(color="black"),            axis.line.y = element_line(color="black"), panel.grid.major = element_line(size=0.5, linetype='solid', color="white"),    panel.grid.minor = element_line(size=0.25, linetype='solid', color="white")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width=30)) +
  scale_y_continuous(limits = c(0,150), expand = c(0, 0))+
  theme(panel.background = element_rect(fill = "white")) +
  geom_text(
    aes(label = Number_of_replicates),
    vjust = 5,
    hjust = 0.9,
    angle = 0,
    position = position_dodge2(width =0.9, preserve = "single"),
    size = 2.5) +
  guides(shape = guide_legend(override.aes = list(size = 10)), color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10),
        legend.key.size = unit(5, "mm"))

#with standard devs for comparison
ggplot(df_replication_2nd_try, aes(x=PID, y=mean_Total_sequences_per_ug, fill=as.character(Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA))) +
  geom_bar (data = subset(df_replication_2nd_try, PID %in% c("PIC1373","PIC1709","PIC1858")), stat="identity", position=position_dodge2(preserve="single"), width = NULL) +
  ylab("Mean HIV Transcripts/Ug")+
  geom_errorbar(data = subset(df_replication_2nd_try, PID %in% c("PIC1373","PIC1709","PIC1858")), aes(ymin=(mean_Total_sequences_per_ug-Stdev_Total_sequences_per_ug), ymax=(mean_Total_sequences_per_ug+Stdev_Total_sequences_per_ug)), position=position_dodge2(preserve="single"))+
  scale_fill_manual(name="Ug RNA Input", labels = c("0.13","0.25","0.5","1","1.9","2","2.9","3"),      values = c("#F8766D","#00BE67","#00BFC4","#00A9FF","#C77CFF","#C77CFF","#FF61CC","#FF61CC"))+
  theme(axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), axis.line.x = element_line(color="black"),            axis.line.y = element_line(color="black"), panel.grid.major = element_line(size=0.5, linetype='solid', color="grey"),    panel.grid.minor = element_line(size=0.25, linetype='solid', color="grey")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width=30)) +
  scale_y_continuous()+
  theme(panel.background = element_rect(fill = "white")) +
  geom_text(
    aes(label = Number_of_replicates),
    vjust = 5,
    hjust = 0.9,
    angle = 0,
    position = position_dodge2(width =0.9, preserve = "single"),
    size = 2.5) +
  guides(shape = guide_legend(override.aes = list(size = 10)), color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10),
        legend.key.size = unit(5, "mm"))

library(gridExtra)
library(grid)

grid.newpage()
grid.draw(cbind(ggplotGrob(g1),ggplotGrob(g2),ggplotGrob(g3), size="last"))

grid.arrange(g1,g2,g3,newpage="TRUE")


# oligo(dT)s vs R6 All PIC ------------------

library(dplyr)
library(tidyverse)

#To do
#INHIBITION CURVE SUGGESTS THAT OPTIMAL RNA INPUT IS BETWEEN 0.5 AND 1UG. ONLY SELECT DATA WITHIN THOSE BOUNDARIES FOR THIS ANALYSIS. (Go up to 1.5ug for 1373 because there is no choice--too few seqs.)

#See if I can also use qPCR data (especially for 1352)
#See if I can only include experiments that have technical replicates and add those error bars.

#df_megasheet_all_observations_summarized3
#agg_tbl_megasheet <-df_megasheet_all_observations_summarized3 %>%
  #filter(Ug_RNA_input_Capture > 0.45)
#df_megasheet_all_observations_summarized_cDNA <- agg_tbl_megasheet %>% as.data.frame()

#df_megasheet_all_observations_summarized_cDNA
#agg_tbl_megasheet <-df_megasheet_all_observations_summarized_cDNA %>%
  #filter(Ug_RNA_input_Capture < 1.5)
#df_megasheet_all_observations_summarized_cDNA <- agg_tbl_megasheet %>% as.data.frame()

#df_megasheet_all_observations_summarized_cDNA 
#agg_tbl_megasheet <-df_megasheet_all_observations_summarized_cDNA %>%
  #filter(Num_of_replicates > 1)
#df_megasheet_all_observations_summarized_cDNA <- agg_tbl_megasheet %>% as.data.frame()

#df_megasheet_all_observations_summarized_cDNA 
#agg_tbl_megasheet <-df_megasheet_all_observations_summarized_cDNA %>%
  #filter(!PIC_ID == 1709)
#df_megasheet_all_observations_summarized_cDNA <- agg_tbl_megasheet %>% as.data.frame()

#WITH concentrations factored in:
df_megasheet_all_observations_summarized3
agg_tbl_megasheet <-df_megasheet_all_observations_summarized3 %>%
  select(PIC_ID, Ug_RNA_input_Capture, cDNA_primers, Sum_total_tx_each_replicate, Replicate)
df_megasheet_all_observations_summarized_cDNA <- agg_tbl_megasheet %>% as.data.frame()

df_megasheet_all_observations_summarized_cDNA
agg_tbl_megasheet <-df_megasheet_all_observations_summarized_cDNA %>%
  group_by(PIC_ID, Ug_RNA_input_Capture, cDNA_primers, Replicate)%>%
  mutate(Sum_cDNA_primer_per_replicate = sum(Sum_total_tx_each_replicate))
df_megasheet_all_observations_summarized_cDNA1 <- agg_tbl_megasheet %>% as.data.frame()

df_megasheet_all_observations_summarized_cDNA1
agg_tbl_megasheet <-df_megasheet_all_observations_summarized_cDNA1 %>%
  group_by(PIC_ID, Ug_RNA_input_Capture, cDNA_primers, Replicate)%>%
  mutate(Sum_total_tx_per_Ug_cDNA_primer_per_replicate = sum(Sum_cDNA_primer_per_replicate)/sum(Ug_RNA_input_Capture))
df_megasheet_all_observations_summarized_cDNA1 <- agg_tbl_megasheet %>% as.data.frame()

df_megasheet_all_observations_summarized_cDNA1
agg_tbl_megasheet <-df_megasheet_all_observations_summarized_cDNA1 %>%
  group_by(PIC_ID, Ug_RNA_input_Capture, cDNA_primers)%>%
  mutate(Sum_cDNA_primer_each_conc = sum(Sum_total_tx_each_replicate))
df_megasheet_all_observations_summarized_cDNA1 <- agg_tbl_megasheet %>% as.data.frame()

df_megasheet_all_observations_summarized_cDNA1
agg_tbl_megasheet <-df_megasheet_all_observations_summarized_cDNA1 %>%
  select(c(PIC_ID, cDNA_primers, Ug_RNA_input_Capture, Sum_total_tx_per_Ug_cDNA_primer_per_replicate, Sum_cDNA_primer_per_replicate, Sum_cDNA_primer_each_conc, Replicate)) 
df_megasheet_all_observations_summarized_cDNA2 <- agg_tbl_megasheet %>% as.data.frame()

#Keep only dates with BOTH R6 and dTs tested
#foo <- split(df_megasheet_all_observations_summarized_cDNA2, f=df_megasheet_all_observations_summarized_cDNA2$Pulldown_date)
#lapply(foo, FUN=function(x){unique(x$cDNA_primers)})

#lapply(foo, FUN=function(x){length(unique(x$cDNA_primers))})

#lapply(foo, FUN=function(x){if(length(unique(x$cDNA_primers))==2){x}else{NA} })

#foo2 <- lapply(foo, FUN=function(x){if(length(unique(x$cDNA_primers))==2){x}else{NA} })

#foo3 <- foo2[!is.na(foo2)]

#df_megasheet_all_observations_summarized_cDNA3 <- do.call(rbind, foo3)


#Get rid of duplicates again
df_megasheet_all_observations_summarized_cDNA3 <- df_megasheet_all_observations_summarized_cDNA2 %>% distinct()


setDT(df_megasheet_all_observations_summarized_cDNA3)
df_megasheet_all_observations_summarized_cDNA3[, Mean_total_tx_per_Ug_cDNA_primer :=mean(Sum_total_tx_per_Ug_cDNA_primer_per_replicate), by = c('PIC_ID','Ug_RNA_input_Capture','cDNA_primers')]

setDT(df_megasheet_all_observations_summarized_cDNA3)
df_megasheet_all_observations_summarized_cDNA3[, SE_total_tx_per_Ug_cDNA_primer :=std.error(Sum_total_tx_per_Ug_cDNA_primer_per_replicate), by = c('PIC_ID','Ug_RNA_input_Capture','cDNA_primers')]

#Select only columns needed
df_megasheet_all_observations_summarized_cDNA
agg_tbl_megasheet <-df_megasheet_all_observations_summarized_cDNA3 %>%
  select(c(PIC_ID, Ug_RNA_input_Capture, cDNA_primers, Sum_cDNA_primer_each_conc,  Mean_total_tx_per_Ug_cDNA_primer, SE_total_tx_per_Ug_cDNA_primer)) 
df_megasheet_all_observations_summarized_cDNA4 <- agg_tbl_megasheet %>% as.data.frame()

#Get rid of duplicates again
df_megasheet_all_observations_summarized_cDNA4 <- df_megasheet_all_observations_summarized_cDNA4 %>% distinct()

#Graph it
ggplot(df_megasheet_all_observations_summarized_cDNA4, aes(x=as.character(Ug_RNA_input_Capture), y=Mean_total_tx_per_Ug_cDNA_primer, fill = cDNA_primers)) +
  geom_bar (stat="identity", position=position_dodge2(preserve="single")) +
  ylab("Counts_per_Ug") +
  xlab("Ug RNA input (Capture)")+
  theme(axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), axis.line.x = element_line(color="black"),            axis.line.y = element_line(color="black"), panel.grid.major = element_line(size=0.5, linetype='solid', color="grey"),    panel.grid.minor = element_line(size=0.25, linetype='solid', color="grey")) +
  
  geom_errorbar(aes(ymin=Mean_total_tx_per_Ug_cDNA_primer-SE_total_tx_per_Ug_cDNA_primer, ymax=Mean_total_tx_per_Ug_cDNA_primer+SE_total_tx_per_Ug_cDNA_primer), position=position_dodge2(preserve="single")) +
geom_text(aes(label= Sum_cDNA_primer_each_conc), vjust = 0.5, hjust = -1, angle = 90, position = position_dodge2(width =0.9, preserve = "single"), size = 2.5) +
  
  scale_x_discrete(guide=guide_axis(angle=90), labels = function(x) str_wrap(x, width=20))+
  
  scale_y_continuous()+
  facet_grid(~ Gene_Tx) +
  theme(panel.background = element_rect(size = 5, fill = "white"), axis.ticks.length = unit(0.5, "cm")) +
  guides(shape = guide_legend(override.aes = list(size = 10)), color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10),
        legend.key.size = unit(5, "mm"))+

  
facet_grid(~ PIC_ID, scales = "free_x", space="free_x") +
  theme(panel.background = element_rect(size = 5, fill = "white"), axis.ticks.length = unit(0.5, "cm")) +
  guides(shape = guide_legend(override.aes = list(size = 10)), colour = guide_legend(nrow=2)) +
  theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10),
        legend.key.size = unit(5, "mm"))

#With NO concentrations factored in
df_megasheet_all_observations_summarized3
agg_tbl_megasheet <-df_megasheet_all_observations_summarized3 %>%
  select(PIC_ID, cDNA_primers, Ug_RNA_input_Capture, Sum_total_tx_each_replicate, Replicate)
df_megasheet_all_observations_summarized_cDNA <- agg_tbl_megasheet %>% as.data.frame()

df_megasheet_all_observations_summarized_cDNA
agg_tbl_megasheet <-df_megasheet_all_observations_summarized_cDNA %>%
  group_by(PIC_ID, cDNA_primers, Replicate)%>%
  mutate(Sum_cDNA_primer_per_replicate = sum(Sum_total_tx_each_replicate))
df_megasheet_all_observations_summarized_cDNA1 <- agg_tbl_megasheet %>% as.data.frame()

df_megasheet_all_observations_summarized_cDNA1
agg_tbl_megasheet <-df_megasheet_all_observations_summarized_cDNA1 %>%
  group_by(PIC_ID, cDNA_primers, Replicate)%>%
  mutate(Sum_total_tx_per_Ug_cDNA_primer_per_replicate = sum(Sum_cDNA_primer_per_replicate)/sum(Ug_RNA_input_Capture))
df_megasheet_all_observations_summarized_cDNA1 <- agg_tbl_megasheet %>% as.data.frame()

df_megasheet_all_observations_summarized_cDNA1
agg_tbl_megasheet <-df_megasheet_all_observations_summarized_cDNA1 %>%
  group_by(PIC_ID, cDNA_primers)%>%
  mutate(Sum_cDNA_primer_each_conc = sum(Sum_total_tx_each_replicate))
df_megasheet_all_observations_summarized_cDNA1 <- agg_tbl_megasheet %>% as.data.frame()

df_megasheet_all_observations_summarized_cDNA1
agg_tbl_megasheet <-df_megasheet_all_observations_summarized_cDNA1 %>%
  select(c(PIC_ID, cDNA_primers, Sum_total_tx_per_Ug_cDNA_primer_per_replicate, Sum_cDNA_primer_per_replicate, Replicate)) 
df_megasheet_all_observations_summarized_cDNA2 <- agg_tbl_megasheet %>% as.data.frame()

#Keep only dates with BOTH R6 and dTs tested
#foo <- split(df_megasheet_all_observations_summarized_cDNA2, f=df_megasheet_all_observations_summarized_cDNA2$Pulldown_date)
#lapply(foo, FUN=function(x){unique(x$cDNA_primers)})

#lapply(foo, FUN=function(x){length(unique(x$cDNA_primers))})

#lapply(foo, FUN=function(x){if(length(unique(x$cDNA_primers))==2){x}else{NA} })

#foo2 <- lapply(foo, FUN=function(x){if(length(unique(x$cDNA_primers))==2){x}else{NA} })

#foo3 <- foo2[!is.na(foo2)]

#df_megasheet_all_observations_summarized_cDNA3 <- do.call(rbind, foo3)


#Get rid of duplicates again
df_megasheet_all_observations_summarized_cDNA3 <- df_megasheet_all_observations_summarized_cDNA2 %>% distinct()


setDT(df_megasheet_all_observations_summarized_cDNA3)
df_megasheet_all_observations_summarized_cDNA3[, Mean_total_tx_per_Ug_cDNA_primer :=mean(Sum_total_tx_per_Ug_cDNA_primer_per_replicate), by = c('PIC_ID','cDNA_primers')]

setDT(df_megasheet_all_observations_summarized_cDNA3)
df_megasheet_all_observations_summarized_cDNA3[, SE_total_tx_per_Ug_cDNA_primer :=std.error(Sum_total_tx_per_Ug_cDNA_primer_per_replicate), by = c('PIC_ID','cDNA_primers')]

#Select only columns needed
df_megasheet_all_observations_summarized_cDNA
agg_tbl_megasheet <-df_megasheet_all_observations_summarized_cDNA3 %>%
  select(c(PIC_ID, cDNA_primers,  Mean_total_tx_per_Ug_cDNA_primer, SE_total_tx_per_Ug_cDNA_primer)) 
df_megasheet_all_observations_summarized_cDNA4 <- agg_tbl_megasheet %>% as.data.frame()

#Get rid of duplicates again
df_megasheet_all_observations_summarized_cDNA4 <- df_megasheet_all_observations_summarized_cDNA4 %>% distinct()

#Graph it
ggplot(df_megasheet_all_observations_summarized_cDNA4, aes(x=as.character(cDNA_primers), y=Mean_total_tx_per_Ug_cDNA_primer, fill = cDNA_primers)) +
  geom_bar (stat="identity", position=position_dodge2(preserve="single")) +
  ylab("Counts_per_Ug") +
  xlab("Ug RNA input (Capture)")+
  theme(axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), axis.line.x = element_line(color="black"),            axis.line.y = element_line(color="black"), panel.grid.major = element_line(size=0.5, linetype='solid', color="grey"),    panel.grid.minor = element_line(size=0.25, linetype='solid', color="grey")) +
  
  geom_errorbar(aes(ymin=Mean_total_tx_per_Ug_cDNA_primer-SE_total_tx_per_Ug_cDNA_primer, ymax=Mean_total_tx_per_Ug_cDNA_primer+SE_total_tx_per_Ug_cDNA_primer), position=position_dodge2(preserve="single")) +
  #geom_text(aes(label= Sum_cDNA_primer), vjust = 0.5, hjust = -1, angle = 90, position = position_dodge2(width =0.9, preserve = "single"), size = 2.5) +
  
  scale_x_discrete(guide=guide_axis(angle=90), labels = function(x) str_wrap(x, width=20))+
  
  scale_y_continuous()+
  facet_grid(~ Gene_Tx) +
  theme(panel.background = element_rect(size = 5, fill = "white"), axis.ticks.length = unit(0.5, "cm")) +
  guides(shape = guide_legend(override.aes = list(size = 10)), color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10),
        legend.key.size = unit(5, "mm"))+
  
  
  facet_grid(~ PIC_ID, scales = "free_x", space="free_x") +
  theme(panel.background = element_rect(size = 5, fill = "white"), axis.ticks.length = unit(0.5, "cm")) +
  guides(shape = guide_legend(override.aes = list(size = 10)), colour = guide_legend(nrow=2)) +
  theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10),
        legend.key.size = unit(5, "mm"))

# #Adding rows with zero value for replicates with no observation ---------

df_megasheet <- read.csv("MegaTxSheet11Jan2024.csv")

agg_tbl_megasheet <- df_megasheet %>%
filter(Normal_or_repetitive == "Normal")
df2_megasheet <- agg_tbl_megasheet %>% as.data.frame()
 
agg_tbl_megasheet_2replicates <- subset(df2_megasheet, Num_of_replicates =='2')

#Apparently in order for .drop = false to work I need to have ONLY 1 column be a factor (the Replicate column)

agg_tbl_megasheet_2replicates$PIC_ID <- as.character(agg_tbl_megasheet_2replicates$PIC_ID)
agg_tbl_megasheet_2replicates$Pulldown_date <- as.character(agg_tbl_megasheet_2replicates$Pulldown_date)
agg_tbl_megasheet_2replicates$cDNA_primers <- as.character(agg_tbl_megasheet_2replicates$cDNA_primers)
agg_tbl_megasheet_2replicates$Ug_RNA_input_Capture <- as.character(agg_tbl_megasheet_2replicates$Ug_RNA_input_Capture)
agg_tbl_megasheet_2replicates$Gene_Tx <- as.character(agg_tbl_megasheet_2replicates$Gene_Tx)
agg_tbl_megasheet_2replicates$Replicate <- as.factor(agg_tbl_megasheet_2replicates$Replicate)
agg_tbl_megasheet_2replicates$Num_of_replicates <- as.character(agg_tbl_megasheet_2replicates$Num_of_replicates)
agg_tbl_megasheet_2replicates_rows <-  agg_tbl_megasheet_2replicates %>%
group_by(PIC_ID, Pulldown_date, cDNA_primers, Ug_RNA_input_Capture, Gene_Tx, Num_of_replicates, Replicate, .drop = FALSE)  %>%
summarize(Sum_total_tx_each_replicate = length(Gene_Tx))
ungroup(agg_tbl_megasheet_2replicates_rows)

#3 replicates
agg_tbl_megasheet_3replicates <- subset(df2_megasheet, Num_of_replicates =='3')

#Apparently in order for .drop = false to work I need to have ONLY 1 column be a factor (the Replicate column)
  
agg_tbl_megasheet_3replicates$PIC_ID <- as.character(agg_tbl_megasheet_3replicates$PIC_ID)
agg_tbl_megasheet_3replicates$Pulldown_date <- as.character(agg_tbl_megasheet_3replicates$Pulldown_date)
agg_tbl_megasheet_3replicates$cDNA_primers <- as.character(agg_tbl_megasheet_3replicates$cDNA_primers)
agg_tbl_megasheet_3replicates$Ug_RNA_input_Capture <- as.character(agg_tbl_megasheet_3replicates$Ug_RNA_input_Capture)
agg_tbl_megasheet_3replicates$Gene_Tx <- as.character(agg_tbl_megasheet_3replicates$Gene_Tx)
agg_tbl_megasheet_3replicates$Replicate <- as.factor(agg_tbl_megasheet_3replicates$Replicate)
agg_tbl_megasheet_3replicates$Num_of_replicates <- as.character(agg_tbl_megasheet_3replicates$Num_of_replicates)

  #add rows
agg_tbl_megasheet_3replicates_rows <-  agg_tbl_megasheet_3replicates %>%
group_by(PIC_ID, Pulldown_date, cDNA_primers, Ug_RNA_input_Capture, Gene_Tx, Num_of_replicates, Replicate, .drop = FALSE)  %>%
#perform counts
summarize(Sum_total_tx_each_replicate = length(Gene_Tx))
ungroup(agg_tbl_megasheet_3replicates_rows)

#4 replicates
agg_tbl_megasheet_4replicates <- subset(df2_megasheet, Num_of_replicates =='4')

#Apparently in order for .drop = false to work I need to have ONLY 1 column be a factor (the Replicate column)

agg_tbl_megasheet_4replicates$PIC_ID <- as.character(agg_tbl_megasheet_4replicates$PIC_ID)
agg_tbl_megasheet_4replicates$Pulldown_date <- as.character(agg_tbl_megasheet_4replicates$Pulldown_date)
agg_tbl_megasheet_4replicates$cDNA_primers <- as.character(agg_tbl_megasheet_4replicates$cDNA_primers)
agg_tbl_megasheet_4replicates$Ug_RNA_input_Capture <- as.character(agg_tbl_megasheet_4replicates$Ug_RNA_input_Capture)
agg_tbl_megasheet_4replicates$Gene_Tx <- as.character(agg_tbl_megasheet_4replicates$Gene_Tx)
agg_tbl_megasheet_4replicates$Replicate <- as.factor(agg_tbl_megasheet_4replicates$Replicate)
agg_tbl_megasheet_4replicates$Num_of_replicates <- as.character(agg_tbl_megasheet_4replicates$Num_of_replicates)

#add rows
agg_tbl_megasheet_4replicates_rows <-  agg_tbl_megasheet_4replicates %>%
group_by(PIC_ID, Pulldown_date, cDNA_primers, Ug_RNA_input_Capture, Gene_Tx, Num_of_replicates, Replicate, .drop = FALSE)  %>%
#perform counts
summarize(Sum_total_tx_each_replicate = length(Gene_Tx))
ungroup(agg_tbl_megasheet_4replicates_rows)

#1 replicate (no need to add rows)
agg_tbl_megasheet_1replicates <- subset(df2_megasheet, Num_of_replicates =='1')

#Apparently in order for .drop = false to work I need to have ONLY 1 column be a factor (the Replicate column). I don't need that here but I'm just keeping things consistent with the other sheets
agg_tbl_megasheet_1replicates$PIC_ID <- as.character(agg_tbl_megasheet_1replicates$PIC_ID)
agg_tbl_megasheet_1replicates$Pulldown_date <- as.character(agg_tbl_megasheet_1replicates$Pulldown_date)
agg_tbl_megasheet_1replicates$cDNA_primers <- as.character(agg_tbl_megasheet_1replicates$cDNA_primers)
agg_tbl_megasheet_1replicates$Ug_RNA_input_Capture <- as.character(agg_tbl_megasheet_1replicates$Ug_RNA_input_Capture)
agg_tbl_megasheet_1replicates$Gene_Tx <- as.character(agg_tbl_megasheet_1replicates$Gene_Tx)
agg_tbl_megasheet_1replicates$Replicate <- as.factor(agg_tbl_megasheet_1replicates$Replicate)
agg_tbl_megasheet_1replicates$Num_of_replicates <- as.character(agg_tbl_megasheet_1replicates$Num_of_replicates)


#perform counts (no need to add rows)
agg_tbl_megasheet_1replicates_rows <-  agg_tbl_megasheet_1replicates %>%
group_by(PIC_ID, Pulldown_date, cDNA_primers, Ug_RNA_input_Capture, Gene_Tx, Num_of_replicates, Replicate)  %>%
summarize(Sum_total_tx_each_replicate = length(Gene_Tx))
ungroup(agg_tbl_megasheet_1replicates_rows)


#sew the sheets together
df_megasheet_all_observations_summarized <- rbind(agg_tbl_megasheet_1replicates_rows, agg_tbl_megasheet_2replicates_rows, agg_tbl_megasheet_3replicates_rows, agg_tbl_megasheet_4replicates_rows)

#Replace 0 with 0.00001
df_megasheet_all_observations_summarized$Sum_total_tx_each_replicate[df_megasheet_all_observations_summarized$Sum_total_tx_each_replicate == "0"] <- 0.0001


#Not sure why but I need to ungroup this now?
df_megasheet_all_observations_summarized %>%
  ungroup()  
df_megasheet_all_observations_summarized <- df_megasheet_all_observations_summarized %>% as.data.frame()

# #Catagorize transcripts by transcript type (GeneType) ------------------------------


#df_megasheet_all_observations_summarized$Replicate <- as.character(df_megasheet_all_observations_summarized$Replicate)
#df_megasheet_all_observations_summarized$Gene_Tx <- as.factor(df_megasheet_all_observations_summarized$Gene_Tx)

#df_megasheet_all_observations_summarized <- update(df_megasheet_all_observations_summarized, na.action = na.exclude)

#Add gene_type
agg_tbl_megasheet <- df_megasheet_all_observations_summarized %>%
  mutate(df_megasheet_all_observations_summarized, Gene_Type=ifelse(Gene_Tx == 'HIV', "HIV-only", ifelse(grepl("5LTR", df_megasheet_all_observations_summarized$Gene_Tx), "5LTR_MSD","3LTR")))
df_megasheet_all_observations_summarized2 <- agg_tbl_megasheet %>% as.data.frame()

agg_tbl_megasheet <-df_megasheet_all_observations_summarized2 %>%
select(c(PIC_ID, Gene_Tx, Gene_Type, Ug_RNA_input_Capture, Replicate, Sum_total_tx_each_replicate, Pulldown_date))
df_megasheet_all_observations_summarized3 <- agg_tbl_megasheet %>% as.data.frame()

# Normalize counts by counts/ug, graph (doesn't make sense to include graph here but has code I guess?) ----------------------------------------

#Normalize by counts/ug
df_megasheet_all_observations_summarized2$Ug_RNA_input_Capture <- as.numeric(df_megasheet_all_observations_summarized2$Ug_RNA_input_Capture)

#is.character(df_megasheet_all_observations_summarized3$Sum_total_tx_each_replicate)
#is.factor(df_megasheet_all_observations_summarized3$Gene_Type)
#is.factor(df_megasheet_all_observations_summarized3$PIC_ID)


#str(df_megasheet_all_observations_summarized2)

df_megasheet_all_observations_summarized2 %>%
  ungroup()  
df_megasheet_all_observations_summarized2 <- df_megasheet_all_observations_summarized2 %>% as.data.frame()

library(dplyr)
df_megasheet_all_observations_summarized2
agg_tbl_megasheet <-df_megasheet_all_observations_summarized2 %>%
  group_by(PIC_ID, cDNA_primers, Gene_Tx, Pulldown_date, Replicate)%>%
  mutate(Sum_total_tx_each_replicate_per_Ug = sum(Sum_total_tx_each_replicate)/sum(Ug_RNA_input_Capture))
df_megasheet_all_observations_summarized3 <- agg_tbl_megasheet %>% as.data.frame()

#Normalize by gene_type/ug
df_megasheet_all_observations_summarized3
agg_tbl_megasheet <-df_megasheet_all_observations_summarized3 %>%
  group_by(PIC_ID,Gene_Type)%>%
  mutate(Sum_Gene_type = n())
df_megasheet_all_observations_summarized3 <- agg_tbl_megasheet %>% as.data.frame()

df_megasheet_all_observations_summarized3
agg_tbl_megasheet <-df_megasheet_all_observations_summarized3 %>%
mutate(Sum_total_Genetype_each_replicate_per_Ug = sum(Sum_Gene_type)/sum(Ug_RNA_input_Capture))
df_megasheet_all_observations_summarized3 <- agg_tbl_megasheet %>% as.data.frame()


#write.csv(df6_1373_3, "df6_1373_3.csv")

#Select only columns I need
library(dplyr)
library(tidyverse)
library(plotrix)

agg_tbl_megasheet <-df_megasheet_all_observations_summarized3  %>%
  select(c(PIC_ID, Gene_Tx, Ug_RNA_input_Capture, Sum_total_tx_each_Gene_Type_per_Ug, Gene_Type, Replicate)) 
df_megasheet_all_observations_summarized4 <- agg_tbl_megasheet %>% as.data.frame()

df_megasheet_all_observations_summarized5 <- df_megasheet_all_observations_summarized4 %>% distinct()

library(data.table)
setDT(df_megasheet_all_observations_summarized5)
df_megasheet_all_observations_summarized5 [, Mean_total_each_tx_type_per_ug :=mean(Sum_total_tx_each_Gene_Type_per_Ug), by = c('PIC_ID', 'Gene_Type')]

setDT(df_megasheet_all_observations_summarized5)
df_megasheet_all_observations_summarized5[, SE_total_tx_each_tx_type_per_ug :=std.error(Sum_total_tx_each_Gene_Type_per_Ug), by = c('PIC_ID', 'Gene_Type')]

agg_tbl_megasheet <-df_megasheet_all_observations_summarized5 %>%
  select(c(PIC_ID, Gene_Type, Mean_total_each_tx_type_per_ug, SE_total_tx_each_tx_type_per_ug))
df_megasheet_all_observations_summarized6 <- agg_tbl_megasheet %>% as.data.frame()

df_megasheet_all_observations_summarized6 <- df_megasheet_all_observations_summarized6 %>% distinct()

library(ggplot2)
library(scales)
library(stringr)

ggplot(df_megasheet_all_observations_summarized6, aes(x=Gene_Type, y=Mean_total_each_tx_type_per_ug, fill=Gene_Type)) +
  geom_bar (stat="identity", position=position_dodge2(preserve="single"), width = NULL) +
  ylab("HIV-cDNA sequences/ug CD4 RNA input") +
  theme(axis.text.x = element_text(angle=90, size=8, hjust = 1, vjust=0.5), axis.line.x = element_line(color="black"),            axis.line.y = element_line(color="black"), panel.grid.major = element_line(size=0.5, linetype='solid', color="grey"),    panel.grid.minor = element_line(size=0.25, linetype='solid', color="grey")) +
  scale_fill_discrete(name="HIV Transcript Type")+
  scale_x_discrete(labels = function(x) str_wrap(x, width=30), name="Gene_Type") +
  scale_y_log10()+
  theme(panel.background = element_rect(fill = "white")) +
  #geom_text(aes(label = Count_Total2), vjust = -1, colour = "black", angle = 0, size = 3) +
  geom_errorbar(aes(ymin=Mean_total_each_tx_type_per_ug-SE_total_tx_each_tx_type_per_ug, ymax=Mean_total_each_tx_type_per_ug+SE_total_tx_each_tx_type_per_ug), position=position_dodge2(preserve="single")) +
  guides(shape = guide_legend(override.aes = list(size = 10)), color = guide_legend(override.aes = list(size = 10))) +
  guides(shape = guide_legend(override.aes = list(size = 10)), color = guide_legend(override.aes = list(size = 10))) +
  theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10),
        legend.key.size = unit(6, "mm"))+
facet_wrap('PIC_ID')
  
  


# #Percent distribution of transcripts in EACH PID from a megasheet --------

#df_megasheet <- read.csv("MegaTxSheet02Oct2023.csv")
library(dplyr)
#df_megasheet
#agg_tbl_megasheet <- df_megasheet %>%
  #filter(Normal_or_repetitive == "Normal")
#df2_megasheet <- agg_tbl_megasheet %>% as.data.frame()

#Catagorize transcripts by transcript type
#agg_tbl_megasheet <- df2_megasheet %>%
  #mutate(df2_megasheet, Gene_Type=ifelse(grepl("3LTR", df2_megasheet$ID), "3LTR",
                                         #ifelse(grepl("5LTR", df2_megasheet$ID), "5LTR_MSD","HIV-only")))
#df2_megasheet <- agg_tbl_megasheet %>% as.data.frame()

#Count transcripts that showed up by Gene_loc, date, input RNA, and replicate
library(dplyr)
#df2_megasheet
#agg_tbl_megasheet <- df2_megasheet %>%
  #group_by(PIC_ID, Gene_Tx, Pulldown_date, Replicate, Ug_RNA_input_Capture) %>%
  #mutate(Count_Total = n())
#df3_megasheet <- agg_tbl_megasheet %>% as.data.frame()

#Count transcripts by broad transcript type
library(dplyr)
df_megasheet_all_observations_summarized3
agg_tbl_megasheet <- df_megasheet_all_observations_summarized3 %>%
  group_by(PIC_ID, Gene_Type, Pulldown_date, Replicate) %>%
  mutate(Count_Total_Gene_Type_per_date_per_replicate = sum(Sum_total_tx_each_replicate))
df_megasheet_all_observations_summarized_genetype_count <- agg_tbl_megasheet %>% as.data.frame()

setDT(df_megasheet_all_observations_summarized_genetype_count)
df_megasheet_all_observations_summarized_genetype_count[, Mean_total_tx_per_Gene_type :=mean(Count_Total_Gene_Type_per_date_per_replicate), by = c('PIC_ID',"Gene_Type", "Pulldown_date")]

setDT(df_megasheet_all_observations_summarized_genetype_count)
df_megasheet_all_observations_summarized_genetype_count[, SE_total_tx_per_Gene_type :=std.error(Count_Total_Gene_Type_per_date_per_replicate), by = c('PIC_ID',"Gene_Type", "Pulldown_date")]


#Get select only the columns needed for broad transcript types and get rid of duplicates
df_megasheet_all_observations_summarized_genetype_count_2 <- df_megasheet_all_observations_summarized_genetype_count %>% select("PIC_ID", "Gene_Type","Pulldown_date", "Mean_total_tx_per_Gene_type", "SE_total_tx_per_Gene_type")
df_megasheet_all_observations_summarized_genetype_count_2 <- df_megasheet_all_observations_summarized_genetype_count_2 %>% distinct()

#Get percentages FIX THIS PART SO EACH PID IS 100% (will probably have to solve tomorrow AM but at least graph looks ok so far)--DONE

library(dplyr)
agg_tbl_megasheet <- df_megasheet_all_observations_summarized_genetype_count_2 %>%
  group_by(PIC_ID, Pulldown_date) %>%
  mutate(Percentage_Mean = round(Mean_total_tx_per_Gene_type/sum(Mean_total_tx_per_Gene_type)*100, 1))
df_megasheet_all_observations_summarized_genetype_count_3 <- agg_tbl_megasheet %>% as.data.frame()

agg_tbl_megasheet <- df_megasheet_all_observations_summarized_genetype_count_3 %>%
  group_by(PIC_ID, Pulldown_date) %>%
  mutate(Percentage_SE = round(SE_total_tx_per_Gene_type/sum(Mean_total_tx_per_Gene_type)*100, 1))
df_megasheet_all_observations_summarized_genetype_count_3 <- agg_tbl_megasheet %>% as.data.frame()

#Graph
ggplot(df_megasheet_all_observations_summarized_genetype_count_3, aes(x=Pulldown_date, y=Percentage_Mean, fill=Gene_Type)) +
  geom_bar (stat="identity", position=position_dodge2(preserve="single"), width = NULL) +
  ylab("% Total transcript distribution") +
  theme(axis.text.x = element_text(angle=90, size=8, hjust = 1, vjust=0.5), panel.grid.major = element_line(size=0.5, linetype='solid', color="white"),panel.grid.minor = element_line(size=0.25, linetype='solid', color="white"), panel.spacing = unit(0, "lines"),panel.border = element_rect(color = "black", fill = NA, size = 0.25),strip.background = element_rect(color = "black", size = 0.25)) +
  geom_errorbar(aes(ymin=Percentage_Mean-Percentage_SE, ymax=Percentage_Mean+Percentage_SE), position=position_dodge2(preserve="single"))+
  scale_fill_discrete(name="HIV Transcript Type")+
  scale_x_discrete(labels = function(x) str_wrap(x, width=30), name="Experiment Date") +
  scale_y_continuous(limits = c(0,80), expand = c(0, 0))+
  theme(panel.background = element_rect(fill = "white")) +
  guides(shape = guide_legend(override.aes = list(size = 10)), color = guide_legend(override.aes = list(size = 10))) +
  guides(shape = guide_legend(override.aes = list(size = 10)), color = guide_legend(override.aes = list(size = 10))) +
  theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10),
        legend.key.size = unit(6, "mm"))+
facet_grid(~ PIC_ID, scales = "free_x", space="free_x")
#geom_text(aes(label = Count_Total_Gene_Type), vjust = -.5, colour = "black", angle = 0, size = 2,  position = position_dodge(width = 1)) 


library(ggplot2)

ggplot(df_megasheet_all_observations_summarized_genetype_count_unique2, aes(fill=Gene_Type, y=Percentage, x=as.character(PIC_ID))) + 
  geom_bar(position='stack', stat='identity')+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.25),strip.background = element_rect(color = "black", size = 0.25), panel.background = element_rect(fill = "white"))+
  scale_y_continuous(limits = c(0,100), expand = c(0, 0))+
  ylab("PIC_ID") +
  geom_text(aes(label = Count_Total_Gene_Type))

# Percent distribution of high frequency transcripts in EACH PID from a megasheet BY DATE)--------

agg_tbl_megasheet <- df_megasheet_all_observations_summarized3 %>%
  group_by(PIC_ID, Gene_Tx) %>%
  mutate(Count_Total_Pulldown_dates_With_This_Gene_tx = n_distinct(Pulldown_date))
df_megasheet_all_observations_summarized_3_genetype_date <- agg_tbl_megasheet %>% as.data.frame()

#Next get count means and SE by date 
setDT(df_megasheet_all_observations_summarized_3_genetype_date)
df_megasheet_all_observations_summarized_3_genetype_date[, Mean_total_tx_per_Gene_Tx :=mean("Sum_total_tx_each_replicate"), by = c('PIC_ID',"Gene_Tx", "Pulldown_date")]

setDT(df_megasheet_all_observations_summarized_3_genetype_date)
df_megasheet_all_observations_summarized_3_genetype_date[, SE_total_tx_per_Gene_Tx :=std.error("Sum_total_tx_each_replicate"), by = c('PIC_ID',"Gene_Tx", "Pulldown_date")]

#Next get percentages by dates
df_megasheet_all_observations_summarized_3_genetype <- 
  df_megasheet_all_observations_summarized_3_genetype %>%
  group_by(PIC_ID) %>%
  mutate(Percentage = round(Count_Total_Gene_Type / sum(Count_Total_Gene_Type)*100, 1))


#Then filter by those showing up on more than 2 dates
agg_tbl_megasheet <- df_megasheet_all_observations_summarized_3_genetype_date %>%
filter(Count_Total_Pulldown_dates_With_This_Gene_tx>2)
df_megasheet_all_observations_summarized_4_genetype_date <- agg_tbl_megasheet %>% as.data.frame()


#I want to group all the high frequency transcripts together.
#Change Gene_Tx names so that those with only 1 count are Unique


library(dplyr)

#Graph
ggplot(df_megasheet_all_observations_summarized_3_genetype, aes(x=Gene_Tx3, y=Percentage, fill=Pulldown_date)) +
  geom_bar (stat="identity", position=position_dodge2(preserve="single"), width = NULL) +
  ylab("% Total transcript distribution") +
  theme(axis.text.x = element_text(angle=90, size=8, hjust = 1, vjust=0.5), panel.grid.major = element_line(size=0.5, linetype='solid', color="white"),panel.grid.minor = element_line(size=0.25, linetype='solid', color="white"), panel.spacing = unit(0, "lines"),panel.border = element_rect(color = "black", fill = NA, size = 0.25),strip.background = element_rect(color = "black", size = 0.25)) +
  scale_fill_discrete(name="HIV Transcript Type")+
  scale_x_discrete(labels = function(x) str_wrap(x, width=30), name="Gene_location") +
  scale_y_continuous(limits = c(0,70), expand = c(0, 0))+
  theme(panel.background = element_rect(fill = "white")) +
  guides(shape = guide_legend(override.aes = list(size = 10)), color = guide_legend(override.aes = list(size = 10))) +
  guides(shape = guide_legend(override.aes = list(size = 10)), color = guide_legend(override.aes = list(size = 10))) +
  theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10),
        legend.key.size = unit(6, "mm"))+
  facet_grid(~ PIC_ID, scales = "free_x", space="free_x")+
  geom_text(aes(label = Count_Total_Gene_Type), vjust = -.5, colour = "black", angle = 0, size = 2,  position = position_dodge(width = 1)) 

library(ggplot2)

ggplot(df_megasheet_all_observations_summarized_genetype_count_unique2, aes(fill=Gene_Type, y=Percentage, x=as.character(PIC_ID))) + 
  geom_bar(position='stack', stat='identity')+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.25),strip.background = element_rect(color = "black", size = 0.25), panel.background = element_rect(fill = "white"))+
  scale_y_continuous(limits = c(0,100), expand = c(0, 0))+
  ylab("PIC_ID") +
  geom_text(aes(label = Count_Total_Gene_Type))

# Percent distribution of high frequency transcripts in EACH PID from a megasheet (in progress, this is dumb)--------


df_megasheet_all_observations_summarized3
agg_tbl_megasheet <- df_megasheet_all_observations_summarized3 %>%
  group_by(PIC_ID, Gene_Type) %>%
  mutate(Count_Total_Gene_Type = sum(Sum_total_tx_each_replicate))
df_megasheet_all_observations_summarized_3_genetype <- agg_tbl_megasheet %>% as.data.frame()


#I want to group all the "unique" transcripts together.
#Change Gene_Tx names so that those with only 1 count are Unique

df_megasheet_all_observations_summarized_3_genetype$Gene_Tx2 <- "Hybrid Transcripts with counts of <3/replicate"
rowsIwant <- df_megasheet_all_observations_summarized_3_genetype$Sum_total_tx_each_replicate>2
df_megasheet_all_observations_summarized_3_genetype$Gene_Tx2[rowsIwant] <- as.character(df_megasheet_all_observations_summarized_3_genetype$Gene_Tx[rowsIwant])

agg_tbl_megasheet <- df_megasheet_all_observations_summarized_3_genetype %>%
  mutate(Gene_Tx3 = if_else(str_detect(Gene_Type,"HIV-only"),"HIV-only",Gene_Tx2))
df_megasheet_all_observations_summarized_3_genetype <- agg_tbl_megasheet %>% as.data.frame()

library(dplyr)
df_megasheet_all_observations_summarized_3_genetype <- 
  df_megasheet_all_observations_summarized_3_genetype %>%
  group_by(PIC_ID) %>%
  mutate(Percentage = round(Count_Total_Gene_Type / sum(Count_Total_Gene_Type)*100, 1))

#Graph
ggplot(df_megasheet_all_observations_summarized_3_genetype, aes(x=Gene_Tx3, y=Percentage, fill=Pulldown_date)) +
  geom_bar (stat="identity", position=position_dodge2(preserve="single"), width = NULL) +
  ylab("% Total transcript distribution") +
  theme(axis.text.x = element_text(angle=90, size=8, hjust = 1, vjust=0.5), panel.grid.major = element_line(size=0.5, linetype='solid', color="white"),panel.grid.minor = element_line(size=0.25, linetype='solid', color="white"), panel.spacing = unit(0, "lines"),panel.border = element_rect(color = "black", fill = NA, size = 0.25),strip.background = element_rect(color = "black", size = 0.25)) +
  scale_fill_discrete(name="HIV Transcript Type")+
  scale_x_discrete(labels = function(x) str_wrap(x, width=30), name="Gene_location") +
  scale_y_continuous(limits = c(0,70), expand = c(0, 0))+
  theme(panel.background = element_rect(fill = "white")) +
  guides(shape = guide_legend(override.aes = list(size = 10)), color = guide_legend(override.aes = list(size = 10))) +
  guides(shape = guide_legend(override.aes = list(size = 10)), color = guide_legend(override.aes = list(size = 10))) +
  theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10),
        legend.key.size = unit(6, "mm"))+
  facet_grid(~ PIC_ID, scales = "free_x", space="free_x")+
  geom_text(aes(label = Count_Total_Gene_Type), vjust = -.5, colour = "black", angle = 0, size = 2,  position = position_dodge(width = 1)) 

library(ggplot2)

ggplot(df_megasheet_all_observations_summarized_genetype_count_unique2, aes(fill=Gene_Type, y=Percentage, x=as.character(PIC_ID))) + 
  geom_bar(position='stack', stat='identity')+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.25),strip.background = element_rect(color = "black", size = 0.25), panel.background = element_rect(fill = "white"))+
  scale_y_continuous(limits = c(0,100), expand = c(0, 0))+
  ylab("PIC_ID") +
  geom_text(aes(label = Count_Total_Gene_Type))


# Reproducibility of assay in all PID ---------------------------------------

#Only have random hexamers in this graph
df_megasheet_all_observations_summarized3
agg_tbl_megasheet <- df_megasheet_all_observations_summarized3 %>%
  filter(cDNA_primers == "R6")
df_megasheet_all_observations_R6_only <- agg_tbl_megasheet %>% as.data.frame()

#Normalize by counts/ug
library(dplyr)
df_megasheet_all_observations_R6_only
agg_tbl_megasheet <-df_megasheet_all_observations_R6_only %>%
  group_by(Sum_total_tx_each_replicate, Ug_RNA_input_Capture) %>%
  mutate(Sum_total_tx_each_replicate_per_Ug = Sum_total_tx_each_replicate/Ug_RNA_input_Capture)
df_megasheet_all_observations_R6_only2 <- agg_tbl_megasheet %>% as.data.frame()

#write.csv(df6_1373_3, "df6_1373_3.csv")

#Select only columns I need
library(dplyr)
library(tidyverse)

#df_megasheet_all_observations_R6_only2
#agg_tbl_megasheet <-df_megasheet_all_observations_R6_only2 %>%
  #select(c(Gene_Tx, Pulldown_date, Sum_total_tx_each_replicate_per_Ug, Gene_Tx, Replicate)) 
#df_megasheet_all_observations_R6_only3 <- agg_tbl_megasheet %>% as.data.frame()

#df_megasheet_all_observations_R6_only3 <- df_megasheet_all_observations_R6_only3 %>% distinct()

setDT(df_megasheet_all_observations_R6_only2)
df_megasheet_all_observations_R6_only2[, Mean_total_tx_each_replicate_per_ug :=mean(Sum_total_tx_each_replicate_per_Ug), by = c('Gene_Tx', 'Pulldown_date')]

setDT(df_megasheet_all_observations_R6_only2)
df_megasheet_all_observations_R6_only2[, SE_total_tx_each_replicate_per_ug :=std.error(Sum_total_tx_each_replicate_per_Ug), by = c('Gene_Tx', 'Pulldown_date')]

library(dplyr)

#Select all columns except "replicates" and "counts"
df_megasheet_all_observations_R6_only2
agg_tbl_megasheet <-df_megasheet_all_observations_R6_only2 %>%
  select(c(PIC_ID, Gene_Tx, Pulldown_date, Sum_total_tx_each_replicate, Sum_total_tx_each_replicate_per_Ug, Replicate, Mean_total_tx_each_replicate_per_ug, SE_total_tx_each_replicate_per_ug)) 
df_megasheet_all_observations_R6_only3 <- agg_tbl_megasheet %>% as.data.frame()

#Get rid of duplicates again
df_megasheet_all_observations_R6_only3 <- df_megasheet_all_observations_R6_only3 %>% distinct()


#other way of selecting, filtering by frequency it shows up between experiment dates
#gene_frequency <- table(df_megasheet_all_observations_R6_only4$Gene_Tx)
#duplicate_genes <- names(gene_frequency)[gene_frequency > 1]
#print(gene_frequency)
#add column to main sheet, then select frequency > 3 and graph this shit.

#Select only columns I need
library(dplyr)
library(tidyverse)

df_megasheet_all_observations_R6_only3
agg_tbl_megasheet <-df_megasheet_all_observations_R6_only3 %>%
  select(c(PIC_ID, Gene_Tx, Pulldown_date, Mean_total_tx_each_replicate_per_ug, SE_total_tx_each_replicate_per_ug, Gene_Tx)) 
df_megasheet_all_observations_R6_only4 <- agg_tbl_megasheet %>% as.data.frame()

df_megasheet_all_observations_R6_only5 <- df_megasheet_all_observations_R6_only4 %>% distinct()

agg_tbl_megasheet <-df_megasheet_all_observations_R6_only5 %>%
  group_by(Gene_Tx, Mean_total_tx_each_replicate_per_ug) %>% 
  filter(PIC_ID == 1352, Mean_total_tx_each_replicate_per_ug >5) 
df_megasheet_all_observations_R6_1352 <- agg_tbl_megasheet %>% as.data.frame()

agg_tbl_megasheet <-df_megasheet_all_observations_R6_only5 %>%
  group_by(Gene_Tx, Mean_total_tx_each_replicate_per_ug) %>% 
  filter(PIC_ID == 1709, Mean_total_tx_each_replicate_per_ug >10) 
df_megasheet_all_observations_R6_1709 <- agg_tbl_megasheet %>% as.data.frame()

agg_tbl_megasheet <-df_megasheet_all_observations_R6_only5 %>%
  group_by(Gene_Tx, Mean_total_tx_each_replicate_per_ug) %>% 
  filter(PIC_ID == 1858, Mean_total_tx_each_replicate_per_ug >2.5) 
df_megasheet_all_observations_R6_1858 <- agg_tbl_megasheet %>% as.data.frame()

agg_tbl_megasheet <-df_megasheet_all_observations_R6_only5 %>%
  group_by(Gene_Tx, Mean_total_tx_each_replicate_per_ug) %>% 
  filter(PIC_ID == 1373, Mean_total_tx_each_replicate_per_ug >2) 
df_megasheet_all_observations_R6_1373 <- agg_tbl_megasheet %>% as.data.frame()   

df_megasheet_all_observations_R6_only6 <- rbind(df_megasheet_all_observations_R6_1352, df_megasheet_all_observations_R6_1709, df_megasheet_all_observations_R6_1373, df_megasheet_all_observations_R6_1858)

#Select the genes I want
#To try: Count at least 2 observations
#Count at least 1 observation on 2 days

df_megasheet_all_observations_R6_only5 <- df_megasheet_all_observations_R6_only6 %>% group_by(Gene_Tx) %>% filter(n()>1)

#DO NOT DELETE BELOW UNTIL I KNOW I NO LONGER NEED IT
df_megasheet_all_observations_R6_1373 <- df_megasheet_all_observations_R6_only5 %>% filter(PIC_ID == 1373, grepl("5LTR_MSE_STAT5B_4223", Gene_Tx) | grepl("FANCA", Gene_Tx) | grepl("EIF", Gene_Tx) | grepl("STAT5B_4227", Gene_Tx) | grepl("CDK13", Gene_Tx) | grepl("NDC80", Gene_Tx))

df_megasheet_all_observations_R6_1858 <- df_megasheet_all_observations_R6_only5 %>% filter(PIC_ID == 1858, grepl("5LTR_MSD_STAT5B", Gene_Tx) | grepl("FLT", Gene_Tx)| grepl("SIPA", Gene_Tx))

df_megasheet_all_observations_R6_1709 <- df_megasheet_all_observations_R6_only5 %>% filter(PIC_ID == 1709, grepl("5LTR_MSD_STAT5B_4223", Gene_Tx) | grepl("MAP4K1_38611", Gene_Tx))

df_megasheet_all_observations_R6_1352 <- df_megasheet_all_observations_R6_only5 %>% filter(PIC_ID == 1352, grepl("ARIH_489", Gene_Tx) | grepl("234411545", Gene_Tx) | grepl("BRE_2791", Gene_Tx) | grepl("RPTOR_807", Gene_Tx) | grepl("98727975", Gene_Tx) | grepl("YWHAZ_10093", Gene_Tx))

#df_megasheet_all_observations_R6_1709 <- df_megasheet_all_observations_R6_only5[!grepl('38611132', df_megasheet_all_observations_R6_only5$Gene_Tx),]
#df_megasheet_all_observations_R6_1709 <- df_megasheet_all_observations_R6_only5[!grepl('42257672', df_megasheet_all_observations_R6_only5$Gene_Tx),]
#df_megasheet_all_observations_R6_1709 <- df_megasheet_all_observations_R6_only5[!grepl('HIVE', df_megasheet_all_observations_R6_only5$Gene_Tx),]
#df_megasheet_all_observations_R6_1709 <- df_megasheet_all_observations_R6_only5[!grepl('38610025', df_megasheet_all_observations_R6_only5$Gene_Tx),]

df_megasheet_all_observations_R6_only6 <- rbind(df_megasheet_all_observations_R6_1352, df_megasheet_all_observations_R6_1709, df_megasheet_all_observations_R6_1373, df_megasheet_all_observations_R6_1858)
#END PART I NEED TO KEEP

#Get rid of HIV
df_megasheet_all_observations_R6_only6 <- df_megasheet_all_observations_R6_only6 %>% filter(!grepl('HIV', Gene_Tx))

#Graph it
ggplot(df_megasheet_all_observations_R6_only6, aes(x=Gene_Tx, y=Mean_total_tx_each_replicate_per_ug, fill = Pulldown_date))+
  geom_bar (stat="identity", position=position_dodge2(preserve="single")) +
  ylab("Counts_per_Ug") +
  xlab("Ug RNA input (capture)")+
  theme(axis.text.x = element_text(angle=90, size=8, hjust = 1, vjust=0.5), axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"), panel.grid.major = element_line(size=0.5, linetype='solid', color="grey"),panel.grid.minor = element_line(size=0.25, linetype='solid', color="grey"),panel.spacing = unit(0, "lines"),panel.border = element_rect(color = "black", fill = NA, size = 1),strip.background = element_rect(color = "black", size = 1)) +
  
geom_errorbar(aes(ymin=Mean_total_tx_each_replicate_per_ug-SE_total_tx_each_replicate_per_ug, ymax=Mean_total_tx_each_replicate_per_ug+SE_total_tx_each_replicate_per_ug), position=position_dodge2(preserve="single")) +
  
  #geom_text(aes(label = Gene_Tx), vjust = 0.5, hjust = 4, angle = 90, position = position_dodge2(width =0.9, preserve = "single"), size = 2.5) +
  
  scale_x_discrete(guide=guide_axis(angle=90), labels = function(x) str_wrap(x, width=20))+
  scale_y_continuous(trans = 'pseudo_log',
                     labels = scales::number_format(accuracy=0.01))+
  facet_grid(~ PIC_ID, scales = "free_x", space="free_x") +
  theme(panel.background = element_rect(size = 5, fill = "white"), axis.ticks.length = unit(0.5, "cm")) +
  guides(shape = guide_legend(override.aes = list(size = 10)), colour = guide_legend(nrow=2)) +
  theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10),
        legend.key.size = unit(5, "mm"))


# Git repository set-up


# as.character/as.factor references ---------------------------------------
df_megasheet_all_observations_summarized_3_genetype_date$Replicate <- as.character(df_megasheet_all_observations_summarized_3_genetype_date$Replicate)

df_megasheet_all_observations_summarized_3_genetype_date$Pulldown_date <- as.character(df_megasheet_all_observations_summarized_3_genetype_date$Pulldown_date)

df_megasheet_all_observations_summarized_3_genetype_date$PIC_ID <- as.factor(df_megasheet_all_observations_summarized_3_genetype_date$PIC_ID)

df_megasheet_all_observations_summarized_3_genetype_date$Gene_Tx <- as.factor(df_megasheet_all_observations_summarized_3_genetype_date$Gene_Tx)

df_megasheet_all_observations_summarized_3_genetype_date$Pulldown_date <- as.numeric(factor((df_megasheet_all_observations_summarized_3_genetype_date$Pulldown_date)))

library(lubridate)

df_megasheet_all_observations_summarized_3_genetype_date$Pulldown_date <- as.Date(df_megasheet_all_observations_summarized_3_genetype_date$Pulldown_date, "%m/%d/%Y")


# #I can probably delete the following ------------------------------------


agg_tbl_megasheet <- df3_megasheet %>%
  select(c(PIC_ID, Pulldown_date, cDNA_primers, Ug_RNA_input_Capture, Gene_Tx, Replicate, Num_of_replicates, Sum_total_tx_each_replicate))
df4_megasheet <- agg_tbl_megasheet %>% as.data.frame() 

df4_megasheet$PIC_ID <- as.factor(df4_megasheet$PIC_ID)
df4_megasheet$Pulldown_date <- as.factor(df4_megasheet$Pulldown_date)
df4_megasheet$cDNA_primers <- as.factor(df4_megasheet$cDNA_primers)
df4_megasheet$Ug_RNA_input_Capture <- as.factor(df4_megasheet$Ug_RNA_input_Capture)
df4_megasheet$Gene_Tx <- as.factor(df4_megasheet$Gene_Tx)
df4_megasheet$Replicate <- as.factor(df4_megasheet$Replicate)
df4_megasheet$Num_of_replicates <- as.factor(df4_megasheet$Num_of_replicates)
df4_megasheet$Sum_total_tx_each_replicate <- as.numeric(df4_megasheet$Sum_total_tx_each_replicate)

agg_tbl_megasheet <- df4_megasheet %>%
group_by(PIC_ID, Pulldown_date, cDNA_primers, Ug_RNA_input_Capture, Gene_Tx, Replicate, Num_of_replicates, .drop=FALSE) %>% summarize(Sum_total_tx_each_replicate=sum(Sum_total_tx_each_replicate))
df5_megasheet <- agg_tbl_megasheet %>% as.data.frame() 

is.factor(agg_tbl_megasheet$Num_of_replicates)


write.csv(df3_megasheet, "/home/mshap1/hybrid_tx_validation/df3_megasheet.csv")

#code from Yeji Bae

df3_megasheet_test <- df3_megasheet %>%
  select(c(PIC_ID, Pulldown_date, cDNA_primers, Ug_RNA_input_Capture, Gene_Tx, Replicate, Num_of_replicates, Sum_total_tx_each_replicate))

length(unique(df3_megasheet_test[!duplicated(df3_megasheet_test),]$Gene_Tx)) # 1837
df3_megasheet_test_unique <- df3_megasheet_test[!duplicated(df3_megasheet_test),]

# df3_megasheet_test_unique[df3_megasheet_test_unique$Gene_Tx == "TCF3_1616348",]
dim(df3_megasheet_test_unique[df3_megasheet_test_unique$Gene_Tx == "TCF3_1616348", ])[1]
df3_megasheet_test_unique[df3_megasheet_test_unique$Gene_Tx == "TCF3_1616348", ]

for (i in df3_megasheet_test_unique$Gene_Tx){
  if (dim(df3_megasheet_test_unique[df3_megasheet_test_unique$Gene_Tx == i, ])[1] != df3_megasheet_test_unique[df3_megasheet_test_unique$Gene_Tx == i, ])
    
    #Yeji Bae: add rows (same for all the other columns (PIC_ID, Pulldown_date, cDNA_primers, Gene_Tx, Replicate, Num_of_replicates, Sum_total_tx_each_replicate) but 0 for Ug_RNA_input_Capture and all the Replicates base on the number of replicates)
    
    #Madelyn: Use something like this code?
    
    {complete(Gene_Tx, Replicate, fill = list(Sum_total_tx_each_replicate=0))
    
    #Yeji Bae: remove the redundant rows where all the other columns (PIC_ID, Pulldown_date, cDNA_primers, Gene_Tx, Replicate, Num_of_replicates, Sum_total_tx_each_replicate) are same but if Ug_RNA_input_Capture are different, remove rows where Ug_RNA_input_Capture == 0, 
   
    #Madelyn: use  something like this?
    unique(agg_tbl_megasheet)
     
  }
}

getwd()

agg_tbl_megasheet <-df3_megasheet %>%
group_by(PIC_ID, Pulldown_date, Ug_RNA_input_Capture, Gene_Tx) %>%
complete(fill = list(Sum_total_tx_each_replicate=0))
df4_megasheet <- agg_tbl_megasheet %>% as.data.frame() 

#closer...but this is making up replicates that aren't there...OS user comment also suggested using .drop = FALSE inside ddplyr. I tried but I got an error so I think I'm not using it in the right syntax.

df4_megasheet <- merge(expand.grid(lapply(df3_megasheet["Gene_Tx","Replicate"], unique)), df4_megasheet, all.x=TRUE)
is.na(df4_megasheet) <- df3_megasheet==0 

df3_megasheet %>%
group_by(PIC_ID, Pulldown_date, Gene_Tx, Replicate, .drop=FALSE) %>%
summarise(mean = mean(Sum_total_tx_each_replicate)) 
df4_megasheet <- df3_megasheet %>% as.data.frame()

#First get count, try to add rows for 0 observations in the other replicates (still working on getting the .drop=false rows to work):


# #probably necessary code to complete my task with the megasheets and statistics ----------------------------



library(dplyr)

#Only have random hexamers in this graph
df_megasheet_all_observations_summarized3
agg_tbl_megasheet <- df3_megasheet %>%
  filter(cDNA_primers == "R6")
df_megasheet_all_observations_summarized_R6 <- agg_tbl_megasheet %>% as.data.frame()

#Normalize by counts/ug, tried adding .drop=FALSE rows here too but it didn't work. :(
library(dplyr)
df_megasheet_all_observations_summarized_R6
agg_tbl_megasheet <-df_megasheet_all_observations_summarized_R6 %>%
group_by(Sum_total_tx_each_replicate, Ug_RNA_input_Capture) %>%
mutate(Sum_total_tx_each_replicate_per_Ug = Sum_total_tx_each_replicate/Ug_RNA_input_Capture)
df_megasheet_all_observations_summarized_R6_2 <- agg_tbl_megasheet %>% as.data.frame()

#crap. I made the Ug_RNA_input_Capture a factor so it won't compute.

#I can get means as below but how do I do standard error if counts show up in 1 replicate but not the other two? It's not like we can make 0 = 0.000001 because there isn't even a list for them, though maybe there is a way I can get the formula to assume that if only one replicate shows counts (when there are 3 replicates), then I can assume the other two replicates are 0.00001?

setDT(df_megasheet_all_observations_summarized_R6_2)
df_megasheet_all_observations_summarized_R6_2[, Mean_total_tx_each_replicate_per_ug :=mean(Sum_total_tx_each_replicate_per_Ug), by = c('Gene_Tx', 'Pulldown_date')]

setDT(df_megasheet_all_observations_summarized_R6_2)
df_megasheet_all_observations_summarized_R6_2[, SE_total_tx_each_replicate_per_ug :=std.error(Sum_total_tx_each_replicate_per_Ug), by = c('Gene_Tx', 'Pulldown_date')]

#write.csv(df6_megasheet, "df6_megasheet.csv")

#Select only columns I need
library(dplyr)
library(tidyverse)

df_megasheet_all_observations_summarized_R6_2
agg_tbl_megasheet <-df_megasheet_all_observations_summarized_R6_2 %>%
select(c(PIC_ID, Gene_Tx, Pulldown_date, Sum_total_tx_each_replicate, Sum_total_tx_each_replicate_per_Ug, Gene_Tx, Replicate)) 
df_megasheet_all_observations_summarized_R6_3 <- agg_tbl_megasheet %>% as.data.frame()

df_megasheet_all_observations_summarized_R6_4 <- df_megasheet_all_observations_summarized_R6_3 %>% distinct()

agg_tbl_megasheet <-df_megasheet_all_observations_summarized_R6_4 %>%
  group_by(PIC_ID, Gene_Tx, Pulldown_date, Sum_total_tx_each_replicate) %>% 
  mutate(Sum_tx_of_both_replicates =sum(Sum_total_tx_each_replicate))
df_megasheet_all_observations_summarized_R6_5 <- agg_tbl_megasheet %>% as.data.frame()

agg_tbl_megasheet <-df_megasheet_all_observations_summarized_R6_5 %>%
  group_by(Gene_Tx, Sum_total_tx_each_replicate) %>% 
  filter(n() >3)
df_megasheet_all_observations_summarized_R6_5 <- agg_tbl_megasheet %>% as.data.frame()

#Get rid of HIV
df_megasheet_all_observations_summarized_R6_5 <- df_megasheet_all_observations_summarized_R6_5 %>% filter(!grepl('HIV', Gene_Tx))

#Take just genes I want
#df8_megasheet <- df7_megasheet %>% filter(grepl("5LTR_MSD_STAT5B_4223", Gene_Tx) | grepl("FANCA", Gene_Tx) | grepl("EIF", Gene_Tx) | grepl("STAT5B_4227", Gene_Tx) | grepl("CDK13", Gene_Tx) | grepl("NDC80", Gene_Tx) | grepl("MAP4K1", Gene_Tx),)

#df9_1709_3 <- df9_1709_3[!grepl('38611132', df9_1709_3$Gene_Tx),]
#df9_1709_3 <- df9_1709_3[!grepl('42257672', df9_1709_3$Gene_Tx),]
#df9_1709_3 <- df9_1709_3[!grepl('HIVE', df9_1709_3$Gene_Tx),]
#df9_1709_3 <- df9_1709_3[!grepl('38610025', df9_1709_3$Gene_Tx),]


#other way of selecting, filtering by frequency it shows up between experiment dates
#gene_frequency <- table(df8_megasheet$Gene_Tx)
#duplicate_genes <- names(gene_frequency)[gene_frequency > 1]

#Graph it
ggplot(df_megasheet_all_observations_summarized_R6_5, aes(x=Gene_Tx, y=Sum_total_tx_each_replicate_per_Ug, fill = Pulldown_date))+
  geom_dotplot (binaxis = "y", stackdir = "center", binpositions="all", binwidth = .5, dotsize =  2.5, stackratio = 1.5)+
  ylab("Counts_per_Ug")+
  xlab("Gene_Location")+
  theme(axis.text.x = element_text(angle=90, size=8, hjust = 1, vjust=0.5), axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"), panel.grid.major = element_line(size=0.5, linetype='solid', color="grey"),panel.grid.minor = element_line(size=0.25, linetype='solid', color="grey"),panel.spacing = unit(0, "lines"),panel.border = element_rect(color = "black", fill = NA, size = 1),strip.background = element_rect(color = "black", size = 1)) +
 
   #geom_errorbar(aes(ymin=Mean_total_tx_each_replicate_per_ug-SE_total_tx_each_replicate_per_ug, ymax=Mean_total_tx_each_replicate_per_ug+SE_total_tx_each_replicate_per_ug), position=position_dodge2(preserve="single")) +
  
  #geom_text(aes(label = Gene_Tx), vjust = 0.5, hjust = 4, angle = 90, position = position_dodge2(width =0.9, preserve = "single"), size = 2.5) +
  
  scale_x_discrete(guide=guide_axis(angle=90), labels = function(x) str_wrap(x, width=20))+
  scale_y_continuous()+
  facet_grid(~ PIC_ID, scales = "free_x", space="free_x") +
  theme(panel.background = element_rect(size = 5, fill = "white"), axis.ticks.length = unit(0.5, "cm")) +
  guides(shape = guide_legend(override.aes = list(size = 10)), colour = guide_legend(nrow=2)) +
  theme(legend.title = element_text(size = 10), 
              legend.text  = element_text(size = 10),
              legend.key.size = unit(5, "mm"))

# Various codes either not working or not being used. ---------------------
  
  df_replication2 <- df_replication %>% select("PID", "Number_of_replicates", "Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA", "Approx_cDNA_HIV_copy_input_per_replicate_determined_via_LTR_qPCR_into_sonication", "Mean_#_3LTR-human_sequences_per_Replicate", "Stdev_#_3LTR-human_sequences_per_Replicate","Total_#_3LTR-human_sequences","")
  
  geom_text(aes(label = Ug_RNA_input_Capture, position = position_dodge(width = .9))+
  
  geom_text(aes(label = Ug_RNA_input_Capture, position = position_dodge(width = .9)))+
              
   scale_x_discrete(labels = function(x) str_wrap(x, width=20)) +
   
   flush_ticks <- function(gg, flush=XY, plot=TRUE, cat=TRUE) +
 
 geom_text(aes(label = Ug_RNA_input_Capture), vjust = 0.5, colour = "black", angle = 0, size = 3) +

 geom_text(aes(label = Ug_RNA_input_Capture), vjust = 0.5, colour = "black", angle = 0, size = 3) +
 
   coord_flip() +
   
  coord_cartesian(clip = "off")+


          
  my_wildcard <- "*"

sheets <- excel_sheets("combined_all_sheets.xlsx")
sheets <- sheets[grep("HIV", sheets)]

excel_data <- lapply(sheets, read_excel, path = "combined_all_sheets.xlsx", range = cell_cols("A:M"))

df13 <- rbindlist(excel_data)

pwd( )

write.csv(df13, file = "HIV_only_tx.csv')

write.csv(df5_1373_3, file = "df5_1373_3.csv")

library(ggplot2)
# Basic scatter plot
# Change the point size, and shape

ggplot(df_replication2, aes(y = Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA, x = Total_HIV_sequences_detected_including_hybrid_transcript/Number_of_replicates, color = PID)) +
  ylab("Total Ug RNA Tested per Replicate") +
  xlab("HIV sequences detected per replicate")+
  geom_point(data = subset(df_replication2, PID %in% c("PIC1373","PIC1709","PIC1858")))+
  geom_smooth(method = "lm", se=FALSE)
  geom_abline(intercept=-5.5, slope=1/2, color="red", linetype="dashed") +
  geom_abline(intercept=-33, slope=1/2, color="blue", linetype="dashed")

coef(lm(subset(df_replication2, PID %in% c("PIC1373","PIC1709","PIC1858"))$Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA)~(subset(df_replication2, PID %in% c("PIC1373","PIC1709","PIC1858")$Total_HIV_sequences_detected_including_hybrid_transcript)/(subset(df_replication2, PID %in% c("PIC1373","PIC1709","PIC1858")$Number_of_replicates))))[1]

ggplot(subset(df_replication2, PID %in% "PIC1373"), aes(Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA, Total_HIV_sequences_detected_per_ug_RNA_input)) +
  geom_point(color = "blue", size = 3) +
ggplot(subset(df_replication2, PID %in% "PIC1709"), aes(Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA, Total_HIV_sequences_detected_per_ug_RNA_input)) +
  geom_point(color = "red", size = 3)
# Change point shapes, colors and sizes
ggplot(subset(df_replication2, PID %in% "PIC1858"), aes(Amount_of_RNA_used_per_replicate_in_1st_strand_cDNA, Total_HIV_sequences_detected_per_ug_RNA_input)) +
  geom_point(color = "green", size = 3)




# Percent transcript types in 1373 ----------------------------------------





#Get rid of repetitive regions.
library(dplyr)
df_1373_3
agg_tbl_1373_3 <- df_1373_3 %>%
  filter(Normal_or_repetitive == "Normal")
df2_1373_3 <- agg_tbl_1373_3 %>% as.data.frame()

#Count transcripts that showed up by Gene_loc, date, input RNA, and replicate
library(dplyr)
df2_1373_3
agg_tbl_1373_3 <- df2_1373_3 %>%
  group_by(Gene_Tx, Pulldown_date, Replicate, Ug_RNA_input_Capture) %>%
  mutate(Count_Total = n())
df3_1373_3 <- agg_tbl_1373_3 %>% as.data.frame()

#Catagorize transcripts by transcript type
agg_tbl_1373_3 <- df3_1373_3 %>%
  mutate(df3_1373_3, Gene_Type=ifelse(grepl("3LTR", df3_1373_3$ID), "3LTR",
                                      ifelse(grepl("5LTR", df3_1373_3$ID), "5LTR_MSD","HIV-only")))
df3_1373_3 <- agg_tbl_1373_3 %>% as.data.frame()

#I want to group all the "unique" transcripts together.
#Change Gene_Tx names so that those with only 1 count are Unique

df3_1373_3$Gene_Tx2 <- "Unique"
rowsIwant <- df3_1373_3$Count_Total>1
df3_1373_3$Gene_Tx2[rowsIwant] <- as.character(df3_1373_3$Gene_Tx[rowsIwant]) 

#Get overall transcript distribution. Percentage. Log10 on Y axis.
library(scales)

#New counts (with Unique values)
df3_1373_3
agg_tbl_1373_3 <- df3_1373_3 %>%
  group_by(Gene_Tx2, Gene_Type) %>%
  mutate(Count_Total2 = n())
df3c_1373_3 <- agg_tbl_1373_3 %>% as.data.frame()

#Get select only the columns needed and get rid of duplicates
df3c_1373_3 <- df3c_1373_3 %>% select("Gene_Tx2", "Count_Total2","Gene_Type")
df3d_1373_3 <- df3c_1373_3 %>% distinct()

#Get percentages
df3d_1373_3 <- df3d_1373_3 %>% mutate(Percentage = Count_Total2/sum(Count_Total2))
df3d_1373_3 <- df3d_1373_3[!(df3d_1373_3$Gene_Type=="HIV-only" & df3d_1373_3$Gene_Tx2=="Unique"),]

#Graph
ggplot(df3d_1373_3, aes(x=Gene_Tx2, y=Percentage, fill=Gene_Type)) +
  geom_bar (stat="identity", position=position_dodge2(preserve="single"), width = NULL) +
  ylab("% Total transcript distribution") +
  theme(axis.text.x = element_text(angle=90, size=8, hjust = 1, vjust=0.5), axis.line.x = element_line(color="black"),            axis.line.y = element_line(color="black"), panel.grid.major = element_line(size=0.5, linetype='solid', color="grey"),    panel.grid.minor = element_line(size=0.25, linetype='solid', color="grey")) +
  scale_fill_discrete(name="HIV Transcript Type")+
  scale_x_discrete(labels = function(x) str_wrap(x, width=30), name="Gene_location") +
  scale_y_continuous(labels = scales::percent_format(accuracy=1))+
  theme(panel.background = element_rect(fill = "white")) +
  geom_text(aes(label = Count_Total2), vjust = -1, colour = "black", angle = 0, size = 3) +
  guides(shape = guide_legend(override.aes = list(size = 10)), color = guide_legend(override.aes = list(size = 10))) +
  guides(shape = guide_legend(override.aes = list(size = 10)), color = guide_legend(override.aes = list(size = 10))) +
  theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10),
        legend.key.size = unit(6, "mm"))




# # replicability of assay in 1373 ------------------------------- --------



#Now determine replicability of assay by testing to see if genes with the highest counts are spread out over multiple experiments.

#10/3/2023 to do: 
#1. Make dot plots. Done 
#2. Get rid of input as a factor. We discovered that there isn't a ton of inhibition. Done.
#3. Replot this with the megasheet and do facet wrapping to display from all 4 PIC.

library(plotrix)

#First get counts
library(dplyr)
df3_1373_3
agg_tbl_1373_3 <-df3_1373_3 %>%
group_by(Pulldown_date, Ug_RNA_input_Capture, Gene_Tx, Replicate) %>%
mutate(Sum_total_tx_each_replicate =n())
df4_1373_3 <- agg_tbl_1373_3 %>% as.data.frame()

#Only have random hexamers in this graph
df4_1373_3
agg_tbl_1373_3 <- df4_1373_3 %>%
  filter(cDNA_primers == "R6")
df4_1373_3 <- agg_tbl_1373_3 %>% as.data.frame()

#Normalize by counts/ug
library(dplyr)
df4_1373_3
agg_tbl_1373_3 <-df4_1373_3 %>%
group_by(Sum_total_tx_each_replicate, Ug_RNA_input_Capture) %>%
mutate(Sum_total_tx_each_replicate_per_Ug = Sum_total_tx_each_replicate/Ug_RNA_input_Capture)
df5_1373_3 <- agg_tbl_1373_3 %>% as.data.frame()

#write.csv(df6_1373_3, "df6_1373_3.csv")

#Select only columns I need
library(dplyr)
library(tidyverse)

df5_1373_3
agg_tbl_1373_3 <-df5_1373_3 %>%
select(c(Gene_Tx, Pulldown_date, Sum_total_tx_each_replicate_per_Ug, Gene_Tx, Replicate, Count_Total)) 
df6_1373_3 <- agg_tbl_1373_3 %>% as.data.frame()

df6_1373_3 <- df6_1373_3 %>% distinct()

setDT(df6_1373_3)
df6_1373_3[, Mean_total_tx_each_replicate_per_ug :=mean(Sum_total_tx_each_replicate_per_Ug), by = c('Gene_Tx', 'Pulldown_date')]

setDT(df6_1373_3)
df6_1373_3[, SE_total_tx_each_replicate_per_ug :=std.error(Sum_total_tx_each_replicate_per_Ug), by = c('Gene_Tx', 'Pulldown_date')]

library(dplyr)
#df8_1373_3 <- df6_1373_3 %>% filter(!is.na(SE_total_tx_each_replicate_per_ug))

#Select all columns except "replicates" and "counts"
df6_1373_3
agg_tbl_1373_3 <-df6_1373_3 %>%
  select(c(Gene_Tx, Pulldown_date, Sum_total_tx_each_replicate_per_Ug, Replicate, Count_Total, Mean_total_tx_each_replicate_per_ug, SE_total_tx_each_replicate_per_ug)) 
df8_1373_3 <- agg_tbl_1373_3 %>% as.data.frame()

#Get rid of duplicates again
df8_1373_3 <- df8_1373_3 %>% distinct()

#Take just genes I want
df8_1373_3 <- df8_1373_3 %>% filter(grepl("5LTR_MSE_STAT5B_4223", Gene_Tx) | grepl("FANCA", Gene_Tx) | grepl("EIF", Gene_Tx) | grepl("STAT5B_4227", Gene_Tx) | grepl("CDK13", Gene_Tx) | grepl("NDC80", Gene_Tx))


#other way of selecting, filtering by frequency it shows up between experiment dates
#gene_frequency <- table(df8_1373_3$Gene_Tx)
#duplicate_genes <- names(gene_frequency)[gene_frequency > 1]

#Graph it
ggplot(df8_1373_3, aes(x=Gene_Tx, y=Sum_total_tx_each_replicate_per_Ug, fill = Pulldown_date)) +
  geom_dotplot (binaxis = "y", stackdir = "center", binpositions="all", binwidth = 0.75, dotsize =  0.5)+
  ylab("Counts_per_Ug") +
  xlab("Gene_Location")+
  theme(axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), axis.line.x = element_line(color="black"),            axis.line.y = element_line(color="black"), panel.grid.major = element_line(size=0.5, linetype='solid', color="grey"),    panel.grid.minor = element_line(size=0.25, linetype='solid', color="grey")) +
 
   #geom_errorbar(aes(ymin=Mean_total_tx_each_replicate_per_ug-SE_total_tx_each_replicate_per_ug, ymax=Mean_total_tx_each_replicate_per_ug+SE_total_tx_each_replicate_per_ug), position=position_dodge2(preserve="single")) +
  
  #geom_text(aes(label = Gene_Tx), vjust = 0.5, hjust = 4, angle = 90, position = position_dodge2(width =0.9, preserve = "single"), size = 2.5) +
  
  scale_x_discrete(guide=guide_axis(angle=90), labels = function(x) str_wrap(x, width=20))
  
  scale_y_continuous()+
  facet_grid(~ Gene_Tx) +
  theme(panel.background = element_rect(size = 5, fill = "white"), axis.ticks.length = unit(0.5, "cm")) +
  guides(shape = guide_legend(override.aes = list(size = 10)), colour = guide_legend(nrow=2)) +
  theme(legend.title = element_text(size = 10), 
              legend.text  = element_text(size = 10),
              legend.key.size = unit(5, "mm"))



# Percent distribution  of transcripts in  1709 ---------------------------
#1709 hybrid transcript analysis
df_1709_3 <- read.csv("hybrid_tx_11Aug2023_3LTR_1709.csv")


#Get rid of repetitive regions
library(dplyr)
df_1373_3
agg_tbl_1373_3 <- df_1373_3 %>%
  filter(Normal_or_repetitive == "Normal")
df2_1373_3 <- agg_tbl_1373_3 %>% as.data.frame()

#Count transcripts that showed up by Gene_loc, date, input RNA, and replicate
library(dplyr)
df2_1709_3
agg_tbl_1709_3 <- df2_1709_3 %>%
  group_by(Gene_Tx, Pulldown_date, Replicate, Ug_RNA_input_Capture) %>%
  mutate(Count_Total = n())
df3_1709_3 <- agg_tbl_1709_3 %>% as.data.frame()


agg_tbl_1709_3 <- df3_1709_3 %>%
  mutate(df3_1709_3, Gene_Type=ifelse(grepl("3LTR", df3_1709_3$ID), "3LTR",
                                      ifelse(grepl("5LTR", df3_1709_3$ID), "5LTR_MSD","HIV-only")))
df3_1709_3 <- agg_tbl_1709_3 %>% as.data.frame()

#I want to group all the "unique" transcripts together.
#Change Gene_Tx names so that those with only 1 count are Unique

df3_1709_3$Gene_Tx2 <- "Unique"
rowsIwant <- df3_1709_3$Count_Total>1
df3_1709_3$Gene_Tx2[rowsIwant] <- as.character(df3_1709_3$Gene_Tx[rowsIwant]) 

#Get overall transcript distribution. Percentage. Log10 on Y axis.
library(scales)

#New counts (with Unique values)
df3_1709_3
agg_tbl_1709_3 <- df3_1709_3 %>%
  group_by(Gene_Tx2, Gene_Type) %>%
  mutate(Count_Total2 = n())
df3c_1709_3 <- agg_tbl_1709_3 %>% as.data.frame()

#Get select only the columns needed and get rid of duplicates
df3c_1709_3 <- df3c_1709_3 %>% select("Gene_Tx2", "Count_Total2","Gene_Type")
df3d_1709_3 <- df3c_1709_3 %>% distinct()

#Get percentages
df3d_1709_3 <- df3d_1709_3 %>% mutate(Percentage = Count_Total2/sum(Count_Total2))

#Got rid of 1 count that was messing up my graph and not contributing much data
#df3d_1709_3 <- df3d_1709_3[!(df3d_1709_3$Gene_Type=="HIV-only" & df3d_1709_3$Gene_Tx2=="Unique"),]

#Graph
g2 <- ggplot(df3d_1709_3, aes(x=Gene_Tx2, y=Percentage, fill=Gene_Type)) +
  geom_bar (stat="identity", position=position_dodge2(preserve="single"), width = NULL) +
  ylab("% Total transcript distribution") +
  theme(axis.text.x = element_text(angle=90, size=8, hjust = 1, vjust=0.5), axis.line.x = element_line(color="black"),            axis.line.y = element_line(color="black"), panel.grid.major = element_line(size=0.5, linetype='solid', color="grey"),    panel.grid.minor = element_line(size=0.25, linetype='solid', color="grey")) +
  scale_fill_discrete(name="HIV Transcript Type")+
  scale_x_discrete(labels = function(x) str_wrap(x, width=30), name="Gene_location") +
  scale_y_continuous(labels = scales::percent_format(accuracy=1))+
  theme(panel.background = element_rect(fill = "white")) +
  geom_text(aes(label = Count_Total2), vjust = -1, colour = "black", angle = 0, size = 3) +
  guides(shape = guide_legend(override.aes = list(size = 10)), color = guide_legend(override.aes = list(size = 10))) +
  guides(shape = guide_legend(override.aes = list(size = 10)), color = guide_legend(override.aes = list(size = 10))) +
  theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10),
        legend.key.size = unit(6, "mm"))


# replicability of assay in 1709 ----------------------------------------



#Now determine replicability of assay by testing to see if genes with the highest counts are spread out over multiple experiments.

#First get counts
df5_1709_3
agg_tbl_1709_3 <-df5_1709_3 %>%
  select(c(Gene_Tx, Pulldown_date, Ug_RNA_input_Capture, Sum_total_tx_each_replicate_per_Ug, Gene_Tx, Replicate, Count_Total)) 
df6_1709_3 <- agg_tbl_1709_3 %>% as.data.frame()

df6_1709_3 <- df6_1709_3 %>% distinct()

setDT(df6_1709_3)
df6_1709_3[, Mean_total_tx_each_replicate_per_ug :=mean(Sum_total_tx_each_replicate_per_Ug), by = c('Gene_Tx', 'Pulldown_date','Ug_RNA_input_Capture')]

setDT(df6_1709_3)
df6_1709_3[, SE_total_tx_each_replicate_per_ug :=std.error(Sum_total_tx_each_replicate_per_Ug), by = c('Gene_Tx', 'Pulldown_date','Ug_RNA_input_Capture')]


library(dplyr)
#df8_1709_3 <- df6_1709_3 %>% filter(!is.na(SD_total_tx_each_replicate_per_ug))

#Select all columns except "replicates" and "counts"
df6_1709_3
agg_tbl_1709_3 <-df6_1709_3 %>%
  select(c(Gene_Tx, Pulldown_date, Ug_RNA_input_Capture, Gene_Tx, Mean_total_tx_each_replicate_per_ug, SE_total_tx_each_replicate_per_ug)) 
df8_1709_3 <- agg_tbl_1709_3 %>% as.data.frame()

#Get rid of duplicates again
df8_1709_3 <- df8_1709_3 %>% distinct()


#Take just genes I want
df9_1709_3 <- df8_1709_3 %>% filter(grepl("5LTR_MSD_STAT5B", Gene_Tx) | grepl("MAP4K1", Gene_Tx),)
df9_1709_3 <- df9_1709_3[!grepl('38611132', df9_1709_3$Gene_Tx),]
df9_1709_3 <- df9_1709_3[!grepl('42257672', df9_1709_3$Gene_Tx),]
df9_1709_3 <- df9_1709_3[!grepl('HIVE', df9_1709_3$Gene_Tx),]
df9_1709_3 <- df9_1709_3[!grepl('38610025', df9_1709_3$Gene_Tx),]

#df8_1709_3 <- df8_1709_3[!grepl('1-Feb-23', df8_1709_3$Pulldown_date),]

#other way of selecting, filtering by frequency it shows up between experiment dates
#gene_frequency <- table(df8_1709_3$Gene_Tx)
#duplicate_genes <- names(gene_frequency)[gene_frequency > 1]

#Graph it
ggplot(df9_1709_3, aes(x=as.character(Ug_RNA_input_Capture), y=Mean_total_tx_each_replicate_per_ug, fill = Pulldown_date)) +
  geom_bar (stat="identity", position=position_dodge2(preserve="single")) +
  ylab("Counts_per_Ug") +
  xlab("Ug RNA input (capture)")+
  theme(axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), axis.line.x = element_line(color="black"),            axis.line.y = element_line(color="black"), panel.grid.major = element_line(size=0.5, linetype='solid', color="grey"),    panel.grid.minor = element_line(size=0.25, linetype='solid', color="grey")) +
  
  geom_errorbar(aes(ymin=Mean_total_tx_each_replicate_per_ug-SE_total_tx_each_replicate_per_ug, ymax=Mean_total_tx_each_replicate_per_ug+SE_total_tx_each_replicate_per_ug), position=position_dodge2(preserve="single")) +
  
  #geom_text(aes(label = Gene_Tx), vjust = 0.5, hjust = 4, angle = 90, position = position_dodge2(width =0.9, preserve = "single"), size = 2.5) +
  
  scale_x_discrete(guide=guide_axis(angle=90), labels = function(x) str_wrap(x, width=20))+
  
  scale_y_continuous()+
  facet_grid(~ Gene_Tx) +
  theme(panel.background = element_rect(size = 5, fill = "white"), axis.ticks.length = unit(0.5, "cm")) +
  guides(shape = guide_legend(override.aes = list(size = 10)), color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10),
        legend.key.size = unit(5, "mm"))


# Percent distribution of transcripts in 1858 -----------------------------


#1858 hybrid transcript analysis
df_1858_3 <- read.csv("hybrid_tx_11Aug2023_3LTR_1858.csv")
#Get rid of repetitive regions
library(dplyr)
df_1373_3
agg_tbl_1373_3 <- df_1373_3 %>%
  filter(Normal_or_repetitive == "Normal")
df2_1373_3 <- agg_tbl_1373_3 %>% as.data.frame()

#Count transcripts that showed up by Gene_loc, date, input RNA, and replicate 
library(dplyr)
df2_1858_3
agg_tbl_1858_3 <- df2_1858_3 %>%
  group_by(Gene_Tx, Pulldown_date, Replicate, Ug_RNA_input_Capture) %>%
  mutate(Count_Total = n())
df3_1858_3 <- agg_tbl_1858_3 %>% as.data.frame()


agg_tbl_1858_3 <- df3_1858_3 %>%
  mutate(df3_1858_3, Gene_Type=ifelse(grepl("3LTR", df3_1858_3$ID), "3LTR",
                                      ifelse(grepl("5LTR", df3_1858_3$ID), "5LTR_MSD","HIV-only")))
df3_1858_3 <- agg_tbl_1858_3 %>% as.data.frame()

#I want to group all the "unique" transcripts together.
#Change Gene_Tx names so that those with only 1 count are Unique

df3_1858_3$Gene_Tx2 <- "Unique"
rowsIwant <- df3_1858_3$Count_Total>1
df3_1858_3$Gene_Tx2[rowsIwant] <- as.character(df3_1858_3$Gene_Tx[rowsIwant]) 

#Get overall transcript distribution. Percentage. Log10 on Y axis.
library(scales)

#New counts (with Unique values)
df3_1858_3
agg_tbl_1858_3 <- df3_1858_3 %>%
  group_by(Gene_Tx2, Gene_Type) %>%
  mutate(Count_Total2 = n())
df3c_1858_3 <- agg_tbl_1858_3 %>% as.data.frame()

#Get select only the columns needed and get rid of duplicates
df3c_1858_3 <- df3c_1858_3 %>% select("Gene_Tx2", "Count_Total2","Gene_Type")
df3d_1858_3 <- df3c_1858_3 %>% distinct()

#Get percentages
df3d_1858_3 <- df3d_1858_3 %>% mutate(Percentage = Count_Total2/sum(Count_Total2))

#Got rid of 1 count that was messing up my graph and not contributing much data
#df3d_1858_3 <- df3d_1858_3[!(df3d_1858_3$Gene_Type=="HIV-only" & df3d_1858_3$Gene_Tx2=="Unique"),]

#Graph
g3 <- ggplot(df3d_1858_3, aes(x=Gene_Tx2, y=Percentage, fill=Gene_Type)) +
  geom_bar (stat="identity", position=position_dodge2(preserve="single"), width = NULL) +
  ylab("% Total transcript distribution") +
  theme(axis.text.x = element_text(angle=90, size=8, hjust = 1, vjust=0.5), axis.line.x = element_line(color="black"),            axis.line.y = element_line(color="black"), panel.grid.major = element_line(size=0.5, linetype='solid', color="grey"),    panel.grid.minor = element_line(size=0.25, linetype='solid', color="grey")) +
  scale_fill_discrete(name="HIV Transcript Type")+
  scale_x_discrete(labels = function(x) str_wrap(x, width=30), name="Gene_location") +
  scale_y_continuous(labels = scales::percent_format(accuracy=1))+
  theme(panel.background = element_rect(fill = "white")) +
  geom_text(aes(label = Count_Total2), vjust = -1, colour = "black", angle = 0, size = 3) +
  guides(shape = guide_legend(override.aes = list(size = 10)), color = guide_legend(override.aes = list(size = 10))) +
  guides(shape = guide_legend(override.aes = list(size = 10)), color = guide_legend(override.aes = list(size = 10))) +
  theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10),
        legend.key.size = unit(6, "mm"))

g3 <- ggplot(df3d_1858_3, aes(x=Gene_Tx2, y=Percentage, fill=Gene_Type)) +
  geom_bar (stat="identity", position=position_dodge2(preserve="single"), width = NULL)

geom_text(aes(label = Count_Total2), vjust = -1, colour = "black", angle = 0, size = 3)


# replicability of assay in 1858 ------------------------------------------


#Now determine replicability of assay by testing to see if genes with the highest counts are spread out over multiple experiments.

#First get counts
library(dplyr)
df3_1858_3
agg_tbl_1858_3 <-df3_1858_3 %>%
  group_by(Pulldown_date, Ug_RNA_input_Capture, Gene_Tx, Replicate) %>%
  mutate(Sum_total_tx_each_replicate =n())
df4_1858_3 <- agg_tbl_1858_3 %>% as.data.frame()

#Only have random hexamers in this graph
df4_1858_3
agg_tbl_1858_3 <- df4_1858_3 %>%
  filter(cDNA_primers == "R6")
df4_1858_3 <- agg_tbl_1858_3 %>% as.data.frame()

#Normalize by counts/ug
library(dplyr)
df4_1858_3
agg_tbl_1858_3 <-df4_1858_3 %>%
  group_by(Sum_total_tx_each_replicate, Ug_RNA_input_Capture) %>%
  mutate(Sum_total_tx_each_replicate_per_Ug = Sum_total_tx_each_replicate/Ug_RNA_input_Capture)
df5_1858_3 <- agg_tbl_1858_3 %>% as.data.frame()

#write.csv(df6_1858_3, "df6_1858_3.csv")

#Select only columns I need
library(dplyr)
library(tidyverse)

df5_1858_3
agg_tbl_1858_3 <-df5_1858_3 %>%
  select(c(Gene_Tx, Pulldown_date, Ug_RNA_input_Capture, Sum_total_tx_each_replicate_per_Ug, Gene_Tx, Replicate, Count_Total)) 
df6_1858_3 <- agg_tbl_1858_3 %>% as.data.frame()

df6_1858_3 <- df6_1858_3 %>% distinct()

setDT(df6_1858_3)
df6_1858_3[, Mean_total_tx_each_replicate_per_ug :=mean(Sum_total_tx_each_replicate_per_Ug), by = c('Gene_Tx', 'Pulldown_date','Ug_RNA_input_Capture')]

setDT(df6_1858_3)
df6_1858_3[, SE_total_tx_each_replicate_per_ug :=std.error(Sum_total_tx_each_replicate_per_Ug), by = c('Gene_Tx', 'Pulldown_date','Ug_RNA_input_Capture')]

library(dplyr)
#df8_1858_3 <- df6_1858_3 %>% filter(!is.na(SD_total_tx_each_replicate_per_ug))

#Select all columns except "replicates" and "counts"
df6_1858_3
agg_tbl_1858_3 <-df6_1858_3 %>%
  select(c(Gene_Tx, Pulldown_date, Ug_RNA_input_Capture, Gene_Tx, Mean_total_tx_each_replicate_per_ug, SE_total_tx_each_replicate_per_ug)) 
df8_1858_3 <- agg_tbl_1858_3 %>% as.data.frame()

#Get rid of duplicates again
df8_1858_3 <- df8_1858_3 %>% distinct()

#Take just genes I want
df8_1858_3 <- df8_1858_3 %>% filter(grepl("5LTR_MSD_ST", Gene_Tx) | grepl("FLT", Gene_Tx)| grepl("VPR", Gene_Tx))


#other way of selecting, filtering by frequency it shows up between experiment dates
#gene_frequency <- table(df8_1858_3$Gene_Tx)
#duplicate_genes <- names(gene_frequency)[gene_frequency > 1]

#Graph it
ggplot(df8_1858_3, aes(x=as.character(Ug_RNA_input_Capture), y=Mean_total_tx_each_replicate_per_ug, fill = Pulldown_date)) +
  geom_bar (stat="identity", position=position_dodge2(preserve="single")) +
  ylab("Counts_per_Ug") +
  xlab("Ug RNA input (capture)")+
  theme(axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), axis.line.x = element_line(color="black"),            axis.line.y = element_line(color="black"), panel.grid.major = element_line(size=0.5, linetype='solid', color="grey"),    panel.grid.minor = element_line(size=0.25, linetype='solid', color="grey")) +
  geom_errorbar(aes(ymin=Mean_total_tx_each_replicate_per_ug-SE_total_tx_each_replicate_per_ug, ymax=Mean_total_tx_each_replicate_per_ug+SE_total_tx_each_replicate_per_ug), position=position_dodge2(preserve="single")) +
  
  #geom_text(aes(label = Gene_Tx), vjust = 0.5, hjust = 4, angle = 90, position = position_dodge2(width =0.9, preserve = "single"), size = 2.5) +
  
  scale_x_discrete(guide=guide_axis(angle=90), labels = function(x) str_wrap(x, width=20))+
  
  scale_y_continuous()+
  facet_grid(~ Gene_Tx) +
  theme(panel.background = element_rect(size = 5, fill = "white"), axis.ticks.length = unit(0.5, "cm")) +
  guides(shape = guide_legend(override.aes = list(size = 10)), color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10),
        legend.key.size = unit(5, "mm"))

library(RColorBrewer)
library(matrixStats)


# oligo(dT)s vs R6 1373 ------------------
#First get counts
library(dplyr)
df3_1373_3
agg_tbl_1373_3 <-df3_1373_3 %>%
  group_by(Pulldown_date, Ug_RNA_input_Capture, cDNA_primers, Gene_Tx, Replicate) %>%
  mutate(Sum_total_tx_each_replicate =n())
df4_1373_3 <- agg_tbl_1373_3 %>% as.data.frame()


library(dplyr)
df4_1373_3
agg_tbl_1373_3 <-df4_1373_3 %>%
  group_by(Sum_total_tx_each_replicate, Ug_RNA_input_Capture) %>%
  mutate(Sum_total_tx_each_replicate_per_Ug = Sum_total_tx_each_replicate/Ug_RNA_input_Capture)
df5_1373_3 <- agg_tbl_1373_3 %>% as.data.frame()

#write.csv(df6_1373_3, "df6_1373_3.csv")

#Select only columns I need
library(dplyr)
library(tidyverse)

df5_1373_3
agg_tbl_1373_3 <-df5_1373_3 %>%
  select(c(Gene_Tx, Pulldown_date, Ug_RNA_input_Capture, Sum_total_tx_each_replicate_per_Ug, cDNA_primers, Gene_Tx, Replicate, Count_Total)) 
df6_1373_3 <- agg_tbl_1373_3 %>% as.data.frame()

df6_1373_3 <- df6_1373_3 %>% distinct()

setDT(df6_1373_3)
df6_1373_3[, Mean_total_tx_each_replicate_per_ug :=mean(Sum_total_tx_each_replicate_per_Ug), by = c('Gene_Tx', 'Pulldown_date','Ug_RNA_input_Capture','cDNA_primers')]

setDT(df6_1373_3)
df6_1373_3[, SE_total_tx_each_replicate_per_ug :=std.error(Sum_total_tx_each_replicate_per_Ug), by = c('Gene_Tx', 'Pulldown_date','Ug_RNA_input_Capture','cDNA_primers')]

library(dplyr)
#df8_1373_3 <- df6_1373_3 %>% filter(!is.na(SE_total_tx_each_replicate_per_ug))

#Select all columns except "replicates" and "counts"
df6_1373_3
agg_tbl_1373_3 <-df6_1373_3 %>%
  select(c(Gene_Tx, Ug_RNA_input_Capture, Pulldown_date, cDNA_primers, Gene_Tx, Mean_total_tx_each_replicate_per_ug, SE_total_tx_each_replicate_per_ug)) 
df8_1373_3 <- agg_tbl_1373_3 %>% as.data.frame()

#Get rid of duplicates again
df8_1373_3 <- df8_1373_3 %>% distinct()

#Keep only dates with BOTH R6 and dTs tested
foo <- split(df8_1373_3, f=df8_1373_3$Pulldown_date)
lapply(foo, FUN=function(x){unique(x$cDNA_primers)})

lapply(foo, FUN=function(x){length(unique(x$cDNA_primers))})

lapply(foo, FUN=function(x){if(length(unique(x$cDNA_primers))==2){x}else{NA} })

foo2 <- lapply(foo, FUN=function(x){if(length(unique(x$cDNA_primers))==2){x}else{NA} })

foo3 <- foo2[!is.na(foo2)]

df8_1373_3 <- do.call(rbind, foo3)

#Take just genes I want
df8_1373_3 <- df8_1373_3 %>% filter(grepl("5LTR_MSE_STAT5B_4223", Gene_Tx) | grepl("FANCA", Gene_Tx) | grepl("EIF", Gene_Tx) | grepl("STAT5B_4227", Gene_Tx) | grepl("CDK13", Gene_Tx) | grepl("NDC80", Gene_Tx))


#other way of selecting, filtering by frequency it shows up between experiment dates
#gene_frequency <- table(df8_1373_3$Gene_Tx)
#duplicate_genes <- names(gene_frequency)[gene_frequency > 1]

#Graph it
ggplot(df8_1373_3, aes(x=as.character(Ug_RNA_input_Capture), y=Mean_total_tx_each_replicate_per_ug, fill = cDNA_primers)) +
  geom_bar (stat="identity", position=position_dodge2(preserve="single")) +
  ylab("Counts_per_Ug") +
  xlab("Ug RNA input (capture)")+
  theme(axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), axis.line.x = element_line(color="black"),            axis.line.y = element_line(color="black"), panel.grid.major = element_line(size=0.5, linetype='solid', color="grey"),    panel.grid.minor = element_line(size=0.25, linetype='solid', color="grey")) +
  
  geom_errorbar(aes(ymin=Mean_total_tx_each_replicate_per_ug-SE_total_tx_each_replicate_per_ug, ymax=Mean_total_tx_each_replicate_per_ug+SE_total_tx_each_replicate_per_ug), position=position_dodge2(preserve="single")) +
  
  #geom_text(aes(label = Gene_Tx), vjust = 0.5, hjust = 4, angle = 90, position = position_dodge2(width =0.9, preserve = "single"), size = 2.5) +
  
  scale_x_discrete(guide=guide_axis(angle=90), labels = function(x) str_wrap(x, width=20))+
  
  scale_y_continuous()+
  facet_grid(~ Gene_Tx) +
  theme(panel.background = element_rect(size = 5, fill = "white"), axis.ticks.length = unit(0.5, "cm")) +
  guides(shape = guide_legend(override.aes = list(size = 10)), color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10),
        legend.key.size = unit(5, "mm"))


# oligo(dt)s vs r6 1709 ---------------------------------------------------


#First get counts
df5_1709_3
agg_tbl_1709_3 <-df5_1709_3 %>%
  group_by(Pulldown_date, Ug_RNA_input_Capture, cDNA_primers, Gene_Tx, Replicate) %>%
  mutate(Sum_total_tx_each_replicate =n())
df6_1709_3 <- agg_tbl_1709_3 %>% as.data.frame()

df6_1709_3 <- df6_1709_3 %>% distinct()

setDT(df6_1709_3)
df6_1709_3[, Mean_total_tx_each_replicate_per_ug :=mean(Sum_total_tx_each_replicate_per_Ug), by = c('Gene_Tx', 'cDNA_primers','Ug_RNA_input_Capture')]

setDT(df6_1709_3)
df6_1709_3[, SE_total_tx_each_replicate_per_ug :=std.error(Sum_total_tx_each_replicate_per_Ug), by = c('Gene_Tx', 'cDNA_primers','Ug_RNA_input_Capture')]


library(dplyr)
#df8_1709_3 <- df6_1709_3 %>% filter(!is.na(SD_total_tx_each_replicate_per_ug))

#Select all columns except "replicates" and "counts"
df6_1709_3
agg_tbl_1709_3 <-df6_1709_3 %>%
  select(c(Gene_Tx, Pulldown_date, cDNA_primers, Ug_RNA_input_Capture, Gene_Tx, Mean_total_tx_each_replicate_per_ug, SE_total_tx_each_replicate_per_ug)) 
df8_1709_3 <- agg_tbl_1709_3 %>% as.data.frame()

#Get rid of duplicates again
df8_1709_3 <- df8_1709_3 %>% distinct()

#Keep only dates with BOTH R6 and dTs tested
foo <- split(df8_1709_3, f=df8_1709_3$Pulldown_date)
lapply(foo, FUN=function(x){unique(x$cDNA_primers)})

lapply(foo, FUN=function(x){length(unique(x$cDNA_primers))})

lapply(foo, FUN=function(x){if(length(unique(x$cDNA_primers))==2){x}else{NA} })

foo2 <- lapply(foo, FUN=function(x){if(length(unique(x$cDNA_primers))==2){x}else{NA} })

foo3 <- foo2[!is.na(foo2)]

df8_1709_3 <- do.call(rbind, foo3)

#Take just genes I want
df9_1709_3 <- df8_1709_3 %>% filter(grepl("HIV", Gene_Tx) | grepl("5LTR_MSD_STAT5B", Gene_Tx) | grepl("MAP4K1", Gene_Tx),)
df9_1709_3 <- df9_1709_3[!grepl('38611132', df9_1709_3$Gene_Tx),]
df9_1709_3 <- df9_1709_3[!grepl('42257672', df9_1709_3$Gene_Tx),]
df9_1709_3 <- df9_1709_3[!grepl('HIVE', df9_1709_3$Gene_Tx),]
df9_1709_3 <- df9_1709_3[!grepl('38610025', df9_1709_3$Gene_Tx),]

#Get rid of duplicates again
df9_1709_3 <- df9_1709_3 %>% distinct()

#df8_1709_3 <- df8_1709_3[!grepl('1-Feb-23', df8_1709_3$Pulldown_date),]

#other way of selecting, filtering by frequency it shows up between experiment dates
#gene_frequency <- table(df8_1709_3$Gene_Tx)
#duplicate_genes <- names(gene_frequency)[gene_frequency > 1]

#Graph it
ggplot(df9_1709_3, aes(x=as.character(Ug_RNA_input_Capture), y=Mean_total_tx_each_replicate_per_ug, fill = cDNA_primers)) +
  geom_bar (stat="identity", position=position_dodge2(preserve="single")) +
  ylab("Counts_per_Ug") +
  xlab("Ug RNA input (capture)")+
  theme(axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), axis.line.x = element_line(color="black"),            axis.line.y = element_line(color="black"), panel.grid.major = element_line(size=0.5, linetype='solid', color="grey"),    panel.grid.minor = element_line(size=0.25, linetype='solid', color="grey")) +
  
  geom_errorbar(aes(ymin=Mean_total_tx_each_replicate_per_ug-SE_total_tx_each_replicate_per_ug, ymax=Mean_total_tx_each_replicate_per_ug+SE_total_tx_each_replicate_per_ug), position=position_dodge2(preserve="single")) +
  
  #geom_text(aes(label = Gene_Tx), vjust = 0.5, hjust = 4, angle = 90, position = position_dodge2(width =0.9, preserve = "single"), size = 2.5) +
  
  scale_x_discrete(guide=guide_axis(angle=90), labels = function(x) str_wrap(x, width=20))+
  
  scale_y_continuous()+
  facet_grid(~ Gene_Tx) +
  theme(panel.background = element_rect(size = 5, fill = "white"), axis.ticks.length = unit(0.5, "cm")) +
  guides(shape = guide_legend(override.aes = list(size = 10)), color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10),
        legend.key.size = unit(5, "mm"))

# oligo(dt)s vs r6 1858 ---------------------------------------------------

library(dplyr)
df3_1858_3
agg_tbl_1858_3 <-df3_1858_3 %>%
  group_by(Pulldown_date, Ug_RNA_input_Capture, Gene_Tx, cDNA_primers, Replicate) %>%
  mutate(Sum_total_tx_each_replicate =n())
df4_1858_3 <- agg_tbl_1858_3 %>% as.data.frame()

#Normalize by counts/ug
library(dplyr)
df4_1858_3
agg_tbl_1858_3 <-df4_1858_3 %>%
  group_by(Sum_total_tx_each_replicate, Ug_RNA_input_Capture) %>%
  mutate(Sum_total_tx_each_replicate_per_Ug = Sum_total_tx_each_replicate/Ug_RNA_input_Capture)
df5_1858_3 <- agg_tbl_1858_3 %>% as.data.frame()

#write.csv(df6_1858_3, "df6_1858_3.csv")

#Select only columns I need
library(dplyr)
library(tidyverse)

df5_1858_3
agg_tbl_1858_3 <-df5_1858_3 %>%
  select(c(Gene_Tx, Pulldown_date, cDNA_primers, Ug_RNA_input_Capture, Sum_total_tx_each_replicate_per_Ug, Gene_Tx, Replicate, Count_Total)) 
df6_1858_3 <- agg_tbl_1858_3 %>% as.data.frame()

df6_1858_3 <- df6_1858_3 %>% distinct()

setDT(df6_1858_3)
df6_1858_3[, Mean_total_tx_each_replicate_per_ug :=mean(Sum_total_tx_each_replicate_per_Ug), by = c('Gene_Tx', 'cDNA_primers','Ug_RNA_input_Capture')]

setDT(df6_1858_3)
df6_1858_3[, SE_total_tx_each_replicate_per_ug :=std.error(Sum_total_tx_each_replicate_per_Ug), by = c('Gene_Tx', 'cDNA_primers','Ug_RNA_input_Capture')]

library(dplyr)
#df8_1858_3 <- df6_1858_3 %>% filter(!is.na(SD_total_tx_each_replicate_per_ug))

#Select all columns except "replicates" and "counts"
df6_1858_3
agg_tbl_1858_3 <-df6_1858_3 %>%
  select(c(Gene_Tx, cDNA_primers, Ug_RNA_input_Capture, Gene_Tx, Mean_total_tx_each_replicate_per_ug, SE_total_tx_each_replicate_per_ug)) 
df8_1858_3 <- agg_tbl_1858_3 %>% as.data.frame()

#Get rid of duplicates again
df8_1858_3 <- df8_1858_3 %>% distinct()

#Take just genes I want
df8_1858_3 <- df8_1858_3 %>% filter(grepl("HIV", Gene_Tx) | grepl("5LTR_MSD_ST", Gene_Tx) | grepl("FLT", Gene_Tx)| grepl("VPR", Gene_Tx))


#other way of selecting, filtering by frequency it shows up between experiment dates
#gene_frequency <- table(df8_1858_3$Gene_Tx)
#duplicate_genes <- names(gene_frequency)[gene_frequency > 1]

#Graph it
ggplot(df8_1858_3, aes(x=as.character(Ug_RNA_input_Capture), y=Mean_total_tx_each_replicate_per_ug, fill = cDNA_primers)) +
  geom_bar (stat="identity", position=position_dodge2(preserve="single")) +
  ylab("Counts_per_Ug") +
  xlab("Ug RNA input (capture)")+
  theme(axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), axis.line.x = element_line(color="black"),            axis.line.y = element_line(color="black"), panel.grid.major = element_line(size=0.5, linetype='solid', color="grey"),    panel.grid.minor = element_line(size=0.25, linetype='solid', color="grey")) +
  geom_errorbar(aes(ymin=Mean_total_tx_each_replicate_per_ug-SE_total_tx_each_replicate_per_ug, ymax=Mean_total_tx_each_replicate_per_ug+SE_total_tx_each_replicate_per_ug), position=position_dodge2(preserve="single")) +
  
  #geom_text(aes(label = Gene_Tx), vjust = 0.5, hjust = 4, angle = 90, position = position_dodge2(width =0.9, preserve = "single"), size = 2.5) +
  
  scale_x_discrete(guide=guide_axis(angle=90), labels = function(x) str_wrap(x, width=20))+
  
  scale_y_continuous()+
  facet_grid(~ Gene_Tx) +
  theme(panel.background = element_rect(size = 5, fill = "white"), axis.ticks.length = unit(0.5, "cm")) +
  guides(shape = guide_legend(override.aes = list(size = 10)), color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10),
        legend.key.size = unit(5, "mm"))

-----------------------------------------------------
  

#Select only columns needed

#agg_tbl_megasheet <-df3_megasheet %>%
 # select(c(PIC_ID, Gene_Tx, Gene_Type, Ug_RNA_input_Capture, Replicate, Count_Total, Pulldown_date))
#df3_megasheet <- agg_tbl_megasheet %>% as.data.frame()

#df3_megasheet <- df3_megasheet %>% distinct()


#agg_tbl_megasheet <-df3_megasheet %>%
  #group_by(PIC_ID, Gene_Type) %>% 
  #mutate(Sum_Gene_Type_Counts =sum(Count_Total))
#df3_megasheet <- agg_tbl_megasheet %>% as.data.frame()

#agg_tbl_megasheet <-df3_megasheet %>%
 # group_by(PIC_ID, Gene_Type, Pulldown_date) %>%
#mutate(Sum_Ug_RNA_input_Capture =sum(Ug_RNA_input_Capture))
#df3_megasheet <- agg_tbl_megasheet %>% as.data.frame()

#agg_tbl_megasheet <-df3_megasheet %>%
#  select(c(PIC_ID, Gene_Tx, Gene_Type, Ug_RNA_input_Capture, Replicate, Count_Total))
#df3_megasheet <- agg_tbl_megasheet %>% as.data.frame()

#df3_megasheet <- df3_megasheet %>% distinct()

#agg_tbl_megasheet <-df3_megasheet %>%
 # group_by(PIC_ID, Gene_Type) %>% 
  #mutate(Sum_Gene_Type_Counts =sum(Count_Total))
#df3_megasheet <- agg_tbl_megasheet %>% as.data.frame()

#df3_megasheet <- df3_megasheet %>% distinct()

#Ok, right here I cannot perform simple R functions on my sheet and I cannot figure out why. 

#df_megasheet_all_observations_summarized3 <- data.frame("PIC_ID"=unlist("PIC_ID"), "Gene_Tx"=unlist("Gene_Tx)"), "Gene_Type"=unlist("Gene_Type"), "Ug_RNA_input_Capture"=unlist("Ug_RNA_input_Capture"), "Replicate"=unlist("Replicate"), "Sum_total_tx_each_replicate"=unlist("Sum_total_tx_each_replicate"), "Pulldown_date"=unlist("Pulldown_date"))


#df_megasheet_all_observations_summarized3$Sum_total_tx_each_replicate <- as.numeric(df_megasheet_all_observations_summarized3$Sum_total_tx_each_replicate)
#df_megasheet_all_observations_summarized3$PIC_ID <- as.factor(df_megasheet_all_observations_summarized3$PIC_ID)
#df_megasheet_all_observations_summarized3$Gene_Type <- as.factor(df_megasheet_all_observations_summarized3$Gene_Type)
#df_megasheet_all_observations_summarized3$Pulldown_date <- as.factor(df_megasheet_all_observations_summarized3$Pulldown_date)

#df_megasheet_all_observations_summarized3$Pulldown_date <- as.factor(df_megasheet_all_observations_summarized3$Pulldown_date)
