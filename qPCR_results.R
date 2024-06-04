##---Plotting results of qPCR analysis---##
library("ggplot2")
library("tidyverse")

#Set working directory
directory <- "../OneDrive - UConn Health Center/Megabase_deletion_lines_resource_paper/"
setwd(directory)

#Import results
allRes_CellType <- read.csv("qPCR_Results_byCellType.csv", header = TRUE)
allRes_toESC <- read.csv("qPCR_Results_toWTESCs.csv", header = TRUE)

#Make lists of genes to parse results files with
pluri <- c("NANOG", "SOX2", "OCT3/4", "FGF4", "LIN28A", "ZFP42")
imprint_ESC <- c("MKRN3", "MAGEL2", "NDN", "SNRPN U4/ex2", "SNRPN ex1/2", "PWAR5", "116HGG1", "116HGG2", "116HGG3", "IPW")
biallelic_ESC <- c("TUBGCP5", "CYFIP1", "UBE3A", "GABRB3", "GABRA5", "GABRG3", "HERC2", "CHRNA7")
imprint_neu <- c("MKRN3", "MAGEL2", "NDN", "SNRPN U4/ex2", "SNRPN ex1/2", "PWAR5", "116HGG1", "116HGG2", "116HGG3", "IPW", "SNORD115", "UBE3A-ATS", "UBE3A")
biallelic_neu <- c("TUBGCP5", "CYFIP1", "GABRB3", "GABRA5", "GABRG3", "HERC2", "CHRNA7")
neu <- c("PAX6", "MAP2", "DLG4")

#Make lists of samples to parse results files with
ESCs <- c("WT_ESC", "PWS_ESC", "AS1_ESC", "AS2_ESC")
neurons <- c("WT_neu", "PWS_neu", "AS1_neu", "AS2_neu")
neu2ESC <- c("WT_ESC", "WT_neu", "PWS_neu", "AS1_neu", "AS2_neu")

#Create column in results files for average RQ
allRes_CellType$RQ_mean <- rowMeans(allRes_CellType[3:5])
allRes_toESC$RQ_mean <- rowMeans(allRes_toESC[3:5])
#Create column in results files for ordering bars later on
allRes_CellType$Order <- "A"
allRes_CellType <- allRes_CellType %>% mutate(Order = ifelse(CellLine %in% c("PWS_ESC", "PWS_neu"), "B", Order))
allRes_CellType <- allRes_CellType %>% mutate(Order = ifelse(CellLine %in% c("AS1_ESC", "AS1_neu"), "C", Order))
allRes_CellType <- allRes_CellType %>% mutate(Order = ifelse(CellLine %in% c("AS2_ESC", "AS2_neu"), "D", Order))
allRes_toESC$Order <- "A"
allRes_toESC <- allRes_toESC %>% mutate(Order = ifelse(CellLine %in% c("WT_neu"), "B", Order))
allRes_toESC <- allRes_toESC %>% mutate(Order = ifelse(CellLine %in% c("PWS_neu"), "C", Order))
allRes_toESC <- allRes_toESC %>% mutate(Order = ifelse(CellLine %in% c("AS1_neu"), "D", Order))
allRes_toESC <- allRes_toESC %>% mutate(Order = ifelse(CellLine %in% c("AS2_neu"), "E", Order))

#Parse results tables
df_pluriESCs <- allRes_CellType[(allRes_CellType$Gene) %in% pluri,]
df_pluriESCs <- df_pluriESCs[(df_pluriESCs$CellLine) %in% ESCs,]
df_imprintESCs <- allRes_CellType[(allRes_CellType$Gene) %in% imprint_ESC,]
df_imprintESCs <- df_imprintESCs[(df_imprintESCs$CellLine) %in% ESCs,]
df_biallelicESCs <- allRes_CellType[(allRes_CellType$Gene) %in% biallelic_ESC,]
df_biallelicESCs <- df_biallelicESCs[(df_biallelicESCs$CellLine) %in% ESCs,]
df_neuNeu <- allRes_CellType[(allRes_CellType$Gene) %in% neu,]
df_neuNeu <- df_neuNeu[(df_neuNeu$CellLine) %in% neurons,]
df_imprintNeu <- allRes_CellType[(allRes_CellType$Gene) %in% imprint_neu,]
df_imprintNeu <- df_imprintNeu[(df_imprintNeu$CellLine) %in% neurons,]
df_biallelicNeu <- allRes_CellType[(allRes_CellType$Gene) %in% biallelic_neu,]
df_biallelicNeu <- df_biallelicNeu[(df_biallelicNeu$CellLine) %in% neurons,]
df_pluriNeu2ESC <- allRes_toESC[(allRes_toESC$Gene) %in% pluri,]
df_pluriNeu2ESC <- df_pluriNeu2ESC[(df_pluriNeu2ESC$CellLine) %in% neu2ESC,]
df_neuNeu2ESC <- allRes_toESC[(allRes_toESC$Gene) %in% neu,]
df_neuNeu2ESC <- df_neuNeu2ESC[(df_neuNeu2ESC$CellLine) %in% neu2ESC,]

#Make plots
pluri_ESCs <- ggplot(df_pluriESCs, aes(x=factor(Gene, level = pluri), y=RQ_mean, fill=Order)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_point(aes(y=RQ1), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ2), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ3), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_errorbar(aes(ymin=RQmin, ymax=RQmax), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values = c(A = "#484848", B="#002AD0", C="#D41159", D="#FB3C9B")) + theme_minimal() + 
  theme(axis.text.x = element_text(size=10, face="italic")) +
  labs(title="Pluipotency Genes in hESCs", x ="", y = "Relative Gene Expression")
pluri_ESCs
imprint_ESCs <- ggplot(df_imprintESCs, aes(x=factor(Gene, level = imprint_ESC), y=RQ_mean, fill=Order)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_point(aes(y=RQ1), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ2), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ3), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_errorbar(aes(ymin=RQmin, ymax=RQmax), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values = c(A = "#484848", B="#002AD0", C="#D41159", D="#FB3C9B")) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=10, face="italic")) +
  labs(title="Imprinted 15q11-q13 Genes in hESCs", x ="", y = "Relative Gene Expression")
imprint_ESCs
biallelic_ESCs <- ggplot(df_biallelicESCs, aes(x=factor(Gene, level = biallelic_ESC), y=RQ_mean, fill=Order)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_point(aes(y=RQ1), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ2), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ3), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_errorbar(aes(ymin=RQmin, ymax=RQmax), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values = c(A = "#484848", B="#002AD0", C="#D41159", D="#FB3C9B")) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=10, face="italic")) +
  labs(title="Biallelic 15q11-q13 Genes in hESCs", x ="", y = "Relative Gene Expression")
biallelic_ESCs
neuronal_neus <- ggplot(df_neuNeu, aes(x=factor(Gene, level = neu), y=RQ_mean, fill=Order)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_point(aes(y=RQ1), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ2), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ3), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_errorbar(aes(ymin=RQmin, ymax=RQmax), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values = c(A = "#484848", B="#002AD0", C="#D41159", D="#FB3C9B")) + theme_minimal() + 
  theme(axis.text.x = element_text(size=10, face="italic")) +
  labs(title="Neuronal Genes in 10-week Neurons", x ="", y = "Relative Gene Expression")
neuronal_neus
imprint_neus <- ggplot(df_imprintNeu, aes(x=factor(Gene, level = imprint_neu), y=RQ_mean, fill=Order)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_point(aes(y=RQ1), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ2), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ3), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_errorbar(aes(ymin=RQmin, ymax=RQmax), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values = c(A = "#484848", B="#002AD0", C="#D41159", D="#FB3C9B")) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=10, face="italic")) +
  labs(title="Imprinted 15q11-q13 Genes in 10-week Neurons", x ="", y = "Relative Gene Expression")
imprint_neus
biallelic_neus <- ggplot(df_biallelicNeu, aes(x=factor(Gene, level = biallelic_neu), y=RQ_mean, fill=Order)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_point(aes(y=RQ1), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ2), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ3), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_errorbar(aes(ymin=RQmin, ymax=RQmax), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values = c(A = "#484848", B="#002AD0", C="#D41159", D="#FB3C9B")) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=10, face="italic")) +
  labs(title="Biallelic 15q11-q13 Genes in 10-week Neurons", x ="", y = "Relative Gene Expression")
biallelic_neus
pluri_Neu2ESC <- ggplot(df_pluriNeu2ESC, aes(x=factor(Gene, level = pluri), y=RQ_mean, fill=Order)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_point(aes(y=RQ1), position = position_jitterdodge(jitter.width = 0.05)) +
  geom_point(aes(y=RQ2), position = position_jitterdodge(jitter.width = 0.05)) +
  geom_point(aes(y=RQ3), position = position_jitterdodge(jitter.width = 0.05)) +
  geom_errorbar(aes(ymin=RQmin, ymax=RQmax), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values = c(A = "#484848", B="#9A9A9A", C="#002AD0", D= "#D41159", E="#FB3C9B")) + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 60, hjust=1, size=10, face="italic")) +
  labs(title="Pluripotency Genes in 10-week Neurons vs WT ESCs", x ="", y = "Relative Gene Expression")
pluri_Neu2ESC
neu_Neu2ESC <- ggplot(df_neuNeu2ESC, aes(x=factor(Gene, level = neu), y=RQ_mean, fill=Order)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_point(aes(y=RQ1), position = position_jitterdodge(jitter.width = 0.05)) +
  geom_point(aes(y=RQ2), position = position_jitterdodge(jitter.width = 0.05)) +
  geom_point(aes(y=RQ3), position = position_jitterdodge(jitter.width = 0.05)) +
  geom_errorbar(aes(ymin=RQmin, ymax=RQmax), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values = c(A = "#484848", B="#9A9A9A", C="#002AD0", D= "#D41159", E="#FB3C9B")) + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 60, hjust=1, size=10, face="italic")) +
  labs(title="Neuronal Genes in 10-week Neurons vs WT ESCs", x ="", y = "Relative Gene Expression")
neu_Neu2ESC

#Save plots
pdf("Plots.pdf", width = 10, height = 4)
pluri_ESCs
imprint_ESCs
biallelic_ESCs
neuronal_neus
imprint_neus
biallelic_neus
pluri_Neu2ESC
neu_Neu2ESC
dev.off()
