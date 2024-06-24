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
pluri_all <- c("NANOG", "SOX2", "OCT3/4", "FGF4", "KLF4", "LIN28A", "ZFP42")
pluri_sig <- c("NANOG", "OCT3/4", "LIN28A", "ZFP42")
pluri_NS <- c("SOX2", "FGF4", "KLF4")
imprint_ESC <- c("MKRN3", "MAGEL2", "NDN", "SNRPN U4/ex2", "SNRPN ex1/2", "PWAR5", "116HGG1", "116HGG2", "116HGG3", "IPW")
biallelic_ESC <- c("TUBGCP5", "CYFIP1", "UBE3A", "GABRB3", "GABRA5", "HERC2", "CHRNA7")
imprint_neu <- c("MKRN3", "MAGEL2", "NDN", "SNRPN U4/ex2", "SNRPN ex1/2", "PWAR5", "116HGG1", "116HGG2", "116HGG3", "IPW", "SNORD115", "UBE3A-ATS", "UBE3A")
biallelic_neu <- c("TUBGCP5", "CYFIP1", "GABRB3", "GABRA5", "GABRG3", "HERC2", "CHRNA7")
neu <- c("PAX6", "MAP2", "DLG4")

#Make lists of samples to parse results files with
ESCs <- c("WT_ESC", "PWS_ESC", "AS1_ESC", "AS2_ESC")
neurons <- c("WT_neu", "PWS_neu", "AS1_neu")
AS2 <- c("WT_neu", "AS2_neu")
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
df_pluriESCs <- allRes_CellType[(allRes_CellType$Gene) %in% pluri_all,]
df_pluriESCs <- df_pluriESCs[(df_pluriESCs$CellLine) %in% ESCs,]
df_imprintESCs <- allRes_CellType[(allRes_CellType$Gene) %in% imprint_ESC,]
df_imprintESCs <- df_imprintESCs[(df_imprintESCs$CellLine) %in% ESCs,]
df_biallelicESCs <- allRes_CellType[(allRes_CellType$Gene) %in% biallelic_ESC,]
df_biallelicESCs <- df_biallelicESCs[(df_biallelicESCs$CellLine) %in% ESCs,]
df_GABRG3 <- allRes_CellType[(allRes_CellType$Gene) %in% "GABRG3",]
df_GABRG3 <- df_GABRG3[(df_GABRG3$CellLine) %in% ESCs,]
df_neuNeu <- allRes_CellType[(allRes_CellType$Gene) %in% neu,]
df_neuNeu_main <- df_neuNeu[(df_neuNeu$CellLine) %in% neurons,]
df_neuAS2 <- df_neuNeu[(df_neuNeu$CellLine) %in% AS2,]
df_imprintNeu <- allRes_CellType[(allRes_CellType$Gene) %in% imprint_neu,]
df_imprintNeu_main <- df_imprintNeu[(df_imprintNeu$CellLine) %in% neurons,]
df_imprintAS2 <- df_imprintNeu[(df_imprintNeu$CellLine) %in% AS2,]
df_biallelicNeu <- allRes_CellType[(allRes_CellType$Gene) %in% biallelic_neu,]
df_biallelicNeu_main <- df_biallelicNeu[(df_biallelicNeu$CellLine) %in% neurons,]
df_biallelicNeuAS2 <- df_biallelicNeu[(df_biallelicNeu$CellLine) %in% AS2,]
df_pluriNeu2ESC <- allRes_toESC[(allRes_toESC$CellLine) %in% neu2ESC,]
df_pluriNeu2ESC_sig <- df_pluriNeu2ESC[(df_pluriNeu2ESC$Gene) %in% pluri_sig,]
df_pluriNeu2ESC_NS <- df_pluriNeu2ESC[(df_pluriNeu2ESC$Gene) %in% pluri_NS,]
df_neuNeu2ESC <- allRes_toESC[(allRes_toESC$CellLine) %in% neu2ESC,]
df_PAX6_Neu2ESC <- df_neuNeu2ESC[(df_neuNeu2ESC$Gene) %in% "PAX6",]
df_MAP2_Neu2ESC <- df_neuNeu2ESC[(df_neuNeu2ESC$Gene) %in% "MAP2",]
df_DLG4_Neu2ESC <- df_neuNeu2ESC[(df_neuNeu2ESC$Gene) %in% "DLG4",]

#Make plots
pluri_ESCs <- ggplot(df_pluriESCs, aes(x=factor(Gene, level = pluri_all), y=RQ_mean, fill=Order)) + 
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
GABRG3_ESCs <- ggplot(df_GABRG3, aes(x=Gene, y=RQ_mean, fill=Order)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_point(aes(y=RQ1), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ2), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ3), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_errorbar(aes(ymin=RQmin, ymax=RQmax), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values = c(A = "#484848", B="#002AD0", C="#D41159", D="#FB3C9B")) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=10, face="italic")) +
  labs(title="GABRG3 in hESCs", x ="", y = "Relative Gene Expression")
GABRG3_ESCs
neuronal_neus <- ggplot(df_neuNeu_main, aes(x=factor(Gene, level = neu), y=RQ_mean, fill=Order)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_point(aes(y=RQ1), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ2), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ3), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_errorbar(aes(ymin=RQmin, ymax=RQmax), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values = c(A = "#484848", B="#002AD0", C="#D41159")) + theme_minimal() + 
  theme(axis.text.x = element_text(size=10, face="italic")) +
  labs(title="Neuronal Genes in 10-week Neurons", x ="", y = "Relative Gene Expression")
neuronal_neus
imprint_neus <- ggplot(df_imprintNeu_main, aes(x=factor(Gene, level = imprint_neu), y=RQ_mean, fill=Order)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_point(aes(y=RQ1), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ2), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ3), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_errorbar(aes(ymin=RQmin, ymax=RQmax), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values = c(A = "#484848", B="#002AD0", C="#D41159")) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=10, face="italic")) +
  labs(title="Imprinted 15q11-q13 Genes in 10-week Neurons", x ="", y = "Relative Gene Expression")
imprint_neus
biallelic_neus <- ggplot(df_biallelicNeu_main, aes(x=factor(Gene, level = biallelic_neu), y=RQ_mean, fill=Order)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_point(aes(y=RQ1), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ2), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ3), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_errorbar(aes(ymin=RQmin, ymax=RQmax), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values = c(A = "#484848", B="#002AD0", C="#D41159")) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=10, face="italic")) +
  labs(title="Biallelic 15q11-q13 Genes in 10-week Neurons", x ="", y = "Relative Gene Expression")
biallelic_neus
neuronal_AS2 <- ggplot(df_neuAS2, aes(x=factor(Gene, level = neu), y=RQ_mean, fill=Order)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_point(aes(y=RQ1), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ2), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ3), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_errorbar(aes(ymin=RQmin, ymax=RQmax), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values = c(A = "#484848", D="#FB3C9B")) + theme_minimal() + 
  theme(axis.text.x = element_text(size=10, face="italic")) +
  labs(title="Neuronal Genes in 10-week Neurons", x ="", y = "Relative Gene Expression")
neuronal_AS2
imprint_AS2 <- ggplot(df_imprintAS2, aes(x=factor(Gene, level = imprint_neu), y=RQ_mean, fill=Order)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_point(aes(y=RQ1), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ2), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ3), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_errorbar(aes(ymin=RQmin, ymax=RQmax), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values = c(A = "#484848", D="#FB3C9B")) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=10, face="italic")) +
  labs(title="Imprinted 15q11-q13 Genes in 10-week Neurons", x ="", y = "Relative Gene Expression")
imprint_AS2
biallelic_AS2 <- ggplot(df_biallelicNeuAS2, aes(x=factor(Gene, level = biallelic_neu), y=RQ_mean, fill=Order)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_point(aes(y=RQ1), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ2), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_point(aes(y=RQ3), position = position_jitterdodge(jitter.width = 0.07)) +
  geom_errorbar(aes(ymin=RQmin, ymax=RQmax), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values = c(A = "#484848", B="#002AD0", C="#D41159", D="#FB3C9B")) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=10, face="italic")) +
  labs(title="Biallelic 15q11-q13 Genes in 10-week Neurons", x ="", y = "Relative Gene Expression")
biallelic_AS2
pluri_Neu2ESC_sig <- ggplot(df_pluriNeu2ESC_sig, aes(x=factor(Gene, level = pluri_sig), y=RQ_mean, fill=Order)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_point(aes(y=RQ1), position = position_jitterdodge(jitter.width = 0.05)) +
  geom_point(aes(y=RQ2), position = position_jitterdodge(jitter.width = 0.05)) +
  geom_point(aes(y=RQ3), position = position_jitterdodge(jitter.width = 0.05)) +
  geom_errorbar(aes(ymin=RQmin, ymax=RQmax), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values = c(A = "#484848", B="#9A9A9A", C="#002AD0", D= "#D41159", E="#FB3C9B")) + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 60, hjust=1, size=10, face="italic")) + ylim(0, 1.5) +
  labs(title="Pluripotency Genes in 10-week Neurons vs WT ESCs", x ="", y = "Relative Gene Expression")
pluri_Neu2ESC_sig
pluri_Neu2ESC_NS <- ggplot(df_pluriNeu2ESC_NS, aes(x=factor(Gene, level = pluri_NS), y=RQ_mean, fill=Order)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_point(aes(y=RQ1), position = position_jitterdodge(jitter.width = 0.05)) +
  geom_point(aes(y=RQ2), position = position_jitterdodge(jitter.width = 0.05)) +
  geom_point(aes(y=RQ3), position = position_jitterdodge(jitter.width = 0.05)) +
  geom_errorbar(aes(ymin=RQmin, ymax=RQmax), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values = c(A = "#484848", B="#9A9A9A", C="#002AD0", D= "#D41159", E="#FB3C9B")) + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 60, hjust=1, size=10, face="italic")) +
  labs(title="Pluripotency Genes in 10-week Neurons vs WT ESCs", x ="", y = "Relative Gene Expression")
pluri_Neu2ESC_NS
PAX6_Neu2ESC <- ggplot(df_PAX6_Neu2ESC, aes(x=Gene, y=RQ_mean, fill=Order)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_point(aes(y=RQ1), position = position_jitterdodge(jitter.width = 0.05)) +
  geom_point(aes(y=RQ2), position = position_jitterdodge(jitter.width = 0.05)) +
  geom_point(aes(y=RQ3), position = position_jitterdodge(jitter.width = 0.05)) +
  geom_errorbar(aes(ymin=RQmin, ymax=RQmax), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values = c(A = "#484848", B="#9A9A9A", C="#002AD0", D= "#D41159", E="#FB3C9B")) + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 60, hjust=1, size=10, face="italic")) +
  labs(title="Neuronal Genes in 10-week Neurons vs WT ESCs", x ="", y = "Relative Gene Expression")
PAX6_Neu2ESC
MAP2_Neu2ESC <- ggplot(df_MAP2_Neu2ESC, aes(x=Gene, y=RQ_mean, fill=Order)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_point(aes(y=RQ1), position = position_jitterdodge(jitter.width = 0.05)) +
  geom_point(aes(y=RQ2), position = position_jitterdodge(jitter.width = 0.05)) +
  geom_point(aes(y=RQ3), position = position_jitterdodge(jitter.width = 0.05)) +
  geom_errorbar(aes(ymin=RQmin, ymax=RQmax), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values = c(A = "#484848", B="#9A9A9A", C="#002AD0", D= "#D41159", E="#FB3C9B")) + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 60, hjust=1, size=10, face="italic")) +
  labs(title="Neuronal Genes in 10-week Neurons vs WT ESCs", x ="", y = "Relative Gene Expression")
MAP2_Neu2ESC
DLG4_Neu2ESC <- ggplot(df_DLG4_Neu2ESC, aes(x=Gene, y=RQ_mean, fill=Order)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_point(aes(y=RQ1), position = position_jitterdodge(jitter.width = 0.05)) +
  geom_point(aes(y=RQ2), position = position_jitterdodge(jitter.width = 0.05)) +
  geom_point(aes(y=RQ3), position = position_jitterdodge(jitter.width = 0.05)) +
  geom_errorbar(aes(ymin=RQmin, ymax=RQmax), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values = c(A = "#484848", B="#9A9A9A", C="#002AD0", D= "#D41159", E="#FB3C9B")) + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 60, hjust=1, size=10, face="italic")) +
  labs(title="Neuronal Genes in 10-week Neurons vs WT ESCs", x ="", y = "Relative Gene Expression")
DLG4_Neu2ESC

#Save plots
pdf("MultiGene_Plots.pdf", width = 10, height = 4)
pluri_ESCs
imprint_ESCs
biallelic_ESCs
imprint_neus
biallelic_neus
imprint_AS2
biallelic_AS2
dev.off()

pdf("NeuGene_Plots.pdf", width = 5, height = 4)
neuronal_neus
neuronal_AS2
dev.off()

pdf("Neu2ESC_Pluri_Plots.pdf", width = 5, height = 4)
pluri_Neu2ESC_sig
pluri_Neu2ESC_NS
dev.off()

pdf("IndivGene_Plots.pdf", width = 4, height = 4)
GABRG3_ESCs
PAX6_Neu2ESC
MAP2_Neu2ESC
DLG4_Neu2ESC
dev.off()
