
# Binning results metabinner Gabon and Guinea, separately 

library("readxl")
library("tidyverse")
library("dplyr")
library("writexl")
library("grid")
library("futile.logger")
library("VennDiagram")
library("ggVennDiagram")
library("ggplot2")

# Ok now reading in the tsv file
GA_mb <- read_tsv("metabinner_result_GA.tsv")
GU_mb <- read_tsv("metabinner_result_GU.tsv")
colnames(GA_mb) <- c("GA_Contigs", "GA_Bin")
colnames(GU_mb) <- c("GU_Contigs", "GU_Bin")
head(GA_mb)
head(GU_mb)
nrow(GA_mb) # 4823
nrow(GU_mb) # 1413

Table_GA_bins <- table(GA_mb) # Make it into a table 
Table_GU_bins <- table(GU_mb) 
head(Table_GA_bins)
head(Table_GU_bins)

# GA 
Counted_GA_bins <- margin.table(Table_GA_bins, 2) 
print(Counted_GA_bins) # 4823 total

GA_Bins_1 <- GA_mb[GA_mb$GA_Bin == "1", c("GA_Contigs", "GA_Bin")]
print(GA_Bins_1)
GA_Bins_2 <- GA_mb[GA_mb$GA_Bin == "2", c("GA_Contigs", "GA_Bin")]
print(GA_Bins_2)
GA_Bins_3 <- GA_mb[GA_mb$GA_Bin == "3", c("GA_Contigs", "GA_Bin")]
print(GA_Bins_3)
GA_Bins_4 <- GA_mb[GA_mb$GA_Bin == "4", c("GA_Contigs", "GA_Bin")]
print(GA_Bins_4)

# GU
Counted_GU_bins <- margin.table(Table_GU_bins, 2) 
print(Counted_GU_bins) # 4823 total

GU_Bins_1 <- GU_mb[GU_mb$GU_Bin == "1", c("GU_Contigs", "GU_Bin")]
print(GU_Bins_1)
GU_Bins_2 <- GU_mb[GU_mb$GU_Bin == "2", c("GU_Contigs", "GU_Bin")]
print(GU_Bins_2)
GU_Bins_3 <- GU_mb[GU_mb$GU_Bin == "3", c("GU_Contigs", "GU_Bin")]
print(GU_Bins_3)

if (any(GA_mb == " ")) {
  blank_count_GA <- sum(GA_mb == " ")
  print(paste("Total number of contigs WITHOUT bins in Gabon", blank_count_GA))
} else {
  print("All contigs accounted for in Gabon dataset")
}

if (any(GU_mb == " ")) {
  blank_count_GU <- sum(GU_mb == " ")
  print(paste("Total number of contigs WITHOUT bins in Guinea", blank_count_GU))
} else {
  print("All contigs accounted for in Guinea dataset")
}

genes_summary <- data.frame(
  Data_Overview = c("Total Genes in GA", "Total Genes in GU", "Shared Genes", "Different Genes"),
  Count = c(length(GA_genes), length(GU_genes), length(Intersect_genes), length(Diff_genes))
)

Binning_sum <- data.frame(
  Data_Overview = c("Total Contigs in GA", "Gabon Bin 1", "Gabon Bin 2", "Gabon Bin 3", "Gabon Bin 4", 
                    "Total Contings in GU", "Guinea Bin 1", "Guinea Bin 2", "Guinea Bin 3"),
  Count = c(nrow(GA_mb), nrow(GA_Bins_1), nrow(GA_Bins_2), nrow(GA_Bins_3), nrow(GA_Bins_4), 
            nrow(GU_mb), nrow(GU_Bins_1), nrow(GU_Bins_2), nrow(GU_Bins_3))
)

write_xlsx(list(
  "Binning_Summary" = Binning_sum,
  "Counted_GA_Bins" = data.frame(Counted_GA_bins),
  "GA_Bin_1" = data.frame(GA_Bins_1),
  "GA_Bin_2" = data.frame(GA_Bins_2),
  "GA_Bin_3" = data.frame(GA_Bins_3),
  "GA_Bin_4" = data.frame(GA_Bins_4),
  "Counted_GU_Bins" = data.frame(Counted_GU_bins),
  "GU_Bin_1" = data.frame(GU_Bins_1),
  "GU_Bin_2" = data.frame(GU_Bins_2),
  "GU_Bin_3" = data.frame(GU_Bins_3)
), "Binning_Sum_GAGU.xlsx")