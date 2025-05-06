
# Binning results metabinner SM, NBM and PBM separately 

library("readxl")
library("tidyverse")
library("dplyr")
library("writexl")


# Ok now reading in the tsv file
SM_mb <- read_tsv("metabinner_result_sm.tsv")
NBM_mb <- read_tsv("metabinner_result_nbm.tsv")
PBM_mb <- read_tsv("metabinner_result_pbm.tsv")
colnames(SM_mb) <- c("SM_Contigs", "SM_Bin")
colnames(NBM_mb) <- c("NBM_Contigs", "NBM_Bin")
colnames(PBM_mb) <- c("PBM_Contigs", "PBM_Bin")
head(SM_mb)
head(NBM_mb)
head(PBM_mb)
nrow(SM_mb) 
nrow(NBM_mb) 
nrow(PBM_mb) 

Table_SM_bins <- table(SM_mb) # Make it into a table 
Table_NBM_bins <- table(NBM_mb) 
Table_PBM_bins <- table(PBM_mb) 
head(Table_SM_bins)
head(Table_NBM_bins)
head(Table_PBM_bins)

# SM
Counted_SM_bins <- margin.table(Table_SM_bins, 2) 
print(Counted_SM_bins) # 1:308 - 2:732

SM_Bins_1 <- SM_mb[SM_mb$SM_Bin == "1", c("SM_Contigs", "SM_Bin")]
print(SM_Bins_1)
SM_Bins_2 <- SM_mb[SM_mb$SM_Bin == "2", c("SM_Contigs", "SM_Bin")]
print(SM_Bins_2)


# NBM
Counted_NBM_bins <- margin.table(Table_NBM_bins, 2) 
print(Counted_NBM_bins) # 1:271 - 2:1309

NBM_Bins_1 <- NBM_mb[NBM_mb$NBM_Bin == "1", c("NBM_Contigs", "NBM_Bin")]
print(NBM_Bins_1)
NBM_Bins_2 <- NBM_mb[NBM_mb$NBM_Bin == "2", c("NBM_Contigs", "NBM_Bin")]
print(NBM_Bins_2)

# PBM  1   2   3   4   5   6   7 
#     122 257 156  40 144 996 727 
Counted_PBM_bins <- margin.table(Table_PBM_bins, 2) 
print(Counted_PBM_bins) # 

PBM_Bins_1 <- PBM_mb[PBM_mb$PBM_Bin == "1", c("PBM_Contigs", "PBM_Bin")]
print(PBM_Bins_1)
PBM_Bins_2 <- PBM_mb[PBM_mb$PBM_Bin == "2", c("PBM_Contigs", "PBM_Bin")]
print(PBM_Bins_2)
PBM_Bins_3 <- PBM_mb[PBM_mb$PBM_Bin == "3", c("PBM_Contigs", "PBM_Bin")]
print(PBM_Bins_3)
PBM_Bins_4 <- PBM_mb[PBM_mb$PBM_Bin == "4", c("PBM_Contigs", "PBM_Bin")]
print(PBM_Bins_4)
PBM_Bins_5 <- PBM_mb[PBM_mb$PBM_Bin == "5", c("PBM_Contigs", "PBM_Bin")]
print(PBM_Bins_5)
PBM_Bins_6 <- PBM_mb[PBM_mb$PBM_Bin == "6", c("PBM_Contigs", "PBM_Bin")]
print(PBM_Bins_6)
PBM_Bins_7 <- PBM_mb[PBM_mb$PBM_Bin == "7", c("PBM_Contigs", "PBM_Bin")]
print(PBM_Bins_7)

if (any(SM_mb == " ")) {
  blank_count_SM <- sum(SM_mb == " ")
  print(paste("Total number of contigs WITHOUT bins in SM", blank_count_SM))
} else {
  print("All contigs accounted for in SM dataset")
}

if (any(NBM_mb == " ")) {
  blank_count_NBM <- sum(NBM_mb == " ")
  print(paste("Total number of contigs WITHOUT bins in NBM", blank_count_NBM))
} else {
  print("All contigs accounted for in NBM dataset")
}

if (any(PBM_mb == " ")) {
  blank_count_PBM <- sum(PBM_mb == " ")
  print(paste("Total number of contigs WITHOUT bins in PBM", blank_count_PBM))
} else {
  print("All contigs accounted for in PBM dataset")
}


Binning_sum_SM <- data.frame(
  Data_Overview = c("Total Contigs in SM", "SM Bin 1", "SM Bin 2",
                    "Total Contings in NBM", "NBM Bin 1", "NBM Bin 2", 
                    "Total Contings in PBM", "PBM Bin 1", "PBM Bin 2", 
                    "PBM Bin 3", "PBM Bin 4","PBM Bin 5", "PBM Bin 6", 
                    "PBM Bin 7"),
  Count = c(nrow(SM_mb), nrow(SM_Bins_1), nrow(SM_Bins_2), 
            nrow(NBM_mb), nrow(NBM_Bins_1), nrow(NBM_Bins_2),
            nrow(PBM_mb), nrow(PBM_Bins_1), nrow(PBM_Bins_2), 
            nrow(PBM_Bins_3), nrow(PBM_Bins_4), nrow(PBM_Bins_5), nrow(PBM_Bins_6), 
            nrow(PBM_Bins_7))
)

write_xlsx(list(
  "Binning_Summary_SM" = Binning_sum_SM,
  "Counted_SM_Bins" = data.frame(Counted_SM_bins),
  "SM_Bin_1" = data.frame(SM_Bins_1),
  "SM_Bin_2" = data.frame(SM_Bins_2),
  "Counted_NBM_Bins" = data.frame(Counted_NBM_bins),
  "NBM_Bin_1" = data.frame(NBM_Bins_1),
  "NBM_Bin_2" = data.frame(NBM_Bins_2),
  "Counted_PBM_Bins" = data.frame(Counted_PBM_bins),
  "PBM_Bin_1" = data.frame(PBM_Bins_1),
  "PBM_Bin_2" = data.frame(PBM_Bins_2),
  "PBM_Bin_3" = data.frame(PBM_Bins_3),
  "PBM_Bin_4" = data.frame(PBM_Bins_4),
  "PBM_Bin_5" = data.frame(PBM_Bins_5),
  "PBM_Bin_6" = data.frame(PBM_Bins_6),
  "PBM_Bin_7" = data.frame(PBM_Bins_7)),"Binning_Sum_SM.xlsx")