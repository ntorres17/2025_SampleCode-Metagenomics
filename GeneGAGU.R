
# Genes only

library("readxl")
library("tidyverse")
library("dplyr")
library("writexl")
library("grid")
library("futile.logger")
library("VennDiagram")
library("ggVennDiagram")
library("ggplot2")

setwd("/Users/nancytorres/ResearchProject/Script/GAGU_DS")
# sheet names are Gabon DIAmond and Guinea Diamond
# Specify the sheet by name, since its named 

# Set up
sheet_names <- excel_sheets("MAGs-DIAMOND.xlsx")
print(sheet_names)

GA <-read_excel("MAGs-DIAMOND.xlsx", sheet = "Gabon DIAmond")
head(GA)

GU<-read_excel("MAGs-DIAMOND.xlsx", sheet = "Guinea Diamond")
head(GU)

# Run through Venn Diagram - letting it self compute
# Results shared with intersect and setdiff
# Genes 2284 shared across both

# Approach 1 - short and sweet
# Create List of Genes 
Genes <-list(
  Gabon = GA$gene...13,
  Guinea = GU$`gene name`
)

Venn_genes <- ggVennDiagram(Genes, label = "both") + 
  ggtitle("Venn Diagram Shared Genes: Gabon & Guinea") + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold") 
  )

ggsave("Venn_genes.png", plot = Venn_genes, width = 10, height = 8)

# GENES SECTION - Approach 2, double checking work
# Comparing genes GA column M vs. GU column L and obtaining common elements shared between.
# 2284 intersect 
# Gabon, GENE column M aka column 13

GA_genes <- GA$gene...13
# gives total number of rows
print(length(GA_genes)) #8911

# Guinea, gene column L - aka column 14
GU_genes <- GU$`gene name`
# gives total number of rows
print(length(GU_genes)) #18765

Intersect_genes <- intersect(GA_genes, GU_genes)

if (length(Intersect_genes) > 0) {
  print(paste("Total number of common GENES between GA col. M and GU col. L:", length(Intersect_genes)))
} else {
  print("No genes shared between GA col. M v. GU col. L")
}

typeof(Intersect_genes) # Character
tail(Intersect_genes) # only contains names 
length(Intersect_genes) # 2284

# Set Diff
# Comparing genes GA column M vs. GU column L and obtaining different elements to each each.
# 3971 different
Diff_genes <- setdiff(GA_genes, GU_genes)

if (length(Diff_genes) > 0) {
  print(paste("Total number of different GENES between GA col. M and GU col. L", length(Diff_genes)))
} else {
  print("No different GENES between GA col. M and v. GU col. L")
}

# Both approach 1 and 2 are consistent in their results
# Summary of genes
genes_summary <- data.frame(
  Data_Overview = c("Total Genes in GA", "Total Genes in GU", "Shared Genes", "Different Genes"),
  Count = c(length(GA_genes), length(GU_genes), length(Intersect_genes), length(Diff_genes))
)

# Further exploring, actual breakdown
# First GA, raw counts 
# Converting into data frame
# Gabon, GENES
# Now getting specifics and frequencies in each df.
GA_genes <- as.data.frame(table(GA_genes))
head(GA_genes)
colnames(GA_genes) <- c("GA_Genes" , "GA_Freq.")
print(sum(GA_genes$GA_Freq.))
GA_genes <- GA_genes[order(GA_genes$GA_Freq.,decreasing = TRUE),]
head(GA_genes)
GA_genes_percent <- (GA_genes$GA_Freq./ sum(GA_genes$GA_Freq.) * 100)
head(GA_genes_percent)
GA_genes <- cbind(GA_genes, GA_genes_percent)
head(GA_genes) # Now have genes, frequency and percent 
(sum(GA_genes_percent)) # 100%
(sum(GA_genes$GA_Freq.)) # 8911
typeof(GA_genes) # list

# Guinea, genes
GU_genes <- as.data.frame(table(GU_genes))
head(GU_genes)
colnames(GU_genes) <- c("GU_Genes" , "GU_Freq.")
print(sum(GU_genes$GU_Freq.))
head(GU_genes)
GU_genes <- GU_genes[order(GU_genes$GU_Freq.,decreasing = TRUE),]
head(GU_genes)
GU_genes_percent <- (GU_genes$GU_Freq./ sum(GU_genes$GU_Freq.) * 100)
head(GU_genes_percent)
GU_genes <- cbind(GU_genes, GU_genes_percent)
head(GU_genes) # Now have genes, frequency and percent for Guinea
(sum(GU_genes_percent)) # 100%
(sum(GU_genes$GU_Freq)) # 18765

# Filter shared genes from GA and GU 
GAGU_genes <- list(GA_genes, GU_genes)
head(n = 3, GAGU_genes) # Mass list with all 6 columns NOT filtered re: Shared
typeof(GAGU_genes) # list
sum(GA_genes$GA_genes_percent) # 100
sum(GU_genes$GU_genes_percent) # 100 

# Filter shared genes from GA and GU
GA_Shared_Genes <- GA_genes %>% filter(GA_Genes %in% Intersect_genes)
GU_Shared_Genes <- GU_genes %>% filter(GU_Genes %in% Intersect_genes)

print(sum(GA_Shared_Genes$Freq.)) #4426
print(sum(GA_Shared_Genes$GA_genes_percent)) # only adds up to 50%
head(n = 3, GA_Shared_Genes) # colnames Genes, GA_Freq., GA_genes_percent
head(n = 3, GU_Shared_Genes) # colnames Genes, GU_Freq., GU_genes_percent
nrow(GA_Shared_Genes) # 2248
nrow(GU_Shared_Genes) # 2248
print(sum(GU_Shared_Genes$GU_Freq.)) # 7982
print(sum(GU_Shared_Genes$GU_genes_percent)) # adds up to 43%

# Merge SHARED genes data between GA and GU; narrows down
Shared_Freq_Genes <- merge(GA_Shared_Genes, GU_Shared_Genes, by.x = "GA_Genes", by.y = "GU_Genes")
head(Shared_Freq_Genes)
Shared_Freq_Genes <- Shared_Freq_Genes %>% rename(Genes = GA_Genes) #r rename since same 2284
head(Shared_Freq_Genes)

# GA to GU in order; GA higher
Shared_Freq_Genes_GA <- Shared_Freq_Genes[order(Shared_Freq_Genes$GA_Freq., decreasing = TRUE),]
head(Shared_Freq_Genes_GA)

# GU to GA in order; GU higher 
Shared_Freq_Genes_GU <- Shared_Freq_Genes[order(Shared_Freq_Genes$GU_Freq., decreasing = TRUE),]
head(Shared_Freq_Genes_GU)

# Sum the frequencies for both GA and GU
sum(Shared_Freq_Genes_GA$GA_Freq.) # 4426
sum(Shared_Freq_Genes_GU$GU_Freq.) # 7982

# Now for pie charts, top 10 using % 
Top_10_Shared_GA_Genes <- head(n=10, Shared_Freq_Genes_GA)
print(Top_10_Shared_GA_Genes)

Shared_GA_Genes <- Top_10_Shared_GA_Genes$GA_genes_percent
labels_Shared_GA_Genes <- paste0(Top_10_Shared_GA_Genes$Genes, " (", round(Top_10_Shared_GA_Genes$GA_genes_percent, 2), "%)")
n <- length(Top_10_Shared_GA_Genes$Genes)
col_pie_Genes <- scale_color_poke(pokemon = "Croconaw")$palette(n)

png("Top 10 Shared Genes in Gabon.png", width = 900, height = 600)

pie(
  Shared_GA_Genes, 
  label = labels_Shared_GA_Genes, 
  col = col_pie_Genes,
  main = "Top 10 Shared Genes in Gabon"
)

dev.off()

# Top 10 respective to Guinea using %
Top_10_Shared_GU_Genes<- head(n=10, Shared_Freq_Genes_GU)
print(Top_10_Shared_GU_Genes)

Shared_GU_Genes <- Top_10_Shared_GU_Genes$GU_genes_percent
labels_Shared_GU_Genes <- paste0(Top_10_Shared_GU_Genes$Genes, " (", round(Top_10_Shared_GU_Genes$GU_genes_percent, 2), "%)")
n <- length(Top_10_Shared_GU_Genes$Genes)
col_pie_genes <- scale_color_poke(pokemon = "Weedle")$palette(n)

png("Top 10 Shared Genes in Guinea.png", width = 900, height = 600)

pie(
  Shared_GU_Genes, 
  label = labels_Shared_GU_Genes, 
  col = col_pie_genes,
  main = "Top 10 Shared Genes in Guinea"
)

dev.off()

write_xlsx(list(
  "Genes_Summary_SM" = genes_summary,
  "Shared_Genes_GA" = data.frame(Shared_Freq_Genes_GA),
  "Shared_Genes_GU" = data.frame(Shared_Freq_Genes_GU),
  "GA_Freq_Genes" = data.frame(GA_genes),
  "GU_Freq_Genes" = data.frame(GU_genes),
  "Different_Genes" = data.frame(Diff_genes)
), "Genes_Summary.xlsx")

