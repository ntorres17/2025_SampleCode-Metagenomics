
# Taxon only

library("readxl")
library("tidyverse")
library("dplyr")
library("writexl")
library("grid")
library("futile.logger")
library("VennDiagram")
library("ggVennDiagram")
library("ggplot2")
library("palettetown")

setwd("PATH")
# sheet names are Gabon DIAmond and Guinea Diamond
# Specify the sheet by name, since its named 
sheet_names <- excel_sheets("EXCEL_SHEET")
print(sheet_names)

GA <-read_excel("MAGs-DIAMOND.xlsx", sheet = "Gabon DIAmond")
head(GA)

GU <-read_excel("MAGs-DIAMOND.xlsx", sheet = "Guinea Diamond")
head(GU)

GA_taxon <- GA$Taxon...14
head(GA_taxon)
print(length(GA_taxon)) #8911

GU_taxon <- GU$Taxon...13
print(length(GU_taxon))
head(GU_taxon) # with ]

# Exclude ]
GU_taxon <- gsub("\\]", "", GU_taxon)
print(length(GU_taxon)) #18765
head(GU_taxon) # without ]

# Approach 1, Venn diagram, minimal data manipulation
# Good for overview but no specifics
Taxon_new <- list(
  Gabon = GA_taxon,
  Guinea = GU_taxon
)
# Run through Venn Diagram - letting it self compute
# Taxon, 475 shared across both
Venn_taxon <- ggVennDiagram(Taxon_new, label = "both") +
  ggtitle("Venn Diagram Shared Taxon: Gabon & Guinea") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold") 
  )
ggsave("Venn_taxon.png", plot = Venn_taxon, width = 10, height = 8)

# Approach 2 - double checking using another method - minimal 
# Comparing TAXON GA column N vs. GU column M and obtaining common elements shared between.
# 475 shared 
Intersect_taxonNM <- intersect(GA_taxon, GU_taxon)

if (length(Intersect_taxonNM) > 0) {
  print(paste("Total number of common TAXON between GA col. N and GU col. M:", length(Intersect_taxonNM)))
} else {
  print("No TAXON shared between GA col. N v. GU col. M!: ")
}
typeof(Intersect_taxonNM) # Character
tail(Intersect_taxonNM) # only contains names 
length(Intersect_taxonNM) # 475

# Set Diff
# Comparing TAXON GA column N vs. GU column M and obtaining different elements to each each.
# 174 different
Diff_taxonNM<-setdiff(GA_taxon, GU_taxon)

if (length(Diff_taxonNM) > 0) {
  print(paste("Total number of different TAXON between GA col. N and GU col. M", length(Diff_taxonNM)))
} else {
  print("No different TAXON between GA col. N and v. GU col. M")
}

# Summary of taxon
taxon_summary <- data.frame(
  Data_Overview = c("Total Rows in GA", "Total Rows in GU", "Shared Taxon", "Different Taxon"),
  Count = c(length(GA_taxon), length(GU_taxon), length(Intersect_taxonNM), length(Diff_taxonNM))
)

# Converting into data frame
# Gabon, TAXON column N, col. 14
# Now getting specifics and frequencies in each df.
GA_taxon <- as.data.frame(table(GA_taxon))
head(GA_taxon)
colnames(GA_taxon) <- c("GA_Taxon" , "GA_Freq.")
print(sum(GA_taxon$Freq.))
GA_taxon <- GA_taxon[order(GA_taxon$GA_Freq.,decreasing = TRUE),]
head(GA_taxon)
GA_taxon_percent <- (GA_taxon$GA_Freq./ sum(GA_taxon$GA_Freq.) * 100)
head(GA_taxon_percent)
GA_taxon <- cbind(GA_taxon, GA_taxon_percent)
head(GA_taxon) # Now have GA_Taxon, GA_Freq. and GA_taxon_percent 
(sum(GA_taxon_percent)) # 100%
(sum(GA_taxon$GA_Freq.)) # 8911
typeof(GA_taxon) # list

# Guinea, taxon column M, col. 13
GU_taxon <- as.data.frame(table(GU_taxon))
head(GU_taxon)
colnames(GU_taxon) <- c("GU_Taxon" , "GU_Freq.")
print(sum(GU_taxon$GU_Freq.))
head(GU_taxon)
GU_taxon <- GU_taxon[order(GU_taxon$GU_Freq.,decreasing = TRUE),]
head(GU_taxon)
GU_taxon_percent <- (GU_taxon$GU_Freq./ sum(GU_taxon$GU_Freq.) * 100)
head(GU_taxon_percent)
GU_taxon <- cbind(GU_taxon, GU_taxon_percent)
head(GU_taxon) # Now have taxon, frequency and percent for Guinea
(sum(GU_taxon_percent)) # 100%
(sum(GU_taxon$GU_Freq)) # 18765

# Filter shared species from GA and GU 
GAGU_taxon <- list(GA_taxon, GU_taxon)
head(n = 3, GAGU_taxon) # Mass list with all 6 columns NOT filtered re: Shared
typeof(GAGU_taxon) # list
sum(GA_taxon$GA_taxon_percent) # 100
sum(GU_taxon$GU_taxon_percent) # 100 

# Filter shared species from GA and GU
GA_Shared <- GA_taxon %>% filter(GA_Taxon %in% Intersect_taxonNM)
GU_Shared <- GU_taxon %>% filter(GU_Taxon %in% Intersect_taxonNM)

print(sum(GA_Shared$GA_Freq.)) #8128
print(sum(GA_Shared$GA_taxon_percent)) # only adds up to 91%
head(n = 3, GA_Shared) # colnames GA_taxon, GA_Freq., GA_taxon_percent
head(n = 3, GU_Shared) # colnames GU_taxon, GU_Freq., GU_taxon_percent
nrow(GA_Shared) # 475
nrow(GU_Shared) # 475 
print(sum(GU_Shared$GU_Freq.)) # 18632
print(sum(GU_Shared$GU_taxon_percent)) # adds up to 99.29 %

# Merge SHARED species data between GA and GU; narrows down
Shared_Freq_Taxon <- merge(GA_Shared, GU_Shared, by.x = "GA_Taxon", by.y = "GU_Taxon")
head(Shared_Freq_Taxon)
Shared_Freq_Taxon <- Shared_Freq_Taxon %>% rename(Taxon = GA_Taxon) #r rename since same 475
head(Shared_Freq_Taxon)
# GA to GU in order; GA higher
Shared_Freq_Taxon_GA <- Shared_Freq_Taxon[order(Shared_Freq_Taxon$GA_Freq., decreasing = TRUE),]
head(Shared_Freq_Taxon_GA)

# GU to GA in order; GU higher 
Shared_Freq_Taxon_GU <- Shared_Freq_Taxon[order(Shared_Freq_Taxon$GU_Freq., decreasing = TRUE),]
head(Shared_Freq_Taxon_GU)

# Sum the frequencies for both GA and GU
sum(Shared_Freq_Taxon_GA$GA_Freq.) # 8128
sum(Shared_Freq_Taxon_GU$GU_Freq.) # 18632

# Now for pie charts, top 10 using % 
Top_10_Shared_GA <- head(n=10, Shared_Freq_Taxon_GA)
print(Top_10_Shared_GA)

Shared_GA <- Top_10_Shared_GA$GA_taxon_percent
labels_Shared_GA <- paste0(Top_10_Shared_GA$Taxon, " (", round(Top_10_Shared_GA$GA_taxon_percent, 2), "%)")
n <- length(Top_10_Shared_GA$Taxon)
col_pie <- scale_color_poke(pokemon = "Croconaw")$palette(n)

png("Top 10 Shared Species in Gabon.png", width = 900, height = 600)

pie(
  Shared_GA, 
  label = labels_Shared_GA, 
  col = col_pie,
  main = "Top 10 Shared Species in Gabon"
)

dev.off()

# Top 10 respective to Guinea using %
Top_10_Shared_GU<- head(n=10, Shared_Freq_Taxon_GU)
print(Top_10_Shared_GU)

Shared_GU <- Top_10_Shared_GU$GU_taxon_percent
labels_Shared_GU <- paste0(Top_10_Shared_GU$Taxon, " (", round(Top_10_Shared_GU$GU_taxon_percent, 2), "%)")
n <- length(Top_10_Shared_GU$Taxon)
col_pie <- scale_color_poke(pokemon = "Weedle")$palette(n)

png("Top 10 Shared Species in Guinea.png", width = 900, height = 600)

pie(
  Shared_GU, 
  label = labels_Shared_GU, 
  col = col_pie,
  main = "Top 10 Shared Species in Guinea"
)

dev.off()

write_xlsx(list(
  "Taxon_Summary" = taxon_summary,
  "Shared_Taxon_GA" = data.frame(Shared_Freq_Taxon_GA),
  "Shared_Taxon_GU" = data.frame(Shared_Freq_Taxon_GU),
  "GA_Freq_Taxon" = data.frame(GA_taxon),
  "GU_Freq_Taxon" = data.frame(GU_taxon),
  "Different_Taxon" = data.frame(Diff_taxonNM)
), "Taxon_Summary.xlsx")

