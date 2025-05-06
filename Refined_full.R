
# Refined FULL Input MicrobeAnnotator

library("readxl")
library("writexl")
library("data.table")
library("stringr")
library("Biostrings")
library("ggVennDiagram")
library("ggplot2")

setwd("/Users/nancytorres/ResearchProject/Refined_MB_Full")
# sheet names are Gabon, Guinea, PBM, NBM and SM

# Read in excel
Excel_path <- "/Users/nancytorres/ResearchProject/Refined_MB_Full/Refined_Input_Full.xlsx"
sheet_names <- excel_sheets(Excel_path)

df_list <- list()
for (sheet in sheet_names) {
  df_list[[sheet]] <- read_excel(Excel_path, sheet = sheet)
}

# Read in fasta files 
Fasta_path <- "/Users/nancytorres/ResearchProject/aa_files"
Fasta_files <- list.files(Fasta_path, pattern = "_aa\\.fa$", full.names = TRUE)

aa_list <- list()
for (file in Fasta_files) {
  name <- tools::file_path_sans_ext(basename(file))
  name <- sub("_aa$", "", name)
  aa_list[[name]] <- readAAStringSet(file)
}

# Pull AA names only 
aa_names <- list()
for (name in names(aa_list)) {
  aa_names[[name]] <- names(aa_list[[name]])
}

# Convert to data.table 
for (sheet_name in names(df_list)) {
  df_list[[sheet_name]] <- setDT(df_list[[sheet_name]])
}

# Pull unique annotations only
Unique_annotations <- list()
for (sheet_name in names(df_list)){
  Unique_annotations[[sheet_name]] <- unique(df_list[[sheet_name]]$query_id)
}
print(Unique_annotations[[sheet]])

# Databases annotated 
Refined_FullDB_Anno <- list()

for (sheet_name in names(df_list)) {
  matched_kofam <- sum(grepl("kofam", df_list[[sheet_name]]$database))
  matched_swiss <- sum(grepl("swissprot", df_list[[sheet_name]]$database))
  matched_refseq <- sum(grepl("refseq", df_list[[sheet_name]]$database))
  matched_trembl <- sum(grepl("trembl", df_list[[sheet_name]]$database))
  no_matches <- sum(grepl("NA", df_list[[sheet_name]]$database))
  Refined_FullDB_Anno[[sheet_name]] <- c(matched_kofam, matched_swiss, matched_refseq, matched_trembl, no_matches)
  Refined_FullDB_Anno <- as.data.frame(Refined_FullDB_Anno)
}

Refined_FullDB_Anno <- t(Refined_FullDB_Anno)
colnames(Refined_FullDB_Anno) <- c("Kofam", "Swissprot","Refseq","Trembl", "NA")
Refined_FullDB_Anno <- as.data.frame(Refined_FullDB_Anno)

# Pulling taxonomy 
Taxon_DB <- list()
Taxon_DB_final <- list()

for (sheet_name in names(df_list)){
  Taxon_DB[[sheet_name]] <- df_list[[sheet_name]]$taxonomy
  Taxon_DB[[sheet_name]] <- strsplit(Taxon_DB[[sheet_name]], split = "[;,]")
  Taxon_DB[[sheet_name]] <- lapply(Taxon_DB[[sheet_name]], function(x) gsub("^ ", "", x))
  Taxon_DB[[sheet_name]] <- table(unlist(Taxon_DB[[sheet_name]]))
  Taxon_DB[[sheet_name]] <- sort(Taxon_DB[[sheet_name]], decreasing = TRUE)
  Taxon_DB[[sheet_name]] <- subset(Taxon_DB[[sheet_name]], names(Taxon_DB[[sheet_name]]) != "NA")
  
  Taxon_DB_final[[sheet_name]] <- data.frame(Taxon_DB[[sheet_name]])
  colnames(Taxon_DB_final[[sheet_name]]) <- c("Taxon", "Freq.")
  colnames(Taxon_DB_final[[sheet_name]])
  print(Taxon_DB_final[[sheet_name]])
  Taxon_DB_final[[sheet_name]] <- data.frame(Taxon_DB_final[[sheet_name]])
}

tail(Taxon_DB[[sheet_name]])
head(Taxon_DB[[sheet_name]])
head(Taxon_DB[["Gabon"]])
head(Taxon_DB[["Guinea"]])
head(Taxon_DB[["PBM"]])
head(Taxon_DB[["SM"]])

# Pie chart with percents 
for (sheet_name in names(Taxon_DB)) {
  top_taxa <- head(Taxon_DB[[sheet_name]], 5)
  total <- sum(top_taxa)
  pct <- round(100 * top_taxa / total, 1)
  lbls <- paste0(names(top_taxa), " (", pct, "%)")
  
  png(filename = paste0("Refined_FullTop_5_Taxa_", sheet_name, ".png"), width = 1325, height = 1000, res = 200, bg = "gray57")
  
  par(col = "white")
  
  pie(
    top_taxa,
    labels = lbls,
    col = terrain.colors(length(top_taxa)),
    cex.main = 1.5, 

  )

  title(main = paste("Refined Full Top 5 Taxa (Overall) - ", sheet_name),
        col.main = "white", cex.main = 1.5, line = -1)
  
  dev.off()
}

# Full taxon rows
Taxon_DB_rows <- list()
Taxon_DB_complete <- list()

for (sheet_name in names(df_list)){
  Taxon_DB_rows[[sheet_name]] <- df_list[[sheet_name]]$taxonomy
  Taxon_DB_rows[[sheet_name]] <- table(unlist(Taxon_DB_rows[[sheet_name]]))
  Taxon_DB_rows[[sheet_name]] <- sort(Taxon_DB_rows[[sheet_name]], decreasing = TRUE)
  Taxon_DB_rows[[sheet_name]] <- subset(Taxon_DB_rows[[sheet_name]], names(Taxon_DB_rows[[sheet_name]]) != "NA")
  
  Taxon_DB_complete[[sheet_name]] <- data.frame(Taxon_DB_rows[[sheet_name]])
  colnames(Taxon_DB_complete[[sheet_name]]) <- c("Taxon", "Freq.")
  colnames(Taxon_DB_complete[[sheet_name]])
  print(Taxon_DB_complete[[sheet_name]])
  Taxon_DB_complete[[sheet_name]] <- data.frame(Taxon_DB_complete[[sheet_name]])
}

# Breaking apart Contigs / MAGs fr CDS using regexp. 
Separated_contigs <- list()
for (sheet_name in names(df_list)) {
  Separated_contigs[[sheet_name]] <- str_extract(df_list[[sheet_name]]$query_id, "_contig_[0-9]+|_MAG_[0-9]+")
}

# identifying duplicated CONTIGs
Dup_contigs <- list ()
Dup_contig_rows <- list()
Dup_contig_values <- list()

for (sheet_name in names(Separated_contigs)){
  contigs <- Separated_contigs[[sheet_name]]
  
  Dup_contigs[[sheet_name]] <- table(duplicated(contigs))
  Dup_values <- contigs[duplicated(contigs)]
  Dup_contig_values[[sheet_name]] <- Dup_values 
  Dup_contig_rows[[sheet_name]] <- contigs[contigs %in% Dup_values]
}

# identifying duplicates 
Dup_counts <- list()
Dup_rows <- list()

for (sheet_name in names(df_list)) {
  data <- df_list[[sheet_name]]
  Dup_counts[[sheet_name]] <- table(duplicated(data$query_id))
  Dup_ids <- data$query_id[duplicated(data$query_id)]
  Dup_rows[[sheet_name]] <- data[data$query_id %in% Dup_ids, ]
}

test <- unique(df_list[[sheet]]$database)
print(test)

df_list[[sheet_name]]$database <- trimws(df_list[[sheet_name]]$database)

# Back-Up Venn diagram
for (sheet_name in names(df_list)) {
  df <- df_list[[sheet_name]]  
  
  kofam_ids <- df$query_id[df$database == "kofam"]
  swissprot_ids <- df$query_id[df$database == "swissprot"]
  refseq_ids <- df$query_id[df$database == "refseq"]
  trembl_ids <- df$query_id[df$database == "trembl"]
  no_hit_ids <- df$query_id[df$database == "NA"]
  
  venn_input <- list(
    KOfam = kofam_ids,
    SwissProt = swissprot_ids,
    RefSeq = refseq_ids,
    TrEMBL = trembl_ids,
    NoMatch = no_hit_ids
  )
  
  venn_plot <- ggVennDiagram(venn_input, label = "both") +
    ggtitle(paste("Back-up Refined Full Annotation Databases Distribution: ", sheet_name)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(
    filename = paste0("Back_Up_Annotation_Refined_Full_Venn_", sheet_name, ".png"),
    plot = venn_plot,
    width = 8,
    height = 8, 
    bg = "gray57"
  )
}

# Clean - Venn diagram
for (sheet_name in names(df_list)) {
  clean_df <- df_list[[sheet_name]]  
  
  clean_kofam_ids <- clean_df$query_id[clean_df$database == "kofam"]
  clean_swissprot_ids <- clean_df$query_id[clean_df$database == "swissprot"]
  clean_refseq_ids <- clean_df$query_id[clean_df$database == "refseq"]
  
  clean_venn_input <- list(
    KOfam = clean_kofam_ids,
    SwissProt = clean_swissprot_ids,
    RefSeq = clean_refseq_ids
  )
  
  venn_plot <- ggVennDiagram(clean_venn_input, label = "both") +
    ggtitle(paste("Refined Full Annotation Databases Distribution: ", sheet_name)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(
    filename = paste0("Clean_Refined_Full_Annotation_Venn_", sheet_name, ".png"),
    plot = venn_plot,
    width = 8,
    height = 8,
    bg = "gray57"
  )
}

# Finding unclassified taxonomy 
Unclass_rows_all <- list()
Unclass_counts <- list()
total_unclass_rows <- 0
Unclass_names <- list()
Unclass_taxa <- list()
Split_taxa <- list()
Taxa_counts <- list()

for (sheet_name in names(df_list)) {
  sheet_data <- df_list[[sheet_name]]
  matches <- grepl("unclass", sheet_data$taxonomy, ignore.case = TRUE)
  
  if (any(matches)) {
    Unclass_rows_all[[sheet_name]] <- sheet_data[matches, ]
    Unclass_counts[[sheet_name]] <- sum(matches)
    total_unclass_rows <- total_unclass_rows + sum(matches)
    
    Unclass_names <- Unclass_rows_all[[sheet_name]]$taxonomy
    Unclass_taxa[[sheet_name]] <- unlist(Unclass_names)
    Split_taxa[[sheet_name]] <- Unclass_taxa[[sheet_name]][grepl("unclass", Unclass_taxa[[sheet_name]], ignore.case = TRUE)]
    Split_taxa[[sheet_name]] <- data.frame(Split_taxa[[sheet_name]])
    Taxa_counts[[sheet_name]] <- table(Split_taxa[[sheet_name]])
    Taxa_counts[[sheet_name]] <-data.frame(Taxa_counts[[sheet_name]])
    colnames(Taxa_counts[[sheet_name]]) <- c("Unclassified Taxa", "Count")
  } else {
    Unclass_rows_all[[sheet_name]] <- "None"
    Unclass_counts[[sheet_name]] <- 0
  }
}

Unclass_rows_all[["Gabon"]]
Unclass_taxa[["Gabon"]]

# Bar graph with Unclass Taxa counts 
New_taxa <- list()

for(sheet_name in names(Unclass_taxa)){
  taxa_names <- Unclass_taxa[[sheet_name]]
  split_taxa_names <- strsplit(taxa_names, ",")
  unlist_split <- unlist(split_taxa_names)
  unlist_names <- unlist_split[grepl("unclassified", unlist_split, ignore.case = TRUE)]
  taxa_table <- sort(table(unlist_names), decreasing = TRUE)
  New_taxa[[sheet_name]] <- data.frame(
    New_taxa = names(taxa_table),
    Freq. = as.numeric(taxa_table)
  )

  New_taxa_bar_data <- ggplot(New_taxa[[sheet_name]], aes(x = reorder(New_taxa, Freq.), y = Freq.)) +
    geom_col(fill = "gold") +
    coord_flip() +
    geom_text(aes(label = Freq.), 
              position = position_dodge(width = 0.8), 
              vjust = 0.5, hjust = 0, angle = 0, color = "white", size = 3.5) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_minimal() +
    theme(
          axis.text.y =  element_text(face = "italic", color = "white"), 
          plot.title = element_text(hjust = 1, face = "bold", color = "white"),
          axis.text.x = element_text(color = "white"),
          axis.title = element_text(color = "white"),
          panel.grid.major = element_line(color = "gray50", size = 0.3),
          panel.grid.minor = element_blank(),
          plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    ) +
    labs(title = paste0("Refined Full: Unclassified Species ", sheet_name))
  
  ggsave(paste0("Refined Full Unclassified Species ", sheet_name, ".png"), plot = New_taxa_bar_data,
         width = 8, height = 5, units = "in", bg = "gray57", dpi = 300)
}

head(New_taxa[["Gabon"]])

# View results
total_unclass_rows
Unclass_rows_all[["Gabon"]]
nrow(Unclass_rows_all[["Gabon"]])
Unclass_counts <- data.frame(Unclass_counts)
Unclass_counts <- t(Unclass_counts)
colnames(Unclass_counts) <- c("Unclassified")
colnames(Taxa_counts[[sheet_name]])

for (sheet in names(Taxa_counts)) {
  Taxa_counts[[sheet]] <- Taxa_counts[[sheet]][order(-Taxa_counts[[sheet]]$Count), ]
}

# Isolating all Known Col. C
Known_ColC_DB <- list()

for (sheet_name in names(df_list)) {
  ColC_list <- df_list[[sheet_name]]$product
  ColC_table <- table(ColC_list)
  ColC_sorted <- sort(ColC_table, decreasing = TRUE)
  ColC_filtered <- subset(ColC_sorted, names(ColC_sorted) != "NA")
  
  Known_ColC_DB[[sheet_name]] <- unlist(ColC_filtered)
}

# KNOWN: Top 5 Products Swissprot and refseq, Column C.
Known_Products_ColC <- list()

for (sheet_name in names(df_list)){
  Known_Products_ColC[[sheet_name]] <- head(Known_ColC_DB[[sheet_name]])
  Known_Products_ColC[[sheet_name]] <- data.frame(
    Known_Product = names(Known_Products_ColC[[sheet_name]]),
    Freq. = as.numeric(Known_Products_ColC[[sheet_name]])
  )
  
  ColC_bar_data <- ggplot(Known_Products_ColC[[sheet_name]], aes(x = reorder(Known_Product, Freq.), y = Freq.)) +
    geom_col(fill = "blue") +
    coord_flip() +
    geom_text(aes(label = Freq.), 
              position = position_dodge(width = 0.8), 
              vjust = 0.5, hjust = 0, angle = 0, color = "white", size = 3.5) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(color = "white"),
      axis.text.x = element_text(color = "white"),
      axis.title = element_text(color = "white"),
      plot.title = element_text(hjust = 1, face = "bold", color = "white"),
      panel.grid.major = element_line(color = "gray50", size = 0.3),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
      ) + 
    labs(title = paste0("Col. C: Refined Full Top 5 Protein Products Counts - ", sheet_name))
  
  ggsave(paste0("Col. C: Refined_Full_Top_5_Protein_Products", sheet_name, ".png"), plot = ColC_bar_data,
         width = 8, height = 6, units = "in", bg = "gray57", dpi = 300)
}

print(Known_Products_ColC[["Gabon"]])

# Column C: UNknowns in products swissprot column
# What are the top proteins listed and those uncharacterized or hypothetical 
All_proteins <- list()
Top_proteins <- list()
Hypo_counts <- list()
Unchar_rows_all <- list()
Unchar_counts <- list()
total_unchar_rows <- 0

for (sheet_name in names(df_list)) {
  sheet_data <- df_list[[sheet_name]]
  
  All_proteins[[sheet_name]] <- sheet_data$product
  Top_proteins[[sheet_name]] <- sort(table(All_proteins[[sheet_name]]), decreasing = TRUE)
  
  is_unchar <- grepl("hypothetical|uncharacterized", sheet_data$product, ignore.case = TRUE)
  Both <- sheet_data$product[is_unchar]
  Hypo_counts[[sheet_name]] <- sort(table(Both), decreasing = TRUE)
  
  if (length(Both) > 0) {
    Unchar_rows_all[[sheet_name]] <- sheet_data[is_unchar, ]
    Unchar_counts[[sheet_name]] <- sum(is_unchar)
    total_unchar_rows <- total_unchar_rows + sum(is_unchar)
  } else {
    Unchar_rows_all[[sheet_name]] <- "None"
    Unchar_counts[[sheet_name]] <- 0
  }
}

# Bar graph of Unknown Hypothetical proteins, Column C 
# Pulling Top 5 Hypothetical proteins and Uncharacterized
Final_Hypo <- list()

for(sheet_name in names(df_list)){
  Final_Hypo[[sheet_name]] <- head(Hypo_counts[[sheet_name]])
  Final_Hypo[[sheet_name]] <- data.frame(
    Hypo_protein = names(Final_Hypo[[sheet_name]]),
    Freq. = as.numeric(Final_Hypo[[sheet_name]])
  )
  
  Hypo_bar_data <- ggplot(Final_Hypo[[sheet_name]], aes(x = reorder(Hypo_protein, Freq.), y = Freq.)) +
    geom_col(fill = "skyblue1") +
    coord_flip() +
    geom_text(aes(label = Freq.), 
              position = position_dodge(width = 0.8), 
              vjust = 0.5, hjust = 0, angle = 0, color = "white", size = 3.5) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(color = "white"),
          axis.text.x = element_text(color = "white"),
          axis.title = element_text(color = "white"),
          plot.title = element_text(hjust = 1, face = "bold", color = "white"),
          panel.grid.major = element_line(color = "gray50", size = 0.3),
          panel.grid.minor = element_blank(),
          plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
      ) +
    labs(title = paste0("Col. C: Refined Full Top 5 Hypothetical/Uncharacteristic Proteins - ", sheet_name))
  
  ggsave(paste0("Col. C: Refined_Full_Protein_Unchar_Product_", sheet_name, ".png"), plot = Hypo_bar_data,
         width = 8, height = 6, units = "in", bg = "gray57", dpi = 300)
}

# Pulling entire KO Products Kofam; Column E
KO_DB <- list()

for (sheet_name in names(df_list)) {
  ko_list <- df_list[[sheet_name]]$ko_product
  ko_table <- table(ko_list)
  ko_sorted <- sort(ko_table, decreasing = TRUE)
  ko_filtered <- subset(ko_sorted, names(ko_sorted) != "NA")
  
  KO_DB[[sheet_name]] <- unlist(ko_filtered)
}

# KNOWN: Pulling top 5 from Ko Product Kofam Database and some RefSeq, Column E
Final_KO <- list()

for (sheet_name in names(df_list)){
  Final_KO[[sheet_name]] <- head(KO_DB[[sheet_name]])
  Final_KO[[sheet_name]] <- data.frame(
    KO_Product = names(Final_KO[[sheet_name]]),
    Freq. = as.numeric(Final_KO[[sheet_name]])
  )
  
  KO_bar_data <- ggplot(Final_KO[[sheet_name]], aes(x = reorder(KO_Product, Freq.), y = Freq.)) +
    geom_col(fill = "maroon") +
    coord_flip() +
    geom_text(aes(label = Freq.), 
              position = position_dodge(width = 0.8), 
              vjust = 0.5, hjust = 0, angle = 0, color = "white", size = 3.5) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(color = "white"),
      axis.text.x = element_text(color = "white"),
      axis.title = element_text(color = "white"),
      plot.title = element_text(hjust = 1, face = "bold", color = "white"),
      panel.grid.major = element_line(color = "gray50", size = 0.3),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)) + 
    labs(title = paste0("Col. E: Refined Full Top 5 Known Products - ", sheet_name))
  
  ggsave(paste0("Col. E: Refined Full Top 5 Known Products ", sheet_name, ".png"), plot = KO_bar_data,
         width = 8, height = 6, units = "in", bg = "gray57", dpi = 300)
}

print(Final_KO[["Gabon"]])

# Bar graph of what consists of Hypothetical proteins; Column E Ko Product 
# Pulling Top 5 Hypothetical proteins and Uncharacterized
# What are the top proteins listed and those uncharacterized or hypothetical 
All_proteins_ColE <- list()
Top_proteins_ColE <- list()
Hypo_counts_ColE <- list()
Unchar_rows_all_ColE <- list()
Unchar_counts_ColE <- list()
total_unchar_rows_ColE <- 0

for (sheet_name in names(df_list)) {
  sheet_data <- df_list[[sheet_name]]
  
  All_proteins_ColE[[sheet_name]] <- sheet_data$ko_product
  Top_proteins_ColE[[sheet_name]] <- sort(table(All_proteins_ColE[[sheet_name]]), decreasing = TRUE)
  
  is_unchar_ColE <- grepl("hypothetical|uncharacterized", sheet_data$ko_product, ignore.case = TRUE)
  Both_ColE <- sheet_data$ko_product[is_unchar_ColE]
  Hypo_counts_ColE[[sheet_name]] <- sort(table(Both_ColE), decreasing = TRUE)
  
  if (length(Both_ColE) > 0) {
    Unchar_rows_all_ColE[[sheet_name]] <- sheet_data[is_unchar_ColE, ]
    Unchar_counts_ColE[[sheet_name]] <- sum(is_unchar_ColE)
    total_unchar_rows_ColE <- total_unchar_rows_ColE + sum(is_unchar_ColE)
  } else {
    Unchar_rows_all_ColE[[sheet_name]] <- "None"
    Unchar_counts_ColE[[sheet_name]] <- 0
  }
}

# UNKNOWN: Pulling top 5 uncharacterized/hypothetical proteins from Kofam Product, Column E
Final_Col_E <- list()

for(sheet_name in names(df_list)){
  Final_Col_E[[sheet_name]] <- head(Hypo_counts_ColE[[sheet_name]])
  Final_Col_E[[sheet_name]] <- data.frame(
    Hypo_protein_ColE = names(Final_Col_E[[sheet_name]]),
    Freq. = as.numeric(Final_Col_E[[sheet_name]])
  )
  
  
  Hypo_bar_data_ColE <- ggplot(Final_Col_E[[sheet_name]], aes(x = reorder(Hypo_protein_ColE, Freq.), y = Freq.)) +
    geom_col(fill = "pink") +
    coord_flip() +
    geom_text(aes(label = Freq.), 
              position = position_dodge(width = 0.8), 
              vjust = 0.5, hjust = 0, angle = 0, color = "white", size = 3.5) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(color = "white"),
      axis.text.x = element_text(color = "white"),
      axis.title = element_text(color = "white"),
      plot.title = element_text(hjust = 1, face = "bold", color = "white"),
      panel.grid.major = element_line(color = "gray50", size = 0.3),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    ) +
    labs(title = paste0("Col. E: Refined Full Hypothetical/Uncharacteristic KO Proteins - ", sheet_name))
  
  ggsave(paste0("Col. E: Refined Full Top 5 Unknown Products ", sheet_name, ".png"), plot = Hypo_bar_data_ColE,
         width = 8, height = 6, units = "in", bg = "gray57", dpi = 300)
}

# Summary 
MicrobeAnnotator_Refined_Fullsummary <- data.table(
  Dataset = names(df_list),
  Un = Unclass_counts,
  Hypothetical = sapply(Unchar_rows_all, nrow),
  Original_Rows = sapply(df_list, nrow),
  Unique_Rows = sapply(Unique_annotations, length),
  Duplicates_plus = sapply(Dup_rows, nrow),
  Duplicated_contigs = sapply(Dup_contig_values, length),
  Annotation_counts = Refined_FullDB_Anno
)

excel_summary <- list("Full Summary" = MicrobeAnnotator_Refined_Fullsummary)

for (sheet in names(df_list)){
  excel_summary[[paste(sheet, "Originals")]] <- df_list[[sheet]]
  excel_summary[[paste(sheet, "Unclass_Taxa_Rows")]] <- data.frame(Unclass_rows_all[[sheet]])
  excel_summary[[paste(sheet, "Unclass_Taxa_only")]] <- data.frame(Taxa_counts[[sheet]])
  excel_summary[[paste(sheet, "All_Taxonomy")]] <- data.frame(Taxon_DB_complete[[sheet]])
  excel_summary[[paste(sheet, "Uncharacterized_Proteins")]] <- data.frame(Unchar_rows_all[[sheet]])
  excel_summary[[paste(sheet, "Unique")]] <- data.frame(Unique_annotations[[sheet]])
  excel_summary[[paste(sheet, "Duplicates+")]] <- Dup_rows[[sheet]][order(query_id)]
  excel_summary[[paste(sheet, "Dup Contigs")]] <- data.frame(Dup_contig_values[[sheet]])
}

write_xlsx(excel_summary, "Refined_Full_Summary.xlsx")

