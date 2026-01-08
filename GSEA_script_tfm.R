###############################################################


# 1. Load libraries

library(dplyr)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# Set your working directory (The structure is always ~/results_selected/disease_results/Common_results) just change disease_results.
setwd("C:/Users/jmper/OneDrive/Escritorio/tfm_bioinformatica/results_selected/Ovx_results/Common_results")


## 2. Load and Transform Data

# Loading the table with proper Ensembl rownames
df <- read.table("hunter_results_table.txt", sep = "\t", header = TRUE, row.names = 1)

# 1. Remove rows tagged as "FILTERED_OUT" by ExpHunter Suite.
# 2. Clean Ensembl IDs (strips version suffixes like .1)
# 3. Create a robust Ranking Score (using meanlogFC)
df_clean <- df %>%
  filter(genes_tag != "FILTERED_OUT") %>%
  mutate(
    ensembl_id = gsub("\\..*$", "", rownames(.)),
    # Add a small value to combined_FDR to avoid log(0)
    score = mean_logFCs
  )

# Create the named numeric vector for GSEA
rnk_vector <- df_clean$score
names(rnk_vector) <- df_clean$ensembl_id
rnk_vector <- sort(rnk_vector, decreasing = TRUE)

# Verification check
message(paste("Successfully processed", length(rnk_vector), "genes for analysis."))


# 3. Retrieve Gene Sets (Mus musculus)

# Fetching directly from MSigDB in 'TERM2GENE' format

# Fetch all C5 (GO) terms first
all_go <- msigdbr(species = "Mus musculus", collection = "C5")

# Split into the three main ontologies
bp_t2g <- all_go %>% filter(gs_subcollection == "GO:BP") %>% dplyr::select(gs_name, ensembl_gene)
mf_t2g <- all_go %>% filter(gs_subcollection == "GO:MF") %>% dplyr::select(gs_name, ensembl_gene)
cc_t2g <- all_go %>% filter(gs_subcollection == "GO:CC") %>% dplyr::select(gs_name, ensembl_gene)

# Then KEGG and Reactome
kegg_t2g     <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:KEGG_LEGACY") %>% dplyr::select(gs_name, ensembl_gene)
reactome_t2g <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:REACTOME") %>% dplyr::select(gs_name, ensembl_gene)


# 4. Run GSEA for each database

set.seed(123)

run_gsea_wrapper <- function(t2g, rnk) {
  GSEA(geneList = rnk, 
       TERM2GENE = t2g,
       minGSSize = 10,
       maxGSSize = 500,
       pvalueCutoff = 0.05,
       pAdjustMethod = "BH",
       eps = 0)# eps=0 for accurate p-value estimation
}

message("Analyzing GO: Biological Process...")
res_bp <- run_gsea_wrapper(bp_t2g, rnk_vector)

message("Analyzing GO: Molecular Function...")
res_mf <- run_gsea_wrapper(mf_t2g, rnk_vector)

message("Analyzing GO: Cellular Component...")
res_cc <- run_gsea_wrapper(cc_t2g, rnk_vector)

message("Analyzing KEGG...")
res_kegg <- run_gsea_wrapper(kegg_t2g, rnk_vector)

message("Analyzing Reactome...")
res_reactome <- run_gsea_wrapper(reactome_t2g, rnk_vector)


# 5. Visualizations & Export


export_results <- function(gsea_obj, label) {
  if(!is.null(gsea_obj) && nrow(as.data.frame(gsea_obj)) > 0) {
    # 1. Save results table
    write.csv(as.data.frame(gsea_obj), paste0("GSEA_", label, "_results.csv"), row.names = FALSE)
    
    # 2. Generate Dotplot (Top 15 categories)
    p_dot <- dotplot(gsea_obj, showCategory = 15) + 
      ggtitle(paste(label, "Pathway Enrichment"))
    ggsave(paste0("GSEA_", label, "_dotplot.png"), p_dot, width = 10, height = 14) 
    
    # 3. Generate Ridgeplot (Distribution of scores)
    p_ridge <- ridgeplot(gsea_obj) + xlab("LogFC-based ranking")
    ggsave(paste0("GSEA_", label, "_ridgeplot.png"), p_ridge, width = 20, height = 15)
    
    message(paste("Saved results for", label))
  } else {
    message(paste("No significant pathways found for", label))
  }
}

export_results(res_bp, "GO_BP")
export_results(res_mf, "GO_MF")
export_results(res_cc, "GO_CC")
export_results(res_kegg, "KEGG")
export_results(res_reactome, "Reactome")


message("Done! Check your working directory for .csv and .png files.")
