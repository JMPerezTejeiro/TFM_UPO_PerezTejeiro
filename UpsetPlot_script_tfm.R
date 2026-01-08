################################################################################
#     Script to convert ENSEMBL IDs into GENE SYMBOLS AND VENN DIAGRAM         #
################################################################################

#Libraries:

library(ComplexHeatmap) # Must be cited in the paper.


# Set working directory:
setwd("C:/Users/jmper/OneDrive/Escritorio/tfm_bioinformatica/prevalent_genes")

# Read each file containing the prevalent genes for each dataset and transform them into vectors.
aged_male <- read.table("Aged_male_prevalent.txt", stringsAsFactors = FALSE)[[1]]
aged_female <- read.table("Aging_female_prevalent.txt", stringsAsFactors = FALSE)[[1]]
aged_HLU <- read.table("Aged_AgedDisuse_male_prevalent.txt", stringsAsFactors = FALSE)[[1]]
male_HLU_old <- read.table("Young_AgedHLU_male_prevalent.txt", stringsAsFactors = FALSE)[[1]]
male_HLU_young <- read.table("Young_YoungHLU_male_prevalent.txt", stringsAsFactors = FALSE)[[1]]
giop <- read.table("GIOP_f1_prevalent.txt", stringsAsFactors = FALSE)[[1]]
fem_HLU_young <- read.table("Disuse_young_female_prevalent.txt", stringsAsFactors = FALSE)[[1]]
hfd <- read.table("HFD_prevalent.txt", stringsAsFactors = FALSE)[[1]]
ovx <- read.table("Ovx_prevalent.txt", stringsAsFactors = FALSE)[[1]]
piezo1 <- read.table("Piezo1_prevalent.txt", stringsAsFactors = FALSE)[[1]]
t1d <- read.table("T1D_prevalent.txt", stringsAsFactors = FALSE)[[1]]
tbi <- read.table("TBI_prevalent.txt", stringsAsFactors = FALSE)[[1]]
hlu_progression <- read.table("YoungHLU_AgedHLU_male_prevalent.txt", stringsAsFactors = FALSE)[[1]]

# Create a named list of gene sets
upset_list <- list(
  Aging_male = aged_male, 
  Aging_fem = aged_female, 
  Aged_HLU = aged_HLU, 
  Male_HLU_old = male_HLU_old,
  Male_HLU_young = male_HLU_young,
  GIOP = giop,
  Fem_HLU_young = fem_HLU_young,
  HFD = hfd,
  OVX = ovx,
  Piezo1 = piezo1,
  T1D = t1d,
  TBI = tbi,
  HLU_progression = hlu_progression
)


# Create a binary matrix from our list, using the function make_comb_mat and the mode distinct to find shared genes:

binary_matrix <- make_comb_mat(upset_list, mode = "distinct")

# Explore combinations and sizes:
binary_matrix
comb_size(binary_matrix)
set_size(binary_matrix)

# --- Compute colors based on the number of datasets per intersection ---
num_datasets_dist <- comb_degree(binary_matrix)  # number of datasets per intersection

# --- Define discrete colors for each number of shared datasets ---
# You can choose your own palette here
degree_colors <- c(
  "1" = "blue3",
  "2" = "green3",
  "3" = "red3",
  "4" = "pink3",
  "5" = "purple3",
  "6" = "yellow3"
)

# --- Assign a color to each intersection based on degree ---
comb_colored_dist <- degree_colors[as.character(num_datasets_dist)]

# Plot the UpSet diagram
UpSet(
  binary_matrix,
  set_order = order(set_size(binary_matrix), decreasing = TRUE),
  comb_order = order(comb_size(binary_matrix), decreasing = TRUE),
  comb_col = comb_colored_dist,
  bg_pt_col = "ivory2",
  bg_col="ivory",
  pt_size = unit(14, "point"), lwd = 5,
  top_annotation = HeatmapAnnotation(
    "Genes in each \n intersection" = anno_barplot(comb_size(binary_matrix), 
                                       border = FALSE,
                                 height = unit(7, "cm"),
                                 bar_width = 1,
                                 add_numbers = TRUE,
                                       gp = gpar(fill = "azure4")),
    annotation_name_side = "left"
  ),
  right_annotation = rowAnnotation(
    "Number of \nprevalent genes" = anno_barplot(set_size(binary_matrix), 
                              border = FALSE, 
                              add_numbers = TRUE,
                              height = unit(15, "cm"),
                              width = unit(4, "cm"),
                              bar_width = 1,
                              axis_param = list(at = c(0, 1000, 2000), labels_rot = 45),
                              row_names_gp = gpar(fontsize = 12),
                              gp = gpar(fill = "deepskyblue")),
    annotation_name_side = "bottom"
  )
)

# Extract all shared genes across conditions:

all_combinations <- comb_name(binary_matrix)  # List all combination names

all_genes <- lapply(all_combinations, function(nm) extract_comb(binary_matrix, nm))

# Combine into a named list
names(all_genes) <- all_combinations

# Save the results
for (nm in names(all_genes)) {
  write.table(all_genes[[nm]],
              file = paste0("intersection_", nm, ".txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Extract only genes from combinations with 4 or more shared genes:

more4_combinations <- comb_name(binary_matrix)[comb_degree(binary_matrix) >= 4]
more4_genes <- lapply(more4_combinations, function(nm) extract_comb(binary_matrix, nm))

# Combine into a named list
names(more4_genes) <- more4_combinations

# Save the results
for (nm in names(more4_genes)) {
  write.table(all_genes[[nm]],
              file = paste0("intersection_", nm, ".txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}


################################################################################
#--------------------------- INERSECT MODE ------------------------------------#
################################################################################

# Create a binary matrix from our list, using the function make_comb_mat and the mode intersect to find shared genes:

binary_matrix_int <- make_comb_mat(upset_list, mode = "intersect")

# Explore combinations and sizes:
binary_matrix_int
comb_size(binary_matrix_int)
set_size(binary_matrix_int)

# --- Compute colors based on the number of datasets per intersection ---
num_datasets_int <- comb_degree(binary_matrix_int)  # number of datasets per intersection

# --- Define discrete colors for each number of shared datasets ---
# You can choose your own palette here
degree_colors_int <- c(
  "1" = "blue2",
  "2" = "green2",
  "3" = "red2",
  "4" = "pink2",
  "5" = "purple2",
  "6" = "yellow2"
)

# --- Assign a color to each intersection based on degree ---
comb_colored_int <- degree_colors_int[as.character(num_datasets_int)]

# Plot the UpSet diagram
UpSet(
  binary_matrix_int,
  set_order = order(set_size(binary_matrix_int), decreasing = TRUE),
  comb_order = order(comb_size(binary_matrix_int), decreasing = TRUE),
  comb_col = comb_colored_int,
  bg_pt_col = "ivory2",
  bg_col="ivory",
  pt_size = unit(6, "point"), lwd = 3,
  top_annotation = HeatmapAnnotation(
    "Genes in each \n intersection" = anno_barplot(comb_size(binary_matrix_int), 
                                                   border = FALSE,
                                                   height = unit(7, "cm"),
                                                   bar_width = 1,
                                                   add_numbers = FALSE,
                                                   gp = gpar(fill = "azure4")),
    annotation_name_side = "left"
  ),
  right_annotation = rowAnnotation(
    "Number of \nprevalent genes" = anno_barplot(set_size(binary_matrix_int), 
                                                 border = FALSE, 
                                                 add_numbers = TRUE,
                                                 height = unit(15, "cm"),
                                                 width = unit(4, "cm"),
                                                 bar_width = 1,
                                                 axis_param = list(at = c(0, 1000, 2000), labels_rot = 45),
                                                 row_names_gp = gpar(fontsize = 12),
                                                 gp = gpar(fill = "deepskyblue")),
    annotation_name_side = "bottom"
  )
)

# Extract all shared genes across conditions:
all_combinations_int <- comb_name(binary_matrix_int)  # List all combination names

all_genes_int <- lapply(all_combinations_int, function(nm) extract_comb(binary_matrix_int, nm))

# Combine into a named list
names(all_genes_int) <- all_combinations_int

# Save the results
for (nm in names(all_genes_int)) {
  write.table(all_genes[[nm]],
              file = paste0("intersection_", nm, ".txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Extract only genes from combinations with 4 or more shared genes:

more5_combinations <- comb_name(binary_matrix_int)[comb_degree(binary_matrix_int) >= 5]
more5_genes <- lapply(more5_combinations, function(nm) extract_comb(binary_matrix_int, nm))

# Combine into a named list
names(more5_genes) <- more5_combinations

# Save the results
for (nm in names(more5_genes)) {
  write.table(all_genes_int[[nm]],
              file = paste0("modeintersect_", nm, ".txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}
