################################################################################
######### SCRIPT TO MAKE PCA and correct PCA ############################
################################################################################

#Libraries required:
library(ggplot2)
library(RColorBrewer)
library(DGEobj)
library(DGEobj.utils)
library(limma)
library(edgeR)
library(devtools)
library(ggfortify)
library(gdata)
library(stats)
library(dplyr)
library(pheatmap)
library(DESeq2) 
library(cowplot)
library(pheatmap)

#Set working directory:
setwd("C:/Users/jmper/OneDrive/Escritorio/tfm_bioinformatica/pca_corrmatrix/final_matrixes")

#Load the data, counts matrix obtained after executing do_counts script in the samples obtained by featureCounts.

counts <- read.table(file = "osteoporosis_bonemarrow_counts.txt",
                     sep = "\t",
                     header = TRUE,
                     row.names = 1)


#Now, execute PCA.
factors <- colnames(counts)
factors <- gsub("_[0-9]+(\\.[0-9]+)?$", "", factors) # Eliminate _1, _2, _n tags from samples, leaving unique name and grouping them into factors.


#3. Define groups.
group <- factor(factors)

#4. Normalize.
d <- DGEList(counts = counts, group = group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)

ncounts <-  cpm (d, normalized.lib.sizes = TRUE, log = TRUE) # We use log for better PCA results.

#5. Represent the PCA.
pca <- prcomp(t(ncounts))

#Percentage of each component.

pca.stat <- stats:::summary.prcomp(pca)$importance

#6. Vector of colours.
# Create vector making each factor unique.
uniq_factors <- unique(factors)

# Color palette.
n<- length(uniq_factors)
col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <-  unlist(mapply(brewer.pal, col_pals$maxcolors, rownames(col_pals)))
col_vector <- unique(col_vector)

# Assign one color to one factor.
colors_factors <- as.vector(factors)

#Generate colours for plots.

ggcolors <- c()

# Assign colours
for (i in 1:length(uniq_factors))
{
  colors_factors <- replace(colors_factors,
                            colors_factors == uniq_factors[i],
                            col_vector[i])
}

ggcolors <- unique(colors_factors)
names(ggcolors) <- uniq_factors



#8. Plotting with ggplot2.

axis_pca <- pca$sdev^2 / sum(pca$sdev^2)

pca_x <- as.data.frame(pca$x)

pca_x <- pca_x[, c(1,2)] # We only need first 2 dimensions.

pca_x$muestra <- factors

gg_pca <- ggplot(data = pca_x,
                 aes(x = PC1, y = PC2, color = muestra)) +
  geom_point(size = 3.5) +
  geom_text(
    label=c(group), 
    stat="identity",
    color = "black",
    check_overlap = TRUE,
    size = 2.5) +
  
  #scale_color_manual(values = ggcolors) + ## ggcolors en un vector con nombre en el que defino los colores
  labs(x = paste0("PC1 (",round(axis_pca[1]*100), "%)"),
       y = paste0("PC2 (",round(axis_pca[2]*100), "%)"),
       color = "") + # nombre leyenda
  guides(color = guide_legend(nrow = 4, # número filas leyenda
                              override.aes = list(size = 4))) + # tamaño puntos leyenda
  scale_shape_manual(values = 1:20)+
  
  ggtitle("PCA datasets") +
  
  theme_bw() +
  theme(legend.position="bottom",
        legend.title = element_text(size = 0), # tamaño titulo leyenda
        legend.text = element_text(size = 8), # tamaño texto leyenda
        axis.text = element_text(size = 7), # tamaño texto ejes
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
  )

gg_pca # Para imprimirlo en pantalla


################################################################################
######################## CORRELATION MATRIX ####################################
################################################################################

# Create a mapping of samples by extracting the sample name before the underscore
sample_names <- sub("_[0-9]+(\\.[0-9]+)?$", "", colnames(counts))  

# Aggregate counts by sample (sum or mean)
merged_counts <- as.data.frame(t(sapply(unique(sample_names), function(sample) {
  rowMeans(counts[, sample_names == sample])
})))

# Normalize the counts (optional)
log_counts <- log2(merged_counts + 1)  # log1p avoids log(0) issues

# Compute Spearman correlation matrix
spearman_cor <- cor(t(log_counts), method = "spearman") 


# Save the correlation matrix
write.table(spearman_cor, "spearman_correlation_matrix.txt", sep = "\t", quote = FALSE)

# Visualize as a heatmap (optional)
spearman_hmap <- pheatmap(
  spearman_cor,
  fontsize_row = 8,  # Adjust row label size
  fontsize_col = 8,  # Adjust column label size
  angle_col = 45     # Rotate column labels to avoid cutting
  
)




