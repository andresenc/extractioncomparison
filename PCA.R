library(dplyr)
library(magrittr)
library(ggplot2)
library(tidyr)

# Load input
load(file.path("input.RData"))

## Method names
method_df <- data.frame(Abb = c("1", "2", "3", "4", "5", 
                                "6", "7a", "7b", "8", "9"),
                        Name = c("100/30 IPA", "100 IPA", "MeOH/ACN/H2O", 
                                 "MeOH/ACN", "EtOH/PP", "100/20 MeOH", 
                                 "75 EtOH/MTBE A", "75 EtOH/MTBE B", 
                                 "MeOH/MTBE/MeOH", "MeOH/ChCl3"))

# Filter metabolites which are LOD in all samples (first rep NA by LOD)
lod_m <- as.matrix(flag_df)
lod_m[which(is.na(lod_m))] <- "LOD"
lod_num_m <- apply(lod_m, MARGIN = c(1,2), FUN = function(x) ifelse(x == "LOD", 1, 0))
lod_vec <- which(colMeans(lod_num_m) == 1)
raw_df <- raw_df[, -lod_vec]

# Filter metabolites which are 0 in all samples - already covered by lod filtering
# raw_df <- raw_df[, -which(colSums(raw_df) == 0)]

# Replace NA with zero - already covered by lod filtering
input_m <- as.matrix(raw_df)

# Zero imputation
imputed_m <- apply(input_m, MARGIN = 2, FUN = function(x){
  if(min(x, na.rm = T) == 0){
    temp <- x[-which(x == 0)]
    out <- x + (min(temp, na.rm = T)*0.2)
  }else{
    out <- x
  }
  return(out)
})

# Log transformation
plot_m <- log2(imputed_m)

# Pareto scaling
plot_m_scaled <- apply(plot_m, MARGIN = 1, function(x) (x - mean(x, na.rm = T)) / sqrt(sd(x, na.rm = T)))

# PCA 
pca_object <- prcomp(t(plot_m_scaled), center = TRUE)
pca_var <- as.character(round((pca_object$sdev^2 / sum(pca_object$sdev^2)*100), digits = 0))

replicate <- gsub("_.+$", "", gsub("[[:alnum:]]+\\.", "", colnames(plot_m_scaled)))
tissue <- gsub(".+_", "", colnames(plot_m_scaled))
methods <- gsub("\\..+", "", colnames(plot_m_scaled))

pca_df <- data.frame(PC1 = pca_object$x[,1],
                     PC2 = pca_object$x[,2],
                     PC3 = pca_object$x[,3],
                     PC4 = pca_object$x[,4],
                     Replicate = replicate,
                     Tissue = tissue,
                     Methods = methods)
pca_df$Methods <- factor(pca_df$Methods, levels = method_df$Name)
pca_df$Tissue <- factor(pca_df$Tissue, levels = c("Liver", "HEK", "HL60", "Bone Marrow"))

# Plot
ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(aes(color = Methods, shape = Tissue), size = 2) +
  scale_color_manual("Methods", values = col_vector_methods) +
  xlab(paste0("PC1: ", pca_var[1], "% variance")) +
  ylab(paste0("PC2: ", pca_var[2], "% variance")) +
  theme_classic()

# Only Liver #
################################################################################
raw_df <- raw_df[grep("Liver", rownames(raw_df)), ]
flag_df <- flag_df[grep("Liver", rownames(flag_df)), ]

# Filter metabolites which are LOD in all samples (first rep NA by LOD)
lod_m <- as.matrix(flag_df)
lod_num_m <- apply(lod_m, MARGIN = c(1,2), FUN = function(x) ifelse(x == "LOD", 1, 0))
lod_vec <- which(colMeans(lod_num_m) == 1)
raw_df <- raw_df[, -lod_vec]

# Filter metabolites which are 0 in all samples - already covered by lod filtering
# raw_df <- raw_df[, -which(colSums(raw_df) == 0)]

# Replace NA with zero - already covered by lod filtering
input_m <- as.matrix(raw_df)

# Zero imputation
imputed_m <- apply(input_m, MARGIN = 2, FUN = function(x){
  if(min(x, na.rm = T) == 0){
    temp <- x[-which(x == 0)]
    out <- x + (min(temp, na.rm = T)*0.2)
  }else{
    out <- x
  }
  return(out)
})

# Log transformation
plot_m <- log2(imputed_m)

# Pareto scaling
plot_m_scaled <- apply(plot_m, MARGIN = 1, function(x) (x - mean(x, na.rm = T)) / sqrt(sd(x, na.rm = T)))

# PCA
pca_object <- prcomp(t(plot_m_scaled), center = TRUE)
pca_var <- as.character(round((pca_object$sdev^2 / sum(pca_object$sdev^2)*100), digits = 0))

replicate <- gsub("_.+$", "", gsub("[[:alnum:]]+\\.", "", colnames(plot_m_scaled)))
tissue <- gsub(".+_", "", colnames(plot_m_scaled))
methods <- gsub("\\..+", "", colnames(plot_m_scaled))

pca_df <- data.frame(PC1 = pca_object$x[,1],
                     PC2 = pca_object$x[,2],
                     PC3 = pca_object$x[,3],
                     PC4 = pca_object$x[,4],
                     Replicate = replicate,
                     Tissue = tissue,
                     Methods = methods)
pca_df$Methods <- factor(pca_df$Methods, levels = method_df$Name)
pca_df$Tissue <- factor(pca_df$Tissue, levels = c("Liver", "HEK", "HL60", "Bone Marrow"))

# Plot
gg_obj <- ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(aes(color = Methods), size = 2) +
  scale_color_manual("Methods", values = col_vector_methods) +
  xlab(paste0("PC1: ", pca_var[1], "% variance")) +
  ylab(paste0("PC2: ", pca_var[2], "% variance")) +
  theme_classic()

gg_obj

