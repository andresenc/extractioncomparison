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

# Chemical KW
rna.pca.nine <- data.frame(PC1 = pca_object$x[,1],
                           PC2 = pca_object$x[,2],
                           PC3 = pca_object$x[,3],
                           PC4 = pca_object$x[,4],
                           PC5 = pca_object$x[,5],
                           PC6 = pca_object$x[,6],
                           PC7 = pca_object$x[,7],
                           PC8 = pca_object$x[,8],
                           PC9 = pca_object$x[,9],
                           PC9 = pca_object$x[,10])

# set explained variance as colnames
colnames(rna.pca.nine) <- paste0(as.character(round((pca_object$sdev^2 / 
                                                       sum(pca_object$sdev^2)*100), 
                                                    digits = 0)), " %")[c(1:10)]
chemical_df <- select(pca_df, c(Methods))
chemical_df$ChCl3 <- ifelse(grepl("ChCl3", chemical_df$Methods), "yes", "no")
chemical_df$PP <- ifelse(grepl("PP", chemical_df$Methods), "yes", "no")
chemical_df$MTBE <- ifelse(grepl("MTBE", chemical_df$Methods), "yes", "no")
chemical_df$EtOH <- ifelse(grepl("EtOH", chemical_df$Methods), "yes", "no")
chemical_df$ACN <- ifelse(grepl("ACN", chemical_df$Methods), "yes", "no")
chemical_df$MeOH <- ifelse(grepl("MeOH", chemical_df$Methods), "yes", "no")
chemical_df$IPA <- ifelse(grepl("IPA", chemical_df$Methods), "yes", "no")

chemical_df <- select(chemical_df, -Methods)

non_cont_test_list <- apply(chemical_df, MARGIN = 2, FUN = function(feat){
  na.ind <- which(is.na(feat))
  if(is.integer(na.ind) && length(na.ind) == 0){
    feat <- feat
  }else{
    feat <- feat[-na.ind] 
  }
  test_sig <- apply(pca_object$x[, c(1:10)], MARGIN = 2,
                    FUN = function(sig){
                      if(is.integer(na.ind) && length(na.ind) == 0){
                        sig <- sig
                      }else{
                        sig <- sig[-na.ind]
                      }
                      kru.test <- kruskal.test(sig, factor(feat))
                      temp_out <- data.frame(kru.test$p.value)
                      return(temp_out)
                    })
  sig.test <- do.call(rbind, test_sig)
  rownames(sig.test) <- paste0("PC ", c(1:10), " (", colnames(rna.pca.nine), ")")
  colnames(sig.test) <- "p.value"
  return(sig.test)
})
kru_p.value_df <- do.call(cbind, non_cont_test_list)
colnames(kru_p.value_df) <- names(non_cont_test_list)

# Visualization
rownames(kru_p.value_df) <- gsub("s", "expS", rownames(kru_p.value_df))
names(kru_p.value_df) <- gsub("clin", "D", names(kru_p.value_df))
p.abs.cat.m <- t(as.matrix(kru_p.value_df))
rownames(p.abs.cat.m)

draw(Heatmap(-log10(p.abs.cat.m),
             name = "-log10(pval)",
             col = colorRamp2(breaks = c(0, -log10(0.05), 10),
                              c("white", "white", muted("red"))),
             cluster_rows = FALSE, cluster_columns = FALSE,
             column_names_side = c("top"),
             rect_gp = gpar(col = "black"),
             heatmap_legend_param = list(direction = "horizontal")), 
     heatmap_legend_side = "bottom")

# Loadings plot
loadings_df <- data.frame(pca_object$rotation[, c(1:2)])
loadings_list <- apply(loadings_df, MARGIN = 2, FUN = function(x){
  in_df <- data.frame(x)
  temp_df <- rownames_to_column(in_df, var = "Metabolite") %>% 
    gather(PC, Loadings, -Metabolite) %>%
    arrange(desc(abs(Loadings)))
  out_df <- head(temp_df, 15)
  return(out_df)
})
for(i in c(1:length(loadings_list))){
  loadings_list[[i]]$PC <- paste0("PC", i)
}
plot_df <- do.call(rbind, loadings_list)

ggplot(plot_df, aes(x = PC, y = Loadings)) +
  geom_point(size = 5, aes(color = Loadings)) +
  geom_point(shape = 1, size = 5, colour = "black") +
  ylab("Loadings") +
  xlab("") +
  scale_color_gradient2() +
  theme_bw() +
  geom_text_repel(aes(x = PC, y = Loadings, label = Metabolite), color = "black", 
                  max.iter = 15000, segment.size = 0.20, size = 2.3, 
                  segment.colour = "grey47", nudge_x = 0.5, max.overlaps = 1000)

pc_plot_list <- split(plot_df, f = plot_df$PC)

plot_list <- lapply(pc_plot_list, FUN = function(x){
  ggplot(x, aes(y = reorder(Metabolite, Loadings), x = Loadings)) +
    geom_bar(stat = "identity", aes(fill = Loadings)) +
    # geom_bar(colour = "black", stat = "identity") +
    facet_wrap(~ PC, scales = "free_x") +
    xlab("Loadings") +
    ylab("") +
    scale_fill_gradient2(low = muted("blue"), high = muted("red")) +
    theme_bw() +
    theme(legend.position = "none")
})

plot_list

  
