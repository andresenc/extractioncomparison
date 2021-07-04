library(dplyr)
library(stringr)
library(agricolae)
library(tidyr)
library(tibble)

source("readMetIDQ.R")

## load raw data ## - output from MetIDQ
outListBM <- readMetIDQ("raw/2020-11-17_BoneMarrow_normalized to pmol_106cells_allvalues.xlsx",
                        colSampleName = 6)
outListHEK_HL60 <- readMetIDQ("raw/2020-11-17_HEK_HL_normalized to pmol_106cells_allvalues.xlsx",
                              colSampleName = 6)
outListHEK <- list(outListHEK_HL60[[1]][c(1:27, 29, 30),], # remove HEK 28
                   outListHEK_HL60[[2]][c(1:27, 29, 30),],
                   outListHEK_HL60[[3]])
outListHL60 <- list(outListHEK_HL60[[1]][31:60,],
                    outListHEK_HL60[[2]][31:60,],
                    outListHEK_HL60[[3]])
outListLiver <- readMetIDQ("raw/2020-11-17_Liver_normalized to pmol_mgTissue_allvalues.xlsx",
                           colSampleName = 6)

## Method names ##
method_vector <- c("M1:100/30 IPA", "M2:100 IPA", "M3:MeOH/ACN/H2O", "M4:MeOH/ACN", "M5:EtOH/PP",
                   "M6:100/20 MeOH", "M7:75 EtOH/MTBE A", "M8:75 EtOH/MTBE B", "M9:MeOH/MTBE/MeOH",
                   "M10:MeOH/ChCl3")
names(method_vector) <- c("1", "2", "3", "4", "5", "6", "7a", "7b", "8", "9")

outListBM <- renameRownames(outListBM, c("BM"), c("Bone Marrow"), startwith = FALSE)
outListHEK <- renameRownames(outListHEK, c("Cells1"), c("HEK"), startwith = FALSE)
outListHL60 <- renameRownames(outListHL60, c("Cells2"), c("HL60"), startwith = FALSE)
outListLiver <- renameRownames(outListLiver, c("Leber"), c("Liver"), startwith = FALSE)

outListBM <- renameRownames(outListBM, names(method_vector), method_vector)
outListHEK <- renameRownames(outListHEK, names(method_vector), method_vector)
outListHL60 <- renameRownames(outListHL60, names(method_vector), method_vector)
outListLiver <- renameRownames(outListLiver, names(method_vector), method_vector)

## Metabolites in Kit ##
class_anno_row_df <- data.frame(Metabolite = outListBM[[3]], 
                          Class = names(outListBM[[3]]))
anno_row_df <- read.csv("biocrates_500_names.csv")
anno_row_df$Class <- class_anno_row_df$Class[match(class_anno_row_df$Metabolite, anno_row_df$Metabolite)]


## Colors ##
col_vector <- c("#a6cee3", "#1f78b4", "#ffff99", "#b2df8a", "#33a02c", "#fb9a99",
                "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#8dd3c7",
                "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69",
                "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f")
names(col_vector) <- unique(anno_row_df$Class)

tissues <- c("Liver", "HEK", "HL60", "Bone Marrow")
col_vector_tissue <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072")
names(col_vector_tissue) <- tissues

col_vector_methods <- c("#2983B1", "#B9DBF4", "#B56478", "#F0D8E9", "#F0DCA9",
                        "#DDD9F3", "#17506C", "#278179", "#C0F0EA", "#F3C2A2")
names(col_vector_methods) <- method_vector

inputBM <- toInputFormat(outListBM)
inputHEK <- toInputFormat(outListHEK)
inputHL60 <- toInputFormat(outListHL60)
inputLiver <- toInputFormat(outListLiver)

input_df <- bind_rows(inputBM[[1]],
                      inputHEK[[1]],
                      inputHL60[[1]],
                      inputLiver[[1]])

raw_df <- bind_rows(outListBM[[1]],
                    outListHEK[[1]],
                    outListHL60[[1]],
                    outListLiver[[1]])

flag_df <- bind_rows(outListBM[[2]],
                     outListHEK[[2]],
                     outListHL60[[2]],
                     outListLiver[[2]])

transformed_df <- bind_rows(inputBM[[2]],
                            inputHEK[[2]],
                            inputHL60[[2]],
                            inputLiver[[2]])

# Remove M prefix
names(col_vector_methods) <- gsub("^M", "", names(col_vector_methods))
input_df$Methods <- gsub("^M", "", input_df$Methods)
input_df$Methods <- factor(input_df$Methods, levels = names(col_vector_methods))
rownames(raw_df) <- gsub("^M", "", rownames(raw_df))
rownames(flag_df) <- gsub("^M", "", rownames(flag_df))

# Save
save(anno_row_df, col_vector, col_vector_tissue, col_vector_methods,
     input_df, raw_df, flag_df,
     file = file.path(output.dir, "input.RData"))
