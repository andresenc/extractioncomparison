library(ggplot2)
library(cowplot)
library(dplyr)
library(magrittr)

# Load input
load(file.path("input.RData"))

input_df$Methods <- gsub(":", "\n", input_df$Methods)
input_df$Methods <- factor(input_df$Methods, levels = unique(input_df$Methods))

cv_df <- filter(input_df, LOD == 1)
unfiltered_df <- input_df

col.vector <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072")
names(col.vector) <- as.character(unique(cv_df$Tissue))

ggplot(unfiltered_df, aes(CV)) +
  geom_histogram(aes(fill = Tissue), alpha=0.5, position = "identity") +
  facet_wrap(~ Methods, ncol = 2, strip.position = "right") +
  # ggtitle("All metabolites") +
  theme_bw() +
  theme(legend.position = "none", strip.text.y = element_text(angle = 270, size = 6)) +
  scale_fill_manual(values = col.vector) +
  ylab("")

ggplot(cv_df, aes(CV)) +
  geom_histogram(aes(fill = Tissue), alpha=0.5, position = "identity") +
  facet_wrap(~ Methods, ncol = 2, strip.position = "right") +
  # ggtitle("LOD filtered") +
  theme_bw() +
  theme(legend.position = "none", strip.text.y = element_text(angle = 270, size = 6)) +
  scale_fill_manual(values = col.vector) +
  ylab("")

