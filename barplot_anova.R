library(ggplot2)
library(cowplot)
library(dplyr)
library(magrittr)

## Set directories
input.dir <- file.path("/icgc/dkfzlsdf/analysis/OE0285_projects/biocrates500/extraction_evaluation/paper_figures")
output.dir <- "/icgc/dkfzlsdf/analysis/OE0285_projects/biocrates500/extraction_evaluation/paper_figures"
load(file.path(input.dir, "input.RData"))

# Remodel data
filtered_input_df <- filter(input_df, LOD != 0, !is.na(Group)) %>%
  droplevels() %>%
  group_by(Methods, Class, Tissue) %>%
  dplyr::summarise(Number = length(grep("A", Group))) %>%
  data.frame()

# Split across tissue types
df_list <- split(filtered_input_df, f = filtered_input_df$Tissue)

# Plot
plot_list <- lapply(df_list, function(x){
  output <- ggplot(x, aes(x = Methods, y = Number, fill = Class)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = col_vector) +
    theme_light() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
    ggtitle(unique(x$Tissue)) +
    ylim(c(0, 520)) +
    xlab("")
  return(output)
})

legend_plot <- ggplot(df_list$Liver, aes(x = Methods, y = Number, fill = Class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col_vector) +
  theme_light() +
  ggtitle(unique(df_list$Liver$Tissue))

legend_graph <- get_legend(legend_plot)

# Combine plots
pdf(file = file.path(output.dir, "Combined_barplot_anova.pdf"), width = 9, height = 7)
plot_grid(plot_grid(plot_list[[1]],
                    plot_list[[2]],
                    plot_list[[3]],
                    plot_list[[4]],
                    ncol = 2,
                    scale = 1),
          plot_grid(legend_graph, scale = 0.3), 
          ncol = 2, 
          rel_widths = c(1, 0.7),
          rel_heights = c(1, 0.5))
dev.off()


