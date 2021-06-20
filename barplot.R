library(ggplot2)
library(cowplot)
library(dplyr)
library(magrittr)

# Load input
load(file.path("input.RData"))

# Remodel data
filtered_input_df <- dplyr::select(input_df, c(Class, Methods, Tissue, LOD)) %>%
  group_by(Class, Tissue, Methods) %>%
  summarize(Number = sum(LOD)) %>% 
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

legend_plot <- ggplot(df_list$`Bone Marrow`, aes(x = Methods, y = Number, fill = Class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col_vector) +
  theme_light() +
  ggtitle(unique(df_list$`Bone Marrow`$Tissue))

legend_graph <- get_legend(legend_plot)

# Combine plots
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

# Return statistics
# stat_list <- lapply(df_list, FUN = function(x){
#   out <- group_by(x, Methods) %>% summarize(Values = sum(Number))
#   out_df <- data.frame(Median = median(out$Values),
#                        Min = min(out$Values),
#                        Max = max(out$Values))
#   return(out_df)
# })
