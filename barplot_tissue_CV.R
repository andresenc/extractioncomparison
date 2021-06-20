library(magrittr)
library(dplyr)
library(tidyr)
library(tibble)

# Load input
load(file.path("input.RData"))

# Filter for LOD and calculate statistics
ablod_df <- filter(input_df, LOD == 1)
statistics_df <- group_by(ablod_df, Tissue, Methods) %>%
  summarise(Number_Metabolites = n(),
            Median_CV = median(CV),
            MAD_CV = mad(CV)) %>%
  data.frame()

# Filter for CV and calculate statistics
statistics_30_df <- filter(ablod_df, CV < 0.3) %>%
  group_by(Tissue, Methods) %>%
  summarise(Sum_CV_30 = n(),
            Median_CV_30 = median(CV),
            MAD_CV_30 = mad(CV)) %>%
  data.frame()

statistics_df$Sum_CV_30 <- statistics_30_df$Sum_CV_30
statistics_df$Median_CV_30 <- statistics_30_df$Median_CV_30
statistics_df$MAD_CV_30 <- statistics_30_df$MAD_CV_30
summary_df <- statistics_df

# Median with error bars
gg_mean <- ggplot(summary_df, aes(Methods, Number_Metabolites)) +
  geom_bar(aes(fill = Tissue), position = "dodge", stat="identity") +
  xlab("") +
  ylab("#") +
  theme_light() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = col_vector_tissue)+
  ylim(0, 500)

gg_cv <- ggplot(summary_df, aes(Methods, Median_CV)) +
  geom_bar(aes(fill = Tissue), position = "dodge", stat="identity") +
    geom_errorbar(aes(ymin=Median_CV, ymax=Median_CV+MAD_CV, fill = Tissue), width=.2,
                  position=position_dodge(.9), col = "black") +
  xlab("") +
  ylab("Median CV + MAD") +
  theme_light() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = col_vector_tissue)

gg_cv_30 <- ggplot(summary_df, aes(Methods, Sum_CV_30)) +
  geom_bar(aes(fill = Tissue), position = "dodge", stat="identity") +
  xlab("") +
  ylab("# (CV < 30%)") +
  theme_light() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = col_vector_tissue)+
  ylim(0, 500)

legend_plot <- ggplot(summary_df, aes(Methods, Median_CV_30)) +
  geom_bar(aes(fill = Tissue), position = "dodge", stat="identity") +
  scale_fill_manual(values = col_vector_tissue) +
  theme(legend.position = "bottom")
legend <- get_legend(legend_plot)

# Plot
plot_grid(gg_mean, 
          gg_cv_30,
          gg_cv,
          legend,
          ncol = 2, rel_heights = c(1,1))

