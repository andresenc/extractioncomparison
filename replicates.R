library(magrittr)
library(dplyr)
library(tidyr)
library(tibble)

# Load input
load(file.path("input.RData"))

# Filter for LOD
ablod_df <- filter(input_df, LOD == 1)
statistics_df <- group_by(ablod_df, Tissue, Methods) %>%
  summarise(Number_Metabolites = n(),
            Median_CV = median(CV),
            MAD_CV = mad(CV)) %>%
  data.frame()

# Calculate sum per replicate
sum_df <- group_by(ablod_df, Tissue, Methods) %>%
  summarise(Sample1 = max(SizeA),
            Sample2 = max(SizeB),
            Sample3 = max(SizeC)) %>%
  data.frame()

statistics_df$Sample1 <- sum_df$Sample1
statistics_df$Sample2 <- sum_df$Sample2
statistics_df$Sample3 <- sum_df$Sample3

# Concentrations over replicates
summary_rep_df <- dplyr::select(statistics_df, c(Methods, Sample1, Sample2, Sample3, Tissue)) %>%
  gather(Replicate, Sum, -c(Methods, Tissue))

# Plot
ggplot(summary_rep_df, aes(x = Tissue, y = Sum, fill = Replicate)) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_wrap(~ Methods, ncol = 2) +
  theme_bw() +
  ylab("Sum of concentrations") +
  scale_fill_manual(values = c("grey40", "grey60", "grey80"))

# Statistics
# stat_df <- summary_rep_df %>%
#   group_by(Tissue, Methods) %>%
#   summarise(Median = median(Sum), MAD = mad(Sum))
