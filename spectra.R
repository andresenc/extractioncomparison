library(ggplot2)
library(cowplot)
library(scales)

# Load input
load(file.path("input.RData"))

cv_df <- input_df
cv_df$Mean <- ifelse(cv_df$LOD == 1, cv_df$Mean, 0)

# Plot liver samples without method one (Figure 3)
liver_temp <- filter(cv_df, Tissue == "Liver" & Method != "1")
ggplot(liver_temp, aes(x = Metabolite, y = Mean, fill = Class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col_vector) +
  facet_wrap(~ Methods, ncol = 1, strip.position = "right") +
  theme_bw() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  ylab("Mean Concentration [pmol/mg]") +
  theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "cm"),
        panel.grid.major= element_line(color = "white"), 
        panel.grid.minor.y = element_blank(),
        strip.text.y = element_text(angle = 270, size = 6)) 

# Plot HEK samples
hek_temp <- filter(cv_df, Tissue == "HEK")
ggplot(hek_temp, aes(x = Metabolite, y = Mean, fill = Class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col_vector) +
  facet_wrap(~ Methods, ncol = 1, strip.position = "right") +
  theme_bw() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  ylab("Mean Concentration [pmol/10^6 cells]") +
  theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "cm"),
        panel.grid.major= element_line(color = "white"), 
        panel.grid.minor.y = element_blank(),
        strip.text.y = element_text(angle = 270, size = 6)) 

# Plot HL60 samples
hl_temp <- filter(cv_df, Tissue == "HL60")
ggplot(hl_temp, aes(x = Metabolite, y = Mean, fill = Class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col_vector) +
  facet_wrap(~ Methods, ncol = 1, strip.position = "right") +
  theme_bw() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  ylab("Mean Concentration [pmol/10^6 cells]") +
  theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "cm"),
        panel.grid.major= element_line(color = "white"), 
        panel.grid.minor.y = element_blank(),
        strip.text.y = element_text(angle = 270, size = 6)) 

# Plot bone marrow samples
bm_temp <- filter(cv_df, Tissue == "Bone Marrow")
ggplot(bm_temp, aes(x = Metabolite, y = Mean, fill = Class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col_vector) +
  facet_wrap(~ Methods, ncol = 1, strip.position = "right") +
  theme_bw() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  ylab("Mean Concentration [pmol/10^6 cells]") +
  theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "cm"),
        panel.grid.major= element_line(color = "white"), 
        panel.grid.minor.y = element_blank(),
        strip.text.y = element_text(angle = 270, size = 6)) 

# Plot method 1 across all tissues
method_temp <- filter(cv_df, Method == "1")
method_temp$Tissue <- factor(method_temp$Tissue, levels = c("Liver", "Bone Marrow", "HEK", "HL60"))
ggplot(method_temp, aes(x = Metabolite, y = Mean, fill = Class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col_vector) +
  facet_wrap(~ Tissue, ncol = 1, strip.position = "top") +
  theme_bw() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  ylab("Mean Concentration [pmol/mg] or [pmol/10^6 cells]") +
  theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "cm"),
        panel.grid.major= element_line(color = "white"), 
        panel.grid.minor.y = element_blank(),
        strip.text.y = element_text(angle = 360, size = 6), legend.position = "none") 
