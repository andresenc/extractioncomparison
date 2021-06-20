library(ggplot2)
library(cowplot)
library(dplyr)
library(magrittr)

# Load input
load(file.path("input.RData"))

# Remodel data
pie_df <- data.frame(Weight = 1, Class = anno_row_df$Class) %>%
  group_by(Class) %>%
  dplyr::summarise(Weight = n()) %>% 
  arrange(desc(Class)) %>% 
  mutate(text_temp = cumsum(Weight),  
         text_y = cumsum(Weight) - Weight/2) %>%
  arrange(Class) %>% 
  data.frame()

# Plot
ggplot(pie_df, aes(x="", y=Weight, fill=Class)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  scale_fill_manual(values = col_vector) + theme_void() + 
  theme(axis.text.x=element_blank()) +
  geom_label_repel(aes(label = Weight, y = text_y), size=4, 
                   show.legend = F, nudge_x = 0.3)



