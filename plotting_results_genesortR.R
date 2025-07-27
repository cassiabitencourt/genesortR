# Loading R packages 
library(ape)
library(tidyverse)
library(tidyr)
library(dplyr)
library(vegan)
library(labdsv)
library(purrr)

# Data preparing
# Loading the main dataset and gene annotation table
table1 <- read.csv("properties_sorted_dataset.csv")
head(table1)
table2 <- read.csv("Table_names_edited.csv")
head(table2)

# Merge both tables on gene names
merged_table <- left_join(table1, table2, by = c("A353.gene" = "a353.genes"))

# Remove rows with missing data after the join
merged_table <- merged_table %>% drop_na()

# Save the cleaned merged table
write.csv(merged_table, "Results/Table_allgenes.csv", row.names = FALSE)

# PCA plotting merged data
# data preparing
gs_all <- read.csv("Results/Table_allgenes.csv", header = TRUE)
head(gs_all)
class(gs_all)

# plotting basic graph to visualise the distribution of axes
ggplot(gs_all, aes(x = PC1, y = PC2, color = categories)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Noise genes" = "darkred", "Signal-bearing genes" = "lightblue")) +
  labs(title = "phylogenetic properties of genesortR",
       color = "pythia") +
  theme_classic() 

ggsave("Rplot_stats_PCA.jpg", dpi = 600, width = 12, height = 6)

# Phylogenetic properties (genesortR) plotting
# Create a long-format dataframe with phylogenetic properties
gs_long <- gs_all %>%
  select(categories, `Root.tip.variance`:`RF.similarity`) %>%
  pivot_longer(
    cols = -categories,
    names_to = "property",
    values_to = "value"
  )

# Run ANOVA for each property
stat_results <- gs_long %>%
  group_by(property) %>%
  group_split() %>%
  map_df(~ {
    fit <- aov(value ~ categories, data = .x)
    p_val <- summary(fit)[[1]][["Pr(>F)"]][1]
    tibble(
      property = unique(.x$property),
      p_value = p_val,
      significance = case_when(
        p_val < 0.001 ~ "***",
        p_val < 0.01 ~ "**",
        p_val < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  })

# Get y-position for significance labels
label_positions <- gs_long %>%
  group_by(property) %>%
  summarise(y_pos = max(value, na.rm = TRUE) + 0.1, .groups = "drop")

# Combine with stat results and plot it
stat_results <- left_join(stat_results, label_positions, by = "property")

ggplot(gs_long, aes(x = categories, y = value, color = categories)) +
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.5) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  facet_wrap(~property, scales = "free") +
  theme_classic() +
  scale_color_manual(values = c("Noise genes" = "darkred", "Signal-bearing genes" = "lightblue")) +
  labs(color = "pythia", x = NULL, y = NULL) +
  geom_text(data = stat_results, aes(x = 1.5, y = y_pos, label = significance), inherit.aes = FALSE) +
  theme(legend.position = "none")

ggsave("Rplot_stats_ANOVA.jpg", dpi = 600, width = 14, height = 8)
