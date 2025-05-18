library(dplyr)
library(data.table)
library(readr)
library(patchwork)

#DESCRIPTION
#The following code was all of the data analysis for my dissertation. The data was from short_hydro_bonds.py

#Full dataset analysis
mega_combined <- read_csv("~/Documents/Protein modeling/Phase 2/mega_combined.csv")

#Count number of rows with sec_type = "both"
mega_combined %>%
  count(sec_type) %>%
  filter(sec_type == "both")

#Create coloum Helix index
mega_combined$helix_num <- NA
ix <- grepl("^HELX_P", mega_combined$`PDB secondary structure`)
mega_combined$helix_num[ix] <- sub("^HELX_P", "", mega_combined$`PDB secondary structure`[ix])

#data frame of Helix lengths and distributions
setDT(mega_combined)
helix_lengths <- read.csv("helix_lengths.csv", stringsAsFactors=FALSE)

helix_lengths <- mega_combined[!is.na(helix_num),
   .(helix_length = .N),
   by = .(pdb_id, model, chain, helix_num)
]

#write.csv(helix_lengths, "helix_lengths.csv", row.names=FALSE)

#Max helix length
max(helix_lengths$helix_length)

#Helix lengths counts and distributions
counts <- helix_lengths %>%
  count(helix_length, sort = TRUE) %>%
  mutate(percentage = round(100 * n / sum(n), 2))

helix_lengths2 <- helix_lengths[
  , .(
    pdb_id,
    model,
    chain,
    helix_num,
    helix_length = fifelse(helix_length > 100, 100, helix_length)
  )
]

#Graph of helix length distributions
hist(pmin(helix_lengths$helix_length,100),breaks=0:101,right=FALSE,xaxt="n",yaxt="n",xlab="Helix length",main="")
axis(1,at=seq(0,100,10),labels=c(as.character(seq(0,90,10)),"100+")); axis(2,at=seq(0,400000,100000),labels=c("0","100000","200000","300000","400000"))

#---------------------------------------------------------------------------------------------------
# Sample 5%, 10%, 20% of unique proteins (based on pdb_id)
sampled_protein_5 <- sample(unique_proteins, size = round(0.05 * length(unique_proteins)))
subset_5 <- new_mega_combined[new_mega_combined$pdb_id %in% sampled_protein_5, ]
write_csv(subset_5, "subset_5.csv")
rm(subset_5)
gc()

sampled_protein_10 <- sample(unique_proteins, size = round(0.10 * length(unique_proteins)))
subset_10 <- new_mega_combined[new_mega_combined$pdb_id %in% sampled_protein_10, ]
write_csv(subset_10, "subset_10.csv")
rm(subset_10)
gc()

sampled_protein_20 <- sample(unique_proteins, size = round(0.20 * length(unique_proteins)))
subset_20 <- new_mega_combined[new_mega_combined$pdb_id %in% sampled_protein_20, ]
write_csv(subset_20, "subset_20.csv")
rm(subset_20)
gc()

# Sample 5%, 10%, 20% random rows
sample_5_row <- new_mega_combined[sample(nrow(new_mega_combined), size = round(0.05 * nrow(new_mega_combined))), ]
write_csv(sample_5_row, "subset_5_row.csv")
rm(sample_5_row)
gc()

sample_10_row <- new_mega_combined[sample(nrow(new_mega_combined), size = round(0.10 * nrow(new_mega_combined))), ]
write_csv(sample_10_row, "subset_10_row.csv")
rm(sample_10_row)
gc()

sample_20_row <- new_mega_combined[sample(nrow(new_mega_combined), size = round(0.20 * nrow(new_mega_combined))), ]
write_csv(sample_20_row, "subset_20_row.csv")
rm(sample_20_row)
gc()

# Read protein-based subsets
subset_5 <- read_csv("subset_5.csv")
subset_10 <- read_csv("subset_10.csv")
subset_20 <- read_csv("subset_20.csv")

# Read row-based subsets
subset_5_row <- read_csv("subset_5_row.csv")
subset_10_row <- read_csv("subset_10_row.csv")
subset_20_row <- read_csv("subset_20_row.csv")

#Processed version of subset_20 with helix data
subset_20 <- read_csv("subset_20_plushelixinfo.csv")

#Add step -3 and -4 -----------------------------------------

subset_20 <- subset_20 %>%
  group_by(pdb_id, model, chain) %>%      
  arrange(index, .by_group = TRUE) %>%   
  mutate(difference_minus_3 = lag(difference_S3, n = 3)) %>%
  ungroup()

subset_20 <- subset_20 %>%
  group_by(pdb_id, model, chain) %>%      
  arrange(index, .by_group = TRUE) %>%    
  mutate(difference_minus_4 = lag(difference_S4, n = 4)) %>%
  ungroup()

subset_20_helix <- subset(subset_20, sec_type == "helix")

#Proportion of helix start end both 
subset_20 %>%
  filter(sec_type == "helix") %>%
  count(helix_region) %>%
  mutate(proportion = n / sum(n))

#Overall SS distribution
subset_20 %>%
  count(sec_type) %>%
  mutate(proportion = n / nrow(subset_20)) %>%
  print()

#End of helices reversed distribution
helix_tail_summary <- subset_20 %>%
  filter(sec_type == "helix", helix_index_reversed <= 20) %>%
  group_by(helix_index_reversed) %>%
  summarise(mean_distance_S4 = mean(distance_S4, na.rm = TRUE), .groups = "drop") %>%
  arrange(helix_index_reversed)
labels <- c("End", paste0("-", 1:19))

helix_tail_summary_s3 <- subset_20 %>%
  filter(sec_type == "helix", helix_index_reversed <= 20) %>%
  group_by(helix_index_reversed) %>%
  summarise(mean_distance_S3 = mean(distance_S3, na.rm = TRUE), .groups = "drop") %>%
  arrange(helix_index_reversed)

#Graphs of helix end distribution
h4_t <- ggplot(helix_tail_summary, aes(x = factor(helix_index_reversed), y = mean_distance_S4)) +
  geom_col(fill = "salmon") +
  scale_x_discrete(
    labels = rev(c("End", paste0("-", 1:19))),
    limits = rev(as.character(1:20))
  ) +
  ylim(0,7.5) +
  labs(x = "Position from Helix End",y = "Mean Distance |O N+4| (Å)") +
  theme_minimal() +
  theme(text = element_text(size = 12))

h3_t <- ggplot(helix_tail_summary_s3, aes(x = factor(helix_index_reversed), y = mean_distance_S3)) +
  geom_col(fill = "salmon") +
  scale_x_discrete(
    labels = rev(c("End", paste0("-", 1:19))),
    limits = rev(as.character(1:20))
  ) +
  ylim(0,7.5) +
  labs(x = "Position from Helix End", y = "Mean Distance |O N+3| (Å)") +
  theme_minimal() +
  theme(text = element_text(size = 12))

(h3_t | h4_t)

# Histogram for Hbond Distance ---------------------------------------------------------------------
h1 <- ggplot(subset_20, aes(x = distance_S3, fill = sec_type)) +
  geom_histogram(binwidth = 0.1, alpha = 0.4, position = "identity") +
  labs(title = "Histogram for Distance |O N+3|", x = "Distance |O N+3|", fill = "Secondary Type") +
  ylim(NA, 1000000) +
  xlim(NA, 10)
h_d1

h2 <- ggplot(subset_20, aes(x = distance_S4, fill = sec_type)) +
  geom_histogram(binwidth = 0.1, alpha = 0.4, position = "identity") +
  labs(title = "Histogram for Distance |O N+4|", x = "Distance |O N+4|", fill = "Secondary Type") +
  ylim(NA, 1000000) +
  xlim(NA, 13)
h2

h3 <- ggplot(subset_20, aes(x = distance_S5, fill = sec_type)) +
  geom_histogram(binwidth = 0.1, alpha = 0.4, position = "identity") +
  labs(title = "Histogram for Distance |O N+5|", x = "Distance |O N+5|", fill = "Secondary Type") +
  ylim(NA, 1000000) +
  xlim(NA, 17)
h3

# Histogram for Hbond Differences ------------------------------------------------------------------
h4 <- ggplot(subset_20, aes(x = difference_S3, fill = sec_type)) +
  geom_histogram(binwidth = 0.02, alpha = 0.4, position = "identity") +
  labs(title = "Histogram for Difference |C N+3| - |O N+3|", x = "Difference |C N+3| - |O N+3|", fill = "Secondary Type") +
  ylim(NA, 400000) +
  xlim(-1.5, 1.5)
h4

h5 <- ggplot(subset_20, aes(x = difference_S4, fill = sec_type)) +
  geom_histogram(binwidth = 0.02, alpha = 0.4, position = "identity") +
  labs(title = "Histogram for Difference |C N+4| - |O N+4|", x = "Difference |C N+4| - |O N+4|", fill = "Secondary Type") +
  ylim(NA, 400000) +
  xlim(-1.5, 1.5)
h5

h6 <- ggplot(subset_20, aes(x = difference_S5, fill = sec_type)) +
  geom_histogram(binwidth = 0.02, alpha = 0.4, position = "identity") +
  labs(title = "Histogram for Difference |C N+5| - |O N+5|", x = "Difference |C N+5| - |O N+5|", fill = "Secondary Type") +
  ylim(NA, 400000) +
  xlim(-1.5, 1.5)
h6

#Subset chains -------------------------------------------------------------------------------------

subset_chain_df1 <- subset(subset_20, pdb_id == "108M" & model == 1 & chain == "A")
subset_chain_df2 <- subset(subset_20, pdb_id == "7W7B" & model == 1 & chain == "C")

subset_chain_df1_helix <- subset_chain_df1[subset_chain_df1$sec_type == "helix", ]

# Subset chain distance S3 and S4
subset_chain_df1_helix <- subset_chain_df1_helix[order(subset_chain_df1_helix$index), ]
group <- cumsum(c(0, diff(subset_chain_df1_helix$index) != 1))
plot(subset_chain_df1_helix$index, subset_chain_df1_helix$distance_S3, type = "n",
     ylim = range(c(subset_chain_df1_helix$distance_S3, subset_chain_df1_helix$distance_S4), na.rm = TRUE),
     xlab = "Residue Index", ylab = "Distance |O N+x| (Å)")
for (g in unique(group)) {
  segment <- subset_chain_df1_helix[group == g, ]
  lines(segment$index, segment$distance_S3, col = "black", lwd = 0.5)
}
points(subset_chain_df1_helix$index, subset_chain_df1_helix$distance_S3, col = "red", pch = 16, cex = 0.6)
for (g in unique(group)) {
  segment <- subset_chain_df1_helix[group == g, ]
  lines(segment$index, segment$distance_S4, col = "gray40", lwd = 0.5)
}
points(subset_chain_df1_helix$index, subset_chain_df1_helix$distance_S4, col = "blue", pch = 16, cex = 0.6)

# Subset chain step negative distance 
subset_chain_df1_helix <- subset_chain_df1_helix[order(subset_chain_df1_helix$index), ]
group <- cumsum(c(0, diff(subset_chain_df1_helix$index) != 1))
plot(subset_chain_df1_helix$index, subset_chain_df1_helix$distance_step_minus_3, type = "n",
     ylim = range(c(subset_chain_df1_helix$distance_step_minus_3, subset_chain_df1_helix$distance_step_minus_4), na.rm = TRUE),
     xlab = "Residue Index", ylab = "Distance |O N-x| (Å)")
for (g in unique(group)) {
  segment <- subset_chain_df1_helix[group == g, ]
  lines(segment$index, segment$distance_step_minus_3, col = "black", lwd = 0.5)
}
points(subset_chain_df1_helix$index, subset_chain_df1_helix$distance_step_minus_3, col = "red", pch = 16, cex = 0.6)
for (g in unique(group)) {
  segment <- subset_chain_df1_helix[group == g, ]
  lines(segment$index, segment$distance_step_minus_4, col = "gray40", lwd = 0.5)
}
points(subset_chain_df1_helix$index, subset_chain_df1_helix$distance_step_minus_4, col = "blue", pch = 16, cex = 0.6)


# Subset chain distance S3 and S4
subset_chain_df1_helix <- subset_chain_df1_helix[order(subset_chain_df1_helix$index), ]
group <- cumsum(c(0, diff(subset_chain_df1_helix$index) != 1))
plot(subset_chain_df1_helix$index, subset_chain_df1_helix$difference_S3, type = "n",
     ylim = range(c(subset_chain_df1_helix$difference_S3, subset_chain_df1_helix$difference_S4), na.rm = TRUE),
     xlab = "Residue Index", ylab = "Difference |C N+x| - |O N+x| (Å)")
for (g in unique(group)) {
  segment <- subset_chain_df1_helix[group == g, ]
  lines(segment$index, segment$difference_S3, col = "black", lwd = 0.5)
}
points(subset_chain_df1_helix$index, subset_chain_df1_helix$difference_S3, col = "red", pch = 16, cex = 0.6)
for (g in unique(group)) {
  segment <- subset_chain_df1_helix[group == g, ]
  lines(segment$index, segment$difference_S4, col = "gray40", lwd = 0.5)
}
points(subset_chain_df1_helix$index, subset_chain_df1_helix$difference_S4, col = "blue", pch = 16, cex = 0.6)


#Scatter plot --------------------------------------------------------------------------------------
s3 <- ggplot(subset_20, aes(x = distance_S3, y = difference_S3, color = sec_type)) +
  geom_point(alpha = 0.01, size = 0.0005) +
  labs(x = "Distance |O N+3| (Å)", y = "Difference |C N+3| - |O N+3| (Å)", color = "Secondary Structure") +
  xlim(2, 9) +
  ylim(-1.25, 1.25) +
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 1, alpha = 1)))

s3_h <- ggplot(subset(subset_20, sec_type == "helix"), aes(x = distance_S3, y = difference_S3)) +
  geom_point(alpha = 0.01, size = 0.0005, color = "salmon") +
  labs(x = "Distance |O N+3| (Å)", y = "Difference |C N+3| - |O N+3| (Å)") +
  xlim(2, 9) +
  ylim(-1.25, 1.25) +
  theme_minimal() +
  theme(text = element_text(size = 18))

s3_s <- ggplot(subset(subset_20, sec_type == "sheet"), aes(x = distance_S3, y = difference_S3)) +
  geom_point(alpha = 0.01, size = 0.0005, color = "steelblue") +
  labs(x = "Distance |O N+3| (Å)", y = "Difference |C N+3| - |O N+3| (Å)") +
  xlim(2, 9) +
  ylim(-1.25, 1.25) +
  theme_minimal() +
  theme(text = element_text(size = 12))

s4 <- ggplot(subset_20, aes(x = distance_S4, y = difference_S4, color = sec_type)) +
  geom_point(alpha = 0.01, size = 0.0005) +
  labs(x = "Distance |O N+4| (Å)", y = "Difference |C N+4| - |O N+4| (Å)", color = "Secondary Structure") +
  xlim(2, 13) +
  ylim(-1.6, 1.5) +
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 1, alpha = 1)))

s4_h <- ggplot(subset(subset_20, sec_type == "helix"), aes(x = distance_S4, y = difference_S4)) +
  geom_point(alpha = 0.01, size = 0.0005, color = "salmon") +
  labs(x = "Distance |O N+4| (Å)", y = "Difference |C N+4| - |O N+4| (Å)") +
  xlim(2, 13) +
  ylim(-1.6, 1.5) +
  theme_minimal() +
  theme(text = element_text(size = 18))

s4_s <- ggplot(subset(subset_20, sec_type == "sheet"), aes(x = distance_S4, y = difference_S4)) +
  geom_point(alpha = 0.01, size = 0.0005, color = "steelblue") +
  labs(x = "Distance |O N+4| (Å)", y = "Difference |C N+4| - |O N+4| (Å)") +
  xlim(2, 13) +
  ylim(-1.6, 1.5) +
  theme_minimal() +
  theme(text = element_text(size = 12))

s5 <- ggplot(subset_20, aes(x = distance_S5, y = difference_S5, color = sec_type)) +
  geom_point(alpha = 0.01, size = 0.0005) +
  labs(x = "Distance |O N+5| (Å)", y = "Difference |C N+5| - |O N+5| (Å)", color = "Secondary Structure") +
  xlim(2, 17) +
  ylim(-1.5, 1.5) +
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 1, alpha = 1)))

s5_h <- ggplot(subset(subset_20, sec_type == "helix"), aes(x = distance_S5, y = difference_S5)) +
  geom_point(alpha = 0.01, size = 0.0005, color = "salmon") +
  labs(x = "Distance |O N+5| (Å)", y = "Difference |C N+5| - |O N+5| (Å)") +
  xlim(2, 17) +
  ylim(-1.5, 1.5) +
  theme_minimal() +
  theme(text = element_text(size = 18))

s5_s <- ggplot(subset(subset_20, sec_type == "sheet"), aes(x = distance_S5, y = difference_S5)) +
  geom_point(alpha = 0.01, size = 0.0005, color = "steelblue") +
  labs(x = "Distance |O N+5| (Å)", y = "Difference |C N+5| - |O N+5| (Å)") +
  xlim(2, 17) +
  ylim(-1.5, 1.5) +
  theme_minimal() +
  theme(text = element_text(size = 12))

rm(s3, s3_h, s3_s, s4, s4_h, s4_s)

#Heatmas 
s3_h <- ggplot(subset_20_helix, aes(x = distance_S3, y = difference_S3)) +
  geom_bin2d(bins = 200) +
  scale_fill_gradientn(colors = rev(rainbow(7)), trans = "log10", name = "Log Count") +
  labs(x = "Distance |O N+3| (Å)", y = "Difference |C N+3| - |O N+3| (Å)") +
  xlim(2, 9) +
  ylim(-1.25, 1.25) +
  theme_minimal() +
  theme(text = element_text(size = 12))
s3_h

s4_h <- ggplot(subset_20_helix, aes(x = distance_S4, y = difference_S4)) +
  geom_bin2d(bins = 200) +
  scale_fill_gradientn(colors = rev(rainbow(7)), trans = "log10", name = "Log Count") +
  labs(x = "Distance |O N+4| (Å)", y = "Difference |C N+4| - |O N+4| (Å)") +
  xlim(2, 13) +
  ylim(-1.6, 1.5) +
  theme_minimal() +
  theme(text = element_text(size = 12))
s4_h

#Set fill amounts for comparison of heatmaps
shared_fill <- scale_fill_gradientn(
  colors = rev(c("black","red","yellow","lightgreen","blue")),
  trans = "log10",
  limits = c(1, 2*(10^5)),  # Raw count scale; log10 applied internally
  name = "Log Count"
)

s3_h_fixed <- ggplot(subset_20_helix, aes(x = distance_S3, y = difference_S3)) +
  geom_bin2d(bins = 200) +
  shared_fill +
  labs(x = "Distance |O N+3| (Å)", y = "Difference |C N+3| - |O N+3| (Å)") +
  xlim(2, 9.5) +
  ylim(-1.3, 1.3) +
  theme_minimal() +
  theme(text = element_text(size = 12))
#s3_h_fixed

s4_h_fixed <- ggplot(subset_20_helix, aes(x = distance_S4, y = difference_S4)) +
  geom_bin2d(bins = 200) +
  shared_fill +
  labs(x = "Distance |O N+4| (Å)", y = "Difference |C N+4| - |O N+4| (Å)") +
  xlim(1.5, 13) +
  ylim(-1.7, 1.5) +
  theme_minimal() +
  theme(text = element_text(size = 12))
#s4_h_fixed

(s3_h_fixed | s4_h_fixed) + plot_layout(guides = "collect")

s3_h_fixed <- ggplot(subset_20_helix, aes(x = distance_S3, y = difference_S3)) +
  geom_bin2d(bins = 200) +
  shared_fill +
  labs(x = "Distance |O N+3| (Å)", y = "Difference |C N+3| - |O N+3| (Å)") +
  xlim(2, 9.5) +
  ylim(-1.3, 1.3) +
  theme_minimal() +
  theme(text = element_text(size = 12))
#s3_h_fixed

s4_h_fixed <- ggplot(subset_20_helix, aes(x = distance_S4, y = difference_S4)) +
  geom_bin2d(bins = 200) +
  shared_fill +
  labs(x = "Distance |O N+4| (Å)", y = "Difference |C N+4| - |O N+4| (Å)") +
  xlim(1.5, 13) +
  ylim(-1.7, 1.5) +
  theme_minimal() +
  theme(text = element_text(size = 12))
#s4_h_fixed

(s3_h_fixed | s4_h_fixed) + plot_layout(guides = "collect")

# Step ± grpahs -----------------------------------------------------------------------------------

#3
ggplot(subset_20_helix, aes(x = distance_S3, y = distance_step_minus_3, color = helix_region)) +
  geom_point(alpha = 0.01, size = 0.0005) +
  labs(x = "Distance |O N+3| (Å)", y = "Distance |O N-3| (Å)", color = "Helix Region") +
  xlim(1.8, 9) + ylim(1.8, 9) +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
  theme(text = element_text(size = 20))

#4
ggplot(subset_20_helix, aes(x = distance_S4, y = distance_step_minus_4, color = helix_region)) +
  geom_point(alpha = 0.01, size = 0.0005) +
  labs(x = "Distance |O N+4| (Å)", y = "Distance |O N-4| (Å)", color = "Helix Region") +
  xlim(2.5, 12.5) + ylim(2.5, 12.5) +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
  theme(text = element_text(size = 20))

region_colors <- c(
  "start" = "#C77CFF",
  "middle" = "#00BFC4",
  "end" = "#00BA38",
  "both" = "#F8766D"
)

# Function to create a plot for a given region
plot_region <- function(region_name) {
  ggplot(filter(subset_20_helix, helix_region == region_name),
         aes(x = distance_S4, y = distance_step_minus_4)) +
    geom_point(color = region_colors[region_name], alpha = 0.01, size = 0.0005) +
    labs(
      x = "Distance |O N+4| (Å)",
      y = "Distance |O N-4| (Å)"
    ) +
    xlim(2.5, 12.5) + ylim(2.5, 12.5) +
    theme_minimal() +
    theme(text = element_text(size = 12))
}

heatmap_pm4 <- function(region_name) {
  ggplot(filter(subset_20_helix, helix_region == region_name), aes(x = distance_S4, y = distance_step_minus_4)) +
    geom_bin2d(bins = 200) +
    scale_fill_gradientn(colors = rev(c("black","red","yellow","lightgreen","blue")), trans = "log10", name = "Log Count") +
    labs(x = "Distance |O N+4| (Å)", y = NULL) +
    xlim(2.5, 12.5) + ylim(2.5, 12.5) +
    theme_minimal() +
    theme(text = element_text(size = 12))
}

# Create plots
s4_start  <- plot_region("start")
s4_middle <- plot_region("middle")
s4_end    <- plot_region("end")
s4_both   <- plot_region("both")

h4_start  <- heatmap_pm4("start")
h4_middle <- heatmap_pm4("middle")
h4_end    <- heatmap_pm4("end")
h4_both   <- heatmap_pm4("both")

(s4_start | h4_start)
(s4_middle | h4_middle)
(s4_end  | h4_end )
(s4_both | h4_both)

# Step -3 scatter heatmap USES FIXED SCALE
sminus3_h <- ggplot(subset_20_helix, aes(x = distance_step_minus_3, y = difference_minus_3)) +
  geom_bin2d(bins = 200) +
  shared_fill +
  labs(x = "Distance |O N-3| (Å)", y = "Difference |C N-3| - |O N-3| (Å)") +
  xlim(2, 9) +
  ylim(-1.25, 1.25) +
  theme_minimal() +
  theme(text = element_text(size = 12))
sminus3_h

# Step -4 scatter heatmap USES FIXED SCALE 
sminus4_h <- ggplot(subset_20_helix, aes(x = distance_step_minus_4, y = difference_minus_4)) +
  geom_bin2d(bins = 200) +
  shared_fill +
  labs(x = "Distance |O N-4| (Å)", y = "Difference |C N-4| - |O N-4| (Å)") +
  xlim(2, 13) + ylim(-1.6, 1.5) +
  theme_minimal() +
  theme(text = element_text(size = 12))
sminus4_h

(sminus3_h | sminus4_h) + plot_layout(guides = "collect")

#Looking at elipse regions ----------------------------------------------------------------

#Subset 20 cluster analysis
subset_20_CA <- subset_20 %>% select(pdb_id, model, chain, index, residue, `PDB secondary structure`, sec_type, distance_S3, difference_S3, distance_S4, difference_S4, distance_step_minus_3, difference_minus_3, distance_step_minus_4, difference_minus_4, helix_num_extracted, helix_num, helix_index, helix_index_reversed)

# 1. Define the ellipse functions
ellipse_df <- function(x0, y0, a, b, theta, n = 200) {
  t  <- seq(0, 2*pi, length.out = n)
  x  <- a * cos(t)
  y  <- b * sin(t)
  xr <- x * cos(theta) - y * sin(theta) + x0
  yr <- x * sin(theta) + y * cos(theta) + y0
  data.frame(x = xr, y = yr)
}

inside_ellipse <- function(x, y, x0, y0, a, b, theta) {
  na_flag <- is.na(x) | is.na(y)
  xp      <- x - x0
  yp      <- y - y0
  x_rot   <- xp * cos(theta) + yp * sin(theta)
  y_rot   <- -xp * sin(theta) + yp * cos(theta)
  inside  <- (x_rot / a)^2 + (y_rot / b)^2 <= 1
  inside[na_flag] <- FALSE
  inside
}

# S3
{
  x0_s3    <- 3.45
  y0_s3    <- 0.45
  a_s3     <- 0.92
  b_s3     <- 0.25
  theta_s3 <- -0.48
ell_S3 <- ellipse_df(x0_s3, y0_s3, a_s3, b_s3, theta_s3)
}
s3_h_fixed +
  geom_path(data = ell_S3, aes(x, y), inherit.aes = FALSE, colour = "black", linewidth = 0.7) +
  coord_cartesian(xlim = c(2, 5), ylim = c(-0.25, 1.5))
sminus3_h + 
  geom_path(data = ell_S3, aes(x, y), inherit.aes = FALSE, colour = "black", linewidth = 0.7) +
  coord_cartesian(xlim = c(2, 5), ylim = c(-0.25, 1.5))

# S4
{
  x0_s4    <- 3.2
  y0_s4    <- 1.04
  a_s4     <- 0.85
  b_s4     <- 0.2
  theta_s4 <- -0.15
ell_S4 <- ellipse_df(x0_s4, y0_s4, a_s4, b_s4, theta_s4)
}
s4_h_fixed +
  geom_path(data = ell_S4, aes(x, y), inherit.aes = FALSE, colour = "black", linewidth = 0.7) +
  coord_cartesian(xlim = c(2, 5), ylim = c(0.7, 1.3))
sminus4_h + 
  geom_path(data = ell_S4, aes(x, y), inherit.aes = FALSE, colour = "black", linewidth = 0.7) +
  coord_cartesian(xlim = c(2, 5), ylim = c(-0.25, 1.5))

#S3 S-3, S4 S-4 side by side heatmaps with helixes
{
  s3_h_zoomed <- s3_h_fixed + 
    geom_path(data = ell_S3, aes(x, y), inherit.aes = FALSE, colour = "black", linewidth = 0.7) +
    coord_cartesian(xlim = c(2, 5), ylim = c(-0.25, 1.5)) +
    guides(fill = "none")
  
  sminus3_h_zoomed <- sminus3_h + 
    geom_path(data = ell_S3, aes(x, y), inherit.aes = FALSE, colour = "black", linewidth = 0.7) +
    coord_cartesian(xlim = c(2, 5), ylim = c(-0.25, 1.5)) +
    guides(fill = "none")
  
  s4_h_zoomed <- s4_h_fixed + 
    geom_path(data = ell_S4, aes(x, y), inherit.aes = FALSE, colour = "black", linewidth = 0.7) +
    coord_cartesian(xlim = c(2, 5), ylim = c(0.7, 1.3)) +
    guides(fill = "none")
  
  sminus4_h_zoomed <- sminus4_h + 
    geom_path(data = ell_S4, aes(x, y), inherit.aes = FALSE, colour = "black", linewidth = 0.7) +
    coord_cartesian(xlim = c(2, 5), ylim = c(0.7, 1.3)) +
    guides(fill = "none")
  
  s3_h_zoomed | sminus3_h_zoomed
  s4_h_zoomed | sminus4_h_zoomed
  }

#Add coloums of in ellipses
subset_20_CA <- subset_20_CA %>%
  mutate(
    in_ellipse_pS3 = inside_ellipse(distance_S3, difference_S3, x0_s3, y0_s3, a_s3, b_s3, theta_s3),
    in_ellipse_mS3 = inside_ellipse(distance_step_minus_3, difference_minus_3, x0_s3, y0_s3, a_s3, b_s3, theta_s3),
    in_ellipse_pS4 = inside_ellipse(distance_S4, difference_S4, x0_s4, y0_s4, a_s4, b_s4, theta_s4),
    in_ellipse_mS4 = inside_ellipse(distance_step_minus_4, difference_minus_4, x0_s4, y0_s4, a_s4, b_s4, theta_s4),
    
    in_p_or_mS4   = in_ellipse_pS4 | in_ellipse_mS4,
    in_p_and_mS4  = in_ellipse_pS4 & in_ellipse_mS4,
  )


ellipse_summary <- subset_20_CA %>%
  summarise(
    inside_pS3          = mean(in_ellipse_pS3),
    outside_pS3.        = mean(!in_ellipse_pS3),
    inside_mS3          = mean(in_ellipse_mS3),
    outside_mS3          = mean(!in_ellipse_mS3),
    inside_pS4          = mean(in_ellipse_pS4), #Start and middle
    outside_pS4          = mean(!in_ellipse_pS4), #Start and middle
    inside_mS4          = mean(in_ellipse_mS4), #Middle and end
    outside_mS4          = mean(!in_ellipse_mS4), #Middle and end
    
    inside_p_or_mS4         = mean(in_p_or_mS4), #Whole helix
    outside_p_or_mS4        = mean(!in_p_or_mS4), #Should be not helix
  )
print(t(ellipse_summary))

ellipse_summary_sec_type <- subset_20_CA %>%
  group_by(sec_type) %>%
  summarise(
    inside_pS3          = sum(in_ellipse_pS3),
    outside_pS3         = sum(!in_ellipse_pS3),
    inside_mS3          = sum(in_ellipse_mS3),
    outside_mS3         = sum(!in_ellipse_mS3),
    inside_pS4          = sum(in_ellipse_pS4),
    outside_pS4         = sum(!in_ellipse_pS4),
    inside_mS4          = sum(in_ellipse_mS4),
    outside_mS4         = sum(!in_ellipse_mS4),
    
    inside_p_or_mS4     = sum(in_p_or_mS4), #Inside helix 
    outside_p_or_mS4    = sum(!in_p_or_mS4),
    
    inside_pS3S4        = sum(in_ellipse_pS3 & in_ellipse_pS4),
    inside_mS3S4        = sum(in_ellipse_mS3 & in_ellipse_mS4),
    
    inside_strong       = sum((in_ellipse_pS3 & in_ellipse_pS4) | (in_ellipse_mS3 & in_ellipse_mS4)), #inside helix strong
    outside_strong      = sum(!((in_ellipse_pS3 & in_ellipse_pS4) | (in_ellipse_mS3 & in_ellipse_mS4))),
    
    .groups = "drop"
  )
print(t(ellipse_summary_sec_type))

# Heatmaps of points outside

# Filter points outside the ±4 region
subset_20_CA_Helix <- subset_20_CA[subset_20_CA$sec_type == "helix", ]
outside_pmS4 <- subset_20_CA_Helix %>% 
  filter(in_p_or_mS4 == FALSE)

# Create heatmap
ggplot(outside_pmS4, aes(x = distance_step_minus_4, y = difference_minus_4)) +
  geom_bin2d(bins = 200) +
  shared_fill +  # assumes shared_fill is already defined
  labs(
    x = "Distance |O N-4| (Å)", 
    y = "Difference |C N-4| - |O N-4| (Å)",
    title = "Step -4 Geometry Outside ±4 Region"
  ) +
  xlim(2, 13) +
  ylim(-1.6, 1.5) +
  theme_minimal() +
  theme(text = element_text(size = 12))

