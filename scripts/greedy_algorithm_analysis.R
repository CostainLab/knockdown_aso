# Set working directory 
setwd("/Volumes/gregorycostain/kd_aso")

# Load libraries
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(scales)

# Load in master table of greedy algorithm output 
master = read.csv("greedy_out/master_greedy_selected_snps_repeat_filtered_29Jun26.tsv", sep = "\t")

# 1. Gene-level summary

# Old tiers:
   # tractability_tier = case_when(
    #  !is.na(snps_to_90) & snps_to_90 <= 10 ~ "High tractability",
    #  !is.na(snps_to_80) & snps_to_80 <= 10 ~ "Moderate tractability",
    #  max_coverage < 0.80 ~ "Low tractability",
    #  TRUE ~ "Intermediate"

gene_summary = master %>% 
  group_by(gene) %>% 
  summarise(
    n_selected_snps = max(rank, na.rm = TRUE),
    max_coverage = max(cumulative_fraction, na.rm = TRUE),
    coverage_after_1_snp = cumulative_fraction[rank == 1][1],
    coverage_after_3_snps = max(cumulative_fraction[rank <= 3], na.rm = TRUE),
    coverage_after_5_snps = max(cumulative_fraction[rank <= 5], na.rm = TRUE),
    snps_to_80 = ifelse(any(cumulative_fraction >= 0.80), min(rank[cumulative_fraction >= 0.80]), NA),
    snps_to_90 = ifelse(any(cumulative_fraction >= 0.90), min(rank[cumulative_fraction >= 0.90]), NA),
    snps_to_95 = ifelse(any(cumulative_fraction >= 0.95), min(rank[cumulative_fraction >= 0.95]), NA),
    .groups = "drop"
  ) %>%
  mutate(
    tractability_tier = case_when(
      !is.na(snps_to_90) & snps_to_90 <= 10 ~ "High Feasibility",
      !is.na(snps_to_80) & snps_to_80 <= 10 ~ "Moderate Feasibility",
      !is.na(snps_to_80) & snps_to_80 > 10 ~ "Limited Feasibility",
      max_coverage < 0.80 ~ "Low Feasibility"
    )
  ) %>%
  arrange(desc(max_coverage), snps_to_90)

write.csv(gene_summary, "analysis/gene_level_tractability_summary_29Jun26.csv", row.names = FALSE)

# -----------------------------
# 2. Overall summary numbers
# -----------------------------

overall_summary = gene_summary %>%
  summarise(
    n_genes = n(),
    
    n_genes_ge_80 = sum(max_coverage >= 0.80, na.rm = TRUE),
    n_genes_ge_90 = sum(max_coverage >= 0.90, na.rm = TRUE),
    n_genes_ge_95 = sum(max_coverage >= 0.95, na.rm = TRUE),
    
    pct_genes_ge_80 = 100 * mean(max_coverage >= 0.80, na.rm = TRUE),
    pct_genes_ge_90 = 100 * mean(max_coverage >= 0.90, na.rm = TRUE),
    pct_genes_ge_95 = 100 * mean(max_coverage >= 0.95, na.rm = TRUE),
    
    median_max_coverage = median(max_coverage, na.rm = TRUE),
    
    median_snps_to_80 = median(snps_to_80, na.rm = TRUE),
    median_snps_to_90 = median(snps_to_90, na.rm = TRUE),
    median_snps_to_95 = median(snps_to_95, na.rm = TRUE),
    
    # Overall SNP count statistics
    median_selected_snps = median(n_selected_snps, na.rm = TRUE),
    min_selected_snps = min(n_selected_snps, na.rm = TRUE),
    max_selected_snps = max(n_selected_snps, na.rm = TRUE),
    
    # Optional formatted range string
    selected_snps_range = paste0(
      min(n_selected_snps, na.rm = TRUE),
      "-",
      max(n_selected_snps, na.rm = TRUE)
    )
  )

write.csv(overall_summary, "analysis/overall_tractability_summary_29Jun26.csv", row.names = FALSE)

# -----------------------------
# 3. Top tractable genes
# -----------------------------

top_genes_by_5_snps = gene_summary %>%
  arrange(desc(coverage_after_5_snps)) %>%
  select(gene, n_selected_snps, max_coverage, coverage_after_1_snp, coverage_after_3_snps, coverage_after_5_snps, snps_to_80, snps_to_90, snps_to_95) %>%
  head(25)

write.csv(top_genes_by_5_snps, "analysis/top_25_genes_by_5_snp_coverage_29Jun26.csv", row.names = FALSE)

# -----------------------------
# 4. Plot: maximum coverage distribution
# -----------------------------

p1 = ggplot(gene_summary, aes(x = max_coverage)) +
  geom_histogram(
    binwidth = 0.05,
    boundary = 0,
    fill = "#3066BE",
    color = "white",
    linewidth = 0.3
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.10),
    labels = scales::number_format(accuracy = 1, scale = 100)
  ) +
  scale_y_continuous(
    limits = c(0, 80),
    breaks = seq(0, 80, 20)
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = "Maximum cumulative population coverage (%)",
    y = "Number of genes"
  )

p1 

ggsave("analysis/max_coverage_distribution_29Jun26.png", p1, width = 7, height = 5, dpi = 300)

# -----------------------------
# 5. Plot: SNPs needed to reach 80% and 90%
# -----------------------------
snps_threshold_long = gene_summary %>%
  select(gene, snps_to_80, snps_to_90, snps_to_95) %>%
  tidyr::pivot_longer(
    cols = starts_with("snps_to_"),
    names_to = "threshold",
    values_to = "n_snps"
  ) %>%
  filter(!is.na(n_snps)) %>%
  mutate(
    threshold = recode(
      threshold,
      snps_to_80 = "80% coverage",
      snps_to_90 = "90% coverage",
      snps_to_95 = "95% coverage"
    ),
    threshold = factor(
      threshold,
      levels = c("80% coverage", "90% coverage", "95% coverage")
    )
  )

p2 = ggplot(snps_threshold_long, aes(x = n_snps)) +
  geom_histogram(
    binwidth = 2,
    boundary = 0,
    fill = "#3066BE",
    color = "white",
    linewidth = 0.3
  ) +
  facet_wrap(~ threshold, nrow = 1) +
  scale_x_continuous(
    breaks = seq(0, 60, 10),
    limits = c(0, 60)
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 20),
    limits = c(0, 100)
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "grey40"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0)
  ) +
  labs(
    x = "Number of selected SNPs required",
    y = "Number of genes"
  )

p2

ggsave("analysis/snps_needed_for_coverage_thresholds_29Jun26.png", p2, width = 7, height = 5, dpi = 300)

# -----------------------------
# 6. Plot: greedy coverage curves for top example genes
# -----------------------------

example_genes = gene_summary %>%
  arrange(desc(coverage_after_5_snps)) %>%
  slice_head(n = 10) %>%
  pull(gene)

curve_data = master %>%
  filter(gene %in% example_genes) %>%
  mutate(gene = factor(gene, levels = example_genes))

p3 = ggplot(curve_data, aes(x = rank, y = cumulative_fraction, group = gene, color = gene)) +
  geom_line(linewidth = 1.1, alpha = 0.9) +
  geom_point(size = 2) +
  scale_y_continuous(
    limits = c(0.45, 1.00),
    breaks = seq(0.5, 1.0, 0.1),
    labels = scales::number_format(accuracy = 1, scale = 100)
  ) +
  scale_x_continuous(
    breaks = 1:5,
    limits = c(1, 5)
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "Number of selected SNPs",
    y = "Cumulative population coverage (%)"
  )


p3

ggsave("analysis/top_gene_greedy_coverage_curves_29Jun26.png", p3, width = 5, height = 5, dpi = 300)

# Plot top 10 genes in high feasibility category

example_genes = gene_summary %>%
  filter(tractability_tier == "High Feasibility") %>%
  arrange(desc(max_coverage)) %>%
  slice_head(n = 10) %>%
  pull(gene)

curve_data = master %>%
  filter(gene %in% example_genes) %>%
  mutate(gene = factor(gene, levels = example_genes))

p3.1 = ggplot(curve_data, aes(x = rank, y = cumulative_fraction, group = gene, color = gene)) +
  geom_line(linewidth = 1.1, alpha = 0.9) +
  geom_point(size = 2) +
  scale_y_continuous(
    limits = c(0.45, 1.00),
    breaks = seq(0.5, 1.0, 0.1),
    labels = scales::number_format(accuracy = 1, scale = 100)
  ) +
  scale_x_continuous(
    breaks = 1:6,
    limits = c(1, 6)
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "Number of selected SNPs",
    y = "Cumulative population coverage (%)"
  )


p3.1

ggsave("analysis/top_gene_greedy_coverage_curves_hfeas_29Jun26.png", p3.1, width = 5, height = 5, dpi = 300)

# -----------------------------
# 7. Tractability tier counts
# -----------------------------

tier_counts = gene_summary %>%
  count(tractability_tier) %>%
  mutate(
    tractability_tier = factor(
      tractability_tier,
      levels = c(
        "Low Feasibility",
        "Limited Feasibility",
        "Moderate Feasibility",
        "High Feasibility"
      )
    )
  ) %>%
  arrange(tractability_tier)

p4 = ggplot(
  tier_counts,
  aes(x = tractability_tier, y = n, fill = tractability_tier)
) +
  geom_col(width = 0.75) +
  geom_text(
    aes(label = n),
    vjust = -0.4,
    size = 3.5
  ) +
  scale_fill_manual(
    values = c(
      "Low Feasibility" = "#963484",
      "Limited Feasibility" = "#3066BE",
      "Moderate Feasibility" = "#60AFFF",
      "High Feasibility" = "#28C2FF"
    )
  ) +
  scale_y_continuous(
    breaks = seq(0, 90, 15),
    limits = c(0,150)
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none"
  ) +
  labs(
    x = "Tractability tier",
    y = "Number of genes",
  )

p4 

ggsave("analysis/tractability_tier_counts_29Jun26.png", p4, width = 7, height = 5, dpi = 300)

# -----------------------------
# 8. Print quick abstract-friendly summary
# -----------------------------

print(overall_summary)
print(tier_counts)
print(head(gene_summary, 20))
