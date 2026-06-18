# set wd
setwd("/Volumes/gregorycostain/kd_aso")

# load libraries
library(readr)
library(dplyr)
library(lme4)
library(broom.mixed)
library(tidyverse)
library(data.table)
library(emmeans)
library(ggrepel)
library(marginaleffects)

#---------MLEM to calculate whole population differences---------#


# read data
cov_df <- read_tsv("coverage/g2p_logofunc_analysis/master_gene_coverage_by_pop_long.tsv")

# setting EUR as reference pop for comparision by reordering
cov_df <- cov_df %>%
  mutate(super_pop = as.factor(super_pop),
         super_pop = relevel(super_pop, ref = "EUR"))

# run mixed-effects logistic regression
model <- glmer(
  cbind(covered, not_covered) ~ super_pop + (1 | gene),
  data = cov_df,
  family = binomial
)

# view raw stats summary
print(summary(model))

# summarize probabilities using marginalized means
# expected coverage of superpopulation by greedy algorithm
marginalized_means <- avg_predictions(model, by = "super_pop")

# compute pairwise comparisons
comparisons <- avg_comparisons(model, variables = list(super_pop = "pairwise"))
# summarize and correct multiple testing
comparisons_summary <- tidy(comparisons, p_adjust = "holm")

# write files
write_tsv(summary, "analysis/population_bias_probabilities.tsv")
write_tsv(pairs_summary, "analysis/population_bias_pairwise_probabilities.tsv")

# plot
p1 <- ggplot(marginalized_means, aes(x = reorder(super_pop, estimate), y = estimate)) +
  geom_col(fill="steelblue") +
  coord_flip() +
  ylab("Probability of coverage") +
  xlab("Superpopulation") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1))

p1

ggsave("analysis/pop_prob_cov.png")

#---------Chi2 test to calculate per gene differences---------#

# calculate chi2 per gene
chi_gene <- cov_df %>%
  group_by(gene) %>%
  group_modify(~ {
    matrix <- as.matrix(.x[, c("covered", "not_covered")])
    test <- chisq.test(matrix)
    data.frame(
      pvalue = test$p.value,
      statistic = as.numeric(test$statistic)
    )
  })

# multiple testing correction
chi_gene <- chi_gene %>%
  mutate(adj_p = p.adjust(pvalue, method = "BH"))

# add in difference in coverage
coverage_diff <- cov_df %>%
  group_by(gene) %>%
  summarise(
    max_cov = max(proportion_covered),
    min_cov = min(proportion_covered),
    max_cov_diff = max_cov - min_cov
  )

chi_gene <- chi_gene %>%
  left_join(coverage_diff, by = "gene")


chi_gene <- chi_gene %>%
  arrange(adj_p, desc(max_cov_diff))

# write file
write_tsv(chi_gene, "analysis/population_bias_per_gene.tsv")


#---------Summarize coverage across populations---------#
cov_wide <- read_tsv("/Volumes/gregorycostain/kd_aso/coverage/g2p_logofunc_analysis/master_gene_coverage_by_pop.tsv")

# Take only coverage per population per gene
# Excluding SIK1 bc there were no SNPs pulled for that gene
ppln_cov <- cov_wide[(cov_wide$gene != "SIK1"),c(1,17:21)]

ppln_mean <- apply(ppln_cov[,2:6], 2, mean)
ppln_median <- apply(ppln_cov[,2:6], 2, median)
ppln_cov_stats <- data.frame(
  super_pop = sub("_.*", "", names(ppln_mean)),
  mean_coverage = ppln_mean,
  median_coverage = ppln_median
)

rownames(ppln_cov_stats) = NULL

write_tsv(ppln_cov_stats, "analysis/ppln_cov_mean_median.tsv")

# Identify genes with moderate/high tractability
gene_summary <- fread("analysis/gene_level_tractability_summary.csv")

eligible_genes_cov <- gene_summary %>%
  filter(tractability_tier %in% c("High tractability", "Moderate tractability")) %>%
  dplyr::select(gene) %>%
  left_join(ppln_cov, by = "gene")

ppln_mean_eligible <- apply(eligible_genes_cov[,2:6], 2, mean)
ppln_median_eligible <- apply(eligible_genes_cov[,2:6], 2, median)
ppln_cov_stats_eligible <- data.frame(
  super_pop = sub("_.*", "", names(ppln_mean_eligible)),
  mean_coverage = ppln_mean_eligible,
  median_coverage = ppln_median_eligible
)

write_tsv(ppln_cov_stats_eligible, "analysis/ppln_cov_mean_median_eligible.tsv")


#---------Summarize population size distribution & coverage correlation---------#
sample_ids <- fread("/Volumes/gregorycostain/kd_aso/metadata/1000G_unrelated_ids.txt", col.names = "sample")

pop_data <- fread("/Volumes/gregorycostain/kd_aso/metadata/20130606_g1k_3202_samples_ped_population.txt")

pop_dist <- pop_data %>%
  filter(SampleID %in% sample_ids$sample) %>%
  count(Superpopulation, name = "n_samples")

write_tsv(pop_dist, "metadata/superpop_dist.tsv")

# summarize gene level coverage
gene_cov <- cov_df %>%
  filter(gene != "SIK1") %>%
  group_by(super_pop, gene) %>%
  summarise(
    prob = sum(covered)/sum(covered + not_covered),
    .groups="drop"
  ) %>%
  merge(pop_dist, by.x = "super_pop", by.y = "Superpopulation")

# group by pop
pop_level <- gene_cov %>%
  group_by(super_pop) %>%
  summarise(
    mean_prob = mean(prob),
    n_samples = first(n_samples),
    .groups = "drop"
  )

cor.test(pop_level$n_samples, pop_level$mean_prob)

p2 <- ggplot(pop_level, aes(x = n_samples, y = mean_prob)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_text_repel(aes(label = super_pop)) +
  xlab("Number of samples") +
  ylab("Mean coverage probability") +
  theme_minimal()

p2

ggsave("analysis/size_coverage_correlation.png")

