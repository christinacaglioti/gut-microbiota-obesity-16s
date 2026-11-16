# ------------------------------------------------------------------
# 05_analysis.R
# Tertiary analysis in R (phyloseq): alpha and beta diversity,
# taxonomic composition and differential abundance (ANCOM-BC2)
# for adults with obesity vs. normal weight (PRJNA962831).
#
# Inputs (exported by scripts/04_qiime2_processing.sh):
#   results/processed/feature-table.biom
#   results/processed/taxonomy.tsv
#   results/processed/tree.nwk
#   metadata/sample-metadata.tsv
#
# Usage:
#   conda activate r-microbiome
#   Rscript scripts/05_analysis.R
#
# Outputs: results/figures/*.png, results/tables/*.tsv
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(phyloseq)
  library(biomformat)
  library(vegan)
  library(picante)
  library(ANCOMBC)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

set.seed(2026)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

GROUP_COLORS <- c(lean = "#2a78d6", obese = "#eb6834")

save_tsv <- function(df, file) {
  write.table(df, file.path("results/tables", file),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# ==================================================================
# 1. Load data and build the phyloseq object
# ==================================================================
message("1. Loading data ...")

counts <- as(biom_data(read_biom("results/processed/feature-table.biom")),
             "matrix")

# Parse SILVA taxonomy strings ("d__Bacteria;k__...;p__...;...;g__...")
tax_raw <- read.delim("results/processed/taxonomy.tsv", row.names = 1,
                      check.names = FALSE)
ranks <- c(d = "Domain", k = "Kingdom", p = "Phylum", c = "Class",
           o = "Order", f = "Family", g = "Genus", s = "Species")

parse_taxon <- function(x) {
  out <- setNames(rep(NA_character_, length(ranks)), ranks)
  for (part in trimws(strsplit(x, ";")[[1]])) {
    prefix <- sub("__.*$", "", part)
    value  <- sub("^[a-z]__", "", part)
    if (prefix %in% names(ranks) && nchar(value) > 0) {
      out[[ranks[[prefix]]]] <- value
    }
  }
  out
}

tax_mat <- do.call(rbind, lapply(tax_raw$Taxon, parse_taxon))
rownames(tax_mat) <- rownames(tax_raw)
tax_mat <- tax_mat[, colSums(!is.na(tax_mat)) > 0, drop = FALSE]

# Unresolved genera ("Incertae_Sedis", "uncultured", missing) are labelled
# with the next higher informative rank, so that taxa from different
# families are not merged just because they share a placeholder name.
placeholder <- function(x) is.na(x) | x %in% c("Incertae_Sedis", "uncultured",
                                               "unidentified")
higher <- ifelse(!placeholder(tax_mat[, "Family"]), tax_mat[, "Family"],
          ifelse(!placeholder(tax_mat[, "Order"]),  tax_mat[, "Order"],
                 "Unknown"))
genus  <- tax_mat[, "Genus"]
unres  <- placeholder(genus)
tax_mat[unres, "Genus"] <- paste0(higher[unres], "_",
                                  ifelse(is.na(genus[unres]), "unclassified",
                                         genus[unres]))

meta <- read.delim("metadata/sample-metadata.tsv", row.names = 1,
                   check.names = FALSE)
meta$group <- factor(meta$group, levels = c("lean", "obese"))

tree <- read_tree("results/processed/tree.nwk")

ps <- phyloseq(otu_table(counts, taxa_are_rows = TRUE),
               tax_table(tax_mat),
               sample_data(meta),
               phy_tree(tree))

# Keep Bacteria/Archaea; remove mitochondria and chloroplasts
ps <- subset_taxa(ps, Domain %in% c("Bacteria", "Archaea"))
ps <- subset_taxa(ps, is.na(Family) | Family != "Mitochondria")
ps <- subset_taxa(ps, is.na(Order)  | Order  != "Chloroplast")
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

message(sprintf("   %d samples, %d ASVs after filtering",
                nsamples(ps), ntaxa(ps)))

# ==================================================================
# 2. Rarefaction (for diversity metrics only)
# ==================================================================
# Depth = smallest library (no sample is lost). The original study
# rarefied to 3,800 reads; our depth is much higher.
depth <- min(sample_sums(ps))
message(sprintf("2. Rarefying to %d reads per sample ...", depth))
ps_rare <- rarefy_even_depth(ps, sample.size = depth, rngseed = 2026,
                             replace = FALSE, verbose = FALSE)

# phyloseq's UniFrac requires a strictly bifurcating tree; after pruning,
# the tree can contain polytomies (warning "data length ... is not a
# sub-multiple ..."). Resolve them into zero-length branches.
phy_tree(ps_rare) <- ape::multi2di(phy_tree(ps_rare))

# ==================================================================
# 3. Alpha diversity
# ==================================================================
# Chao1 is not used: DADA2 removes singletons, which Chao1 relies on.
message("3. Alpha diversity ...")

alpha <- estimate_richness(ps_rare, measures = c("Observed", "Shannon"))
rownames(alpha) <- sample_names(ps_rare)

otu_r <- as(otu_table(ps_rare), "matrix")
if (taxa_are_rows(ps_rare)) otu_r <- t(otu_r)
faith <- pd(otu_r, phy_tree(ps_rare), include.root = TRUE)
alpha$Faith_PD <- faith[rownames(alpha), "PD"]
alpha$group <- sample_data(ps_rare)$group
alpha$sample_id <- rownames(alpha)
save_tsv(alpha, "alpha_diversity.tsv")

metrics <- c("Observed", "Shannon", "Faith_PD")
alpha_tests <- do.call(rbind, lapply(metrics, function(m) {
  w <- wilcox.test(alpha[[m]] ~ alpha$group, exact = FALSE)
  data.frame(metric       = m,
             median_lean  = median(alpha[alpha$group == "lean",  m]),
             median_obese = median(alpha[alpha$group == "obese", m]),
             W            = unname(w$statistic),
             p_value      = w$p.value)
}))
alpha_tests$p_adj_BH <- p.adjust(alpha_tests$p_value, method = "BH")
save_tsv(alpha_tests, "alpha_diversity_tests.tsv")
print(alpha_tests)

alpha_long <- alpha |>
  pivot_longer(all_of(metrics), names_to = "metric", values_to = "value") |>
  mutate(metric = factor(metric, levels = metrics,
                         labels = c("Observed ASVs", "Shannon index",
                                    "Faith's PD")))

p_alpha <- ggplot(alpha_long, aes(group, value, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.35,
               colour = "grey30") +
  geom_jitter(aes(colour = group), width = 0.12, size = 1.8, alpha = 0.8) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = GROUP_COLORS) +
  scale_colour_manual(values = GROUP_COLORS) +
  labs(x = NULL, y = NULL,
       title = "Alpha diversity: obesity vs. normal weight",
       subtitle = sprintf("Rarefied to %s reads per sample; Wilcoxon rank-sum test",
                          format(depth, big.mark = ","))) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", panel.grid.minor = element_blank())
ggsave("results/figures/alpha_diversity.png", p_alpha,
       width = 9, height = 4, dpi = 300)

# ==================================================================
# 4. Beta diversity
# ==================================================================
message("4. Beta diversity ...")

dists <- list(
  Bray_Curtis        = phyloseq::distance(ps_rare, method = "bray"),
  Unweighted_UniFrac = phyloseq::distance(ps_rare, method = "unifrac"),
  Weighted_UniFrac   = phyloseq::distance(ps_rare, method = "wunifrac")
)
md <- as(sample_data(ps_rare), "data.frame")

beta_tests <- do.call(rbind, lapply(names(dists), function(n) {
  d  <- dists[[n]]
  pa <- adonis2(d ~ group, data = md, permutations = 999)
  bd <- betadisper(d, md$group)
  pt <- permutest(bd, permutations = 999)
  data.frame(distance      = n,
             PERMANOVA_R2  = pa$R2[1],
             PERMANOVA_F   = pa$F[1],
             PERMANOVA_p   = pa$`Pr(>F)`[1],
             betadisper_p  = pt$tab$`Pr(>F)`[1])
}))
save_tsv(beta_tests, "beta_diversity_tests.tsv")
print(beta_tests)

for (n in names(dists)) {
  ord <- ordinate(ps_rare, method = "PCoA", distance = dists[[n]])
  bt  <- beta_tests[beta_tests$distance == n, ]
  p_ord <- plot_ordination(ps_rare, ord, color = "group") +
    geom_point(size = 2.5, alpha = 0.85) +
    stat_ellipse(type = "t", linewidth = 0.6) +
    scale_colour_manual(values = GROUP_COLORS) +
    labs(title = paste("PCoA -", gsub("_", " ", n)),
         subtitle = sprintf("PERMANOVA R2 = %.3f, p = %.3f; betadisper p = %.3f",
                            bt$PERMANOVA_R2, bt$PERMANOVA_p, bt$betadisper_p),
         colour = NULL) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  ggsave(sprintf("results/figures/pcoa_%s.png", tolower(n)), p_ord,
         width = 6.5, height = 5, dpi = 300)
}

# ==================================================================
# 5. Taxonomic composition (top genera per sample)
# ==================================================================
message("5. Taxonomic composition ...")

genus_rel <- psmelt(transform_sample_counts(ps, function(x) x / sum(x))) |>
  group_by(Sample, group, Genus) |>
  summarise(Abundance = sum(Abundance), .groups = "drop")

top_genera <- genus_rel |>
  group_by(Genus) |>
  summarise(mean = mean(Abundance), .groups = "drop") |>
  arrange(desc(mean)) |>
  slice_head(n = 8) |>
  pull(Genus)

genus_means <- genus_rel |>
  group_by(group, Genus) |>
  summarise(mean_rel_abundance = mean(Abundance), .groups = "drop") |>
  pivot_wider(names_from = group, values_from = mean_rel_abundance) |>
  arrange(desc(lean + obese))
save_tsv(genus_means, "genus_mean_relative_abundance.tsv")

comp <- genus_rel |>
  mutate(Genus = ifelse(Genus %in% top_genera, Genus, "Other")) |>
  group_by(Sample, group, Genus) |>
  summarise(Abundance = sum(Abundance), .groups = "drop") |>
  mutate(Genus = factor(Genus, levels = c(top_genera, "Other")))

sample_order <- comp |>
  filter(Genus == top_genera[1]) |>
  arrange(desc(Abundance)) |>
  pull(Sample)
comp$Sample <- factor(comp$Sample, levels = sample_order)

genus_palette <- c("#2a78d6", "#eb6834", "#1f9e89", "#d9a21b",
                   "#8e5bd6", "#d6457a", "#5c8a2e", "#7a6f66", "#c9c5bd")
names(genus_palette) <- c(top_genera, "Other")

p_comp <- ggplot(comp, aes(Sample, Abundance, fill = Genus)) +
  geom_col(width = 0.95) +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = genus_palette) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = "Samples", y = "Relative abundance", fill = "Genus",
       title = "Gut microbiota composition (top 8 genera)") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank())
ggsave("results/figures/composition_top_genera.png", p_comp,
       width = 11, height = 5, dpi = 300)

# ==================================================================
# 6. Differential abundance with ANCOM-BC2 (unrarefied counts)
# ==================================================================
# ANCOM-BC2 accounts for compositionality and sampling fraction.
# Taxa present in < 10 % of samples are excluded (prv_cut = 0.10).
# FDR control with Benjamini-Hochberg.
message("6. Differential abundance (ANCOM-BC2) - this can take a while ...")

run_ancombc2 <- function(level) {
  out <- ancombc2(data = ps, tax_level = level,
                  fix_formula = "group", group = "group",
                  p_adj_method = "BH", prv_cut = 0.10, lib_cut = 1000,
                  struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05,
                  n_cl = 1, verbose = FALSE)
  res <- out$res
  res <- res[order(res$q_groupobese), ]
  list(res = res, zero = out$zero_ind)
}

da_genus  <- run_ancombc2("Genus")
da_family <- run_ancombc2("Family")
save_tsv(da_genus$res,  "ancombc2_genus.tsv")
save_tsv(da_family$res, "ancombc2_family.tsv")
if (!is.null(da_genus$zero)) save_tsv(da_genus$zero, "ancombc2_genus_structural_zeros.tsv")

cols <- c("taxon", "lfc_groupobese", "se_groupobese", "p_groupobese",
          "q_groupobese", "diff_groupobese", "passed_ss_groupobese")

sig_genus <- da_genus$res[da_genus$res$diff_groupobese &
                          da_genus$res$passed_ss_groupobese, cols]
message(sprintf("   %d genera differentially abundant (q < 0.05, passed sensitivity analysis)",
                nrow(sig_genus)))
print(sig_genus)

# Taxa from the study hypothesis and the literature
hyp_pattern <- "Akkermansia|Christensenellaceae|Alistipes|Odoribacter|Segatella|Lachnospira|Faecalibacterium"
hypothesis_taxa <- rbind(
  cbind(level = "Genus",  da_genus$res[grepl(hyp_pattern, da_genus$res$taxon), cols]),
  cbind(level = "Family", da_family$res[grepl(hyp_pattern, da_family$res$taxon), cols])
)
save_tsv(hypothesis_taxa, "hypothesis_taxa.tsv")
print(hypothesis_taxa)

# Plot: significant genera, or the 15 lowest q-values if none passes
plot_df <- if (nrow(sig_genus) > 0) sig_genus else head(da_genus$res[, cols], 15)
plot_df <- plot_df |>
  mutate(direction = ifelse(lfc_groupobese > 0, "higher in obese", "lower in obese"),
         taxon = reorder(taxon, lfc_groupobese))

p_da <- ggplot(plot_df, aes(lfc_groupobese, taxon, colour = direction)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_errorbarh(aes(xmin = lfc_groupobese - 1.96 * se_groupobese,
                     xmax = lfc_groupobese + 1.96 * se_groupobese),
                 height = 0.25) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c("higher in obese" = "#eb6834",
                                 "lower in obese"  = "#2a78d6")) +
  labs(x = "Log fold change (obese vs. lean) with 95% CI", y = NULL,
       colour = NULL,
       title = "Differential abundance (ANCOM-BC2, genus level)",
       subtitle = if (nrow(sig_genus) > 0)
         "Genera with q < 0.05 that passed the sensitivity analysis"
       else "No genus reached q < 0.05; showing the 15 lowest q-values") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave("results/figures/ancombc2_genus.png", p_da,
       width = 8, height = max(4, 0.35 * nrow(plot_df) + 2), dpi = 300)

# ==================================================================
# 7. Reproducibility
# ==================================================================
writeLines(capture.output(sessionInfo()), "results/sessionInfo_R.txt")
message("Done. Figures in results/figures/, tables in results/tables/")