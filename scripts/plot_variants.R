## Plots 
# Load required libraries
library(tidyverse)
library(patchwork)

# Set parameters
plasmid_id <- "NC_003277.2"
BIN <- 100L   # variant bins

# Load variants 
variants <- read.table("variants_hq.bed", header = FALSE)
colnames(variants) <- c("chr", "start", "end", "type", "qual")
variants <- as_tibble(variants)

# Load fine-window coverage 
coverage_small <- read.table("coverage_100bp.tsv", header = FALSE)
colnames(coverage_small) <- c("chr", "start", "end", "depth")
coverage_small <- as_tibble(coverage_small)

# Filter variants on plasmid 
v <- variants %>%
  filter(chr == plasmid_id) %>%
  mutate(pos = start + 1) %>%
  mutate(type = case_when(
    type %in% c("SNP", "SNV") ~ "SNP",
    type %in% c("INS", "INSERTION") ~ "INS",
    type %in% c("DEL", "DELETION") ~ "DEL",
    TRUE ~ as.character(type)
  ))

stopifnot(nrow(v) > 0)

# Determine plasmid length 
plasmid_len <- max(
  coverage_small %>% filter(chr == plasmid_id) %>% pull(end),
  v %>% pull(pos),
  na.rm = TRUE
)

# Create variant bins (100 bp) 
bins <- tibble(
  bin = 0:floor((plasmid_len - 1) / BIN),
  bin_mid = (bin * BIN) + BIN / 2
)

vb <- v %>%
  mutate(bin = floor((pos - 1) / BIN)) %>%
  count(bin, type, name = "n") %>%
  right_join(bins, by = "bin") %>%
  mutate(n = replace_na(n, 0))

vb_total <- vb %>%
  group_by(bin, bin_mid) %>%
  summarise(n = sum(n), .groups = "drop")

# Create variant plot with BLACK lollipops
p_var <- ggplot() +
  geom_segment(
    data = vb_total,
    aes(x = bin_mid, xend = bin_mid, y = 0, yend = n),
    linewidth = 0.35, alpha = 0.20, color = "grey35"
  ) +
  geom_point(
    data = vb_total,
    aes(x = bin_mid, y = n),
    size = 0.7, alpha = 0.20, color = "grey35"
  ) +
  geom_segment(
    data = vb %>% filter(n > 0),
    aes(x = bin_mid, xend = bin_mid, y = 0, yend = n, color = type),
    linewidth = 0.55
  ) +
  geom_point(
    data = vb %>% filter(n > 0),
    aes(x = bin_mid, y = n, color = type),
    size = 1.2
  ) +
  scale_color_manual(values = c(SNP = "black")) +
  coord_cartesian(xlim = c(0, plasmid_len), expand = FALSE) +
  labs(
    title = "Variant Density Across Plasmid NC_003277.2",
    subtitle = "Raw Nanopore Reads Aligned to ASM694v2 Reference",
    y = "Variants per bin",
    x = NULL,
    color = "Type"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

# Create coverage heatmap strip with GREEN shades
cov_strip <- coverage_small %>%
  filter(chr == plasmid_id) %>%
  mutate(depth_cap = pmin(depth, as.numeric(quantile(depth, 0.99, na.rm = TRUE))))

p_cov <- ggplot(cov_strip) +
  geom_rect(aes(
    xmin = start, xmax = end,
    ymin = 0, ymax = 1,
    fill = depth_cap
  ), color = NA) +
  scale_fill_gradientn(colors = c("white", "lightgreen", "darkgreen"), name = "Coverage") +
  coord_cartesian(xlim = c(0, plasmid_len), expand = FALSE) +
  labs(x = "Position on plasmid (bp)", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  )

# Combine plots
final_plot <- p_var / p_cov + plot_layout(heights = c(4, 0.9))

# Save plot to results/plots directory
ggsave("results/plots/plasmid_variants_coverage.png", final_plot, width = 12, height = 7, dpi = 300)

# Print coverage summary
cat("\nCoverage summary:\n")
print(summary(cov_strip$depth))
cat("\nCoverage quantiles:\n")
print(quantile(cov_strip$depth, c(.9, .95, .99, 1)))

cat("\nPlot saved as: results/plots/plasmid_variants_coverage.png\n")


##plot2

# Load required libraries
library(tidyverse)
library(patchwork)

# Set parameters
BIN <- 100L   # variant bins

# Load variants 
variants <- read.table("variants_hq.bed", header=FALSE)
colnames(variants) <- c("chr", "start", "end", "type", "qual")
variants <- as_tibble(variants)

# Load fine-window coverage 
coverage_small <- read.table("coverage_100bp.tsv", header=FALSE)
colnames(coverage_small) <- c("chr", "start", "end", "depth")
coverage_small <- as_tibble(coverage_small)

# Function to create variant + coverage plot for a given sequence
create_variant_plot <- function(seq_id, seq_name) {
  
  # Filter variants on this sequence
  v <- variants %>%
    filter(chr == seq_id) %>%
    mutate(pos = start + 1) %>%
    mutate(type = case_when(
      type %in% c("SNP","SNV") ~ "SNP",
      type %in% c("INS","INSERTION") ~ "INS",
      type %in% c("DEL","DELETION") ~ "DEL",
      TRUE ~ as.character(type)
    ))
  
  if(nrow(v) == 0) {
    warning(paste("No variants found for", seq_id))
    return(NULL)
  }
  
  # Determine sequence length 
  seq_len <- max(
    coverage_small %>% filter(chr == seq_id) %>% pull(end),
    v %>% pull(pos),
    na.rm = TRUE
  )
  
  # Create variant bins (100 bp) 
  bins <- tibble(
    bin = 0:floor((seq_len - 1) / BIN),
    bin_mid = (bin * BIN) + BIN/2
  )
  
  vb <- v %>%
    mutate(bin = floor((pos - 1) / BIN)) %>%
    count(bin, type, name = "n") %>%
    right_join(bins, by = "bin") %>%
    mutate(n = replace_na(n, 0))
  
  vb_total <- vb %>%
    group_by(bin, bin_mid) %>%
    summarise(n = sum(n), .groups="drop")
  
  # Create variant plot with BLACK lollipops
  p_var <- ggplot() +
    geom_segment(
      data = vb_total,
      aes(x = bin_mid, xend = bin_mid, y = 0, yend = n),
      linewidth = 0.35, alpha = 0.20, color = "grey35"
    ) +
    geom_point(
      data = vb_total,
      aes(x = bin_mid, y = n),
      size = 0.7, alpha = 0.20, color = "grey35"
    ) +
    geom_segment(
      data = vb %>% filter(n > 0),
      aes(x = bin_mid, xend = bin_mid, y = 0, yend = n, color = type),
      linewidth = 0.55
    ) +
    geom_point(
      data = vb %>% filter(n > 0),
      aes(x = bin_mid, y = n, color = type),
      size = 1.2
    ) +
    scale_color_manual(values = c(SNP="black", INS="black", DEL="black")) +
    coord_cartesian(xlim = c(0, seq_len), expand = FALSE) +
    labs(
      title = paste0("Variant Density Across ", seq_name, " (", seq_id, ")"),
      subtitle = "Raw Nanopore Reads Aligned to ASM694v2 Reference",
      y = "Variants per bin",
      x = NULL,
      color = "Type"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  # Create coverage heatmap strip with GREEN shades
  cov_strip <- coverage_small %>%
    filter(chr == seq_id) %>%
    mutate(depth_cap = pmin(depth, as.numeric(quantile(depth, 0.99, na.rm = TRUE))))
  
  cov_strip <- cov_strip %>%
    mutate(depth_plot = ifelse(depth == 0, NA, depth_cap))
  
  p_cov <- ggplot(cov_strip) +
    geom_rect(aes(
      xmin = start, xmax = end,
      ymin = 0, ymax = 1,
      fill = depth_cap
    ), color = NA) +
    scale_fill_gradientn(colors = c("white", "lightgreen", "darkgreen"), name = "Coverage") +
    coord_cartesian(xlim = c(0, seq_len), expand = FALSE) +
    labs(x = paste("Position on", seq_name, "(bp)"), y = NULL) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank()
    )
  
  # Combine plots
  final_plot <- p_var / p_cov + plot_layout(heights = c(4, 0.9))
  
  # Print coverage summary
  cat("\n=== Coverage summary for", seq_name, "===\n")
  print(summary(cov_strip$depth))
  cat("\nCoverage quantiles:\n")
  print(quantile(cov_strip$depth, c(.9, .95, .99, 1)))
  
  return(final_plot)
}

# Generate plots for both chromosome and plasmid
cat("Generating chromosome plot...\n")
chrom_plot <- create_variant_plot("NC_003197.2", "Chromosome")
if(!is.null(chrom_plot)) {
  ggsave("results/plots/chromosome_variants_coverage.png", chrom_plot, 
         width = 12, height = 7, dpi = 300)
  cat("✓ Chromosome plot saved to results/plots/chromosome_variants_coverage.png\n")
}

cat("\nGenerating plasmid plot...\n")
plasmid_plot <- create_variant_plot("NC_003277.2", "Plasmid")
if(!is.null(plasmid_plot)) {
  ggsave("results/plots/plasmid_variants_coverage.png", plasmid_plot, 
         width = 12, height = 7, dpi = 300)
  cat("✓ Plasmid plot saved to results/plots/plasmid_variants_coverage.png\n")
}

cat("\n=== All plots generated successfully! ===\n")
