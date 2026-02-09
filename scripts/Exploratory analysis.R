
setwd("~/binf6110/assignment1")

library(vcfR)
library(dplyr)
library(ggplot2)

vcf <- read.vcfR("/School/November 1st 2025/Master/Winter 2026/DATA 6110/medaka.af.vcf.gz")

pos_df <- data.frame(
  chrom = vcf@fix[, "CHROM"],
  pos   = as.integer(vcf@fix[, "POS"]),
  ref   = vcf@fix[, "REF"],
  alt   = vcf@fix[, "ALT"],
  af    = as.numeric(vcf@fix[, "INFO"] %>% sub(".*AF=([^;]+).*", "\\1", .))
)

snps <- pos_df %>%
  filter(!is.na(pos),
         chrom == "NC_003197.2",
         nchar(ref) == 1,
         nchar(alt) == 1,
         !is.na(af))

ggplot(snps, aes(x = pos, y = af)) +
  geom_segment(aes(xend = pos, yend = 0), linewidth = 0.2) +
  geom_point(size = 0.7) +
  theme_minimal() +
  labs(
    x = "Genomic position on NC_003197.2",
    y = "Allele frequency",
    title = "SNV allele frequencies called by Medaka vs ASM694v2"
  )

library(Rsamtools)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Paths
reads_bam <- "results/align/reads_vs_ASM694v2.bam"
asm_bam <- "results/align/polished_flye_vs_ASM694v2.bam"
depth_reads <- "results/plots/depth_ASM694v2.txt"  # Existing

# 1. Coverage Histogram (reads)
depth_df <- read.table(depth_reads, col.names=c("ref","pos","depth"))
p1 <- ggplot(depth_df, aes(x=depth)) +
  geom_histogram(bins=100, fill="steelblue", alpha=0.7) +
  labs(title="Reads-to-Ref Coverage (150x mean)", x="Depth", y="Positions") +
  theme_minimal()

ggsave("results/plots/reads_coverage_hist.png", p1, width=8, height=6)

# 2. Depth per Contig (boxplot)
p2 <- ggplot(depth_df %>% filter(ref %in% c("NC_003197.2","NC_003277.2")), 
             aes(x=ref, y=depth, fill=ref)) +
  geom_boxplot() + labs(title="Depth: Chrom vs Plasmid", y="Depth") +
  theme_minimal()

ggsave("results/plots/reads_depth_per_contig.png", p2)

# 3. Assembly Coverage (flagstat summary or contig %)
# Run samtools coverage first
system("samtools coverage results/align/polished_flye_vs_ASM694v2.bam > results/plots/asm_coverage.txt")
asm_cov <- read.table("results/plots/depth_ASM694v2.txt", header=T, skip=1)
print(asm_cov)  # Table for report

# Assembly dotplot proxy: Use nucmer coords if generated

##3
library(Rsamtools)
library(ggplot2)
library(dplyr)



plot_alignment_dotplot <- function(bam_file, title, out_png) {
  # Param with mapq filter (MAPQ>0 primary)
  param <- ScanBamParam(what=c("rname","qname","pos","qwidth","flag","mapq"),
                        flag=scanBamFlag(isSecondaryAlignment=FALSE,
                                         isSupplementaryAlignment=FALSE,
                                         isDuplicate=FALSE),
                        mapqFilter=1)
  
  alns <- scanBam(bam_file, param=param)[[1]]
  
  # Dataframe
  df <- data.frame(
    ref = as.character(alns$rname),
    ref_pos = alns$pos / 1e6,
    query = as.character(alns$qname),
    query_pos = (alns$pos + alns$qwidth) / 1e6,  # End pos approx
    mapq = alns$mapq
  ) %>% filter(!is.na(ref))
  
  p <- ggplot(df, aes(ref_pos, query_pos, color=ref, alpha=mapq/60)) +
    geom_point(size=0.3) +
    facet_wrap(~ref, scales="free", ncol=1) +
    labs(title=title, x="Ref Pos (Mb)", y="Query Pos (Mb)") +
    theme_minimal() + guides(alpha="none")
  
  ggsave(out_png, p, width=10, height=8, dpi=300)
}

# Run
plot_alignment_dotplot("results/align/reads_vs_ASM694v2.bam", 
                       "Ref x Raw Reads Dotplot", "results/plots/reads_dotplot.png")

plot_alignment_dotplot("results/align/polished_flye_vs_ASM694v2.bam", 
                       "Ref x Assembly Dotplot", "results/plots/asm_dotplot.png")

###4

install.packages("RcppRoll")

library(ggplot2); library(dplyr); library(RcppRoll)


depth_reads <- read.table("results/plots/depth_ASM694v2.txt", 
                          col.names=c("contig","pos","depth"))

# Smooth
depth_reads <- depth_reads %>%
  group_by(contig) %>%
  mutate(depth_smooth = roll_mean(depth, 1001, fill=NA)) %>%  # 1kb smooth
  ungroup()

# Plots as before...
p1 <- ggplot(depth_reads) + geom_line(aes(pos/1e6, depth_smooth/1e3, color=contig), size=1) +
  facet_grid(contig ~ ., scales="free") + theme_minimal()
ggsave("results/plots/reads_overview.png", p1, w=10, h=8)

##5 

library(ggplot2); library(dplyr); library(RcppRoll)  # install.packages(c("ggplot2","dplyr")); bioconda r-rcpproll

depth <- read.table("results/plots/depth_ASM694v2.txt", 
                    col.names=c("contig","pos","depth"), stringsAsFactors=F)

depth$depth_smooth <- roll_mean(depth$depth, 1001, fill=0)  # 1kb window

chrom <- filter(depth, contig=="NC_003197.2")
plasm <- filter(depth, contig=="NC_003277.2")

p_chrom <- ggplot(chrom, aes(pos/1e6, depth_smooth)) + geom_line(color="#1f77b4", size=0.8) + 
  theme_minimal(base_size=12) + labs(x="Chromosome Position (Mb)", y="Coverage Depth", 
                                     title="Raw Reads Coverage vs. ASM694v2 (Mean: 150x, Breadth: 96.8%)") + 
  scale_x_continuous(breaks=seq(0,5,1))

p_plasm <- ggplot(plasm, aes(pos/1000, depth_smooth)) + geom_line(color="#ff7f0e", size=0.8) + 
  theme_minimal(base_size=12) + labs(x="Plasmid Position (kb)", y="Coverage Depth", 
                                     title="Raw Reads Coverage: Plasmid (Median Depth ~200x)") + 
  scale_x_continuous(breaks=seq(0,100,20))

ggsave("results/plots/reads_coverage_chrom.png", p_chrom, dpi=300, width=10, height=4)
ggsave("results/plots/reads_coverage_plasm.png", p_plasm, dpi=300, width=8, height=4)

##6
dot_raw <- read.table("results/plots/dotplot.tsv", skip=3, header=F, fill=T)
head(dot_raw); ncol(dot_raw)  # Paste output here next time

# Robust: Assume standard 8 cols, select REF_ST (col4), QRY_ST (col6), %ID (col8)
dot <- data.frame(ref_start=dot_raw[,4], asm_start=dot_raw[,6], perc_id=as.numeric(gsub("%","",as.character(dot_raw[,8]))))

# Filter/plot (drop NA/low ID)
dot <- filter(dot, !is.na(perc_id), perc_id>80, ref_start>0)

ggplot(dot, aes(ref_start/1e6, asm_start/1e6)) + 
  geom_point(aes(color=perc_id), alpha=0.6, size=0.3) + theme_minimal() +
  labs(title="Assembly vs. Ref Dotplot", x="Ref Pos (Mb)", y="Asm Pos (Mb)", color="%ID") +
  scale_color_gradient(low="red", high="darkblue")
ggsave("results/plots/dotplot_fixed.png", width=10, height=6, dpi=300)


##7

library(ggplot2); library(dplyr); library(RcppRoll)

# Reads coverage
depth_reads <- read.table("results/plots/depth_ASM694v2.txt", col.names=c("contig","pos","depth"))
depth_reads$depth_smooth <- roll_mean(depth_reads$depth, 1001, fill=0)

chrom_r <- filter(depth_reads, contig=="NC_003197.2")
ggplot(chrom_r, aes(pos/1e6, depth_smooth)) + geom_line(color="blue") + theme_minimal() +
  labs(title="Raw Reads Alignment: Chromosome Coverage (150x, 96.8% Breadth)", x="Position (Mb)", y="Depth")
ggsave("results/plots/reads_chrom.png", width=10, height=4, dpi=300)

plasm_r <- filter(depth_reads, contig=="NC_003277.2")
ggplot(plasm_r, aes(pos/1000, depth_smooth)) + geom_line(color="blue") + theme_minimal() +
  labs(title="Raw Reads Alignment: Plasmid Coverage", x="Position (kb)", y="Depth")
ggsave("results/plots/reads_plasm.png", width=8, height=4, dpi=300)

# Assembly coverage
depth_asm <- read.table("results/plots/asm_depth.txt", col.names=c("contig","pos","depth"))
depth_asm$depth_smooth <- roll_mean(depth_asm$depth, 1001, fill=0)

chrom_a <- filter(depth_asm, contig=="NC_003197.2")
ggplot(chrom_a, aes(pos/1e6, depth_smooth)) + geom_line(color="green") + theme_minimal() +
  labs(title="Assembly Alignment: Polished Contigs vs Ref (96% Mapped)", x="Position (Mb)", y="Depth")
ggsave("results/plots/asm_chrom.png", width=10, height=4, dpi=300)0

plasm_a <- filter(depth_asm, contig=="NC_003277.2")
ggplot(plasm_a, aes(pos/1000, depth_smooth)) + geom_line(color="green") + theme_minimal() +
  labs(title="Assembly Alignment: Plasmid Contig", x="Position (kb)", y="Depth")
ggsave("results/plots/asm_plasm.png", width=8, height=4, dpi=300)

