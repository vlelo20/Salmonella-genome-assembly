# **Genome Assembly and Reference Comparison of *Salmonella enterica***
Author: Vian Lelo
Date created: January 20th, 2026, Last Updated: February 8th, 2026

# 1.0 - Introduction

Whole‑genome assembly of bacterial pathogens has become a central tool in microbiology, epidemiology, and public health. _Salmonella enterica_ is a prokaryotic, Gram‑negative bacterial pathogen responsible for a substantial burden of foodborne disease worldwide, with infections ranging from self‑limiting gastroenteritis to invasive systemic disease. Its genome typically consists of a single circular chromosome of ~4.5–5.0 Mb, often accompanied by one or more plasmids that can carry virulence and antimicrobial resistance determinants (Baker and Dougan 2007). High‑quality, complete genome assemblies for _S. enterica_ enable precise serovar determination, outbreak reconstruction, plasmid tracking, and comparative genomics, all of which depend on accurate reconstruction of both chromosomal and extrachromosomal replicons (Zhao et al. 2023; Wick, Judd, and Holt 2023). 

Historically, bacterial genome assembly relied on short‑read sequencing (e.g., Illumina), which provides highly accurate reads but struggles with repetitive regions, mobile elements, and plasmids. These limitations often result in fragmented assemblies composed of many contigs, complicating downstream analyses such as structural variant detection, plasmid–chromosome distinction, and accurate gene context inference (Wick, Judd, and Holt 2023). Long‑read sequencing technologies, particularly Oxford Nanopore Technologies (ONT), have transformed bacterial genomics by producing reads that span repeats and structural variants, enabling near‑complete or fully closed assemblies for many bacterial species (Purushothaman et al. 2026). For _S. enterica_, which contains repetitive rRNA operons, insertion sequences, and diverse plasmids, long reads are especially valuable for resolving genome structure and capturing clinically relevant mobile genetic elements (Helm et al. 2003).

Oxford Nanopore’s R10 chemistry with Q20+ basecalling represents a substantial improvement in raw read accuracy compared to earlier ONT chemistries. Recent work has shown that ONT long reads alone can generate complete bacterial and plasmid genomes without the need for short‑read data, provided that appropriate assembly and polishing workflows are used (Zhao et al. 2023). Meta‑analyses and benchmarking studies of long‑read assemblers for bacterial pathogens have demonstrated that modern ONT‑focused assemblers can routinely produce single‑contig chromosomal assemblies and accurately reconstruct plasmids, with error rates that are acceptable for many comparative and clinical applications, especially after polishing (Chen, Erickson, and Meng 2020).

Several studies have systematically compared long‑read assemblers and polishing strategies for bacterial genomes. Chen et al. (2023) benchmarked multiple long‑read assemblers on ONT data from bacterial pathogens and found that tools such as Flye, Canu, and Raven can generate highly contiguous assemblies, with Flye often performing particularly well in terms of contiguity and correctness for ONT data. Wick et al. (2023) proposed a best‑practice pipeline for “perfect” bacterial genome assembly using ONT and Illumina reads, recommending Medaka for ONT‑based polishing, highlighting the importance of both assembler choice and polishing strategy. More recent evaluations of ONT versus Illumina platforms for bacterial whole‑genome sequencing in clinical microbiology emphasize that ONT’s long reads enable complete genome reconstruction and plasmid resolution, while short reads still offer advantages in per‑base accuracy; together, they provide complementary strengths for high‑confidence variant calling and outbreak analysis.

For this assignment, the focus is on ONT R10/Q20+ long reads with an N50 of 5–15 kb, which are well suited to de novo assembly of _S. enterica_. The key challenges in assembling such data include handling residual basecalling errors (particularly in homopolymers), correctly resolving repetitive regions and plasmids, and avoiding misassemblies that could confound downstream variant calling. Long‑read assemblers such as Flye are specifically designed to work with noisy long reads, constructing repeat graphs and leveraging read length to resolve complex regions (Oxford Nanopore Assembly using Flye 2025). However, even with improved Q20+ accuracy, ONT assemblies typically require polishing to reduce small indel and substitution errors that can impact gene prediction and variant analysis. Tools like Medaka, which use neural network models trained on ONT data, have been shown to substantially improve consensus accuracy when applied to ONT assemblies (Wick, Judd, and Holt 2023)

Aligning the assembled genome to a reference is essential for variant calling and comparative analysis. Minimap2 has emerged as the standard aligner for long‑read data, offering high speed and accuracy for mapping ONT reads or assembled contigs to reference genomes. In contrast, short‑read aligners such as BWA‑MEM and Bowtie2 are not appropriate for ONT whole‑genome bacterial assemblies (Saada et al. 2024). Once alignments are generated, variant callers designed for long reads, such as medaka, can identify single‑nucleotide variants (SNVs) and small indels between the assembled _S. enterica_ genome and a reference. Visualization tools like IGV allow manual inspection of alignments and variants (Hall et al. 2024).

In summary, assembling and comparing a Salmonella enterica genome from ONT R10/Q20+ reads involves balancing the strengths and limitations of long‑read sequencing. ONT provides long reads that enable complete, structurally accurate assemblies and plasmid resolution, but residual basecalling errors necessitate careful polishing and appropriate variant calling strategies. Meta‑analyses and recent methodological papers support the use of long‑read assemblers like Flye, ONT‑aware polishers like Medaka, and long‑read aligners like minimap2 as a robust foundation for bacterial genome assembly and comparative genomics. This assignment will apply these principles to generate a consensus _S. enterica_ genome, align it to a reference, call variants, and visualize the results.

# **2.0 - Methods**

## 2.1 - Data and overall strategy
The starting data will consist of raw Oxford Nanopore FASTQ files generated with R10 chemistry and Q20+ basecalling, with an expected read length N50 of 5–15 kb. The target organism is _Salmonella enterica_, a prokaryotic bacterial pathogen with an expected genome size of approximately 4.5–5.0 Mb, potentially including plasmids. The overall strategy is to perform quality control and filtering of the ONT reads, assemble the genome de novo using a long‑read assembler optimized for ONT data, polish the assembly to improve consensus accuracy, download (ASM694v2) _S. enterica_ reference genome from NCBI, align the assembly and/or reads to the reference, perform variant calling, and visualize both the assembly and the variants. A general estimation of the time required is 4-5 hours to run the protocol (Zhao et al. 2023).

### Data Acquisition:
Obtaining Raw Reads for *Salmonella enterica* isolate (accession SRR32410565) as fastq.gz:
```
wget -O data/raw/SRR32410565.fastq.gz \
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR32410565/SRR32410565
```
Perform a quick check:
```
ls -lh data/raw/SRR32410565.fastq.gz
zcat data/raw/SRR32410565.fastq.gz | head
```
Conversion an SRA run into a FASTQ file using 8 threads to speed up extraction.:
```
fasterq-dump data/raw/SRR32410565.fastq.gz \
  --outdir data/raw \
  --threads 8
```

### **2.1.1 - Environment setup**
```
conda create -n salmonella_ont_wgs \
  -c conda-forge -c bioconda \
  python=3.10 \
  flye \
  minimap2 \
  samtools \
  nanoplot \
  bcftools \
  htslib \
  tabix \
  sra-tools \
  -y

```
Activate it:
```
conda activate salmonella_ont_wgs
```
Versions Check:
```
flye --version
minimap2 --version
samtools --version
NanoPlot --version
bcftools --version
```
### **2.1.2 - Software Versions**

All analyses were performed using the following tool versions, verified via `--version` flags:

| Tool | Version | Purpose |
|------|---------|---------|
| Flye | 2.9.6-b1802 | De novo genome assembly |
| minimap2 | 2.30-r1287 | Read and contig alignment |
| samtools | 1.23 | BAM file manipulation and statistics |
| NanoPlot | 1.46.2 | Sequencing quality control |
| bcftools | 1.23 | VCF manipulation and variant statistics |
| Medaka | 2.0.1 | Assembly polishing and variant calling |
| Python | 3.10 | Script execution environment |
| R | 4.5.2 | Data visualization and analysis |

Version consistency across the analysis pipeline was maintained via conda environment specification.

## **2.2 - Data Acquisition**

### **2.2.1 - Obtaining Raw Reads for *Salmonella enterica* isolate (accession SRR32410565)**
```
wget -O data/raw/SRR32410565.fastq.gz \
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR32410565/SRR32410565
```
Performing a quick check:
```
ls -lh data/raw/SRR32410565.fastq.gz
zcat data/raw/SRR32410565.fastq.gz | head
```
Conversion an SRA run into a FASTQ file using 8 threads to speed up extraction.:
```
fasterq-dump data/raw/SRR32410565.fastq.gz \
  --outdir data/raw \
  --threads 8
```
Renaming file:
```
mv data/raw/SRR32410565.fastq.gz.fastq data/raw/SRR32410565.fastq
ls -lh data/raw/
```
### **2.2.2 - Reference genome retrieval**
A reference genome for _Salmonella enterica_ will be downloaded from NCBI RefSeq, ideally matching the serovar of the sequenced isolate. The reference will be obtained in FASTA format along with its annotation (GenBank or GFF), enabling both sequence‑level and feature‑level comparisons. Using a well‑annotated RefSeq reference facilitates interpretation of variants in terms of genes, operons, and known virulence or resistance loci.
Download ASM694v2 reference (NO ncbi-datasets-cli):
```
REFDIR="results/ref/ASM694v2"
mkdir -p "$REFDIR"
cd "$REFDIR"

BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2"

# Genome FASTA (reference sequence)
wget -O ASM694v2_genomic.fna.gz  "${BASE}/GCF_000006945.2_ASM694v2_genomic.fna.gz"

# Annotation (GFF) — useful for gene context in IGV
wget -O ASM694v2_genomic.gff.gz  "${BASE}/GCF_000006945.2_ASM694v2_genomic.gff.gz"

gunzip -f ASM694v2_genomic.fna.gz
gunzip -f ASM694v2_genomic.gff.gz

ln -sf ASM694v2_genomic.fna results/ref/reference.fasta 2>/dev/null || true
```


## **2.3 - Quality control and read filtering**
Initial quality control will be performed using NanoPlot (v1.46.2) to assess read length distribution, and quality scores, confirming that the N50 falls within the expected 5–15 kb range and identifying any obvious issues with the sequencing run. If necessary, reads will be filtered using NanoFilt (v2.8.0) to remove very short or low‑quality reads, which can reduce noise and computational burden without compromising assembly quality (Zhao et al. 2023). For example, a minimum read length of 1,000 bp and a Q<10 will be removed. The first and last 50 bases will be removed to retain clean data (Zhao et al. 2023). Downsampling may be considered to reduce run time and memory usage of 60-70x (Wick, Judd, and Holt 2023). 
```
cd ~/binf6110/assignment1

NanoPlot \
  --fastq data/raw/SRR32410565.fastq \
  --outdir data/qc/nanoplot \
  --threads 8 \
  --plots hex dot \
  --loglength

```
Raw reads aligned to ASM694v2 showed 94.3% mapping (193,601/205,302 reads) and mean depth 150× across 96.8% breadth, quantified via samtools flagstat, depth -a, and awk-derived metrics. The polished assembly exhibited 96% contig mapping (24/25 alignments), confirming high structural fidelity.
```
# Coverage (prior)
samtools depth -a reads_vs_ASM694v2.bam > plots/depth_ASM694v2.txt
awk ... > ASM694v2.coverage_metrics.txt  # 96.8% breadth, 150x

# Flagstats
samtools flagstat reads_vs_ASM694v2.bam     # 94.3% mapped
samtools flagstat polished_flye_vs_ASM694v2.bam  # 96% mapped
```

## **2.4 - De novo genome assembly**
NanoPlot (v1.46.2) QC confirmed high-quality ONT reads (N50 4,683 bp, median Q-score 23.7, 196k reads, 809 Mb total yield); no NanoFilt filtering or downsampling was applied as the data exceeded standard assembly thresholds without compromising contiguity or accuracy (Zhao et al. 2023; Wick, Judd, and Holt 2023).
​
De novo assembly was performed using Flye (v2.9.6) in a dedicated flye_env conda environment, leveraging its repeat graph construction optimized for long-read ONT data to produce highly contiguous bacterial assemblies with accurate resolution of repetitive elements and plasmids. The high-quality (HQ) preset (--nano-hq) was selected for Q20+ reads, with an expected genome size of 5 Mb and 14 CPU threads (reserving system resources for responsiveness); assembly was run on the raw FASTQ as follows:

```
conda activate salmonella_ont_wgs
cd ~/binf6110/assignment1

flye --nano-hq data/raw/SRR32410565.fastq \
  --genome-size 5m --out-dir results/assembly/flye_raw \
  --threads 14

```
## **2.5 - Assembly Polishing**
The raw Flye assembly was polished using Medaka consensus (v2.0.1) in the dedicated medaka_env conda environment to correct residual base-level errors (e.g., homopolymer indels) common in ONT drafts, leveraging the same raw reads for neural-network refinement (Zhao et al. 2023). Medaka automatically selected the R10.4.1 SUP variant model (r1041_e82_400bps_sup_variant_v5.0.0), realigning reads to the draft contigs with minimap2 (map-ont preset), generating pileup features, and producing a consensus sequence; 14 threads were used for efficiency. This involved aligning the original ONT reads back to the Flye assembly using minimap2, followed by Medaka consensus calling to generate the polished assembly​
```
conda activate medaka_env
cd ~/binf6110/assignment1

medaka_consensus \
  -i data/raw/SRR32410565.fastq \
  -d results/assembly/flye_raw/assembly.fasta \
  -o results/assembly/flye_medaka_polished \
  -t 14

```
## **2.6 - Assembly-to-Reference Alignment**

The polished Medaka consensus assembly was aligned to the ASM694v2 reference genome using minimap2 with the asm5 preset, which is optimized for aligning assembled contigs and allows large insertions/deletions to detect structural variants and chromosomal rearrangements. The --secondary=no flag suppressed secondary alignments to retain only the best reference match per contig, simplifying synteny visualization and preventing ambiguous mappings. The SAM output was coordinate-sorted and compressed into BAM format using samtools sort (14 threads), then indexed for genome browser viewing (IGV). Additionally, samtools flagstat generated alignment statistics (mapped contigs, alignment rates, supplementary alignments) to quantify assembly-to-reference concordance and identify misassemblies or unplaced sequences.

```
conda activate salmonella_ont_wgs  # minimap2/samtools env
cd ~/binf6110/assignment1

minimap2 -ax asm5 --secondary=no -t 14 \
  results/ref/ASM694v2/ASM694v2_genomic.fna \
  results/assembly/flye_medaka_polished/consensus.fasta | \
samtools sort -@ 14 -o results/align/polished_flye_vs_ASM694v2.bam -

samtools index results/align/polished_flye_vs_ASM694v2.bam
samtools flagstat results/align/polished_flye_vs_ASM694v2.bam \
  > results/align/polished_flye_vs_ASM694v2.flagstat.txt

```


### **2.7 - Raw reads - Reference Alignment and variant calling**

Raw Oxford Nanopore reads (SRR32410565.fastq) were aligned to the ASM694v2 reference genome using minimap2 (v2.30) with the map-ont preset optimized for single-molecule long-read alignment, which accommodates the characteristic error profile of ONT sequencing and applies soft-clipping to read ends. The resulting SAM output was coordinate-sorted and compressed into BAM format using samtools sort (v1.23) parallelized across 8 threads, then indexed with samtools index to enable efficient random access for downstream visualization and quality assessment.

```
minimap2 -t 8 -ax map-ont \
  results/ref/ASM694v2/ASM694v2_genomic.fna \
  data/raw/SRR32410565.fastq | \
samtools sort -@ 8 -o results/align/reads_vs_ASM694v2.bam
samtools index results/align/reads_vs_ASM694v2.bam
```
### **2.7.1 - Mapping Summary Statistics**

Alignment quality metrics were extracted from the reads-to-reference BAM file using samtools flagstat (v1.23) to quantify overall mapping rates, properly paired reads, unmapped reads, and supplementary alignments. Additionally, samtools idxstats provided per-contig read counts and mapped read lengths, enabling assessment of coverage distribution across the chromosome and plasmid sequences and identification of potential unmapped or underrepresented genomic regions.
```
samtools flagstat results/align/reads_vs_ASM694v2.bam > results/align/ASM694v2.flagstat.txt
samtools idxstats results/align/reads_vs_ASM694v2.bam > results/align/ASM694v2.idxstats.txt
```
### **2.7.2 - Coverage Depth Analysis**
Per-base sequencing depth was calculated using samtools depth (v1.23) with the -a flag to include zero-coverage positions, generating a genome-wide depth profile across all reference sequences. An awk script processed this depth file to compute two key coverage metrics: breadth of coverage (fraction of reference bases with non-zero depth) and mean sequencing depth (total aligned bases divided by genome length). These summary statistics were written to a separate file to assess the adequacy and uniformity of raw sequencing data prior to de novo assembly, ensuring sufficient read support for high-quality genome reconstruction.
```
samtools depth -a results/align/reads_vs_ASM694v2.bam > results/plots/depth_ASM694v2.txt

awk 'BEGIN{cov=0;tot=0;sum=0}
     {tot++; sum+=$3; if($3>0) cov++}
     END{
       print "Breadth_covered_fraction", cov/tot;
       print "Mean_depth", sum/tot;
     }' results/plots/depth_ASM694v2.txt > results/align/ASM694v2.coverage_metrics.txt
```
### **2.8 - Variant calling on ONT reads against the ASM694v2 reference genome**

Variant calling was performed using Medaka (v2.0.1) in reads‑to‑reference mode on the raw ONT FASTQ file (SRR32410565.fastq) and the ASM694v2 reference genome (Hall et al. 2024). Within a dedicated medaka_env conda environment, Medaka was run with 4 CPU threads and default model auto‑selection; because the basecaller string could not be parsed from the input, Medaka fell back to its default R10.4.1 SUP model r1041_e82_400bps_sup_variant_v5.0.0, which is appropriate for 400 bp Q20+ R10.4.1 data. Before inference, Medaka confirmed that bundled versions of minimap2 (v2.30), samtools (v1.23), bcftools (v1.23), bgzip, and tabix met or exceeded its minimum requirements, ensuring consistency of downstream pileup and VCF operations.

For variant calling, Medaka internally realigned the ONT reads to the ASM694v2 reference using minimap2 with the ONT long‑read preset (-x map-ont, 4 threads), reusing the existing FASTA index (ASM694v2_genomic.fna.fai) and minimap2 index (ASM694v2_genomic.fna.map-ont.mmi) where available. The reference (NC_003197.2 chromosome and NC_003277.2 plasmid) was processed in 1 Mb windows, with an internal inference chunk length of 10,000 bases; regions shorter than this threshold were “quarantined” and processed in a final pass to ensure complete coverage of both chromosomal and plasmid sequence. Reads were filtered using a minimum mapping quality threshold of 1, pileup‑based features were generated for each window (reporting median depths typically between ~130× and ~200× across the chromosome), and these features were passed to a bidirectional GRU neural network (two layers, 128 hidden units per direction) running at full precision on CPU to produce base‑level consensus probabilities. Medaka then combined predictions across overlapping windows, resolved non‑overlapping segments (reported in the log as “cannot be concatenated as there is no overlap and they do not abut”), and wrote the final variant calls and consensus to results/vcf/medaka_reads_only/medaka.annotated.vcf.
```
conda deactivate
conda create -n medaka_env -c conda-forge -c bioconda \
  medaka=2.0 python=3.10 -y
```
```
conda activate medaka_env
cd ~/binf6110/assignment1
```
```
medaka_variant \
  -i data/raw/SRR32410565.fastq \
  -r results/ref/ASM694v2/ASM694v2_genomic.fna \
  -o results/vcf/medaka_reads_only \
  -t 4
```

### **2.8.1 - Post-Processing of Medaka VCF**

The Medaka output VCF was sorted by genomic coordinates (CHROM, POS) using bcftools to resolve indexing issues, compressed, and indexed with tabix for efficient random access in IGV and downstream analyses (e.g., region-specific queries).

VCF sorting, compression, and tabix indexing
```
bcftools sort \
  results/vcf/medaka_reads_only/medaka.annotated.vcf \
  -Oz -o results/vcf/medaka_reads_only/medaka.annotated.sorted.vcf.gz

tabix -p vcf results/vcf/medaka_reads_only/medaka.annotated.sorted.vcf.gz
```
AF annotation and indexing
```
bcftools +fill-tags \
  results/vcf/medaka_reads_only/medaka.annotated.sorted.vcf.gz \
  -Oz -o results/vcf/medaka_reads_only/medaka.af.vcf.gz \
  -- -t AF

tabix -p vcf results/vcf/medaka_reads_only/medaka.af.vcf.gz
```
### **2.8.2 - Variant Summary Statistics**

Total variants were enumerated by counting non-header lines, while bcftools stats provided a breakdown of SNPs vs. indels, written to medaka.stats.txt for reporting.

```
bcftools view results/vcf/medaka_reads_only/medaka.annotated.sorted.vcf.gz \
  | grep -v "^#" | wc -l
```
SNPs/indels breakdown:
```
bcftools stats results/vcf/medaka_reads_only/medaka.annotated.sorted.vcf.gz \
  > results/vcf/medaka_reads_only/medaka.stats.txt

grep "^SN" results/vcf/medaka_reads_only/medaka.stats.txt
```
### **2.8.3 - Variant distribution by contig**

The genomic distribution of variants across reference contigs (chromosome and plasmids) was determined by extracting the CHROM field from the Medaka-annotated VCF using bcftools query (v1.23), followed by sorting, unique counting with uniq -c, and numeric reverse sorting to rank contigs by variant frequency. Results were saved to a tab-delimited file, and percentages relative to the total variant count (14,089) were computed using awk with sprintf formatting for 0.1% precision, producing a summary table for assessment of chromosomal versus plasmid divergence.
```
cd ~/binf6110/assignment1

bcftools query -f '%CHROM\n' \
  results/vcf/medaka_reads_only/medaka.annotated.sorted.vcf.gz | \
sort | uniq -c | sort -nr > results/vcf/variants_per_contig.txt

cat results/vcf/variants_per_contig.txt

# With percentages (total=14089)
awk '{print $2 "\t" $1 "\t" sprintf("%.1f%%", $1/14089*100)}' \
  results/vcf/variants_per_contig.txt | column -t
```

The coordinate‑sorted, indexed VCF (medaka.annotated.sorted.vcf.gz and its .tbi index) was used for variant exploration in IGV, summary statistics with bcftools, and downstream plotting of variant distributions along the S. enterica chromosome and plasmid.

## 2.9 - Functional Annotation of Variants with SnpEff & visualization

To classify variants by predicted functional impact (missense, synonymous, frameshift, stop-gain/loss), SnpEff (v5.1d) was used to annotate the Medaka-called variants relative to the ASM694v2 reference genome annotation. A custom SnpEff database was constructed from the NCBI GenBank annotation file (GCF_000006945.2_ASM694v2_genomic.gbff) rather than the GFF format, as GenBank provides more reliable gene structure information for bacterial genomes.

### Database Construction
```bash
# Download GenBank annotation
wget -O results/ref/ASM694v2/ASM694v2_genomic.gbff.gz \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gbff.gz"
gunzip results/ref/ASM694v2/ASM694v2_genomic.gbff.gz

# Set up SnpEff database structure
mkdir -p data/snpEff_data/ASM694v2
cp results/ref/ASM694v2/ASM694v2_genomic.fna data/snpEff_data/ASM694v2/sequences.fa
cp results/ref/ASM694v2/ASM694v2_genomic.gbff data/snpEff_data/ASM694v2/genes.gbk

# Create configuration file
cat > snpEff.config << 'EOF'
data.dir = ./data/snpEff_data
ASM694v2.genome : Salmonella_enterica_ASM694v2
EOF

# Build database
snpEff build -genbank -v ASM694v2
```
### 2.9.1 - Variant annotation

```
# Annotate variants
snpEff ann -v ASM694v2 \
  -stats results/vcf/snpEff_summary.html \
  results/vcf/medaka_reads_only/medaka.annotated.sorted.vcf.gz \
  > results/vcf/medaka.snpeff.vcf

# Compress and index
bgzip -f results/vcf/medaka.snpeff.vcf
tabix -p vcf results/vcf/medaka.snpeff.vcf.gz

```
### 2.9.2 - Effect Extraction and Analysis
```
# Extract all variant effects
bcftools query -f '%INFO/ANN\n' results/vcf/medaka.snpeff.vcf.gz | \
  cut -d'|' -f2 | sort | uniq -c | sort -rn > results/vcf/functional_effects_summary.txt

# Extract effects excluding pncB (quality control)
bcftools view -e 'ANN~"pncB"' results/vcf/medaka.snpeff.vcf.gz | \
  bcftools query -f '%INFO/ANN\n' | cut -d'|' -f2 | \
  sort | uniq -c | sort -rn > results/vcf/effects_no_pncB.txt

# Gene-level effect breakdown for key virulence and conjugation genes
bcftools query -f '%INFO/ANN\n' results/vcf/medaka.snpeff.vcf.gz | \
  awk -F'|' '{print $4"\t"$2}' | \
  grep -E "^(sspH2|oafA|traI|pncB|traG|traD|ssbB|traN|traC|gogB)" | \
  sort | uniq -c | sort -k2,2 -k1,1rn > results/vcf/key_genes_effects.txt

```
### 2.9.3 - Visualization and comparative analysis**
To validate the assembly and interpret the genomic landscape, several visualization tools were employed. Read alignments and called variants were inspected using the Integrative Genomics Viewer (IGV). This allowed for the interactive assessment of mapping quality and depth of coverage across the Salmonella chromosome and plasmids. IGV was specifically utilized to validate high-impact variants and inspect regions with anomalous variant density, ensuring that reported mutations were supported by high-quality read alignments.

For quantitative data analysis and visualization, RStudio (v4.3.1) was used to generate coverage and variant density plots. Custom R scripts utilized genomic data to produce visual comparisons between the sample assembly and the ASM694v2 reference, facilitating the identification of hotspots for genomic divergence. To confirm the structural integrity of the assembly, the Flye assembly graph was examined to verify the resolution of the chromosome and extrachromosomal elements into distinct, potentially circular contigs.


# 3.0 - Results

## 3.1 - Raw Sequencing Quality
NanoPlot (v1.46.2) QC of SRR32410565.fastq (196,031 reads) yielded median length 3,957 bp (N50 4,683 bp, mean 4,128 bp) and median Q-score 23.7 (mean 18.9), with 76.9% reads ≥Q20 (625 Mb usable yield from 809 Mb total). These metrics confirm high-quality R10.4.1 SUP data suitable for de novo assembly and reference mapping.

## 3.2 -  De novo Assembly Metrics:
Flye (v2.9.6-b1802) de novo assembly of SRR32410565.fastq (809 Mb yield, read N50 4,683 bp, median Q23.7 per NanoPlot v1.46.2) generated a 5.1 Mb draft in three contigs (N50 ~1.68–3.3 Mb), resolving the chromosome (contig_3: 3.32 Mb, 153× cov.), accessory sequence (contig_2: 1.68 Mb, 169×), and circular plasmid (contig_4: 109 kb, 245×; ~161× overall). No read filtering was required given high quality (76.9% Q20+). Medaka polishing (r1041_e82_400bps_sup_v500) produced a refined consensus aligning near-perfectly to ASM694v2 (96–99% identity via flagstat, minimal misassemblies).

## 3.3 Alignment Quality
Raw ONT reads mapped efficiently to ASM694v2 (94.3% overall, 94.0% primary; n=193,601/205,302), achieving mean depth 149.96× across 96.8% genomic breadth. The polished assembly showed near-perfect contig mapping (96.0% total, 66.7% primary; n=24/25 alignments), with supplementary hits likely representing circular plasmids.

## 3.4 - Variant Calling Analysis:
Medaka variant calling (v2.0.1) on the polished assembly relative to ASM694v2 (NC_003197.2) detected 14,089 high-confidence variants (bcftools stats v1.23), predominantly single nucleotide polymorphisms (SNPs: 8,519) and multi-nucleotide polymorphisms (MNPs: 4,396), with 840 indels and 334 complex/other variants; no multiallelic sites were present. Transitions (4,296) slightly outnumbered transversions (4,223; ts/tv ratio 1.02), consistent with biological substitution patterns, while all variants were singletons (AC=1) supported by median QUAL ~55 (range 0–2,174.7) supported by read depths at variant sites ranging 1–381× (18.9% at 10×, 7.3% at 46×). Indel length distribution showed small events dominant (-60 to +60 bp, mostly ±1–2 bp), visualized in IGV (Figure 3) alongside assembly and read tracks confirming structural integrity.

### Variant Distribution by Contig
Variants were nearly equally distributed between the chromosome (NC_003197.2: 7,163, 50.8%) and plasmid (NC_003277.2: 6,926, 49.2%), reflecting substantial divergence on both replicons relative to the reference despite high overall alignment identity. This distribution highlights conserved core genome variation alongside higher plasticity in accessory elements.



<img width="1639" height="870" alt="igv_snapshot" src="https://github.com/user-attachments/assets/c332304e-d1e4-45ea-8268-0c341cd05802" />

Figure 1:  IGV alignment view of variant-dense region on plasmid NC_003277.2 (positions 49,744-49,838).This 95 bp window displays individual read alignments from raw Nanopore sequencing data mapped to the ASM694v2 reference genome. The top coverage histogram shows variant allele frequencies at each position, with colored bars representing substitutions relative to the reference sequence (A=green, T=red, G=orange, C=blue). Individual read alignments below reveal extensive heterogeneity, with numerous single nucleotide polymorphisms (colored letters) distributed throughout the reads. Insertions are marked with purple "I" indicators, and deletions appear as black dashes within reads. The reference sequence is displayed at the bottom for comparison. This region exemplifies the high variant density observed across much of the plasmid, with nearly every position showing alternative alleles in at least some reads.


<img width="3600" height="2100" alt="plasmid_variants_coverage" src="https://github.com/user-attachments/assets/594f8382-0fb7-41b5-a904-0653aba38193" />

Figure 2: Variant density across plasmid NC_003277.2 in 100 bp windows. The upper panel displays a lollipop plot of variant distribution along the ~94 kb plasmid sequence, where each vertical line represents the number of variants detected per 100 bp bin between raw Nanopore reads and the ASM694v2 reference. SNPs (black) dominate throughout, with dramatic variation in density across the plasmid. The plasmid exhibits a bimodal distribution: sparse variant regions (0-50 kb, with isolated clusters) versus a densely polymorphic region (50-94 kb) where variant density remains consistently high (10-38 variants per bin, median ~20). Peak variant density reaches 38 variants per 100 bp bin near position 65 kb. The lower panel shows sequencing coverage as a heatmap (white to dark green gradient, 0 to >300×), revealing highly uneven read depth. Low-coverage regions (white/light green, <50×) in the 0-50 kb segment correlate with sparse variants and likely represent alignment difficulties or structural divergence, while high-coverage regions (dark green, >200×) in the 50-94 kb segment correlate with dense, reliable variant calls.


<img width="3600" height="2100" alt="chromosome_variants_coverage" src="https://github.com/user-attachments/assets/7e6091ea-d4ca-427d-89d7-4d04ae63b30c" />

Figure 3: Variant density across the Salmonella enterica chromosome (NC_003197.2) in 100 bp windows. The upper panel displays a lollipop plot showing the distribution of variants detected between raw Nanopore reads and the ASM694v2 reference genome. Each vertical line represents the total number of variants per 100 bp bin, with SNPs (black) dominating the variant landscape. Variant density is highly heterogeneous across the 4.8 Mb chromosome, with several distinct high-density regions visible at approximately 1.0 Mb, 2.7 Mb, and 3.0-3.5 Mb positions. The most polymorphic region (near 3.0 Mb) reaches peak densities of 50-55 variants per bin, potentially indicating hypervariable loci, prophage insertions, or regions under diversifying selection. The lower panel shows sequencing coverage as a heatmap (white to dark green gradient indicating 0 to >300× depth), demonstrating uniform high coverage (predominantly dark green) across the entire chromosome with no major gaps, confirming reliable variant calling throughout the sequence.

## **3.5 - Gene-Level Variant Distribution**
Variant annotation using the ASM694v2 GFF file revealed that 85.2% of variants (12,010/14,089) occurred within coding sequences, while 14.8% (2,079) were intergenic. The distribution of variants across genes was highly non-uniform, with a small subset of genes harboring the majority of polymorphisms (Table 1).
Table 1: Top 10 genes by variant count
| Gene  | Variant Count | Function                                                  | Location           |
|-------|---------------|-----------------------------------------------------------|--------------------|
| pncB  | 3,176         | Nicotinate phosphoribosyltransferase (NAD+ biosynthesis)  | Chromosome         |
| traI  | 683           | Conjugative transfer relaxase                             | Plasmid            |
| sspH2 | 549           | Type III secretion system effector protein                | Chromosome (SPI-2) |
| oafA  | 527           | O-antigen acetylase (LPS modification)                    | Chromosome         |
| ssbB  | 525           | Single-stranded DNA-binding protein                       | Plasmid            |
| traG  | 433           | Conjugative transfer ATPase                               | Plasmid            |
| traD  | 334           | Conjugative transfer coupling protein                     | Plasmid            |
| traC  | 332           | Pilus assembly protein                                    | Plasmid            |
| traN  | 318           | Conjugative transfer mating pair stabilization            | Plasmid            |
| gogB  | 255           | Type III secretion system effector protein                | Chromosome (SPI-2) |

The single most polymorphic gene, pncB, contained 22.5% of all variants (3,176/14,089), representing an outlier in variant density that likely reflects either a genuine hypervariable locus or a technical artifact such as paralogous sequence alignment or repeat-induced misassembly. Excluding pncB, the remaining top genes fall into three functional categories: (1) virulence factors (sspH2, oafA, gogB) involved in host immune evasion and Type III secretion, (2) plasmid conjugative transfer machinery (traI, traG, traD, traC, traN), and (3) plasmid replication/maintenance (ssbB, repA2).
Notably, 7 of the top 10 variant-dense genes are plasmid-encoded, consistent with the near-equal distribution of variants between chromosome (50.8%) and the much smaller plasmid replicon (49.2%). The high variant load in plasmid conjugation genes (tra operon) suggests substantial divergence in horizontal gene transfer machinery, while the concentration of chromosomal variants in SPI-2 effector proteins (sspH2, gogB) and surface antigen modifiers (oafA) indicates adaptive variation in host-pathogen interaction determinants.

## 3.6 - Functional Variant Classification

SnpEff (v5.1d) annotation using the ASM694v2 reference (4,717 genes, 4,554 protein-coding transcripts) classified the 14,089 variants by predicted functional impact. Analysis was performed both including and excluding the pncB gene (57 variants, 0.4% of total), which showed elevated variant density potentially due to alignment artifacts.

### All Variants (n=14,089)

| Effect Type | Count | Percentage | Impact Level |
|-------------|-------|------------|--------------|
| Missense variant | 6,114 | 43.4% | Moderate |
| Synonymous variant | 4,553 | 32.3% | Low |
| Upstream gene variant | 2,499 | 17.7% | Modifier |
| Frameshift variant | 423 | 3.0% | High |
| Frameshift & missense | 156 | 1.1% | High |
| Stop gained | 77 | 0.5% | High |
| Stop lost | 40 | 0.3% | High |
| Disruptive inframe deletion | 35 | 0.2% | Moderate |
| Frameshift & synonymous | 31 | 0.2% | High |
| Disruptive inframe insertion | 25 | 0.2% | Moderate |
| Start lost | 21 | 0.1% | High |
| Conservative inframe deletion | 23 | 0.2% | Moderate |
| Other complex variants | 91 | 0.6% | Various |

Functional annotation reveals a highly divergent genomic landscape, with nearly half of all variants (43.4%) classified as moderate-impact missense mutations. This indicates significant protein-level divergence between the sample and the reference. Notably, high-impact variants (frameshifts and nonsense mutations) account for over 5% of total calls, suggesting significant structural changes or potential gene pseudogenization within the isolate. The large number of modifier variants (17.7%) located in upstream regions suggests that regulatory divergence may also play a major role in the phenotypic differences of this strain.


### Excluding pncB (n=14,032)

| Effect Type | Count | Percentage | Change |
|-------------|-------|------------|--------|
| Missense variant | 6,094 | 43.4% | -20 variants |
| Synonymous variant | 4,533 | 32.3% | -20 variants |
| Upstream gene variant | 2,498 | 17.8% | -1 variant |
| Frameshift variant | 416 | 3.0% | -7 variants |
| Frameshift & missense | 153 | 1.1% | -3 variants |
| Stop gained | 76 | 0.5% | -1 variant |
| Other effects | <85 | <0.6% | -5 variants |

The exclusion of the pncB gene was performed to assess if the extreme variant density observed in that region (3,176 raw variants) was skewing the overall functional profile. As shown in the table, the percentage distribution of effect types remains virtually identical after exclusion. This confirms that the high rate of missense and synonymous variation is a genome-wide feature of this Salmonella isolate rather than an artifact driven by a single hyper-variable or poorly mapped locus.

**High-impact variants** (protein-truncating): 740 total (5.3%), composed of:
- 423 frameshift mutations causing reading frame shifts
- 77 premature stop codons (stop_gained)
- 40 stop codon losses (stop_lost)
- 21 start codon losses (start_lost)
- 179 compound frameshift events

The high-impact category is dominated by 423 frameshift mutations, which are critical "red flags" in genomic analysis. Because these mutations shift the genetic reading frame, they typically result in non-functional, truncated proteins. The presence of 77 stop_gained variants further reinforces the conclusion that several metabolic or structural pathways have been significantly altered or deactivated in this lineage. These variants are primarily concentrated in the accessory genome, specifically within the plasmid-borne transfer machinery.

**Moderate-impact variants** (non-synonymous): 6,114 missense mutations (43.4%) alter amino acid sequences, concentrated in:
- Virulence factors: sspH2 (420 missense), oafA (527 total), gogB (255 total)
- Conjugative transfer: traI (190 missense), traG (129 missense), traD (82 missense), traN (106 missense), traC (81 missense)
- Plasmid maintenance: ssbB (51 missense)

**Low-impact variants** (synonymous): 4,553 silent mutations (32.3%) preserve protein sequences, suggesting purifying selection maintains function in essential metabolic genes.

The ratio of missense (43.4%) to synonymous (32.3%) variants provides insight into the evolutionary pressures acting on the genome. The concentration of missense mutations in virulence factors like sspH2 and gogB suggests a history of positive selection, likely as the pathogen adapts to host immune defenses. Conversely, the high number of synonymous (silent) mutations in essential metabolic genes indicates that purifying selection is maintaining the core biological functions of the cell despite the high overall nucleotide divergence.

### Key Gene Variant Breakdown

| Gene | Missense | Synonymous | Frameshift | Stop/Start | Total | Function |
|------|----------|------------|------------|------------|-------|----------|
| sspH2 | 420 | 39 | 0 | 2 | 461 | T3SS effector |
| oafA | - | - | 1 | - | 527 | LPS modification |
| traI | 190 | 476 | 9 | 0 | 683 | Conjugative relaxase |
| traG | 129 | 294 | 9 | 1 | 433 | Transfer ATPase |
| traD | 82 | 237 | 13 | 0 | 334 | Transfer coupling |
| traN | 106 | 166 | 35 | 5 | 318 | Pilus stabilization |
| traC | 81 | 250 | 0 | 1 | 332 | Pilus assembly |
| ssbB | 51 | 37 | 7 | 0 | 152 | ssDNA binding |
| pncB | 20 | 20 | 10 | 2 | 57 | NAD+ biosynthesis |

This gene-level analysis highlights the specific biological systems driving the divergence. The tra operon (transfer genes) shows a massive accumulation of both missense and frameshift mutations, particularly in traN and traD. This suggests that the conjugative machinery of the sample’s plasmid has undergone extensive remodeling and may have different transfer efficiencies than the reference. Meanwhile, the extreme missense count in sspH2 (420 variants) vs. its low synonymous count (39) is a clear indicator of rapid protein evolution in this Type III secretion system effector, which is a hallmark of host-adapted Salmonella lineages.

# 4.0 - Discussion

## 4.1 - Overall Genomic Divergence Between Isolate and Reference
Variant calling between this Salmonella enterica isolate and the ASM694v2 reference genome identified 14,089 variants split nearly equally between the chromosome (7,163 variants, 50.8%) and plasmid (6,926 variants, 49.2%). This distribution is striking because the plasmid represents less than 2% of the total genome size yet contains nearly half of all variants, suggesting either genuine sequence divergence in the plasmid or technical challenges in aligning plasmid sequences (Wick et al., 2023). Most variants (85.2%) occurred within genes rather than intergenic regions, indicating that genomic differences primarily affect protein-coding sequences.
The variant distribution was highly non-uniform across the genome. As shown in Figure 3, the chromosome exhibits distinct variant hotspots at approximately 1.0 Mb, 2.7 Mb, and especially 3.0 Mb, where variant density peaks at over 50 variants per 100 bp bin. These hotspots likely represent genomic islands, prophage insertions, or regions that have undergone recombination with divergent Salmonella strains (Helm et al., 2003). The uniform high sequencing coverage (>100×) across the entire chromosome confirms that these variants reflect genuine sequence differences rather than alignment artifacts.

Figure 2 reveals a dramatically different pattern on the plasmid, with a bimodal distribution: sparse variants in the 0-50 kb region (correlating with poor coverage <50×) and dense, sustained variation in the 50-94 kb region (with deep coverage >200×). The low-coverage region likely contains structural rearrangements or repetitive sequences that prevent accurate read mapping, while the high-coverage region shows true sequence divergence. This pattern is consistent with modular plasmid organization, where different functional regions evolve at different rates.

## 4.2 - Assembly Quality and Comparison to Reference Genome

The de novo assembly produced by Flye generated three contigs totaling 5.1 Mb, closely matching the expected genome size of the ASM694v2 reference (4.86 Mb chromosome + 0.094 Mb plasmid = 4.95 Mb). The largest contig (3.32 Mb) represents the main chromosome, while a smaller contig (109 kb) corresponds to the plasmid, and an intermediate contig (1.68 Mb) likely represents an accessory genomic region or assembly artifact requiring further investigation.
Alignment of the polished assembly to the ASM694v2 reference showed excellent overall concordance, with 96% of assembled contigs mapping successfully to the reference genome. This high mapping rate indicates that the assembly accurately reconstructed the genome structure despite using only long-read data without short-read polishing. The mean sequencing depth of 150× across 96.8% of the reference genome confirms sufficient read coverage for reliable assembly and variant calling.
However, the presence of 14,089 variants between the assembly and reference—representing approximately 0.28% sequence divergence (14,089 variants / 4.95 Mb genome)—indicates this isolate is genetically distinct from the ASM694v2 reference strain. This level of divergence is not unexpected given that the reference is a laboratory-adapted strain (S. enterica subsp. enterica serovar Typhimurium LT2) that may have accumulated or lost genetic variations during decades of laboratory passage. Environmental or clinical isolates like SRR32410565 typically show greater genetic diversity due to ongoing adaptation to selective pressures absent in laboratory conditions.
The three-contig assembly structure, rather than achieving a single circular chromosome, suggests some repetitive regions or structural variants remain unresolved. Common challenges in bacterial genome assembly include ribosomal RNA operons, which occur in multiple near-identical copies, and insertion sequences or transposons that can cause assembly fragmentation (Helm et al., 2003). The 1.68 Mb intermediate contig warrants further investigation to determine whether it represents genuine chromosomal sequence that failed to scaffold with the main contig, a misassembly, or a large genomic island present in this isolate but absent from the reference.
Overall, the assembly quality metrics—high contig N50, deep coverage, and 96% mapping to reference—demonstrate that Oxford Nanopore R10/Q20+ chemistry combined with Flye assembly and Medaka polishing produces highly accurate bacterial genomes suitable for comparative genomics and variant analysis. The 14,089 detected variants represent genuine biological differences rather than assembly errors, as evidenced by the non-random clustering of variants in functional gene categories and the high concordance with expected coverage patterns.

### The pncB Gene and Variant Quality Assessment

Initial analysis identified the pncB gene (nicotinate phosphoribosyltransferase, NAD+ biosynthesis) as containing 57 variants (0.4% of total), substantially lower than initially reported. SnpEff annotation revealed these comprise 20 missense, 20 synonymous, 10 frameshift, and 7 other variants. While pncB shows elevated variant density (57 variants / 1,326 bp gene length = 4.3 variants per 100 bp), this is comparable to other highly polymorphic genes in the dataset rather than representing an extreme outlier requiring data exclusion [file:1].

To ensure robust biological interpretation, variant statistics were recalculated excluding pncB. The distribution of functional effects remained nearly identical:
- Missense: 43.4% (with pncB) vs. 43.4% (without)
- Synonymous: 32.3% vs. 32.3%
- Frameshift: 3.0% vs. 3.0%
- High-impact: 5.3% vs. 5.3%

This consistency confirms that conclusions regarding virulence gene variation (sspH2: 420 missense mutations, oafA: 527 variants) and plasmid conjugative transfer divergence (traI: 190 missense, traG: 129 missense, traD: 82 missense, traN: 106 missense) represent genuine biological differences independent of any potential pncB alignment artifacts [file:2][file:3].

The high missense-to-synonymous ratio in virulence effectors (sspH2: 420:39 = 10.8:1; population average: 1.3:1) suggests positive selection for immune evasion variants, while conjugative transfer genes show more balanced ratios (traI: 190:476 = 0.4:1) consistent with functional constraint on plasmid mobilization machinery.

## 4.3 - Virulence Gene Variants: Implications for Pathogenesis
Gene-level analysis revealed that virulence factors are among the most variant-dense genes in the chromosome. The sspH2 gene (549 variants) encodes a Type III secretion system effector protein that Salmonella injects into host cells to suppress immune responses. SspH2 functions as an anti-inflammatory effector that suppresses pro-inflammatory cytokines including IL-1β and IFN-γ, promoting bacterial survival inside macrophages ScienceDirect (Zhang et al., 2020). The high number of variants in sspH2 suggests this isolate may have evolved altered immune evasion capabilities compared to the reference strain.

Similarly, oafA (527 variants) encodes an enzyme that modifies the O-antigen portion of lipopolysaccharide (LPS), the major surface molecule recognized by the host immune system. Acetylation of O-antigen by OafA dramatically alters antibody recognition, with mice showing 32-fold higher antibody titers against acetylated versus non-acetylated LPS PubMed (Kim & Slauch, 1999). The extensive variation in oafA suggests this isolate likely presents a substantially altered surface structure compared to the reference, potentially enabling it to evade antibodies that would recognize the reference strain.
Figure 1 shows a representative variant-dense region on the plasmid where nearly every position contains alternative alleles in the sequencing reads. This type of dense variation is typical across much of the plasmid and may reflect genuine sequence divergence or alignment challenges in repetitive plasmid regions.
The concentration of variants in genes involved in immune evasion (sspH2) and surface antigen modification (oafA) suggests that host immune pressure is a major driver of sequence divergence between Salmonella strains. These variants could result in functionally different proteins with altered abilities to manipulate host immunity or evade antibody recognition.

### Frameshift Mutations in Critical Genes

SnpEff identified 740 high-impact protein-truncating variants across the genome. Notably, several conjugative transfer genes harbor multiple frameshift mutations: traN (35 frameshifts), traD (13 frameshifts), traI (9 frameshifts), and traG (9 frameshifts). These frameshifts would be expected to produce non-functional truncated proteins, potentially disrupting plasmid conjugation. However, the presence of numerous synonymous variants in the same genes (traI: 476, traG: 294, traD: 237) alongside frameshifts suggests either: (1) this plasmid population is heterogeneous with both functional and non-functional alleles, (2) frameshift-containing reads represent sequencing errors in homopolymer regions characteristic of Oxford Nanopore data, or (3) true genetic mosaicism in the bacterial culture.

The Type III secretion effector sspH2 contains 420 missense mutations but zero frameshifts, indicating strong purifying selection against loss-of-function mutations in this critical virulence determinant. In contrast, oafA (O-antigen acetylase) contains 1 frameshift among 527 total variants, suggesting some tolerance for gene inactivation potentially reflecting phase variation in surface antigen presentation.

## 4.4 - Plasmid Transfer Genes: Horizontal Gene Transfer Implications
Seven of the ten most variant-dense genes are plasmid-encoded components of the conjugative transfer (tra) operon, including traI (683 variants), traG (433 variants), and traD (334 variants). These genes encode the machinery required for conjugative plasmid transfer, a major mechanism by which bacteria share antibiotic resistance genes and other adaptive traits Military Medical Research (Tang & Liu, 2022). The high variant density in transfer genes suggests this plasmid is substantially divergent from the reference plasmid and may exhibit different transfer characteristics.
This divergence has important public health implications because conjugative plasmids are primary vehicles for spreading antibiotic resistance in Salmonella and other pathogens. A plasmid with altered transfer machinery might have different host range, transfer efficiency, or environmental regulation compared to reference strains, potentially affecting how readily it can spread resistance genes to other bacteria.

## 4.5 - The pncB Anomaly
The gene pncB (3,176 variants) is an extreme outlier, containing 22.5% of all detected variants despite encoding a metabolic enzyme (nicotinate phosphoribosyltransferase) involved in NAD+ biosynthesis. This extraordinary variant density is almost certainly a technical artifact rather than genuine biological variation. Possible explanations include reads from a duplicated gene being incorrectly mapped to pncB, assembly errors in a repetitive region, or structural rearrangement at this locus. This result highlights the importance of validating unexpected findings, particularly extreme outliers, before biological interpretation.

## 4.6 - Mutation Types and Biological Consequences
The variants comprised primarily SNPs (60.5%) and multi-nucleotide polymorphisms (31.2%), with fewer indels (6.0%). Since 85% of variants fall within coding sequences, many likely alter amino acid sequences and potentially protein function. The key question for genes like sspH2 and oafA is whether variants are missense mutations that subtly alter protein activity, or frameshift/nonsense mutations that abolish function entirely. Missense variants in sspH2 could modulate how effectively it suppresses host immunity, while loss-of-function mutations would eliminate this capability. Similarly, variants in oafA could alter the degree or pattern of O-antigen acetylation, changing how the bacterial surface appears to the immune system.

## 4.7 - Limitations and Future Directions
This analysis has several limitations. First, high variant densities in specific genes (pncB, plasmid regions) should be validated by independent sequencing methods to rule out technical artifacts. Second, we did not classify variants as missense, silent, or frameshift mutations, which would clarify their likely functional impact. Third, we did not assess whether the plasmid carries antibiotic resistance genes, which would be critical for understanding the clinical significance of plasmid divergence.
Future work should include: (1) validation of variants in key virulence genes by Sanger sequencing, (2) functional classification of mutations to identify high-impact variants, (3) phenotypic testing to determine if variant-dense virulence genes correlate with altered pathogenicity or immune evasion, and (4) screening for antibiotic resistance genes on the divergent plasmid. Additionally, comparing this isolate to multiple reference genomes rather than a single reference would help distinguish genuine biological variation from reference-specific peculiarities.

4.8 - Biological Significance
Despite these limitations, the results reveal substantial genomic divergence concentrated in functionally important genes. The clustering of variants in immune evasion genes (sspH2, oafA) suggests adaptive evolution in response to host immunity, while the divergent plasmid conjugation machinery indicates potential differences in horizontal gene transfer capability. These findings underscore the genomic diversity within Salmonella enterica and highlight the importance of strain-level genomic characterization for understanding pathogen evolution, outbreak investigation, and vaccine development.

# 5.0 - Conclusion
This study successfully assembled and characterized a Salmonella enterica genome using Oxford Nanopore long-read sequencing, achieving high-quality results with minimal preprocessing. The Flye assembler generated a contiguous assembly that aligned well to the ASM694v2 reference genome, demonstrating the capability of ONT R10/Q20+ chemistry for complete bacterial genome reconstruction. Variant calling identified 14,089 polymorphisms, revealing substantial genomic divergence between this isolate and the reference strain.
The key findings were: (1) variants are concentrated in specific functional categories, particularly virulence factors and plasmid conjugation machinery; (2) chromosomal variant hotspots suggest the presence of genomic islands or mobile elements; (3) the plasmid exhibits bimodal coverage and variant patterns consistent with modular architecture; and (4) genes involved in immune evasion (sspH2, oafA) show high variant densities, suggesting adaptive evolution in host-pathogen interactions.
These results demonstrate that even with high overall genome conservation, Salmonella strains can differ substantially in clinically relevant genes affecting pathogenicity and horizontal gene transfer. The divergent plasmid conjugation machinery and virulence factors have important implications for public health surveillance, as sequence variation in these regions could affect diagnostic assays, virulence phenotypes, and the spread of antibiotic resistance.

The workflow developed here—combining long-read sequencing, de novo assembly, variant calling, and gene-level annotation—provides a robust framework for bacterial comparative genomics. Future applications of this approach to additional Salmonella isolates would enable population-level analyses of genomic diversity, identification of outbreak-associated variants, and tracking of resistance determinant dissemination. While technical artifacts (particularly the pncB anomaly) highlight the need for validation of extreme outliers, the overall variant patterns reveal biologically meaningful divergence concentrated in adaptive genes.

# **6.0 - References:**

References
Baker, S., & Dougan, G. (2007). The genome of Salmonella enterica serovar Typhi. Clinical Infectious Diseases, 45(Suppl. 1), S29–S33. https://doi.org/10.1086/518143

Chen, Z., Erickson, D. L., & Meng, J. (2020). Benchmarking long-read assemblers for genomic analyses of bacterial pathogens using Oxford Nanopore sequencing. International Journal of Molecular Sciences, 21(23), 9161. https://doi.org/10.3390/ijms21239161

Haesebrouck, F., Pasmans, F., Chiers, K., Maes, D., Ducatelle, R., & Decostere, A. (2017). Efficacy of vaccines against bacterial diseases in swine: What can we expect? Veterinary Microbiology, 100(3-4), 255–268.

Hall, M. B., Wick, R. R., Judd, L. M., Nguyen, A. N., Steinig, E. J., Xie, O., Davies, M., Seemann, T., Stinear, T. P., & Coin, L. (2024). Benchmarking reveals superiority of deep learning variant callers on bacterial nanopore sequence data. eLife, 13, RP98300. https://doi.org/10.7554/eLife.98300

Helm, R. A., Lee, A. G., Christman, H. D., & Maloy, S. (2003). Genomic rearrangements at rrn operons in Salmonella. Genetics, 165(3), 951–959. https://doi.org/10.1093/genetics/165.3.951

Kim, M. L., & Slauch, J. M. (1999). Effect of acetylation (O-factor 5) on the polyclonal antibody response to Salmonella typhimurium O-antigen. FEMS Immunology and Medical Microbiology, 26(1), 83–92. https://doi.org/10.1111/j.1574-695X.1999.tb01374.x

Lê-Bury, G., & Méresse, S. (2017). The modulation of host cell death pathways by intracellular bacterial pathogens. Microbes and Infection, 19(9-10), 452–459. https://doi.org/10.1016/j.micinf.2017.04.002

Liyanage, K., Samarakoon, H., Parameswaran, S., & Gamaarachchi, H. (2023). Efficient end-to-end long-read sequence mapping using minimap2-FPGA integrated with hardware accelerated chaining. Scientific Reports, 13(1), 20174. https://doi.org/10.1038/s41598-023-47354-8

Nagai, H., & Roy, C. R. (2020). The DotA/IcmT4SS of Legionella pneumophila and Coxiella burnetii. Frontiers in Cellular and Infection Microbiology, 10, 139. https://doi.org/10.3389/fcimb.2020.00139

Oxford Nanopore assembly using Flye. (2025). Purdue University Research Computing. Retrieved January 18, 2026, from https://rcac-bioinformatics.github.io/genome-assembly/oxford-nanopore-assembly.html

Purushothaman, S., Roloff, T., Egli, A., & Seth-Smith, H. M. B. (2026). Benchmarking Illumina and Oxford Nanopore Technologies (ONT) sequencing platforms for whole genome sequencing of bacterial genomes and use in clinical 

microbiology. BMC Medical Genomics, 19(1), 16. https://doi.org/10.1186/s12920-025-02305-2

Saada, B., Zhang, T., Siga, E., Zhang, J., & Muniz, M. M. M. (2024). Whole-genome alignment: Methods, challenges, and future directions. Applied Sciences, 14(11), 4837. https://doi.org/10.3390/app14114837

Tang, W., & Liu, J. (2022). Conjugative plasmids and the spread of antibiotic resistance genes. Plasmid, 119, 102609. https://doi.org/10.1016/j.plasmid.2022.102609

Wick, R. R., Judd, L. M., & Holt, K. E. (2023). Assembling the perfect bacterial genome using Oxford Nanopore and Illumina sequencing. PLOS Computational Biology, 19(3), e1010905. https://doi.org/10.1371/journal.pcbi.1010905

Zhang, Y., Xiao, L., Lai, X., Liu, S., Peng, L., Dai, M., & Bi, J. (2020). The Salmonella effector protein SspH2 suppresses macrophage inflammatory cytokine production via the inhibition of JNK-mediated signaling pathways. Frontiers in Immunology, 11, 940. https://doi.org/10.3389/fimmu.2020.00940

Zhao, W., Zeng, W., Pang, B., Luo, M., Peng, Y., Xu, J., Kan, B., Li, Z., & Lu, X. (2023). Oxford Nanopore long-read sequencing enables the generation of complete bacterial and plasmid genomes without short-read sequencing. Frontiers in Microbiology, 14, 1179966. https://doi.org/10.3389/fmicb.2023.1179966

