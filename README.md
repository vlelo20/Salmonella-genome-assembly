# **Genome Assembly and Reference Comparison of *Salmonella enterica***
Author: Vian Lelo
Date created: January 20th, 2026, Updated: February 8th, 2026

## **Introduction**

Whole‑genome assembly of bacterial pathogens has become a central tool in microbiology, epidemiology, and public health. _Salmonella enterica_ is a prokaryotic, Gram‑negative bacterial pathogen responsible for a substantial burden of foodborne disease worldwide, with infections ranging from self‑limiting gastroenteritis to invasive systemic disease. Its genome typically consists of a single circular chromosome of ~4.5–5.0 Mb, often accompanied by one or more plasmids that can carry virulence and antimicrobial resistance determinants (Baker and Dougan 2007). High‑quality, complete genome assemblies for _S. enterica_ enable precise serovar determination, outbreak reconstruction, plasmid tracking, and comparative genomics, all of which depend on accurate reconstruction of both chromosomal and extrachromosomal replicons (Zhao et al. 2023; Wick, Judd, and Holt 2023). 

Historically, bacterial genome assembly relied on short‑read sequencing (e.g., Illumina), which provides highly accurate reads but struggles with repetitive regions, mobile elements, and plasmids. These limitations often result in fragmented assemblies composed of many contigs, complicating downstream analyses such as structural variant detection, plasmid–chromosome distinction, and accurate gene context inference (Wick, Judd, and Holt 2023). Long‑read sequencing technologies, particularly Oxford Nanopore Technologies (ONT), have transformed bacterial genomics by producing reads that span repeats and structural variants, enabling near‑complete or fully closed assemblies for many bacterial species (Purushothaman et al. 2026). For _S. enterica_, which contains repetitive rRNA operons, insertion sequences, and diverse plasmids, long reads are especially valuable for resolving genome structure and capturing clinically relevant mobile genetic elements (Helm et al. 2003).

Oxford Nanopore’s R10 chemistry with Q20+ basecalling represents a substantial improvement in raw read accuracy compared to earlier ONT chemistries. Recent work has shown that ONT long reads alone can generate complete bacterial and plasmid genomes without the need for short‑read data, provided that appropriate assembly and polishing workflows are used (Zhao et al. 2023). Meta‑analyses and benchmarking studies of long‑read assemblers for bacterial pathogens have demonstrated that modern ONT‑focused assemblers can routinely produce single‑contig chromosomal assemblies and accurately reconstruct plasmids, with error rates that are acceptable for many comparative and clinical applications, especially after polishing (Chen, Erickson, and Meng 2020).

Several studies have systematically compared long‑read assemblers and polishing strategies for bacterial genomes. Chen et al. (2023) benchmarked multiple long‑read assemblers on ONT data from bacterial pathogens and found that tools such as Flye, Canu, and Raven can generate highly contiguous assemblies, with Flye often performing particularly well in terms of contiguity and correctness for ONT data. Wick et al. (2023) proposed a best‑practice pipeline for “perfect” bacterial genome assembly using ONT and Illumina reads, recommending Medaka for ONT‑based polishing, highlighting the importance of both assembler choice and polishing strategy. More recent evaluations of ONT versus Illumina platforms for bacterial whole‑genome sequencing in clinical microbiology emphasize that ONT’s long reads enable complete genome reconstruction and plasmid resolution, while short reads still offer advantages in per‑base accuracy; together, they provide complementary strengths for high‑confidence variant calling and outbreak analysis.

For this assignment, the focus is on ONT R10/Q20+ long reads with an N50 of 5–15 kb, which are well suited to de novo assembly of _S. enterica_. The key challenges in assembling such data include handling residual basecalling errors (particularly in homopolymers), correctly resolving repetitive regions and plasmids, and avoiding misassemblies that could confound downstream variant calling. Long‑read assemblers such as Flye are specifically designed to work with noisy long reads, constructing repeat graphs and leveraging read length to resolve complex regions (Oxford Nanopore Assembly using Flye 2025). However, even with improved Q20+ accuracy, ONT assemblies typically require polishing to reduce small indel and substitution errors that can impact gene prediction and variant analysis. Tools like Medaka, which use neural network models trained on ONT data, have been shown to substantially improve consensus accuracy when applied to ONT assemblies (Wick, Judd, and Holt 2023)

Aligning the assembled genome to a reference is essential for variant calling and comparative analysis. Minimap2 has emerged as the standard aligner for long‑read data, offering high speed and accuracy for mapping ONT reads or assembled contigs to reference genomes. In contrast, short‑read aligners such as BWA‑MEM and Bowtie2 are not appropriate for ONT whole‑genome bacterial assemblies (Saada et al. 2024). Once alignments are generated, variant callers designed for long reads, such as medaka, can identify single‑nucleotide variants (SNVs) and small indels between the assembled _S. enterica_ genome and a reference. Visualization tools like IGV allow manual inspection of alignments and variants (Hall et al. 2024).

In summary, assembling and comparing a Salmonella enterica genome from ONT R10/Q20+ reads involves balancing the strengths and limitations of long‑read sequencing. ONT provides long reads that enable complete, structurally accurate assemblies and plasmid resolution, but residual basecalling errors necessitate careful polishing and appropriate variant calling strategies. Meta‑analyses and recent methodological papers support the use of long‑read assemblers like Flye, ONT‑aware polishers like Medaka, and long‑read aligners like minimap2 as a robust foundation for bacterial genome assembly and comparative genomics. This assignment will apply these principles to generate a consensus _S. enterica_ genome, align it to a reference, call variants, and visualize the results.

## **2.0 - Methods proposed:**

### 2.1 - Data and overall strategy
The starting data will consist of raw Oxford Nanopore FASTQ files generated with R10 chemistry and Q20+ basecalling, with an expected read length N50 of 5–15 kb. The target organism is _Salmonella enterica_, a prokaryotic bacterial pathogen with an expected genome size of approximately 4.5–5.0 Mb, potentially including plasmids. The overall strategy is to perform quality control and filtering of the ONT reads, assemble the genome de novo using a long‑read assembler optimized for ONT data, polish the assembly to improve consensus accuracy, download (ASM694v2) _S. enterica_ reference genome from NCBI, align the assembly and/or reads to the reference, perform variant calling, and visualize both the assembly and the variants. A general estimation of the time required is 4-5 hours to run the protocol (Zhao et al. 2023).

### Environment setup
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
### 2.2 - Data Acquisition:
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
Renaming file:
```
mv data/raw/SRR32410565.fastq.gz.fastq data/raw/SRR32410565.fastq
ls -lh data/raw/
```
Reference genome retrieval:
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


### 2.3 - Quality control and read filtering
Initial quality control will be performed using NanoPlot (v1.46.2) to assess read length distribution, and quality scores, confirming that the N50 falls within the expected 5–15 kb range and identifying any obvious issues with the sequencing run. If necessary, reads will be filtered using NanoFilt (v2.8.0) to remove very short or low‑quality reads, which can reduce noise and computational burden without compromising assembly quality (Zhao et al. 2023). For example, a minimum read length of 1,000 bp and a Q<10 will be removed. The first and last 50 bases will be removed to retain clean data (Zhao et al. 2023). Downsampling may be considered to reduce run time and memory usage of 60-70x (Wick, Judd, and Holt 2023). 

QC plot for reads (NanoPlot): 
```
cd ~/binf6110/assignment1

NanoPlot \
  --fastq data/raw/SRR32410565.fastq \
  --outdir data/qc/nanoplot \
  --threads 8 \
  --plots hex dot \
  --loglength

```

### De novo genome assembly

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
## Assembly Polishing

The raw Flye assembly was polished using Medaka consensus (v2.0.1) in the dedicated medaka_env conda environment to correct residual base-level errors (e.g., homopolymer indels) common in ONT drafts, leveraging the same raw reads for neural-network refinement. Medaka automatically selected the R10.4.1 SUP variant model (r1041_e82_400bps_sup_variant_v5.0.0), realigning reads to the draft contigs with minimap2 (map-ont), generating pileup features, and producing a consensus sequence; 14 threads were used for efficiency. Polishing was executed as follows:
​
```
conda activate medaka_env
cd ~/binf6110/assignment1

medaka_consensus \
  -i data/raw/SRR32410565.fastq \
  -d results/assembly/flye_raw/assembly.fasta \
  -o results/assembly/flye_medaka_polished \
  -t 14

```
## Assembly-to-Reference Alignment

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


### Raw reads - Reference Alignment and variant calling

Raw Oxford Nanopore reads (SRR32410565.fastq) were aligned to the ASM694v2 reference genome using minimap2 (v2.30) with the map-ont preset optimized for single-molecule long-read alignment, which accommodates the characteristic error profile of ONT sequencing and applies soft-clipping to read ends. The resulting SAM output was coordinate-sorted and compressed into BAM format using samtools sort (v1.23) parallelized across 8 threads, then indexed with samtools index to enable efficient random access for downstream visualization and quality assessment.

```
minimap2 -t 8 -ax map-ont \
  results/ref/ASM694v2/ASM694v2_genomic.fna \
  data/raw/SRR32410565.fastq | \
samtools sort -@ 8 -o results/align/reads_vs_ASM694v2.bam
samtools index results/align/reads_vs_ASM694v2.bam
```
Mapping Summary Statistics:

Alignment quality metrics were extracted from the reads-to-reference BAM file using samtools flagstat (v1.23) to quantify overall mapping rates, properly paired reads, unmapped reads, and supplementary alignments. Additionally, samtools idxstats provided per-contig read counts and mapped read lengths, enabling assessment of coverage distribution across the chromosome and plasmid sequences and identification of potential unmapped or underrepresented genomic regions.
```
samtools flagstat results/align/reads_vs_ASM694v2.bam > results/align/ASM694v2.flagstat.txt
samtools idxstats results/align/reads_vs_ASM694v2.bam > results/align/ASM694v2.idxstats.txt
```
Coverage Depth Analysis
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
### Medaka-based variant calling on ONT reads against the ASM694v2 reference genome.

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

Post-Processing of Medaka VCF

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
Variant Summary Statistics

Total variants were enumerated by counting non-header lines, while bcftools stats provided a breakdown of SNPs vs. indels, written to medaka.stats.txt for reporting.

Finding Total number of Variants:
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

The coordinate‑sorted, indexed VCF (medaka.annotated.sorted.vcf.gz and its .tbi index) was used for variant exploration in IGV, summary statistics with bcftools, and downstream plotting of variant distributions along the S. enterica chromosome and plasmid.

### Assembly polishing
To improve the base‑level accuracy of the Flye assembly, ONT‑specific polishing will be performed using Medaka (v2.2.0), which applies neural network models trained on ONT data to correct systematic errors (Zhao et al. 2023). The appropriate Medaka model is "sup" for R10/Q20+ chemistry with Q20+. Polishing will be performed by aligning the original ONT reads back to the Flye assembly using minimap2, followed by Medaka consensus calling to generate a polished assembly.

### Visualization and comparative analysis
Visualization of read alignments and variants will be performed using the Integrative Genomics Viewer (IGV), which allows interactive inspection of coverage, alignment quality, and specific variant sites. IGV will be used to validate variant calls in regions of interest and to inspect any suspicious regions. To visualize the assembly structure, Bandage may be used to inspect the Flye assembly graph, confirming that the chromosome and plasmids are resolved into complete circular contigs.

## Results

## De novo Assembly

Flye assembled the 809 Mb raw ONT reads (N50 4,683 bp, estimated coverage 161x) into a 5.1 Mb draft genome spanning 3 contigs (N50 3.3 Mb, largest contig 3.3 Mb), with mean coverage 160x and no scaffolding required. The process generated 4 initial disjointigs (overlap coverage 116x, median divergence 1.6%), followed by repeat graph simplification (mean edge coverage 127x, 95.3% reads aligned) and a single polishing iteration (alignment error rate 1.8%), yielding assembly.fasta likely comprising the chromosome and two plasmids based on contig metrics summarized in assembly_info.txt.


## Discussion
biological interpretation.... a gene that has snps..

## Conclusion



## **References:**

Baker, Stephen, and Gordon Dougan. 2007. “The Genome of Salmonella Enterica Serovar Typhi.” Clinical Infectious Diseases 45(Supplement_1):S29–33. doi:10.1086/518143.

Chen, Zhao, David L. Erickson, and Jianghong Meng. 2020. “Benchmarking Long-Read Assemblers for Genomic Analyses of Bacterial Pathogens Using Oxford Nanopore Sequencing.” International Journal of Molecular Sciences 21(23):9161. doi:10.3390/ijms21239161.

Oxford Nanopore Assembly using Flye. 2025. Retrieved January 18, 2026. https://rcac-bioinformatics.github.io/genome-assembly/oxford-nanopore-assembly.html.

Hall, Michael B., Ryan R. Wick, Louise M. Judd, An N. Nguyen, Eike J. Steinig, Ouli Xie, Mark Davies, Torsten Seemann, Timothy P. Stinear, and Lachlan Coin. 2024. “Benchmarking Reveals Superiority of Deep Learning Variant Callers on Bacterial Nanopore Sequence Data.” eLife 13:RP98300. doi:10.7554/eLife.98300.

Helm, R. Allen, Alison G. Lee, Harry D. Christman, and Stanley Maloy. 2003. “Genomic Rearrangements at Rrn Operons in Salmonella.” Genetics 165(3):951–59. doi:10.1093/genetics/165.3.951.


Liyanage, Kisaru, Hiruna Samarakoon, Sri Parameswaran, and Hasindu Gamaarachchi. 2023. “Efficient End-to-End Long-Read Sequence Mapping Using Minimap2-Fpga Integrated with Hardware Accelerated Chaining.” Scientific Reports 13(1):20174. doi:10.1038/s41598-023-47354-8.

Purushothaman, Srinithi, Tim Roloff, Adrian Egli, and Helena MB Seth-Smith. 2026. “Benchmarking Illumina and Oxford Nanopore Technologies (ONT) Sequencing Platforms for Whole Genome Sequencing of Bacterial Genomes and Use in Clinical Microbiology.” BMC Medical Genomics. doi:10.1186/s12920-025-02305-2.

Saada, Bacem, Tianchi Zhang, Estevao Siga, Jing Zhang, and Maria Malane Magalhães Muniz. 2024. “Whole-Genome Alignment: Methods, Challenges, and Future Directions.” Applied Sciences 14(11):4837. doi:10.3390/app14114837.

Wick, Ryan R., Louise M. Judd, and Kathryn E. Holt. 2023. “Assembling the Perfect Bacterial Genome Using Oxford Nanopore and Illumina Sequencing.” PLOS Computational Biology 19(3):e1010905. doi:10.1371/journal.pcbi.1010905.

Zhao, Wenxuan, Wei Zeng, Bo Pang, Ming Luo, Yao Peng, Jialiang Xu, Biao Kan, Zhenpeng Li, and Xin Lu. 2023. “Oxford Nanopore Long-Read Sequencing Enables the Generation of Complete Bacterial and Plasmid Genomes without Short-Read Sequencing.” Frontiers in Microbiology 14. doi:10.3389/fmicb.2023.1179966.

