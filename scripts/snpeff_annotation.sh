#!/bin/bash
# SnpEff Functional Annotation of Variants
# Author: Vian Lelo
# Date: February 8, 2026
# Description: Annotates variants with predicted functional effects using SnpEff

set -e  # Exit on error

cd ~/binf6110/assignment1

echo "=== SnpEff Functional Annotation Pipeline ==="
echo "Start time: $(date)"

# ============================================
# Step 1: Download GenBank annotation
# ============================================
echo ""
echo "[Step 1] Downloading GenBank annotation for ASM694v2..."

cd results/ref/ASM694v2/

# Check if already downloaded
if [ ! -f "ASM694v2_genomic.gbff" ]; then
    wget -O ASM694v2_genomic.gbff.gz \
      "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gbff.gz"
    gunzip ASM694v2_genomic.gbff.gz
    echo "GenBank file downloaded and extracted."
else
    echo "GenBank file already exists, skipping download."
fi

cd ~/binf6110/assignment1

# ============================================
# Step 2: Set up SnpEff custom database
# ============================================
echo ""
echo "[Step 2] Setting up SnpEff custom database..."

# Create directory structure
mkdir -p data/snpEff_data/ASM694v2

# Copy reference files
cp results/ref/ASM694v2/ASM694v2_genomic.fna data/snpEff_data/ASM694v2/sequences.fa
cp results/ref/ASM694v2/ASM694v2_genomic.gbff data/snpEff_data/ASM694v2/genes.gbk

echo "Reference files copied to snpEff data directory."

# Create snpEff config file
cat > snpEff.config << 'EOF'
# snpEff configuration
data.dir = ./data/snpEff_data

# ASM694v2 genome
ASM694v2.genome : Salmonella_enterica_ASM694v2
EOF

echo "snpEff.config created."

# ============================================
# Step 3: Build SnpEff database
# ============================================
echo ""
echo "[Step 3] Building SnpEff database from GenBank annotation..."

snpEff build -genbank -v ASM694v2

echo "Database build complete."

# ============================================
# Step 4: Annotate variants with SnpEff
# ============================================
echo ""
echo "[Step 4] Annotating variants with SnpEff..."

snpEff ann -v ASM694v2 \
  -stats results/vcf/snpEff_summary.html \
  results/vcf/medaka_reads_only/medaka.annotated.sorted.vcf.gz \
  > results/vcf/medaka.snpeff.vcf

echo "Variant annotation complete."

# ============================================
# Step 5: Compress and index annotated VCF
# ============================================
echo ""
echo "[Step 5] Compressing and indexing annotated VCF..."

bgzip -f results/vcf/medaka.snpeff.vcf
tabix -p vcf results/vcf/medaka.snpeff.vcf.gz

echo "VCF compressed and indexed."

# ============================================
# Step 6: Extract functional effect summaries
# ============================================
echo ""
echo "[Step 6] Extracting functional effect summaries..."

# All variants
echo "Extracting all variant effects..."
bcftools query -f '%INFO/ANN\n' results/vcf/medaka.snpeff.vcf.gz | \
  cut -d'|' -f2 | \
  sort | uniq -c | sort -rn > results/vcf/functional_effects_summary.txt

echo "Total variants by effect type:"
cat results/vcf/functional_effects_summary.txt

# Variants excluding pncB
echo ""
echo "Extracting variant effects excluding pncB..."
bcftools view -e 'ANN~"pncB"' results/vcf/medaka.snpeff.vcf.gz | \
  bcftools query -f '%INFO/ANN\n' | \
  cut -d'|' -f2 | \
  sort | uniq -c | sort -rn > results/vcf/effects_no_pncB.txt

echo "Variants excluding pncB:"
head -20 results/vcf/effects_no_pncB.txt

# Key genes effects
echo ""
echo "Extracting effects for key genes..."
bcftools query -f '%INFO/ANN\n' results/vcf/medaka.snpeff.vcf.gz | \
  awk -F'|' '{print $4"\t"$2}' | \
  grep -E "^(sspH2|oafA|traI|pncB|traG|traD|ssbB|traN|traC|gogB)" | \
  sort | uniq -c | sort -k2,2 -k1,1rn > results/vcf/key_genes_effects.txt

echo "Effects in key genes:"
head -30 results/vcf/key_genes_effects.txt

# pncB-specific analysis
echo ""
echo "Analyzing pncB variants..."
PNCB_COUNT=$(bcftools view -i 'ANN~"pncB"' results/vcf/medaka.snpeff.vcf.gz | grep -v "^#" | wc -l)
echo "Total pncB variants: $PNCB_COUNT"

bcftools view -i 'ANN~"pncB"' results/vcf/medaka.snpeff.vcf.gz | \
  bcftools query -f '%INFO/ANN\n' | \
  cut -d'|' -f2 | \
  sort | uniq -c | sort -rn > results/vcf/pncB_effects.txt

echo "pncB effect breakdown:"
cat results/vcf/pncB_effects.txt

# ============================================
# Step 7: Generate summary statistics
# ============================================
echo ""
echo "[Step 7] Generating summary statistics..."

cat > results/vcf/snpeff_analysis_summary.txt << EOF
=== SnpEff Functional Annotation Summary ===
Generated: $(date)

TOTAL VARIANTS: 14,089

HIGH-IMPACT VARIANTS (protein-truncating):
$(bcftools view results/vcf/medaka.snpeff.vcf.gz | bcftools query -f '%INFO/ANN\n' | grep -E "frameshift_variant|stop_gained|stop_lost|start_lost" | wc -l) variants

MODERATE-IMPACT VARIANTS (missense):
$(bcftools view results/vcf/medaka.snpeff.vcf.gz | bcftools query -f '%INFO/ANN\n' | grep "missense_variant" | wc -l) variants

LOW-IMPACT VARIANTS (synonymous):
$(bcftools view results/vcf/medaka.snpeff.vcf.gz | bcftools query -f '%INFO/ANN\n' | grep "synonymous_variant" | wc -l) variants

pncB VARIANTS: $PNCB_COUNT ($(echo "scale=2; $PNCB_COUNT*100/14089" | bc)%)

FILES GENERATED:
- results/vcf/medaka.snpeff.vcf.gz (annotated VCF)
- results/vcf/snpEff_summary.html (HTML report)
- results/vcf/functional_effects_summary.txt
- results/vcf/effects_no_pncB.txt
- results/vcf/key_genes_effects.txt
- results/vcf/pncB_effects.txt

For visualization, open: results/vcf/snpEff_summary.html
EOF

cat results/vcf/snpeff_analysis_summary.txt

echo ""
echo "=== Pipeline Complete ==="
echo "End time: $(date)"
echo ""
echo "View HTML report: results/vcf/snpEff_summary.html"
echo "View summary files in: results/vcf/"
