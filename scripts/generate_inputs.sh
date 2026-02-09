#!/bin/bash
# Script to generate variants_hq.bed and coverage_100bp.tsv

# Set paths based on your tree structure
VCF="results/vcf/variants_vs_ASM694v2.vcf.gz"
BAM="results/align/reads_vs_ASM694v2.bam"
PLASMID="NC_003277.2"

# 1. Generate variants_hq.bed from VCF
echo "Generating variants_hq.bed..."
bcftools view -i 'QUAL>=20' "$VCF" | \
bcftools query -f '%CHROM\t%POS0\t%END\t%TYPE\t%QUAL\n' | \
awk 'BEGIN{OFS="\t"} {if($3=="") $3=$2+1; print}' > variants_hq.bed

echo "Created variants_hq.bed with $(wc -l < variants_hq.bed) variants"

# 2. Generate 100bp coverage windows
echo "Generating coverage_100bp.tsv..."
samtools depth -a "$BAM" | \
awk -v bin=100 -v seqid="$PLASMID" '
BEGIN {OFS="\t"}
$1 == seqid {
    b = int(($2-1)/bin)
    sum[b] += $3
    count[b]++
}
END {
    for (b in sum) {
        start = b * bin
        end = start + bin
        avg = sum[b] / count[b]
        print seqid, start, end, avg
    }
}' | sort -k2,2n > coverage_100bp.tsv

echo "Created coverage_100bp.tsv with $(wc -l < coverage_100bp.tsv) windows"
echo "Done! You can now run your R script."
