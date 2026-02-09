#!/bin/bash
# Script to annotate variants with gene information from GFF

cd ~/binf6110/assignment1

VCF="results/vcf/medaka_reads_only/medaka.annotated.sorted.vcf.gz"
GFF="results/ref/ASM694v2/ASM694v2_genomic.gff"
OUT="results/vcf/variant_gene_annotation.tsv"

echo "Extracting variant positions..."
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE\n' "$VCF" > /tmp/variants.txt

echo "Annotating variants with genes from GFF..."
# Extract gene features from GFF
awk '$3=="gene"' "$GFF" | \
awk '{
    split($9, a, ";");
    for(i in a) {
        if(a[i] ~ /^gene=/) {
            gsub(/gene=/, "", a[i]);
            gene = a[i];
        }
        if(a[i] ~ /locus_tag=/) {
            gsub(/locus_tag=/, "", a[i]);
            locus = a[i];
        }
    }
    print $1"\t"$4"\t"$5"\t"gene"\t"locus;
}' > /tmp/genes.bed

echo "Finding overlaps between variants and genes..."
# Header
echo -e "CHROM\tPOS\tREF\tALT\tTYPE\tGENE\tLOCUS_TAG" > "$OUT"

# For each variant, find overlapping genes
while IFS=$'\t' read -r chrom pos ref alt type; do
    gene_info=$(awk -v chr="$chrom" -v p="$pos" \
        '$1==chr && p>=$2 && p<=$3 {print $4"\t"$5; exit}' /tmp/genes.bed)
    
    if [ -n "$gene_info" ]; then
        echo -e "$chrom\t$pos\t$ref\t$alt\t$type\t$gene_info"
    else
        echo -e "$chrom\t$pos\t$ref\t$alt\t$type\tintergenic\t-"
    fi
done < /tmp/variants.txt >> "$OUT"

echo "Generating summary statistics..."
# Count variants per gene
echo -e "\n=== Top 20 genes by variant count ===" > results/vcf/variant_gene_summary.txt
awk -F'\t' 'NR>1 && $6!="intergenic" {print $6}' "$OUT" | \
    sort | uniq -c | sort -rn | head -20 >> results/vcf/variant_gene_summary.txt

# Count by variant type
echo -e "\n=== Variant types ===" >> results/vcf/variant_gene_summary.txt
awk -F'\t' 'NR>1 {print $5}' "$OUT" | sort | uniq -c | sort -rn >> results/vcf/variant_gene_summary.txt

# Intergenic vs genic
echo -e "\n=== Genic vs Intergenic ===" >> results/vcf/variant_gene_summary.txt
awk -F'\t' 'NR>1 {if($6=="intergenic") print "intergenic"; else print "genic"}' "$OUT" | \
    sort | uniq -c >> results/vcf/variant_gene_summary.txt

echo "Done!"
echo "Results saved to:"
echo "  - $OUT (full annotation)"
echo "  - results/vcf/variant_gene_summary.txt (summary stats)"

# Show preview
echo -e "\n=== Preview of annotated variants ==="
head -20 "$OUT" | column -t

echo -e "\n=== Summary ==="
cat results/vcf/variant_gene_summary.txt
