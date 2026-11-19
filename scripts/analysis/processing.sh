
bcftools query -f "%INFO/SVTYPE" chr2.bcf | sort -u
bcftools view --include "SVTYPE!='CTX' & SVTYPE!='CPX' & SVTYPE!='CNV'" 1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz |
    sed -e 's/SVTYPE=DUP/SVTYPE=INS/' | bcftools view --output 1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.simplesv.vcf.gz
bcftools query -f "%INFO/SVTYPE" 1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.simplesv.vcf.gz | sort -u

vcf_file=./1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.simplesv.vcf.gz

# FILTER is always PASS for this data
# Assuming for now that DP and GQ filters are not necessary
# HWE filtering should not be used for rare variants
# Need to ask about DP and GQ filters
# Filtering QUAL at 125 appears to more conservative than in GTEx
bcftools view -Ou --min-alleles 2 --max-alleles 2 --include "AF<=0.01" ./1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.simplesv.vcf.gz |
    bcftools +fill-tags -Ou -- -t F_MISSING |
    bcftools view -Ou --include "F_MISSING < 0.1 & QUAL >= 125" --output ./1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.simplesv.maxmaf01.qual125.vcf.gz
bcftools view -H ./1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.simplesv.maxmaf01.qual125.vcf.gz | wc -l

bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%ALT\t%INFO/SVLEN' ./1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.simplesv.maxmaf01.qual125.vcf.gz | \
    awk '{print $1, $2, $3, $3-$2, $4, $5, $6}' >! rare_sv_positions

parallel -j16 --eta --progress "bcftools view -Ou -s {} ./1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.simplesv.maxmaf01.qual125.vcf.gz | bcftools query -f \"%CHROM\t%POS\t%END\t%ID\t%ALT\t%SVLEN\t%AF\" --include \"GT='alt'\" > ./individual_svs/{}" < /oak/stanford/groups/smontgom/syang/AFGR/expr_ids

