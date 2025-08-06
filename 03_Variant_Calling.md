# A weaponized phage suppresses competitors in historical and modern metapopulations of pathogenic bacteria
# 3. Variant Calling

Program                  | Location
------------------------ | ----------------------------
*bcftools v.1.11*        | (https://github.com/samtools/bcftools)
*vcftools v.0.1.16*      | (https://github.com/vcftools/vcftools)


Genotype Likelihoods and *de-novo* SNP calls were done using *bcftools mpileup* and *bcftools call*, respectively. We filtered out reads with Mapping Qualities lower than 30.
```bash
reference=p25.C2.fa

bcftools mpileup -f $reference -C 30 -Ov $sample.mapped_to_Pseudomonas.dd.bam | bcftools call -c --ploidy 1 | bgzip > $sample.vcf.gz
```

All the individual VCFs were merged using *bcftools merge* and only biallelic SNPs were kept for downstream analysis. The number 
```bash
ALL_SAMPLES=*.vcf.gz
bcftools merge $ALL_SAMPLES.vcf.gz | bcftools view -v snps -m2 -M2| bgzip > ALL.SNPs.biallelic.vcf.gz
```

The [VCF file can be found here](/data/03_Variant_Calling/)

---
[Main README](/README.md) | [Previous - 02. Historic Authentication](/02_Historic_Authentication.md) | [Next - 04. Phylogenetic Analyses](/04_Phylogenetic_Analyses.md)
