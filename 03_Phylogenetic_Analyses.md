# A weaponized phage suppresses competitors in historical and modern metapopulations of pathogenic bacteria
# 3. Phylogenetic Analyses

Program                  | Location
------------------------ | ----------------------------
*bcftools v.1.11*        | (https://github.com/samtools/bcftools)
*vcftools v.0.1.16*      | (https://github.com/vcftools/vcftools)
*PLINK v.1.9*            | (https://www.cog-genomics.org/plink/1.9)
*BEAST2*
*tped2fasta.py*          | (/scripts/03_Phylogenetic_Analyses/tped2fasta.py)


## 

We kept samples for which at least 60% of the *Pseudomonas* sp reference was covered. Aditionally, we kept positions in which less than 10% of the genoptyes were missing. The correspondant [VCF file can be found here](/data/03_Phylogenetic_Analyses/Pseudomonas_modern_historic.low_coverage_removed.10miss.vcf.gz)

We converted the VCF format into a [transposed ped (tped/tfam) format](/data/03_Phylogenetic_Analyses/Pseudomonas_modern_historic.low_coverage_removed.10miss.tped.gz) and subsequently into a [pseudo-fasta format](/data/03_Phylogenetic_Analyses/Pseudomonas_modern_historic.low_coverage_removed.10miss.fasta.gz).

```bash
plink --chr-set -37 --vcf Pseudomonas_modern_historic.low_coverage_removed.10miss.vcf.gz --recode transpose --out Pseudomonas_modern_historic.low_coverage_removed.10miss

python tped2fasta.py Pseudomonas_modern_historic.low_coverage_removed.10miss > Pseudomonas_modern_historic.low_coverage_removed.10miss.fasta
``` 

---
[Main README](/README.md) | [Previous - 02. Variant Calling](/02_Variant_Calling.md) | [Next - 04. Tailocin Analyses](/04_Tailocin_Analyses.md)
