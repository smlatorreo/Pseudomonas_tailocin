# A weaponized phage suppresses competitors in historical and modern metapopulations of pathogenic bacteria
# 4. Phylogenetic Analyses

Program                  | Location
------------------------ | ----------------------------
*bcftools v.1.11*        | (https://github.com/samtools/bcftools)
*vcftools v.0.1.16*      | (https://github.com/vcftools/vcftools)
*PLINK v.1.9*            | (https://www.cog-genomics.org/plink/1.9)
*IQ-TREE v.2.0.3*        | (http://www.iqtree.org/)
*ClonalFrameML v.1.12*   | (https://github.com/xavierdidelot/clonalframeml)
*mask_positions.py*      | (/scripts/03_Phylogenetic_Analyses/mask_positions.py)
*BEAST2 v.2.4.8*         | (https://github.com/CompEvol/beast2)
*tped2fasta.py*          | (/scripts/03_Phylogenetic_Analyses/tped2fasta.py)


## 

We kept samples for which at least 60% of the *Pseudomonas* sp reference was covered. Aditionally, we kept positions in which less than 10% of the genoptyes were missing. The correspondant [VCF file can be found here](/data/04_Phylogenetic_Analyses/Pseudomonas_modern_historic.low_coverage_removed.10miss.vcf.gz)

We converted the VCF format into a [transposed ped (tped/tfam) format](/data/04_Phylogenetic_Analyses/Pseudomonas_modern_historic.low_coverage_removed.10miss.tped.gz) and subsequently into a [pseudo-fasta format](/data/04_Phylogenetic_Analyses/Pseudomonas_modern_historic.low_coverage_removed.10miss.fasta.gz).

```bash
plink --chr-set -37 --vcf Pseudomonas_modern_historic.low_coverage_removed.10miss.vcf.gz --recode transpose --out Pseudomonas_modern_historic.low_coverage_removed.10miss

python tped2fasta.py Pseudomonas_modern_historic.low_coverage_removed.10miss > Pseudomonas_modern_historic.low_coverage_removed.10miss.fasta
``` 

We computed a full phylogeny based on the concatenated genome-wide SNP set using *IQ-TREE*. We explored the best fitting model (*ModelFinder*) and performed ultrafast bootstrap with 1000 repetitions as well as the SH-aLRT test.
```bash
iqtree -s Pseudomonas_modern_historic.low_coverage_removed.10miss.fasta --alrt 1000 -B 1000 --prefix Pseudomonas_modern_historic -T 20
```

The full [IQtree output](/data/04_Phylogenetic_Analyses/Pseudomonas_modern_historic.iqtree) and the [consensus tree](/data/04_Phylogenetic_Analyses/Pseudomonas_modern_historic.contree) can be examined.

We then used the consensus tree to identify putative SNPs that might arise by recombination, using *ClonalFrameML* with deafult parameters. Then we used *mask_positions.py* to create a new fasta file in which positions targetet as recombination are masked.
```bash
ClonalFrameML Pseudomonas_modern_historic.contree Pseudomonas_modern_historic.low_coverage_removed.10miss.fasta

python mask_positions.py Pseudomonas_modern_historic.low_coverage_removed.10miss.fasta Pseudomonas_herbaria.CFML.importation_status.txt > Pseudomonas_herbaria.CFML.masked.fasta
```

Output files from [ClonalFrameML can be found here](/data/04_Phylogenetic_Analyses/ClonalFrameML/) as well as the generated [fasta file](/data/04_Phylogenetic_Analyses/Pseudomonas_herbaria.CFML.masked.fasta)

We used the recombination-masked file as input for a full Bayesian phylogenetic reconstruction using *BEAST2*. We ran 4 independent chains, each with a length of 10 Million repetitions with parameters that can be found at the [XML configuration file](/data/04_Phylogenetic_Analyses/BEAST2/Pseudomonas_herbaria.xml)
```bash
java -jar beast.jar -threads 20 Pseudomonas_herbaria.xml
```

Each chain was analyzed using *Tracer* from the *BEAST package* and [all chains were combined](/data/04_Phylogenetic_Analyses/BEAST2/Pseudomonas_herbaria.COMBINED.log.gz) using *logcombiner*. Finally, a [Maximum Credibility Tree](/data/04_Phylogenetic_Analyses/BEAST2/Pseudomonas_herbaria.COMBINED.MC.tree) was computed using *treeannotator*.

---
[Main README](/README.md) | [Previous - 03. Variant Calling](/03_Variant_Calling.md) |)
