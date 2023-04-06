# A weaponized phage suppresses competitors in historical and modern metapopulations of pathogenic bacteria
# 1. Processing, Mapping and Variant Calling of Historic Samples

Program                  | Location
------------------------ | ----------------------------
*AdapterRemoval v2*      | (https://github.com/mikkelschubert/adapterremoval)
*bwa  v.0.7.17*          | (https://github.com/lh3/bwa)
*bwa-mem2*               | (https://github.com/bwa-mem2/bwa-mem2)
*samtools v.1.11*        | (https://github.com/samtools/samtools)
*sambamba v0.8.0*        | (https://github.com/biod/sambamba)
*bcftools v.1.11*        | (https://github.com/samtools/bcftools)
*vcftools v.0.1.16*      | (https://github.com/vcftools/vcftools)

## Alignment of short reads to the *Pseudomonas* sp. reference genome

We used a representative strain of the *Pseudomonas* ATUE5 (p25.C2) assembly as the reference genome and indexed it with  *bwa index*.
```bash
bwa index p25.C2.fa
```

Historical and modern samples were preprocessed and mapped to the reference alignment differently.
#### Alignment of short reads from Historical samples
We expect that most of historical-type reads are short and therefore can be merged (collapsed) from each pair of sequenced reads, so we used *AdapterRemoval2* allowing to `--collapse` reads when possible. We then used *bwa aln* with seed deactivation (flag -l 1024) to allow the alignment of of substitutions present at the termini of the reads. Further, we used *bwa samse* to create a SAM file from the aligned format.
```bash
AdapterRemoval --file1 $sample1.R1.fastq.gz --file2 $sample1.R2.fastq.gz --collapse --gzip --basename $sample

reference=p25.C2.fa
bwa aln -l 1024 -f $sample.sai $reference $sample.collapsed.fq.gz
bwa samse -r "@RG\tID:$sample\tSM:$sample" -f $sample.sam $reference $sample.sai $sample.collapsed.fq.gz
```

### Alignment of short reads from Modern samples
In contrast, reads from modern samples are expected to yield longer fragment lengths, so no merging step is required. Furthermore, no substitution enrichment pattern is expected to happen so we used *bwa-mem2* for mapping the trimmed reads to the reference alignment.
```bash
AdapterRemoval --file1 $sample1.R1.fastq.gz --file2 $sample1.R2.fastq.gz --gzip --basename $sample

bwa-mem2 mem -R "@RG\tID:$sample\tSM:$sample" $reference $sample.pair1.truncated.gz $sample.pair2.truncated.gz > $sample.sam
```

At this stage both historical and modern samples are in SAM format. Therefore we continue to treat them without distinction. We used *sambamba* to discard non-mapped reads, sort and mark PCR optical duplicates.
```bash
sambamba view --num-filter /4 -t $threads -f bam -o $sample.uns.bam -S $sample.sam
rm $sample.sam $sample.sai

sambamba sort --tmpdir=$tmpdir -m $maxmem -t $threads -o $sample.mapped_to_Pseudomonas.bam $sample.uns.bam
rm $sample.uns.bam

sambamba markdup --tmpdir=$tmpdir -t $threads $sample.mapped_to_Pseudomonas.bam $sample.mapped_to_Pseudomonas.dd.bam
rm $sample.mapped_to_Pseudomonas.bam*
```

Coverage and depth and statistics were computed using *samtools coverage* and *samtools depth*. We kept only reads with a Mapping Quality value above 30.
```bash
# Genome coverage
samtools coverage -q 30 $sample.mapped_to_Pseudomonas.dd.bam

# Average depth
samtools depth -Q 30 -aa $sample.mapped_to_Pseudomonas.dd.bam | awk '{sum += $3}END{print sum / NR}'
```
A summary of these statistics can be found at this [table](/data/01_Processing_Mapping_VarCalling_Historic_Samples/Coverage_and_Depth_Samples_Pseudomonas.tsv)

## Variant Calling

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



######## HERE!!!!!!! #######







## Variant calling
We used the *HaplotypeCaller* from *GATK* to generate genomic haplotype calls per individual using the duplicate-marked BAM file as input.
```bash
gatk HaplotypeCaller -R 70-15.fa -I sample1_mapped_sorted.dd.bam -O sample1.g.vcf.gz
```

We used *CombineGVCFs*, *GenotypeGVCFs* and *SelectVariants* from *GATK* to combine the individual genomic VCFs, call genotypes and filter SNPs, respectively.
```bash
gatk CombineGVCFs -R 70-15.fa -V sample1.g.vcf.gz -V sample2.g.vcf.gz -V sampleN.g.vcf.gz -O wheat-blast.g.vcf.gz
gatk GenotypeGVCFs -R 70-15.fa -ploidy 1 -V wheat-blast.g.vcf.gz -O wheat-blast.raw.vcf.gz
gatk SelectVariants -select-type SNP -V wheat-blast.raw.vcf.gz -O wheat-blast.raw.snps.vcf.gz
```

We extracted all Quality-by-Depth (QD) values
```bash
bcftools view -H wheat-blast.raw.snps.vcf.gz | cut -f8 |
awk -F "QD=" '{print $2}' | cut -f1 -d ";" | gzip >  wheat-blast.raw.snps.QD.gz
```

Based on the distribution of Quality-by-Depth values, we set filters of one standard deviation around the median value.
```python
# Python
import pandas as pd
QD = pd.read_csv('wheat-blast.raw.snps.QD.gz', header = None, compression = 'gzip')
med = QD.median()
lower = med - QD.std()
upper = med + QD.std()
print(lower, upper)
```

Finally, using the above-mentioned scheme, we filtered SNPs using *GATK VariantFiltration* and created a new VCF file, keeping non-missing positions, using *bcftools*.
```bash
gatk VariantFiltration --filter-name "QD" \
--filter-expression "QD <= $lower || QD >= $upper" \
-V wheat-blast.raw.snps.QD.gz \
-O wheat-blast.snps.filter.vcf.gz

bcftools view -g ^miss wheat-blast.snps.filter.vcf.gz | bgzip > wheat-blast.snps.filtered.vcf.gz
```

The filtered [VCF file can be found here](/data/02_Preprocessing_and_Variant_Calling/wheat.snps.filtered.vcf.gz)

---
[Main README](/README.md) | [Previous - 01. Monsterplex Analyses](/01_Monsteplex_Analyses.md) | [Next - 03. Population Structure Analyses](/03_Population_Structure.md)
