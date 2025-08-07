# A weaponized phage suppresses competitors in historical and modern metapopulations of pathogenic bacteria
# 1. Processing and Mapping of Historic and Modern Samples

Program                  | Location
------------------------ | ----------------------------
*AdapterRemoval v2*      | (https://github.com/mikkelschubert/adapterremoval)
*bwa  v.0.7.17*          | (https://github.com/lh3/bwa)
*bwa-mem2*               | (https://github.com/bwa-mem2/bwa-mem2)
*samtools v.1.11*        | (https://github.com/samtools/samtools)
*sambamba v0.8.0*        | (https://github.com/biod/sambamba)

Short reads were mapped both to *Arabidopsis thaliana* and *Pseudomonas* sp. reference genomes. We describe the pipeline used for th bacterial genome in detail. However, same steps and parameters were used to align to the plat genome.  

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

# Proportion of mapped reads
We summarised the proportion of mapped reads
![Mapped reads proportions](/data/01_Processing_Mapping/MappedReadsProps.png)


---
[Main README](/README.md) | [Next - 02. Historic Authentication](/02_Historic_Authentication.md)
