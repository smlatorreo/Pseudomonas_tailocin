# A weaponized phage suppresses competitors in historical and modern metapopulations of pathogenic bacteria
# 2. Historic Authentication

Program                  | Location
------------------------ | ----------------------------
*mapDamage v2.2.3*       | (https://github.com/ginolhac/mapDamage)


In order to authenticate the historic origin of both plant and pathogen associated DNA molecules, we inspected hallmarks of double-stranded *aDNA*:
- Deamination of Cytosine residues, affecting specially termini of the molecules
- Short DNA frament sizes
- High nick frequency in DNA molecules

We used *mapDamage2* on the mapped reads:

```bash
ref_genome=Pseudomonas.fa.gz
mapped_reads=sample_i.sorted.dedup.bam
mapDamage -i ${mapped_reads} -r ${ref_genome} --merge-reference-sequences -y 0.1 --no-stats

ref_genome=A_thaliana.fa.gz
mapped_reads=sample_i.sorted.dedup.bam
mapDamage -i ${mapped_reads} -r ${ref_genome} --merge-reference-sequences -y 0.1 --no-stats
```

MapDamage outputs can be found [here](/data/02_Historic_Authentication/)

## 2.1. Authentication in the host plant *Arabidopsis thaliana*

![Athaliana auth](/data/02_Historic_Authentication/Athaliana_Auth.png)

## 2.2 Authentication in *Pseudomonas* sp.

![Pseudomonas auth](/data/02_Historic_Authentication/Pseudomonas_Auth.png)

## 2.3. Deamination at first base in both organisms

![Ath Pseudomonas Deam1Base](/data/02_Historic_Authentication/Athaliana_Pseudomonas_Deam1stbase.png)

---
[Main README](/README.md) | [Previous - 01. Processing Mapping](/01_Processing_Mapping.md) | [Next - 03. Variant Calling](/03_Variant_Calling.md)
