# NCTC 11168 reference genome

*Campylobacter jejuni* subsp. *jejuni* NCTC 11168 — the first *C. jejuni*
genome sequenced (Parkhill et al., *Nature*, 2000) and the standard
laboratory reference/model strain for *Campylobacter* flagellar biology.
Well-curated RefSeq annotation, extensively used in the flagellar
literature cited throughout this project (FlhF/FlhG, FliW, flaA/flaB,
PflA, etc.).

NCBI RefSeq assembly: **GCF_000009085.1**
Not one of this project's own 166 genomes (checked explicitly) — used
purely as an *external* reference to seed orthology searches for genes
this project's own annotations never named properly.

## Files
- `GCF_000009085.1.gff` — RefSeq annotation
- `GCF_000009085.1.faa` — RefSeq protein FASTA

## Used by (external-reference orthology seeding)
- `orthology_crossref_maf.py` — maf gene family (flagellin glycosylation)
- `orthology_crossref_flagellin.py` — flaA/flaB (flagellin structural genes;
  fliC/flaA/flaB all had ~zero direct gene= hits in this project's own
  166-genome annotations)
- `orthology_crossref_regulon.py` — comprehensive sweep: every gene in
  this reference whose product description matches
  flagell*/chemotax*/chemorecept*/methyl-accepting/motility (50 proteins,
  47 unique orthogroups), used to catch gaps the incremental gene-by-gene
  approach missed (e.g. fliI, flaG, multiple unnamed methyl-accepting
  chemotaxis paralogs)

Downloaded from NCBI FTP:
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/085/GCF_000009085.1_ASM908v1/
