# microbiome_pipelines


## Key concepts

16s rRNA (taxonomy) vs metagenomic (potential function) vs metatranscriptomic (expression)

Broad vs targeted sequencing

Read-based vs Assembly-based

Sample preparations, redundancy, reference samples

Amplicon Sequence Variant (ASV) vs Operation Taxonomic Unit (OTU)

## File formats

BIOM

fastq

## Steps

QC - FastQC, FaQCs, fastp

Trimming - Trimmomatic

Chimera removal

Filter - Human DNA

Read merging - PEAR

Assembly - MEGAHIT, MetaSPAdes, IDBA-UD

Normalization

## Sequencing

Illumina MiSeq

Ion Torrent PGM

Roche 454 GS FLX Titanium

## Software

Processing:

1. QIIME
2. mothur
3. UPARSE
4. DADA2
5. PICRUst

6. MultiQC

Taxonomy profiling:

7. MetaPhlan
8. Kraken
9. MetaDEGalaxy
10. RDP classifier
11. CLARK
12. Taxy-PRO

ORF-prediction:
13. Prodigal
14. FragGeneScan

ORF/sequence aligner:

15. DIAMOND

Read-based Function profiling:

16. MetaCLADE 
17. HMM-GRASPx 
18. UProC 

Functional Databases:

19. metacyc
20. KEGG
21. NCBI RefSeq
22. UniProt

Taxonomy databases:

23. GreenGenes
24. Silva
25. RDP

Metabolic pathway:

26. MinPath
27. IPath

Amplicon sequence variants:

28. SWARM
29. DADA2
30. MED

Evaluate assembly:

31. MetaQUAST

## Interactive

32. Rhea
33. Calypso
34. phyloseq shiny
35. Galaxy platform

## Statistics

36. DeSeq2
37. EdgeR
38. limma
39. phyloseq
40. vegan
41. metagenomeSeq


## Pipelines

Galaxy platform

42. Prokka
43. EDGE Bioinformatics 
44. MG-RAST
45. Trinity
46. MetaPathways

47. Trinotate

48. MGnify

49. Generally Applicable Gene-Set/Pathway Analysis (GAGE)
 
50. HUMAnN2
51. MetaTrans
52. COMAN
53. FMAP
54. SAMSA2
55. SqueezeMeta
56. IMP
57. MOSCA
58. MEGAN
59. MOCAT (EMBL)

60. Genestack https://genestack-user-tutorials.readthedocs.io/tutorials/Microbiome_analysis/

61. EMBL (European Molecular Biology Lab) Microbiome Tools pipeline https://microbiome-tools.embl.de/

62. MetaDEGalaxy https://f1000researchdata.s3.amazonaws.com/manuscripts/20677/a21689b0-7fc7-4325-b1cf-1de545058127_18866_-_matt_field.pdf?doi=10.12688/f1000research.18866.1&numberOfBrowsableCollections=17&numberOfBrowsableGateways=23

63. MetaWRAP

## Companies

64. Genestack


## Tutorials

65. Stanford https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html

66. Genestack https://genestack-user-tutorials.readthedocs.io/tutorials/Microbiome_analysis/

67. Galaxy https://galaxyproject.github.io/training-material/topics/metagenomics/

## Sources (in order of usefulness and detail)

68. Advances and Challenges in Metatranscriptomic Analysis - has pipeline details
https://www.frontiersin.org/articles/10.3389/fgene.2019.00904/full

69. Analysing Microbial Community Composition through Amplicon Sequencing: From Sampling to Hypothesis Testing
https://www.frontiersin.org/articles/10.3389/fmicb.2017.01561/full
