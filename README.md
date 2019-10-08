# Diatom Metatranscriptomics Workflow

This is a step-by-step walkthrough for the RNAseq workflow developed for the Marchetti lab projects in 2019. The data for the Marchetti lab is publicly available and linked to in their various papers, but for the purposes of this the data included will be dummy samples of data so as to make it as straightforward as possible. 

## 1. Introduction

Diatoms are poorly annotated, but the amount of environmental RNAseq data is increasing in both marine and other relevant environments. It is possible, though harder, to make use of this data without having well-annotated reference genomes/transcriptomes with some computational tools and experimental assumptions. Most importantly from an experimental standpoint it is assumed that without an organisms reference genome/transcriptome the most we can do is annotate contigs and regions. With that in mind, the workflow is as follows:

![ out line ](/images/outline_wide.png)  


It's possible to run this entire workflow on the provided sample data as well as RNAseq data from another experiment, but edits will have to be made for files based on directory locations and other nuances. All work was run on the Longleaf cluster at University of North Carolina Chapel Hill with a SLURM scheduler.

To download this repository use:
```
git clone https://github.com/Lswhiteh/diatom-metatranscriptomics
```
The structure will look as such:
```
Diatom_workflow/
├── Raw_Reads
├── Quality_Control
├── Assembly
├── Annotation
├── Alignment
└── DEA
```









