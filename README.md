# Diatom Metatranscriptomics Workflow

This is a step-by-step walkthrough for the RNAseq workflow developed for the Marchetti lab projects in 2019. The data for the Marchetti lab is publicly available and linked to in their various papers, but for the purposes of this the data included will be dummy samples of data so as to make it as straightforward as possible. 

## 1. Introduction

Diatoms are poorly annotated, but the amount of environmental RNAseq data is increasing in both marine and other relevant environments. It is possible, though harder, to make use of this data without having well-annotated reference genomes/transcriptomes with some computational tools and experimental assumptions. Most importantly from an experimental standpoint it is assumed that without an organisms reference genome/transcriptome the most we can do is annotate contigs and regions. With that in mind, the workflow is as follows:

![ out line ](/images/outline_wide.png)  


It's possible to run this entire workflow on the provided sample data as well as RNAseq data from another experiment, but edits will have to be made for files based on directory locations and other nuances. All work was run on the Longleaf cluster at University of North Carolina Chapel Hill with a SLURM scheduler.

Note that all *.sb* files can be run straight from the command line using bash, and for some things such as QC this may be preferable since they're not resource intensive. For most parts of this workflow, however, you should use the SLURM scheduler and submit the procedures as jobs.

To download this repository use:
```
git clone https://github.com/Lswhiteh/diatom-metatranscriptomics
```
The structure will look as such:
```
Diatom_workflow/
├── Reads
├── QC
├── Assembly
├── Annotation
├── Alignment
├── PhyloDB
└── DEA
```

Each directory contains a subfolder labeled "code" and "data" for easier segregation of the two. 

In the Reads directory you'll find this:

```
Reads/
├── Samp_1
|   ├── sample_1_R1.fq
|   └── sample_1_R2.fq
└── Samp_2
    ├── sample_2_R1.fq
    └── sample_2_R2.fq
```

### General Notes For Users

These are *very* small samples of reads that I took from some real data.

I should note that if you see any programs being used as modules here and your system does not have them they're all easily available with a Google search and download from their respective websites or Git repo. I would suggest using Conda for pacakge management if you're familiar, but that's for advanced users and past the scope of this.

Most (essentially all) programs have a help function or manual (man) page that you can easily call either by just calling the program without input, using `man <tool>` or `<tool> --help`. Otherwise documentation can be found online. 



## 2. Quality Control

Using Trimmomatic to get rid of bad reads and trim adapters if necessary is straightfoward, and we can combine this with a tool called FastQC to get summary statistics about our reads. This creates "before" and "after" folders with FastQC reports in them. It then uses a tool called MultiQC to generate a thorough report of the quality information output by both Trimmomatic and FastQC. (It can also be used for Salmon later on, if you're interested.)

```
#!/bin/bash

#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --time=0-4:00:00
#SBATCH --mem=60G
#SBATCH --ntasks=16
#SBATCH -J clean_align
#SBATCH -o cleanalign.%A.out
#SBATCH -e cleanalign.%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lswhiteh@unc.edu

module load trimmomatic
module load fastqc

cd /pine/scr/l/s/lswhiteh/diatom-metatranscriptomics/QC/code/

fastqc -t 8 --outdir ../data/before_trimming/

for f in `ls /pine/scr/l/s/lswhiteh/diatom-metatranscriptomics/Reads/*/*R1* | cut -f1,2,3 -d'_'`
do
        echo $f
        trimmomatic PE -threads 8 \
                ${f}_R1.fq ${f}_R2.fq \
                ${f}_fwd_paired_trimmed.fq.gz ${f}_fwd_unpaired_trimmed.fq.gz \
                ${f}_rev_paired_trimmed.fq.gz ${f}_rev_unpaired_trimmed.fq.gz \
                #ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads \
                LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:36

done

fastqc -t 8 --outdir ../data/after_trimming/

```

I've read some recent papers suggesting that trimming adapters doesn't matter with a lot of reads nowadays, but it's traditionally a good way to go so I've included it commented out in case they're present in the sample. 

After the fastqc and trimming has been done MultiQC





