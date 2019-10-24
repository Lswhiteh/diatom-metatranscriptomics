# Diatom Metatranscriptomics Workflow

This is a step-by-step walkthrough for the RNAseq workflow developed for the Marchetti lab projects in 2019. The data for the Marchetti lab is publicly available and linked to in their various papers, but for the purposes of this the data included will be dummy samples of data so as to make it as straightforward as possible. 

## 1. Introduction

Diatoms are poorly annotated, but the amount of environmental RNAseq data is increasing in both marine and other relevant environments. It is possible, though harder, to make use of this data without having well-annotated reference genomes/transcriptomes with some computational tools and experimental assumptions. 

Most importantly: from an experimental standpoint it is assumed that without an organisms reference genome/transcriptome the most we can do is annotate contigs and regions using the algorithms built for transcriptome assembly. In other words, nobody has validated the quality of metatranscriptome assemblies using transcriptome assemblers, and they're not explicitly built for that.

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

I've read some recent papers suggesting that trimming adapters doesn't matter with a lot of reads nowadays, or that it's done on the sequencing facility side, but it's traditionally a good way to at least check for so I've included it commented out in case they're present in the sample. 

After the fastqc and trimming has been done MultiQC allows for an easy way to look at all samples in tandem.

`multiqc <reports_dir>` looks through directories recursively and searches for the common strings among report files made by tools, then it compiles them into a nice clean set of graphs on an html page. [Here's an example report](https://github.com/Lswhiteh/diatom-metatranscriptomics/blob/master/QC/multiqc_report.html) of the data provided in this repo. Download it, then open it with your web browser to look at it.


## 3. Assembly

Trinity is the tool of choice for the Marchetti lab workflows. It's built for fast, parallellized de novo transcriptome assembly with an automated 3-step pipeline. The output of Trinity is a fasta file that contains contigs (contiguous sequences of DNA/RNA) it was able to assemble. Here's an example of it's usage:
```
#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=4-00:00:00
#SBATCH --mem=500G
#SBATCH --ntasks=24
#SBATCH -J trinity
#SBATCH -o trinity.%A.out
#SBATCH -e trinity.%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lswhiteh@email.unc.edu

cd /pine/scr/l/s/lswhiteh/p_line/LineP2015/2019_analysis/fastqs

module load trinity

#Easy way to condense filenames is to make variables that contain the path to the files
#That way you can organize better and still easily use them in scripts
samp1path="/pine/scr/l/s/lswhiteh/diatom-metatranscriptomics/Reads/Samp_1"
samp2path="/pine/scr/l/s/lswhiteh/diatom-metatranscriptomics/Reads/Samp_2"

Trinity --seqType fq \
	--max_memory 500G \
	--left ${samp1path}/sample_1_R1.fq,${samp2path}/sample_2_R1.fq \
	--right ${samp1path}/sample_1_R2.fq,${samp2path}/sample_2_R2.fq \
	--CPU 24 \
	--output /pine/scr/l/s/lswhiteh/diatom-metatranscriptomics/Assembly/trinity_assembly 
	#Can use this section to run parallelized, do testing to see how much memory you need per submission
	#--grid_exec "/pine/scr/l/s/lswhiteh/HpcGridRunner-1.0.2/hpc_cmds_GridRunner.pl \
	#		--grid_conf /pine/scr/l/s/lswhiteh/HpcGridRunner-1.0.2/trinity_slurm.conf -c"

```
This example doesn't need much time or resources at all (as with any of the tools in this tutorial), but with real data it can take days, even with parallelization. The nice thing about Trinity is that it seems to be able to pick up where it left off if it fails before completion, so don't worry too much if it doesn't finish the first time. 

To use HPC_gridrunner you need a SLURM configuration file. There's a template provided by HPC_gridrunner, but there's a config included in this repo (`SLRUM_grid.conf`) that I've used on the Longleaf cluster.

It's a good idea to cluster the contigs based on similarity. This reduces reduncancy and improves our alignment qualities by removing the possibility of multi-mapping due to badly-assembled contigs. This is easily done using CD-HIT-EST as such:
```
#!/bin/bash

#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --time=0-10:00:00
#SBATCH --mem=400G
#SBATCH --ntasks=12
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lswhiteh@email.unc.edu
#SBATCH -J cdhitremoveredundancy
#SBATCH -o cdhitremoveredundanc.%A.out
#SBATCH -e cdhitremoveredundancy.%A.err

source activate metagenomics

cd /pine/scr/l/s/lswhiteh/diatom-metatranscriptomics/Assembly/trinity_assembly

cd-hit-est -i Trinity.fasta -o clustered_assembly.fa -c 0.98 -n 10 -d 100 -M 400000 -T 12

```
This clusters contigs based on a 98% similarity score and outputs a new fasta file with reduced redundancy. This new clustered assembly will be the one we use for downstream analysis for the rest of the process.


## 4. Annotation

The next step in the process is annotating our contigs. We do this so we know both what gene and organism the contig is likely to be from as well as the likely function of the gene. For the first part PhyloDB and some scripts I wrote ([found here](https://github.com/Lswhiteh/phylodbannotation)). 

[DIAMOND](https://github.com/bbuchfink/diamond) is a very fast BLAST-like tool that can quickly align many sequences with not a lot of overhead. We can build indexes for and align to both PhyloDB and KEGG in one script as such:
```
#!/bin/bash

#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --time=0-24:00:00
#SBATCH --mem=300G
#SBATCH --ntasks=12
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lswhiteh@email.unc.edu
#SBATCH -J diamond
#SBATCH -o diamond.%A.out
#SBATCH -e diamond.%A.err

module load diamond

cd /pine/scr/l/s/lswhiteh/diatom-metatranscriptomics/Annotation

#Need to blast against phyloDB and KEGG

#diamond makedb --in /proj/marchlab/data/phylodb/phylodb_1.076.pep.fa -d phylodb

#diamond blastx -d /proj/marchlab/data/phylodb/diamond_db/phylodb \
#	-q /pine/scr/l/s/lswhiteh/diatom-metatranscriptomics/Assembly/trinity_assembly/clustered_assembly.fa \
#	-o clustered_assembly_phylo_annotation.m8 \
#	-p 12 -e 0.000001 -k 1

diamond makedb --in /nas/longleaf/data/KEGG/KEGG/genes/fasta/genes.pep.fasta -d keggdb
diamond blastx -d keggdb \
	-q ../Assembly/trinity_assembly/clustered_assembly.fa \
	-o clustered_assembly_kegg_annotation.m8 \
	-p 12 -e 0.000001 -k 1
```

Now there are tabular files that have the alignments and qualities for both PhyloDB and KEGG. With the scripts I wrote to parse this output it's fairly straightforward to get a tsv file containing all PhyloDB information. Note that you need to have downloaded and unzipped PhyloDB to use this tool, as it requires the text files included. If you're not in the Marchetti lab you'll need to edit `fastaannotation.py` to point to the PhyloDB files. 
```
git clone git@github.com:Lswhiteh/phylodbannotation.git
python3 phylodbannotation/fastaannotation.py <clustered_contigs> <diamond_phylodb_output.m8> <output_fasta.fa> <output_tabular.tsv>

```
This takes in the clustered contigs and the output file generated by using DIAMOND on PhyloDB and outputs both a fasta file with the gene and taxonomy information in the headers as well as a tabular file that is easily parsed in Excel or R. We'll come back to this later when we compile all the data, but for now let's move onto the KEGG annotation.

To get annotation data from the KEGG mapping we'll use a tool called [KEGGANNOT](https://github.com/ctberthiaume/keggannot.git). 
Running keggannot is as follows:
```
git clone https://github.com/ctberthiaume/keggannot.git
cd keggannot/
python2 setup.py install --user
keggannot_genes2ko -m <input_kegg_alignment.m8> <path/to/KEGGdb> > annotated_kegg_information.tsv

```
As you can see we pipe these hits to a file called `annotated_kegg_information.tsv` that we can then use similarly to the PhyloDB tabular output. Before we compile, however, we need to align to our contigs to get read quantifications.

## 5. Alignment 

We know know *what* is present in our samples, but the bigger question is *how much?* To do this, we will align our reads back to the contigs.
[Salmon](https://github.com/COMBINE-lab/salmon.git) is a tool specifically built for pseudo-alignment of RNA seq reads. It's extremely fast and works well for purposes. Similarly to DIAMOND Salmon needs to create an index and align, here's an example:
```#!/bin/bash

#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --time=5:00:00
#SBATCH --mem=200G
#SBATCH --ntasks=24
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lswhiteh@email.unc.edu
#SBATCH -J salmon
#SBATCH -o salmonrun.%A.out
#SBATCH -e salmonrun.%A.err


module load salmon

cd /pine/scr/l/s/lswhiteh/diatom-metatranscriptomics/Reads/

#First need to create index, then map
#salmon index -i ../Alignment/assemblyindex --transcripts ../Annotation/annotated_assembly.fa

salmon quant -l A -i ../Alignment/assemblyindex \
	-1 ./Samp_1/sample_1_R1.fq \
	-2 ./Samp_1/sample_1_R2.fq \
	-p 24 \
	-o ../Alignment/Samp1_quants

salmon quant -l A -i ../Alignment/assemblyindex \
	-1 ./Samp_2/sample_2_R1.fq \
	-2 ./Samp_2/sample_2_R2.fq \
	-p 24 \
	-o ../Alignment/Samp2_quants
```

This creates directories with transcript abundances for each sample that can be passed into a differential expression analysis pipeline such as DESeq2, the one we will use. 

## 6. Differential Expression Analysis

http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#input-data
