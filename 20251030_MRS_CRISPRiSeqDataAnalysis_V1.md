#Analysis of CRISPRi-seq data 
Started by Marissa Roghair Stroud, October 30, 2025

---
---

## Description of the data
There were 76 samples submitted for sequencing on the NovaSeq X Plus at SeqCenter with Paired End sequencing. The samples were barcoded with the Illumina Unique Dual indexes, so when demultiplexed, there were 76x2 `.fastq.gz`files. 

Each read should be 140 bp in length, with the central 20 bp being the sgRNA sequence. 


## Downloading and transferring data
The data was shared with us on Globus. Each file needs downloaded individually. I chose to download them to my computer, upload those files to the LSS, then use the data transfer node to move them to the HPC. Downloading and transferring each file using my computer only took 2-5 minutes per file, so it wasn't too slow to do in the background while I worked on other things.

###Logging in to data transfer node

```
ssh mroghair@nova.its.iastate.edu
	OR
ssh mroghair@novadtn.las.iastate.edu

Login is the same for both:
(mroghair@nova.its.iastate.edu) Verification code: Microsoft Authenticator
(mroghair@nova.its.iastate.edu) Password: ISU password

```
```
cd /lustre/hdd/LAS/larryh-lab/MRS_Consensus_Assemblies/
```



###Requesting an interactive session

You log into the head node, but try to only use that for basic tasks and logging in. To open an interactive session run:

```bash
srun --nodes 1 --tasks 32 --partition interactive --time 03:00:00 --pty bash

```
This requests 1 node with 32 threads for 3 hours (feel free to change the time). You can exit a node by calling ```exit```.


###Moving over files

Copy/paste all `*.fastq.gz` files from the LSS to Nova

```
cp /lss/research/larryh-lab/Marissa/CRISPRi-seq_expts/SeqCenter/*.fastq.gz /lustre/hdd/LAS/larryh-lab/MRS_Consensus_Assemblies/

ls /lustre/hdd/LAS/larryh-lab/MRS_Consensus_Assemblies/CRISPRi-seq/rawreads/
```

## Install miniconda3

You can't module load miniconda, it won't work right. So we will install our own. **Only do this once.** _(This should already have been done while doing genome assemblies!)_

```bash
module load gcc
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
```

Now we need to configure miniconda so it won't install the packages into your home directory. Create a ```.condarc``` file by configuring the channels (which we have to do anyways). There will now be a file ```home/user/.condarc``` that you can edit.

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Now we need to edit the ```.condarc``` file using a text editor, like vim. ```vim /home/mroghair/.condarc```

```
refresher vim commands:
	i = insert
	
once done editing:
	esc, then :wq (write and quit)

```

The contents should look like this. 

```bash
channels:
    - conda-forge
    - bioconda
    - defaults

pkgs_dirs:
    - /work/LAS/larryh-lab/MRS_Consensus_Assemblies/pkgs
```

We also need to create the directory we just referenced.

```bash
mkdir /work/LAS/larryh-lab/MRS_Consensus_Assemblies/pkgs
```

Now conda knows that it should install all packages into that folder, and not into the default ```/pkgs``` folder in your home directory. This will keep you from running out of space.



## Merge FASTQ files with PEAR

Install PEAR:

```
conda create --prefix ./pear-env
source activate pear-env/
conda install bioconda::pear
```

PEAR useage: [LINK](https://cme.h-its.org/exelixis/web/software/pear/doc.html)

Another useage: [LINK2](https://learnmetabarcoding.github.io/LearnMetabarcoding/processing/pair_merging.html)


**Pear needs to be run BEFORE cutadapt or it will throw errors!!** Cutadapt will toss reads if they don't meet its standards, and if you have paired end reads and one of the paired reads gets tossed, PEAR won't run at all.


```
pear -f rawreads/InocBA_S41_R1_001.fastq.gz -r rawreads/InocBA_S41_R2_001.fastq.gz -o merged_fastq/InocBA_merged.fastq -m 140 -n 140
```

Ashley wrote a batch script array to perform this for all 76 samples at once :)

```
#!/bin/bash 

#SBATCH --time=0-4:00:00  # max job runtime
#SBATCH --cpus-per-task=32  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --mem=200G  # max memory
#SBATCH -J "PEAR array"  # job name
#SBATCH --mail-user=mroghair@iastate.edu  # email address
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=1-76


# go to the relevant directory
cd /lustre/hdd/LAS/larryh-lab/MRS_Consensus_Assemblies

# activate conda environment, if using conda
eval "$(conda shell hook --shell bash)"
source activate pear-env/
cd CRISPRi-seq/rawreads

# in the slurm set up parameters at the top it specifies array 1-76.
# This will create 76 simultaneous jobs, with a corresponding number that can be referenced as $SLURM_ARRAY_TASK_ID

# ls *S"${SLURM_ARRAY_TASK_ID}"_R1_001.fastq.gz lists the files that fit the pattern (should be just one per $SLURM_ARRAY_TASK_ID) and
# by putting $() I am taking the output of the ls command, which is then saved in the variable READ1 
# so for the example you sent me READ1 would be InocBA_S41_R1_001.fastq.gz
# I did this so I can then extract the first bit of the name to use in the new file name
# "${READ1%%_*}" will keep everything before the first underscore, e.g. InocBA
READ1=$( ls *S"${SLURM_ARRAY_TASK_ID}"_R1_001.fastq.gz )
READ2=$( ls *S"${SLURM_ARRAY_TASK_ID}"_R2_001.fastq.gz )
BASENAME="${READ1%%_*}"

cd /lustre/hdd/LAS/larryh-lab/MRS_Consensus_Assemblies/CRISPRi-seq

#also, I looked at the PEAR software and to set the number of threads you have to use -j
pear -f rawreads/$READ1 -r rawreads/$READ2 -o merged_fastq/"${BASENAME}"_merged.fastq -m 140 -n 140 -j 32
```



## Cutadapt

Installation: [LINK](https://cutadapt.readthedocs.io/en/stable/installation.html)

```
conda config --add channels bioconda  #ALREADY DONE
conda config --add channels conda-forge  #ALREADY DONE
conda config --set channel_priority strict

conda create -n cutadapt-env cutadapt

conda activate cutadapt-env

```

### Sequences to trim off 
Upstream sequence: (Direction = Adapter -> sgRNA) ... AKA the 3' adapter for R1
`TAGAATTATAATTTGGGGACCTAGGCCGCGGCCGCGCGAATTCGAGCTCGGTCTCATCAC`
Downstream, reverse complement ... AKA the 5' adapter for R1
`GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGT`


### Array for trimming

This website showed how to run through an array when your file names don't have numbers in them: [LINK](https://hbctraining.github.io/Training-modules/Accelerate_with_automation/lessons/arrays_in_slurm.html#:~:text=How%20can%20I%20use%20$%7BSLURM_ARRAY_TASK_ID%7D?&text=There%20are%20several%20ways%20we,to%20pull%20the%20file%20names.&text=So%20what%20is%20this%20script,Chan%20Bioinformatics%20Core%20(HBC).)

```
cutadapt -g ^TAGAATTATAATTTGGGGACCTAGGCCGCGGCCGCGCGAATTCGAGCTCGGTCTCATCAC -a GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGT -n 2 -m 20 -M 20 -o merged_trimmed/InocBA_trim.fastq merged_fastq/InocBA_merged.fastq.assembled.fastq

```

- `-g` is the 3' adapter
- `-a` is the 5' adapter
- `-n` is the number of times to cut (2, for the 2 adapters)
- `-m` and `-M` are the minimum and maximum read length post-trimming (I wanted only exactly 20 bp, so they are same)
- `-o` is the output file




## MAGeCK

Install MAGeCK. MAGeCK useage: [LINK](https://sourceforge.net/p/mageck/wiki/Home/)

```
conda create --prefix ./mageck-env
source activate mageck-env/
conda install bioconda::mageck

#MAGeCK won't run with NumPy version 2+
pip install 'numpy<2'
```

### Using MAGeCK
MAGeCK basically has only two commands: `count` and `test`. First, you have MAGeCK create a count table, then you have it test differences between the control and treatment. You can also use count tables generated using other software!

First, obtain read counts. You can have MAGeCK count all samples at once, or split it up into smaller chunks. It is important that all replicates and both the control and treatment are counted in the same file.

I ran this in my terminal (didn't submit a batch script) and it only took like 15 minutes to count all 40ish samples.

`Library_for_MAGeCK_v3.txt` is a file that contains a list of all sgRNAs in the library. Any sequences not matching this file are discarded. I usually get 95% of the reads to map, which seems fine. This file needs to be formatted as below. If the gene is not listed, or if the columns are not in this order, MAGeCK will not run. 

- column 1 = sgRNA id
- column 2 = gene
- column 3 = sgRNA sequence

```
#Count one treatment and the control
mageck count -l Library_for_MAGeCK_v3.txt -n NoneP1 --sample-label NoneP1A,NoneP1B,NoneP1C,NoneP1D,InocPA,InocPB,InocPC,InocPD  --fastq NoneP1A_trim.fastq NoneP1B_trim.fastq NoneP1C_trim.fastq NoneP1D_trim.fastq InocPA_trim.fastq InocPB_trim.fastq InocPC_trim.fastq InocPD_trim.fastq

#Count all treatments and the control at once
mageck count -l Library_for_MAGeCK_v3.txt -n AllCounts --sample-label InocPA,InocPB,InocPC,InocPD,NoneP1A,NoneP1B,NoneP1C,NoneP1D,NoneP2A,NoneP2B,NoneP2C,NoneP2D,NoneP4A,NoneP4B,NoneP4C,NoneP4D,PlusP1A,PlusP1B,PlusP1C,PlusP1D,PlusP2A,PlusP2B,PlusP2C,PlusP2D,PlusP4A,PlusP4B,PlusP4C,PlusP4D,LateP2A,LateP2B,LateP2C,LateP2D,LateP4A,LateP4B,LateP4C,LateP4D,SeedA,SeedB,SeedC,SeedD --fastq InocPA_trim.fastq InocPB_trim.fastq InocPC_trim.fastq InocPD_trim.fastq NoneP1A_trim.fastq NoneP1B_trim.fastq NoneP1C_trim.fastq NoneP1D_trim.fastq NoneP2A_trim.fastq NoneP2B_trim.fastq NoneP2C_trim.fastq NoneP2D_trim.fastq NoneP4A_trim.fastq NoneP4B_trim.fastq NoneP4C_trim.fastq NoneP4D_trim.fastq PlusP1A_trim.fastq PlusP1B_trim.fastq PlusP1C_trim.fastq PlusP1D_trim.fastq PlusP2A_trim.fastq PlusP2B_trim.fastq PlusP2C_trim.fastq PlusP2D_trim.fastq PlusP4A_trim.fastq PlusP4B_trim.fastq PlusP4C_trim.fastq PlusP4D_trim.fastq LateP2A_trim.fastq LateP2B_trim.fastq LateP2C_trim.fastq LateP2D_trim.fastq LateP4A_trim.fastq LateP4B_trim.fastq LateP4C_trim.fastq LateP4D_trim.fastq SeedA_trim.fastq SeedB_trim.fastq SeedC_trim.fastq SeedD_trim.fastq
```

Second, test for differences. This time you only call for the treatment and control you are interested in.

This command runs super fast (less than a minute usually)

```
#With one large counts file
mageck test -k AllCounts.count.txt -t NoneP1A,NoneP1B,NoneP1C,NoneP1D -c InocPA,InocPB,InocPC,InocPD -n NoneP1

#with a smaller counts file (same command!)
mageck test -k NoneP1.count.txt -t NoneP1A,NoneP1B,NoneP1C,NoneP1D -c InocPA,InocPB,InocPC,InocPD -n NoneP1

```

###Description of the MAGeCK commands (from their website)
* `mageck` The main portal of the MAGeCK program
* `count`	A sub-command to ask MAGeCK to generate sgRNA read count table.
* `-l library.txt` The provided sgRNA information, including the sgRNA id, the sequence, and the gene it is targeting. See input files for a detailed explanation.
* `-n demo` The prefix of the output files.
* `--sample-label L1,CTRL`	The labels of the two samples are L1 (test1.fastq) and CTRL (test2.fastq).
* `--fastq test1.fastq test2.fastq` The provided fastq file, separated by space. (Technical replicates of the same sample can also indicated using comma as a separator; for example, "sample1_replicate1.fastq,sample1_replicate2.fastq")
* `test` A sub-command to ask MAGeCK to perform sgRNA and gene ranking based on provided read count tables
* `-k sgrna_count.txt` The provided read count table file. The format of the file is specified here.
* `-t HL60.final,KBM7.final` The treatment samples are defined as HL60.final,KBM7.final (or the 2nd and 3rd sample, starting from 0) in sgrna_count.txt. See input files for a detailed explanation.
* `-c HL60.initial,KBM7.initial` The control samples are defined as HL60.initial,KBM7.initial (or the 0th and 1st sample, starting from 0) in sgrna_count.txt. See input files for a detailed explanation.



### MAGeCK details (from their website)

__How to deal with biological replicates and technical replicates?__

A: Usually you can pool the read counts for technical replicates of the same sample. To do this, use comma (,) to separate the fastq files of the technical replicates from the same sample in the --fastq option. For example, `--fastq sample1_replicate1.fastq,sample1_replicate2.fastq sample2_replicate1.fastq,sample2_replicate2.fastq` indicates two samples with 2 technical replicates for each sample.

For biological replicates, treat them as separate samples and use them together when doing the comparison; so MAGeCK can analyze the variance of these samples. For example in the test command, 
`-t sample1_bio_replicate1,sample1_bio_replicate2 -c sample2_bio_replicate1,sample2_bio_replicate2` compares 2 samples (with 2 biological replicates in each sample).




#### MAGeCK sgRNA_summary.txt file

* `sgrna`	sgRNA ID
* `Gene`	The targeting gene
* `control_count`	Normalized read counts in control samples
* `treatment_count`	Normalized read counts in treatment samples
* `control_mean`	Median read counts in control samples
* `treat_mean`	Median read counts in treatment samples
* `LFC`	The log2 fold change of sgRNA
* `control_var`	The raw variance in control samples
* `adj_var`	The adjusted variance in control samples
* `score`	The score of this sgRNA
* `p.low`	p-value (lower tail)
* `p.hig`h	p-value (higher tail)
* `p.twosided`	p-value (two sided)
* `FDR`	false discovery rate
* `high_in_treatment`	Whether the abundance is higher in treatment samples



## DeSeq2

```
conda create --prefix ./deseq2-env
source activate deseq2-env/
conda install bioconda::bioconductor-deseq2
```

This is how you can install DESeq2, but I actually learned after doing this that it's easier to use DESeq2 in R.



## SeqKit

```
conda create --prefix ./seqkit-env
source activate seqkit-env/
conda install bioconda::seqkit
```
I used SeqKit to learn how large my files were to determine the coverage I had for my libraries. This one takes awhile to run, so I did a batch script. The length of the files is saved in the SLURM log file.

```
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# Job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --nodes=1   # Number of nodes to use
#SBATCH --ntasks-per-node=32   # Use 32 processor cores per node 
#SBATCH --time=0-8:0:0   # Walltime limit (DD-HH:MM:SS)
#SBATCH --mail-user=mroghair@iastate.edu   # Email address
#SBATCH --mail-type=BEGIN,END,FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

cd /lustre/hdd/LAS/larryh-lab/MRS_Consensus_Assemblies/

source activate seqkit-env/

cd CRISPRi-seq/rawreads/

seqkit stats *.fastq.gz 
```



