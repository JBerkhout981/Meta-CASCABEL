# Meta-CASCABEL

Snakemake pipeline for assembly and binning of metagenomics reads.

**Current version:** 5.0

The pipeline creates different output files which allow the user to explore the data and results in a simple way, as well as facilitate downstream analysis based on the generated output files.

* Different quality control steps on the reads.
* Read trimming and filtering
* Assembly
* Binning 
* Bin evaluation
* Taxonomic classification

## Quick start

**Download or clone the repository**

> git clone -b MetaCASCABEL_v5 https://github.com/AlejandroAb/Meta-CASCABEL

**Required input files**

The pipeline is designed to analyze one or more metagenomes.
For each metagenome you should supply the paired end raw reads:

Forward raw reads (fastq or fastq.gz)
Reverse raw reads (fastq or fastq.gz)

When using unzipped reads make sure to set 'gzip_input' to 'F' in config.yaml, or to 'T' when working with zipped reads.
If you want to analyze pre trimmed reads, you can supply these as input and then change 'trimming' in the configfile to 'F'.
In order to only perform the binning, you can also supply a fasta file containing your assembly. In such case, you also need to supply the raw data.

**Configure the input files**

There are two ways to configure the input files, depending on whether you are working with a one sample or multiple samples.

*One sample*

This method can be used when analysing one sample.

Go to the config.yaml file and fill in the following:
* Enter the name of your sample in 'SAMPLES'. 
* Enter the absolute path to the forward reads in 'fw_reads:' between quotes.
* Enter the absolute path to the reverse reads after 'rv_reads:' between quotes.
* Leave 'input_files:' empty.

For example:
```yaml
SAMPLES: ["NIOZ118"]

fw_reads: "/export/lv4/projects/workshop_2023/S10_Assembly/rawdata_1/NIOZ118_R1.fastq.gz"
rv_reads: "/export/lv4/projects/workshop_2023/S10_Assembly/rawdata_1/NIOZ118_R2.fastq.gz"
input_files: ""
```

*Multiple samples*

* Create a 'input.txt' file. The filename can be anything but in this example, we use 'input.txt'.
* Add one sample per line. Each line must contain three tab-separated columns:
1. Sample name
2. Absolute path to the forward reads
3. Absolute path to the reverse reads
* Add the sample names to 'SAMPLES' in config.yaml. The names must be identical to the sample names in the first column of 'input.txt'.
* Leave 'fw_reads:' and 'rv_reads:' empty.
* Enter the name of the '.txt' file after 'input_files:'.

An example of the config file:
```yaml
SAMPLES: ["NIOZ114","NIOZ118","NIOZ130"]

fw_reads: ""
rv_reads: ""
input_files: "input.txt"
```

An example of input.txt:
```
NIOZ114 /export/lv4/projects/workshop_2023/S10_Assembly/rawdata_1/NIOZ114_R1.fastq.gz   /export/lv4/projects/workshop_2023/S10_Assembly/rawdata_1/NIOZ114_R2.fastq.gz
NIOZ118 /export/lv4/projects/workshop_2023/S10_Assembly/rawdata_1/NIOZ118_R1.fastq.gz   /export/lv4/projects/workshop_2023/S10_Assembly/rawdata_1/NIOZ118_R2.fastq.gz
NIOZ130 /export/lv4/projects/workshop_2023/S10_Assembly/rawdata_1/NIOZ130_R1.fastq.gz   /export/lv4/projects/workshop_2023/S10_Assembly/rawdata_1/NIOZ130_R2.fastq.gz
```

**Edit configuration file**

To run the script you need to go through the configuration file (config.yaml). 

Some mandatory options are left empty as default:
* PROJECT
* RUN
* Input configuration (explained above)
* ANALYSIS
* ASSEMBLER
* BINNING

Make sure to go through all these options, otherwise the script won't run. 

IMPORTANT! If you run the pipeline on SLURM set 'interactive' to 'F'

**Run the pipeline using SLURM**

If you open hpc.sh you can set -j (number of jobs) and -c (number of cpu's) according to your needs and available recources
>  sbatch hpc.sh

**Run the pipeline without SLURM**

*Activating environment*

>  module load anaconda/2024.02
> 
>  conda activate snakemake_v7.14.2
> 
>  export GTDBTK_DATA_PATH="/export/lv13/databases/gtdb/release232"

*dry run*

> snakemake --configfile config.yaml -p

*Run*

Set -j (number of jobs) and -c (number of cpu's) according to your needs and available recources
> snakemake --configfile config.yaml  -j2 -c35 --use-conda --conda-frontend conda 

*Generating report file*

You can set the name to anything you want
> snakemake --configfile config.yaml --report report_name.zip

**Output files structure**

Needs to be updated

```
<PROJECT>
├── runs
│   └── <RUN>
│       └── <SAMPLE>_data
│           ├── taxonomy  #Output from taxonomy profiling tool
│           │   └── <TAXONOMY_PROFILING>.taxonomy.report
│           ├── trimmed
│           │   ├── qc  #FastQC result for trimmed reads
│           │   ├── read1_paired.fq  #Trimmed reads
│           │   ├── read1_singles.fq
│           │   ├── read2_paired.fq
│           │   └── read2_singles.fq
│           ├── assembly_<ASSEMBLER> 
│           │   ├── contigs.fasta   # Assembly - contigs
│           │   ├── scaffolds.fasta # Assembly - scaffolds (if available)
│           │   └── quast  # Assembly statistics
│           ├── bwa-mem  #Assembly mapping against raw reads
│           │   ├── <ANALYSIS>_<ASSEMBLER>_depth.txt  # depth coverage
│           │   ├── <ANALYSIS>_<ASSEMBLER>_mapped_against_cross-assembly_sorted.bam # bam file
│           │   └── <ANALYSIS>_<ASSEMBLER>_mapped_against_cross-assembly_sorted.flagstat #stats
│           ├── binning #The location for the bins vary per method 
│           │   ├── abundance.<method>.tsv  #Information about the bin abundance per method
│           │   ├── binsanity
│           │   │   └── <ANALYSIS>_<ASSEMBLER>
│           │   │       └── BinSanity-Final-bins  #BinSanity bins folder
│           │   ├── concoct
│           │   │   └── <ANALYSIS>_<ASSEMBLER>    #Concoct bins
│           │   ├── das
│           │   │   └── <ANALYSIS>_<ASSEMBLER>
│           │   │       └──DasOut_DASTool_bins    #DASTool bins
│           │   ├── maxbin
│           │   │   └── <ANALYSIS>_<ASSEMBLER>    #MaxBin bins
│           │   ├── metabat2
│           │   │   └── <ANALYSIS>_<ASSEMBLER>    #Metabat bins
│           │   ├── checkM_<bin_method>
│           │   │   └── summary.txt
│           │   ├── gtdbtk_<bin_method>
│           │   │   ├── gtdbtk.ar122.summary.tsv -> classify/gtdbtk.ar122.summary.tsv
│           │   │   └── gtdbtk.bac120.summary.tsv -> classify/gtdbtk.bac120.summary.tsv
│           │   ├── FinalBins
│           │   │   ├── contig_coverage.txt
│           │   │   ├── new_names.txt
│           │   │   ├── NIOZ114-1.fna
│           │   │   ├── NIOZ114-2.fna
│           │   │   └── NIOZ114-3.fna
│           │   └── FinalBins.summary.tsv
│           └── unbinned
│               ├── unbinned_contigs_list.txt # List of unbinned contigs
│               └── unbinned.fasta # fasta file with unbinned contigs
└── samples
```


**Environment and dependencies** 

In the directory '/envs' are yaml files of all the different tools with their versions and dependencies.

'snake_env.yaml' shows the tools and dependencies used to run the script. The parameters '--use-conda' and '--conda-frontend conda' allow the pipeline to use different conda environments according to the necessary tools per rule.

Right now the pipeline is configured to use global environments available on the server. If these are not available to you or you want to run the pipeline on your own server, you can change all 'conda:' occurrences in the Snakefile from conda: "sequali_v1.0.2" to conda: "envs/sequali.yaml" and repeat this for every tool. This way it will automatically create all the new environments needed. 

