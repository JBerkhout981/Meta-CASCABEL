# Meta-CASCABEL

Snakemake pipeline for assembly and binning of metagenomics reads.

**Current version:** 5.0

The pipeline creates different output files which allow the user to explore the data and results in a simple way, as well as facilitate downstream analysis based on the generated output files.

* Different quality control steps on the reads.
* Taxonomy assessment at different levels
* Read trimming and filtering
* Assembly
* Gene calling (bins, contigs or scaffolds)
* Binning 
* Bin evaluation

## Quick start


**Required input files**

The pipeline is designed to analyze one or more metagenomes.
For each single metagenome you should supply the paired end raw reads:

Forward raw reads (fastq or fastq.gz)
Reverse raw reads (fastq or fastq.gz)

In order to only perform the binning, you can also supply a fasta file containing your assembly. In such case, you also need to supply the raw data.

**Download or clone the repository**

> git clone -b MetaCASCABEL_v5 https://github.com/AlejandroAb/Meta-CASCABEL

**Initialize directory structure**

There are two ways to configure the input.

The first method can only be used when working with one sample. 

Go to the config.yaml file and fill in the following:
1. Enter the name of your sample in 'SAMPLES', for example: SAMPLES: ["SAMPLE_NAME"]
2. Enter the absolute path to the forward reads after 'fw_reads:' between quotes.
3. Enter the absolute path to the reverse reads after 'rv_reads:' between quotes.
4. Leave 'input_files:' empty.

The second method can be used when working with one or more samples.

1. Create a '.txt' file. The filename can be anything but in this example, we use 'input.txt'.
2. Add one sample per line. Each line must contain three tab-separated columns:
- Sample name
- Absolute path to the forward reads
- Absolute path to the reverse reads
3. Add the sample names to 'SAMPLES' in config.yaml. The names must be identical to the sample names in the first column of `input.txt`.
4. Leave 'fw_reads:' and 'rv_reads:' empty.
5. Enter the name of the '.txt' file after 'input_files:'.

Here is an example of what input.txt should look like:
```
NIOZ118 /export/lv4/projects/workshop_2023/S10_Assembly/rawdata_1/NIOZ118_R1.fastq.gz   /export/lv4/projects/workshop_2023/S10_Assembly/rawdata_1/NIOZ118_R2.fastq.gz
NIOZ130 /export/lv4/projects/workshop_2023/S10_Assembly/rawdata_1/NIOZ130_R1.fastq.gz   /export/lv4/projects/workshop_2023/S10_Assembly/rawdata_1/NIOZ130_R2.fastq.gz
NIOZ114 /export/lv4/projects/workshop_2023/S10_Assembly/rawdata_1/NIOZ114_R1.fastq.gz   /export/lv4/projects/workshop_2023/S10_Assembly/rawdata_1/NIOZ114_R2.fastq.gz
```

**Edit configuration file**

<ins>Project name</ins>

```yaml
#------------------------------------------------------------------------------#
#                             Project Name                                     #
#------------------------------------------------------------------------------#
# The name of the project for which the pipeline will be executed. This should#
# be the same name used as the first parameter on init_sample.sh script        #
#------------------------------------------------------------------------------#
PROJECT: "test_metagenomes"
```

<ins>Samples</ins>

```yaml
#------------------------------------------------------------------------------#
#                               SAMPLES                                        #
#------------------------------------------------------------------------------#
# SAMPLES/Libraries you will like to include on the analysis                   #
# Same sample names used  with init_sample.sh script                           #
# Include all the names between quotes, and comma separated                    #
#------------------------------------------------------------------------------#
SAMPLES: ["sampleA", "sampleB"]
```

Go through the rest of the configuration file and choose your options before running the script. 

**Run the pipeline using SLURM**

>  sbatch hpc.sh

**Run the pipeline without SLURM**

*Activating environment*

>  module load anaconda/2024.02
>  conda activate /export/lv10/user/jberkhout/.conda/envs/snake_env_test
>  export GTDBTK_DATA_PATH="/export/lv13/databases/gtdb/release232"

*dry run*

> snakemake --configfile config.yaml -j2 -c35 --use-conda --conda-frontend conda -np

*Run*

> snakemake --configfile config.yaml  -j2 -c35 --use-conda --conda-frontend conda 

*Generating report file*
>  snakemake --configfile config.yaml --report report_name.zip


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


**Dependencies** 

* Assembly (you don't need to have all the tools installed, only the one for your target analysis)
  * Spades
  * Megahit
  * IDBA
* Mapping back reads to the assembly
  * BWA
* Binning (you don't need to have all the tools installed, only the one for your target analysis)
  * Maxbin
  * Metabat2
  * CONCOCT
  * Bin Sanity
  * DAS Tool
* Binning evaluation
  * CheckM
  * GTDB-Tk

