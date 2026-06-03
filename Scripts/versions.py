import subprocess
import sys

versions = open(snakemake.output[0],"w")
versions.write("Tool\tVersion\tInformation\n")

########################################################
#                Meta Cascabel version                 #
########################################################
version = "ND"
try:
  with open('resources/metacascabel.version') as f:
    version = f.readline().strip('\n')
except FileNotFoundError:
        print("Meta Cascabel's file version missing: ../cascabel.version\nYou can see README file for Cascabel version.")
        version = "see README file"
versions.write("Meta Cascabel\t"+version+"\tThis pipeline\n")

pythonv=str(sys.version_info[0])+"."+str(sys.version_info[1])+"."+str(sys.version_info[2])
versions.write("Python\t"+pythonv+"\tSeveral scripts and snakemake engine.\n")


#--fastq
fqVersion="ND"
if snakemake.config["QC"]["tool"].lower() ==  "fastqc" or snakemake.config["QC"]["tool"].lower() ==  "both":
    fqv = subprocess.run([snakemake.config["fastQC"]["command"], '--version'], stdout=subprocess.PIPE)
    fqVersion = fqv.stdout.decode('utf-8').strip() 
    versions.write("FastQC\t"+fqVersion+"\tQuality on fastq files\n")
sequaliVersion = "ND"
if snakemake.config["QC"]["tool"].lower() ==  "sequali" or snakemake.config["QC"]["tool"].lower() ==  "both":
    sqv = subprocess.run(['sequali', '--version'], stdout=subprocess.PIPE)
    sequaliVersion = sqv.stdout.decode('utf-8').strip()
    versions.write("Sequali\t"+sequaliVersion+"\tQuality on fastq files\n")

 
#spades
assemblyVersion="ND"
if snakemake.config["ASSEMBLER"].lower() ==  "spades" :
    av = subprocess.run(['spades.py', '--version'], stdout=subprocess.PIPE)
    assemblyVersion = av.stdout.decode('utf-8').strip().split('\n')[0] 
    versions.write("Spades\t"+assemblyVersion+"\tRead assembly\n")
elif snakemake.config["ASSEMBLER"].lower() ==  "idba" :
    versions.write("IDBA UD\t1.1.6\tRead assembly\n")
elif snakemake.config["ASSEMBLER"].lower() ==  "megahit" :
    av = subprocess.run(['megahit', '-v'], stdout=subprocess.PIPE)
    assemblyVersion = av.stdout.decode('utf-8').strip() 
    versions.write("Megahit\t"+assemblyVersion+"\tRead assembly\n")
#elif snakemake.config["ASSEMBLER"].lower() ==  "assembled" :
#    versions.write("Metagenome assembled\t1.1.6\tRead assembly\n")


biomVersion="ND"
biomV = subprocess.run([snakemake.config["biom"]["command"], '--version'], stdout=subprocess.PIPE)
biomVersion =  biomV.stdout.decode('utf-8').strip() 
versions.write("Biom\t"+biomVersion+"\tConvert tables to biom format.\n")

dada2Version="ND"
dada2V = subprocess.run([snakemake.config["Rscript"]["command"],'Scripts/dada2Version.R'], stdout=subprocess.PIPE)
dada2Version = dada2V.stdout.decode('utf-8').strip() 
versions.write("Dada2\t"+dada2Version+"\tdada2's QA & Error plots; ASV generation; Tax. asssignation (with RDP)\n")

if  snakemake.config["krona"]["report"].casefold() == "t" or snakemake.config["krona"]["report"].casefold() == "true":
    kVersion="ND"
    try:
        kv=subprocess.run([snakemake.config["krona"]["command"]], stdout=subprocess.PIPE)
        kVersion=kv.stdout.decode('utf-8').split('\n')[1].strip().replace("_","").replace("/","").replace("\\","").split("-")[0].strip()
    except Exception as e:
        kVersion="ND"
    versions.write("Krona Tools\t"+kVersion+"\tKrona interactive report.\n")




