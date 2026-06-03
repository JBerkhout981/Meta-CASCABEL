import sys
import subprocess
#run = sys.argv[0]
#samplesout = snakemake.wildcards.PROJECT + "/runs/" + snakemake.wildcards.run + "/samples.log"
#samplesout2 = snakemake.wildcards.PROJECT + "/runs/" + snakemake.wildcards.run + "/cat_samples.log"
bamlog=snakemake.output[0]
r1s=snakemake.input.r1
#r2=snakemake.input.r2
sample=snakemake.wildcards.sample
#assemblyIdx=snakemake.input.idx

with open(bamlog, "w") as logc:
    for r1 in r1s:
        r2=r1.replace('read1_paired','read2_paired')
        #idx_dir=idx.rsplit('/',1)[0]+"/"
        samp2=r1.split('/')[3].replace("_data","")
        outn=snakemake.wildcards.PROJECT + "/runs/" + snakemake.wildcards.run +"/"+ snakemake.wildcards.sample+"_data/bwa-mem/"+snakemake.config["ANALYSIS"]+"_"+snakemake.config["ASSEMBLER"]+"vs_"+samp2+"_mapped_against_cross-assembly_sorted.bam" 
        subprocess.run(["nice -"+ snakemake.config["bwa"]["nice"] +" bwa mem -t "+ snakemake.config["bwa"]["threads"] +" " + snakemake.params.idx + " " +r1 + " "+r2+" | samtools view --threads "+ snakemake.config["bwa"]["threads"] + " -b - | samtools sort - -o "+ outn + " --threads "+ snakemake.config["bwa"]["threads"]],stdout=subprocess.PIPE, shell=True)
        logc.write("nice -"+ snakemake.config["bwa"]["nice"] +" bwa mem -t "+ snakemake.config["bwa"]["threads"] +" " + snakemake.params.idx + " " +r1 + " "+r2+" | samtools view --threads "+ snakemake.config["bwa"]["threads"] + " -b - | samtools sort - -o "+ outn + " --threads "+ snakemake.config["bwa"]["threads"]+"\n")
#subprocess.run(["nice -"+ snakemake.config[config["bwa"]["nice"] +" bwa mem -t "+ snakemake.config["bwa"]["threads"] +" " + idx_dir + " " +r1 + " "+r2+" | samtools view --threads "+ snakemake.config["bwa"]["threads"] + " -b - | samtools sort - -o "+ outn + " --threads "+ snakemake.config["bwa"]["threads"]],stdout=subprocess.PIPE, shell=True)
      #  "{input.read1_paired} {input.read2_paired} "
      #  "| samtools view --threads {config[bwa][threads]} -b - | samtools sort - -o {output} --threads {config[bwa][threads]}"
#subprocess.run(["cutadapt -g "+ snakemake.config["primers"]["fw_primer"]  + " -G " + snakemake.config["primers"]["rv_primer"]  + " " +extra_params+" -O "+ snakemake.config["primers"]["min_overlap"]+" -m "+ snakemake.config["primers"]["min_length"] +" -o "+snakemake.params[0]+"/primer_removed_fw/"+sample+"_1.fastq.gz -p "+snakemake.params[0]+"/primer_removed_fw/"+sample+"_2.fastq.gz "+discard_untrimmed +" "+ fw_fq + " " +  rv_fq + " >> "+snakemake.params[0]+"/primer_removed_fw/"+sample+".cutadapt.log"],stdout=subprocess.PIPE, shell=True)

