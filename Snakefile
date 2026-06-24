"""
Metagenomics Workflow for NIOZ MMBL.
@Company: NIOZ                 
@Author: Alejandro Abdala and Julia Engelmann
@Version: 4.0     
@Last update: 16/04/2026                
"""

run=config["RUN"]
rule all:
    input:
        expand("{PROJECT}/runs/{run}/{sample}_data/report_f.html", PROJECT=config["PROJECT"],sample=config["SAMPLES"], run=run)

if len(config["SAMPLES"])==1 and len(config["fw_reads"])>0 and len(config["rv_reads"])>0:
    rule init_structure:
        input:
            fw = config["fw_reads"],
            rv = config["rv_reads"]
        output:
            r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
            r2="{PROJECT}/samples/{sample}/rawdata/rv.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz"
        shell:
            "Scripts/init_sample.sh "+config["PROJECT"]+" "+config["SAMPLES"][0]+"  {input.fw} {input.rv}"

elif len(config["input_files"])>2 and len(config["SAMPLES"])>=1:
    rule init_structure:
        input:
            file_list = config["input_files"]
        output:
            r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq"  if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
            r2="{PROJECT}/samples/{sample}/rawdata/rv.fastq"  if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz"
        benchmark:
            "{PROJECT}/samples/{sample}/benchmark/init_structure.benchmark"
        script:
            "Scripts/init_sample.py"

if config["QC"]["onRawReads"].lower() == "t":
    rule sequali:
        """
        Runs QC with sequali on raw reads
        """
        input:
            r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
            r2="{PROJECT}/samples/{sample}/rawdata/rv.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz"
        output:
            o1="{PROJECT}/samples/{sample}/qc/sequali/sequali.html",
            s2="{PROJECT}/samples/{sample}/qc/sequali/sequali.json",
            r=report(directory("{PROJECT}/samples/{sample}/qc/sequali/"), htmlindex="sequali.html", category="2. Quality", subcategory="Raw reads",labels={"sample":"{sample}","Tool":"Sequali", "Read":"fw & rv"},)
        params:
            outdir="{PROJECT}/samples/{sample}/qc/sequali/"
        threads:
            int(config["QC"]["threads"])
        benchmark:
            "{PROJECT}/samples/{sample}/benchmark/sequali.benchmark"
        shell:
            "sequali --outdir {params.outdir} --html  sequali.html --json sequali.json -t {config[QC][threads]}  {config[QC][extra_params]}  {input}"

# run trimmomatic to remove adapter contamination and trim very low quality parts (ends) of the reads.
# trimmomatic-0.35.jar PE -threads 2 $inFile1 $inFile2 $fol/read1_paired.fq $fol/read1_singles.fq $fol/read2_paired.fq $fol/read2_singles.fq
# ILLUMINACLIP:/usr/local/bioinf/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:12:2:TRUE MAXINFO:40:0.6 MINLEN:40
rule trimmomatic:
    input:
        fw="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
        rv="{PROJECT}/samples/{sample}/rawdata/rv.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz",
        tmp_qcr="{PROJECT}/samples/{sample}/qc/sequali/sequali.html",
    output:
        read1_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_paired.fq",
        read1_single="{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_singles.fq",
        read2_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_paired.fq",
        read2_single="{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_singles.fq",
        log="{PROJECT}/runs/{run}/{sample}_data/trimmed/trimmomatic.log"
    benchmark:
        "{PROJECT}/runs/{run}/{sample}_data/trimmed/trimmomatic.benchmark"
    threads:
        int(config["trimm"]["threads"])
    shell:
        #"java -jar /opt/biolinux/Trinity/trinity-plugins/Trimmomatic-0.36/trimmomatic-0.36.jar {config[trimm][mode]} -threads {config[trimm][threads]} {input.fw} {input.rv} "
        "trimmomatic {config[trimm][mode]} -threads {config[trimm][threads]} {input.fw} {input.rv} "
        "{output.read1_paired} {output.read1_single} {output.read2_paired} {output.read2_single} "
        "{config[trimm][clip][type]}:{config[trimm][clip][adapter]}:{config[trimm][clip][seed]}:{config[trimm][clip][palindrome_ct]}:"
        "{config[trimm][clip][simple_ct]}:{config[trimm][clip][minAdpLength]}:{config[trimm][clip][keepBoth]} "
        "{config[trimm][sliding][type]} "
        "{config[trimm][maxinfo][type]}:{config[trimm][maxinfo][targetLength]}:{config[trimm][maxinfo][strictness]} "
        "{config[trimm][minlen][type]}:{config[trimm][minlen][len]} > {output.log} 2>&1"

rule merge_trimmomatic_stats:
    input:
        expand("{PROJECT}/runs/{run}/{sample}_data/trimmed/trimmomatic.log", PROJECT=config["PROJECT"],sample=config["SAMPLES"], run=run)
    output:
        #"{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/stats_assembly.tsv"
        "{PROJECT}/runs/{run}/tables/stats_trimmomatic.tsv"
    benchmark:
        "{PROJECT}/runs/{run}/tables/stats_trimmomatic.benchmark"
    params:
        report_dir="{PROJECT}/runs/{run}/*_data/trimmed/trimmomatic.log"
    shell:
        #" echo -e \"Sample\\tAssembly\\tNum. contigs (>= 0 bp)\\tNum. contigs (>= 1000 bp)\\tNum. contigs (>= 5000 bp)\\tNum. contigs (>= 10000 bp)\\tNum. contigs (>= 25000 bp)\\tNum. contigs (>= 50000 bp)\\tTotal length (>= 0 bp)\\tTotal length (>= 1000 bp)\\tTotal length (>= 5000 bp)\\tTotal length (>= 10000 bp)\\tTotal length (>= 25000 bp)\\tTotal length (>= 50000 bp)\\tNum. contigs\\tLargest contig\\tTotal length\\tGC (%)\\tN50\\tN90\\tauN\\tL50\\tL90\\tNum. N's per 100 kbp\" > {output} ;"
        " echo -e \"Sample\\tInput Read Pairs\\tBoth Surviving\\tBoth %\\tForward Only Surviving\\tForward Only %\\tReverse Only Surviving\\tReverse Only %\\tDropped\\tDropped %\" > {output} ;"
        " for file in {params.report_dir} ; "
        " do "
        "   sample=$(echo $file | awk -F\"/\" '{{gsub(\"_data\",\"\",$4); print $4}}');"
        "   cat $file | grep  \"^Input Read\" | sed 's/(//g ; s/)//g' | awk -v samp=${{sample}} 'BEGIN{{OFS=\"\\t\"}} {{print samp,$4,$7,$8,$12,$13,$17,$18,$20,$21 >> \"{output}\"}}';"
        " done " 

rule create_yaml_trimmomatic_stats_tbl:
    output:
        "{PROJECT}/runs/{run}/tables/trimmomatic.yaml"
    shell:
        "cp resources/datavzrd/trimmomatic.yaml {output}"

rule datavzrd_trimmomatic:
    input:
        config="{PROJECT}/runs/{run}/tables/trimmomatic.yaml",
        table="{PROJECT}/runs/{run}/tables/stats_trimmomatic.tsv"
        #table="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/stats_assembly.tsv"
    output:
        report(
            directory("{PROJECT}/runs/{run}/tables/trimmomatic"),
            #directory("{PROJECT}/runs/{run}/tables/assembly_"+config["ASSEMBLER"]+"/quast/tbl/"),
            htmlindex="index.html",
            #category="1. Project Information",
            category="3. Read trimming",
            #subcategory="Assembler: " + config["ASSEMBLER"],
            labels={"table":"Trimming results"},
        ),
    wrapper:
        "v4.7.2/utils/datavzrd"

if config["QC"]["onTrimmedReads"].lower() == "t":
    rule sequali_trimmed_reads:
        input:
            r1="{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_paired.fq",
            r2="{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_paired.fq"
        output:
            o1="{PROJECT}/runs/{run}/{sample}_data/trimmed/sequali/sequali.html",
            s2="{PROJECT}/runs/{run}/{sample}_data/trimmed/sequali/sequali.json",
            r=report(directory("{PROJECT}/runs/{run}/{sample}_data/trimmed/sequali/"), htmlindex="sequali.html", category="2. Quality", subcategory="Trimmed reads",labels={"sample":"{sample}","Tool":"Sequali", "Read":"fw & rv"},)
        params:
            outdir="{PROJECT}/runs/{run}/{sample}_data/trimmed/sequali/"
        benchmark:
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/sequali.benchmark"
        threads:
            int(config["QC"]["threads"])
        shell:
            "sequali --outdir {params.outdir} --html  sequali.html --json sequali.json -t {config[QC][threads]}  {config[QC][extra_params]}  {input}"

if config["TAXONOMY"]["PROFILING"] == "KRAKEN" or config["TAXONOMY"]["PROFILING"] == "ALL":
    rule kraken:
        """
            Execute Kraken taxonomy profiling
        """
        input:
            fw="{PROJECT}/samples/{sample}/rawdata/fw.fastq"
            if config["TAXONOMY"]["KRAKEN"]["raw_reads"] == "Y" else "{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_paired.fq",
            rv="{PROJECT}/samples/{sample}/rawdata/rv.fastq"
            if config["TAXONOMY"]["KRAKEN"]["raw_reads"] == "Y" else "{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_paired.fq"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kraken.taxonomy.out"
        threads:
            int(config["TAXONOMY"]["KRAKEN"]["threads"])
        shell:
            "kraken --preload --db {config[TAXONOMY][KRAKEN][db]} --paired {input.fw} {input.rv} "
            "--threads {config[TAXONOMY][KRAKEN][threads]} {config[TAXONOMY][KRAKEN][extra_params]} > {output}"
#--gzip-compressed
    rule prepare_kraken_report:
        """
            Prepare the input file for addTaxonNames script
        """
        input:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kraken.taxonomy.out"
        output:
            temp("{PROJECT}/runs/{run}/{sample}_data/taxonomy/kraken.taxonomy.out.tmp")
        shell:
            "cat {input} | cut -f1,2,3 > {output}"
    rule kraken_labels:
        """
            Prepare the input file for addTaxonNames script
        """
        input:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kraken.taxonomy.out.tmp"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kraken.taxonomy.out.labels"
        shell:
            "kaiju-addTaxonNames -t {config[TAXONOMY][KRAKEN][nodes]} -n {config[TAXONOMY][KRAKEN][names]} "
            "-i {input} {config[TAXONOMY][taxonomy_path]}  -o {output}"
    rule kraken_report:
        """
            Prepare the input file for addTaxonNames script
        """
        input:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kraken.taxonomy.out.labels"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kraken.taxonomy.report"
        shell:
            "kaiju2table -t {config[TAXONOMY][KRAKEN][nodes]} -n {config[TAXONOMY][KRAKEN][names]} "
            " {config[TAXONOMY][taxonomy_path]}  -o {output} {input}"

if config["TAXONOMY"]["PROFILING"] == "KAIJU" or config["TAXONOMY"]["PROFILING"] == "ALL":
    rule kaiju:
        """
            Execute Kaiju taxonomy profiling
        """
        input:
            fw="{PROJECT}/samples/{sample}/rawdata/fw.fastq"
            if config["TAXONOMY"]["KRAKEN"]["raw_reads"] == "Y" else "{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_paired.fq",
            rv="{PROJECT}/samples/{sample}/rawdata/rv.fastq"
            if config["TAXONOMY"]["KRAKEN"]["raw_reads"] == "Y" else "{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_paired.fq"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kaiju.taxonomy.out"
        threads:
            int(config["TAXONOMY"]["KAIJU"]["threads"])
        shell:
            "kaiju -i {input.fw} -j {input.rv} "
            " -t {config[TAXONOMY][KAIJU][nodes]}  -f {config[TAXONOMY][KAIJU][db]} "
            "-z {config[TAXONOMY][KAIJU][threads]} {config[TAXONOMY][KAIJU][extra_params]} -o {output}"
    rule kaiju_labels:
        """
            addTaxonLabels for kaiju
        """
        input:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kaiju.taxonomy.out"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kaiju.taxonomy.out.labels"
        shell:
            "kaiju-addTaxonNames -t {config[TAXONOMY][KAIJU][nodes]} -n {config[TAXONOMY][KAIJU][names]} "
            "-i {input} {config[TAXONOMY][taxonomy_path]}  -o {output}"
    rule kaiju_report:
        """
            addTaxonLabels for kaiju
        """
        input:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kaiju.taxonomy.out.labels"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kaiju.taxonomy.out.report"
        shell:
           "kaiju2table -t {config[TAXONOMY][KAIJU][nodes]} -n {config[TAXONOMY][KAIJU][names]} "
            " {config[TAXONOMY][taxonomy_path]}  -o {output} {input}"
           #old kaiju version
           # "kaijuReport -t {config[TAXONOMY][KAIJU][nodes]} -n {config[TAXONOMY][KAIJU][names]} "
           # "-i {input} {config[TAXONOMY][taxonomy_path]}  -o {output}"

if config["TAXONOMY"]["PROFILING"] not in "KRAKEN KAIJU ALL":
    rule create_taxo_out:
        """
            As there is no taxonomy profiling touch one file to generate a "silly" file
        """
        output:
            "{PROJECT}/runs/{run}/{sample}_data/no_tax.txt"
        shell:
            "touch {output}"
elif config["TAXONOMY"]["PROFILING"] == "ALL":
    rule merge_taxonomy_outs:
        """
            Merge the output of the three taxonomic profiler into one single outfile
        """
        input:
            kraken="{PROJECT}/runs/{run}/{sample}_data/taxonomy/kraken.taxonomy.out",
            kaiju="{PROJECT}/runs/{run}/{sample}_data/taxonomy/kaiju.taxonomy.out"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/all.taxonomy.out"
        shell:
            "touch {output}"




if config["ASSEMBLER"] == "SPADES":
    rule fq2fasta:
        input:
            read1_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_paired.fq",
            read2_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_paired.fq",
        output:
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/reads_merged.fasta" if config["gzip_input"] == "F"
            else "{PROJECT}/runs/{run}/{sample}_data/trimmed/reads_merged.fastq.gz"
        benchmark:
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/fq2fasta.benchmark"
        shell:
            "fq2fa --merge {input.read1_paired} {input.read2_paired} {output}"

    rule concat_single_reads:
        input:
            read1_single="{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_singles.fq",
            read2_single="{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_singles.fq",
            tmp_seq="{PROJECT}/runs/{run}/{sample}_data/trimmed/sequali/sequali.html"
              
        output:
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/all_singles.fq"
        benchmark:
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/concat_single_reads.benchmark"
        shell:
            "cat {input.read1_single} {input.read2_single} > {output}"
    #This rule is ment for running spades with already merged reads, and this we need the 
    # merged files created above, but for the moment just using meta_spades rules, that works directly
    #with trimmomatic output. The above procedure is also used for IDBA UD. 
    rule meta_spades_to_fix:
        input:
            reads_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/reads_merged.fasta" if config["gzip_input"] == "F"
            else "{PROJECT}/runs/{run}/{sample}_data/trimmed/reads_merged.fastq.gz",
            read12_singles="{PROJECT}/runs/{run}/{sample}_data/trimmed/all_singles.fq"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.fasta_tmp",
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffolds.fasta_tmp"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/"
        benchmark:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/assembly.benchmark"
        shell:
            "nice -{config[spades][nice]} spades.py --meta -t {config[spades][threads]} -m {config[spades][memory]} "
            "-k {config[spades][kmers]} --12 {input.reads_paired} -s {input.read12_singles} "
            "{config[spades][extra_params]} -o {params}"


    rule meta_spades:
        input:
            read1_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_paired.fq",
            read2_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_paired.fq",
            read12_singles="{PROJECT}/runs/{run}/{sample}_data/trimmed/all_singles.fq"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.fasta",
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffolds.fasta"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/"
        benchmark:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/assembly.benchmark"
        threads:
            int(config["spades"]["threads"])
        shell:
            "nice -{config[spades][nice]} spades.py --meta -t {config[spades][threads]} -m {config[spades][memory]} "
            "-k {config[spades][kmers]} --pe1-1 {input.read1_paired} --pe1-2 {input.read2_paired} --pe1-s {input.read12_singles} "
            "{config[spades][extra_params]} -o {params}"

    rule std_assembly_meta_spades:
        input:
            contigs="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.fasta",
            scaffolds="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffolds.fasta"
        output:
            contigs="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_contigs.fasta",
            scaffolds="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_scaffolds.fasta"
        shell:
            "mv {input.contigs} {output.contigs} && mv {input.scaffolds} {output.scaffolds}"

if config["ASSEMBLER"] == "MEGAHIT":
    rule megahit:
        input:
            read1_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_paired.fq",
            read2_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_paired.fq",
            tmp_seq="{PROJECT}/runs/{run}/{sample}_data/trimmed/sequali/sequali.html"
            
        output:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/final.contigs.fa"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]
        threads:
            int(config["megahit"]["cpus"])
        shell:
            "megahit -1 {input.read1_paired} -2 {input.read2_paired} -f --k-min {config[megahit][kmin]} "
            "--k-max {config[megahit][kmax]} --k-step {config[megahit][kstep]} {config[megahit][extra_params]} "
            "-m {config[megahit][memory]} -t {config[megahit][cpus]} -o {params}"
    rule std_assembly_megahit:
        input:
            contig="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/final.contigs.fa"
        output:
            contigs="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_contigs.fasta",
            scaffolds="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_scaffolds.fasta"
        shell:
            "mv {input.contig} {output.contigs} && ln -sr {output.contigs} {output.scaffolds}"
    #rule std_assembly_megahit_step2:
    #    input:
    #        contig="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.fasta"
    #    output:
    #        scaffolds="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffolds.fasta"
    #    shell:
    #        "ln -sr {input.contig} scaffolds.fasta"

if config["ASSEMBLER"] == "IDBA":
    rule fq2fasta:
        input:
            read1_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_paired.fq",
            read2_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_paired.fq",
            tmp_seq="{PROJECT}/runs/{run}/{sample}_data/trimmed/sequali/sequali.html"
            
        output:
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/reads_merged.fasta"
        benchmark:
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/fq2fasta.benchmark"
        shell:
            "fq2fa --merge {input.read1_paired} {input.read2_paired} {output}"
    #IN order to run idba it is needed to make somechanges into the source code:
    #https://groups.google.com/forum/#!topic/hku-idba/NE2JXqNvTFY and
    #http://seqanswers.com/forums/showthread.php?t=29109
    #change
    #static const uint32_t kMaxShortSequence = 128;
    #by
    #static const uint32_t kMaxShortSequence = 256;
    #For this reason, the pipe line uses my local installation, try to make it for all
    rule idba:
        input:
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/reads_merged.fasta"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contig.fa",
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffold.fa"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/"
        benchmark:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/assembly.benchmark"
        threads:
            int(config["idba"]["threads"])
        shell:
            "idba_ud -r {input} -o {params} "
            "--step {config[idba][step]} --num_threads {config[idba][threads]} {config[idba][extra_params]}"
    #As spades output contigs.fasta and scaffolds.fasta, we standarize those
    #names in order to decress complexity on downstream rules
    rule std_assembly_idba:
        input:
            contig="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contig.fa",
            scaffold="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffold.fa"
        output:
            contigs="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_contigs.fasta",
            scaffolds="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_scaffolds.fasta"
        shell:
            "mv {input.contig} {output.contigs} && mv {input.scaffold} {output.scaffolds}"

if config["ASSEMBLER"] == "ASSEMBLED":
    rule std_assembly:
        input:
            contigs=config["contigs"]
        output:
            contigs="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_contigs.fasta",
            scaffolds="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_scaffolds.fasta"
        shell:
            "ln -s {input.contigs} {output.contigs} && ln -s {input.contigs} {output.scaffolds}"



rule quast_libs:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_contigs.fasta"
    output:
        temp("{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/quast_lib_config.txt")
    shell:
        "export PYTHONPATH=/opt/Downloads/biolinux/quast-4.6.3/ && touch {output}"

rule quast_contigs:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_contigs.fasta"
        #"{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/quast_lib_config.txt"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/contigs/report.txt"
    params:
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/contigs/"
    benchmark:
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/contigs/quast_contigs.benchmark"
    threads:
        int(config["quast"]["threads"])
    shell:
        "quast.py -t {config[quast][threads]} -o {params} {config[quast][extra_params]} {input[0]}"

rule quast_scaffolds:
    input:
        scaffolds="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_scaffolds.fasta"
        #tmp_i="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/quast_lib_config.txt"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/scaffolds/report.txt"
    params:
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/scaffolds/"
    benchmark:
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/scaffolds/quast_scaffolds.benchmark"
    threads:
        int(config["quast"]["threads"])
    shell:
        "quast.py -t {config[quast][threads]} -o {params} -s {config[quast][extra_params]}  {input.scaffolds}"



if config["SPLIT_ASSEMBLY"] == "T":
    """
    This option will split the assembled reads into smaller chuncks.
    In order to do not affect the downstream rules (previously implemented)
    the splitted contig file is renamed as the original contig file and the original
    is renamed as contigs.complete.fasta at the end
    """
    rule split_assembly:
        input:
            contigs="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_contigs.fasta"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_contigs.chunks.fasta"
        shell:
            # "/opt/biolinux/anaconda2.2019.07/bin/cut_up_fasta.py -c {config[SPLIT_SIZE]} "
            # This pathway currently only works for ada, not ada94
            # "/opt/biolinux/anaconda/73/2022.05/envs/metacascabel_c_env/bin/cut_up_fasta.py  -c {config[SPLIT_SIZE]} "
            "cut_up_fasta.py  -c {config[SPLIT_SIZE]} -o 0 -m {input.contigs} > {output}"
    rule std_splitted_assembly:
        input:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_contigs.chunks.fasta"
        output:
            flag="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/split_flag.txt",
            complete="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_contigs.complete.fasta"
        params:
            contigs="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_contigs.fasta"
        shell:
            "mv {params.contigs} {output.complete} && mv {input} {params.contigs} "
            "&& echo \"contigs.fasta has been splitted by cut_up_fasta.py original contig file is: contigs.complete.fasta\" > {output.flag}"
else:
    rule skip_split_assembly:
        output:
            flag="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/split_flag.txt"
        shell:
            "echo \"contigs.fasta has not been splitted\" > {output.flag}"

rule validate_assembly:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/scaffolds/report.txt",
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/contigs/report.txt",
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/split_flag.txt"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/validate_assembly.txt"
    benchmark:
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/validate_assembly.benchmark"
    params:
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/"
    script:
        "Scripts/validateQuast.py"

#TBC: cat Diff_matrix/runs/all_vs_all/NIOZ114_data/assembly_SPADES/quast/scaffolds/report.tsv  | awk -F"\t" 'BEGIN{print "Category\tValue"}{gsub("#","Number of",$0); print $0}'

rule merge_assembly_stats:
    input:
        expand("{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/validate_assembly.txt", PROJECT=config["PROJECT"],sample=config["SAMPLES"], run=run)
    output:
        #"{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/stats_assembly.tsv"
        "{PROJECT}/runs/{run}/tables/stats_assembly.tsv"
    benchmark:
        "{PROJECT}/runs/{run}/tables/stats_assembly.benchmark"
    params:
        report_dir="{PROJECT}/runs/{run}/*_data/assembly_"+config["ASSEMBLER"]+"/quast/*/transposed_report.tsv"
    shell:
        " echo -e \"Sample\\tAssembly\\tNum. contigs (>= 0 bp)\\tNum. contigs (>= 1000 bp)\\tNum. contigs (>= 5000 bp)\\tNum. contigs (>= 10000 bp)\\tNum. contigs (>= 25000 bp)\\tNum. contigs (>= 50000 bp)\\tTotal length (>= 0 bp)\\tTotal length (>= 1000 bp)\\tTotal length (>= 5000 bp)\\tTotal length (>= 10000 bp)\\tTotal length (>= 25000 bp)\\tTotal length (>= 50000 bp)\\tNum. contigs\\tLargest contig\\tTotal length\\tGC (%)\\tN50\\tN90\\tauN\\tL50\\tL90\\tNum. N's per 100 kbp\" > {output} ;"
        " for file in `ls {params.report_dir}`;  "
        "  do "
        "    sample=$(echo $file | awk -F'/' '{{gsub(\"_data\",\"\",$4); print $4}}'); "
        "    cat $file | awk -F\"\\t\" -v s=$sample 'NR>1{{l=s;for(i=1;i<=NF;i++){{l=l\"\\t\"$i}}; print l >> \"{output}\"}}'; "
        " done "
        #" echo -e \"Sample\\tAssembly\\tNum. contigs (>= 0 bp)\\tNum. contigs (>= 1000 bp)\\tNum. contigs (>= 5000 bp)\\tNum. contigs (>= 10000 bp)\\tNum. contigs (>= 25000 bp)\\tNum. contigs (>= 50000 bp)\\tTotal length (>= 0 bp)\\tTotal length (>= 1000 bp)\\tTotal length (>= 5000 bp)\\tTotal length (>= 10000 bp)\\tTotal length (>= 25000 bp)\\tTotal length (>= 50000 bp)\\tNum. contigs\\tLargest contig\\tTotal length\\tGC (%)\\tN50\\tN90\\tauN\\tL50\\tL90\\tNum. N's per 100 kbp\" > {output} ;"
        #" for file in `ls Diff_matrix/runs/all_vs_all/*_data/assembly_SPADES/quast/*/transposed_report.tsv`;  "
        #"  do "
        #"    sample=$(echo $file | awk -F'/' '{{gsub(\"_data\",\"\",$4); print $4}}'); "
        #"    cat $file | awk -F\"\\t\" -v s=$sample 'NR>1{{l=s;for(i=1;i<=NF;i++){{l=l\"\\t\"$i}}; print l >> \"{output}\"}}'; "
        #" done "

rule create_yaml_assembly_stats_tbl:
    output:
        "{PROJECT}/runs/{run}/tables/assembly.yaml"
    shell:
        "cp resources/datavzrd/assembly.yaml {output}"

rule datavzrd_assembly:
    input:
        config="{PROJECT}/runs/{run}/tables/assembly.yaml",
        table="{PROJECT}/runs/{run}/tables/stats_assembly.tsv"
        #table="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/stats_assembly.tsv"
    output:
        report(
            directory("{PROJECT}/runs/{run}/tables/assembly_"+config["ASSEMBLER"]+"/quast/"),
            #directory("{PROJECT}/runs/{run}/tables/assembly_"+config["ASSEMBLER"]+"/quast/tbl/"),
            htmlindex="index.html",
            #category="1. Project Information",
            category="4. Assembly",
            #subcategory="Assembler: " + config["ASSEMBLER"],
            labels={"Assembler": ""+ config["ASSEMBLER"]},
        ),
    wrapper:
        "v4.7.2/utils/datavzrd"

  #cp gm_key_64 ~/.gm_key
rule bwa_index:
    input:
        #tmp_flw="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/validate_assembly.txt",
        tmp_flw="{PROJECT}/runs/{run}/tables/assembly_"+config["ASSEMBLER"]+"/quast/",
        assembly="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_scaffolds.fasta"
        if config["ANALYSIS"] == "SCAFFOLDS" else "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_contigs.fasta"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_assembly.bwt"
    benchmark:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/bwaindex.benchmark"
    params:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_assembly"
    shell:
        "bwa index -p {params} {input.assembly}"
#nice -5 bwa mem -t 4 $folderOut/cross-assembly $fq1 $fq2 \
#   | samtools view -b - \
#   | samtools sort - -o $fol/mapped_against_cross-assembly_sorted.bam
rule bwa_mem:
    input:
        read1_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_paired.fq",
        read2_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_paired.fq",
        tmp_flw="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_assembly.bwt"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_mapped_against_cross-assembly_sorted.bam"
    benchmark:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/bwamem.benchmark"
    params:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_assembly"
    threads:
        int(config["bwa"]["threads"])
    shell:
        "nice -{config[bwa][nice]} bwa mem -t {config[bwa][threads]} {params} "
        "{input.read1_paired} {input.read2_paired} "
        "| samtools view --threads {config[bwa][threads]} -b - | samtools sort - -o {output} --threads {config[bwa][threads]}"

rule bwa_mem_new:
    input:
        r1="{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_paired.fq",
        r2="{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_paired.fq",
        allidx=expand("{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_assembly.bwt", PROJECT=config["PROJECT"],sample=config["SAMPLES"], run=run)
        #tmp_flw="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_assembly.bwt"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_mapped_against_cross-assembly_sorted.bam___"
    benchmark:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/bwamem.benchmark"
    params:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_assembly"
    threads:
        int(config["bwa"]["threads"])
    shell:
        "nice -{config[bwa][nice]} bwa mem -t {config[bwa][threads]} {params} "
        "{input.read1_paired} {input.read2_paired} "
        "| samtools view --threads {config[bwa][threads]} -b - | samtools sort - -o {output} --threads {config[bwa][threads]}"

rule bwa_mem_mtx:
    input:
        assembly="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_scaffolds.fasta"
        if config["ANALYSIS"] == "SCAFFOLDS" else "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_contigs.fasta",
        r1=expand("{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_paired.fq", PROJECT=config["PROJECT"],sample=config["SAMPLES"], run=run),
        r2=expand("{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_paired.fq",  PROJECT=config["PROJECT"],sample=config["SAMPLES"], run=run),
        idx="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_assembly.bwt"
        #tmp_flw="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_assembly.bwt"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/bwa-mem_cmds.log"
    benchmark:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/bwamem.benchmark"
    params:
        idx="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_assembly"
    threads:
        int(config["bwa"]["threads"])
    script:
        "Scripts/bwa_mem.py"

rule sam_flags:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_mapped_against_cross-assembly_sorted.bam"
        if config["bwa"]["differential_coverage_matrix"].lower() == "f" else
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/bwa-mem_cmds.log"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_mapped_against_cross-assembly_sorted.flagstat"
        if config["bwa"]["differential_coverage_matrix"].lower() == "f" else
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/sam-flagstat_cmds.log"
    benchmark:
         "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/samflags.benchmark"
    threads:
        int(config["bwa"]["threads"])
    params:
        idx="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"
    shell:
        "samtools flagstat --threads {config[bwa][threads]} {input} > {output}"
        if config["bwa"]["differential_coverage_matrix"].lower() == "f" else
        "for sample in `ls {params.idx}"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"*_sorted.bam`; do "
        "  out=$(echo $sample | sed 's/sorted.bam/sorted.flagstat/');"
        "  samtools flagstat --threads {config[bwa][threads]} ${{sample}} > $out ;"
        "  echo 'samtools flagstat --threads {config[bwa][threads]} ${{sample}} > $out' >> {output} ;"
        " done"
         #"samtools idxstats ${sample} > ${sample}_mapped.txt ; done"

rule summarize_bam:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_mapped_against_cross-assembly_sorted.bam"
        if config["bwa"]["differential_coverage_matrix"].lower() == "f" else
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/bwa-mem_cmds.log"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth.txt"
    benchmark:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/jgi_summ.benchmark"
    params:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+""
    shell:
        "jgi_summarize_bam_contig_depths --outputDepth {output} {config[jgi_summ][extra_params]} {input}"
        if config["bwa"]["differential_coverage_matrix"].lower() == "f" else
        "jgi_summarize_bam_contig_depths --outputDepth {output} {config[jgi_summ][extra_params]} {params}*sorted.bam"

rule merge_mapping_stats:
    input:
        expand("{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_mapped_against_cross-assembly_sorted.flagstat"
        if config["bwa"]["differential_coverage_matrix"].lower() == "f" else
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/sam-flagstat_cmds.log", PROJECT=config["PROJECT"],sample=config["SAMPLES"], run=run)
    output:
        #"{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/stats_assembly.tsv"
        "{PROJECT}/runs/{run}/tables/stats_bwa.tsv"
    benchmark:
        "{PROJECT}/runs/{run}/tables/stats_bwa.benchmark"
    params:
        report_dir="{PROJECT}/runs/{run}/*_data/bwa-mem/*.flagstat"
    shell:
        #" echo -e \"Sample\\tAssembly\\tNum. contigs (>= 0 bp)\\tNum. contigs (>= 1000 bp)\\tNum. contigs (>= 5000 bp)\\tNum. contigs (>= 10000 bp)\\tNum. contigs (>= 25000 bp)\\tNum. contigs (>= 50000 bp)\\tTotal length (>= 0 bp)\\tTotal length (>= 1000 bp)\\tTotal length (>= 5000 bp)\\tTotal length (>= 10000 bp)\\tTotal length (>= 25000 bp)\\tTotal length (>= 50000 bp)\\tNum. contigs\\tLargest contig\\tTotal length\\tGC (%)\\tN50\\tN90\\tauN\\tL50\\tL90\\tNum. N's per 100 kbp\" > {output} ;"
        " echo -e \"Assembly\\tRaw reads\\tReads\\tMapped reads\\tMapped percentage\\tPaired\\tPaired prc\" > {output} ;"
        " for file in {params.report_dir} ; "
        " do "
        "   sample=$(echo $file | awk -F\"/\" '{{gsub(\"_data\",\"\",$4); split($NF,a,\"_\"); print $4\"\\t\"a[3]}}');  "
        "   cat $file | awk -v s=\"${{sample}}\" 'BEGIN{{OFS=\"\\t\"}} {{if(NR==1){{tot=$1}}else if(NR==7){{map=$1;gsub(/\\(|\\)/,\"\",$0);mapprc=$5}}else if(NR==12){{gsub(/\\(|\\)/,\"\",$0);paired=$1;pairedprc=$6}}}} END{{print s,tot,map,mapprc,paired,pairedprc >> \"{output}\"}}';"
        " done"
       # "   sample=$(echo $file | awk -F\"/\" '{{gsub(\"_data\",\"\",$4); print $4}}');"
       # "   cat $file | grep  \"^Input Read\" | sed 's/(//g ; s/)//g' | awk -v samp=${{sample}} '{{print samp,$4,$7,$8,$12,$13,$17,$18,$20,$21 >> \"{output}\"}}';"
       # " done " 

rule create_yaml_bwa_stats_tbl:
    output:
        "{PROJECT}/runs/{run}/tables/bwa.yaml"
    shell:
        "cp resources/datavzrd/bwa.yaml {output}"

rule datavzrd_bwa:
    input:
        config="{PROJECT}/runs/{run}/tables/bwa.yaml",
        table="{PROJECT}/runs/{run}/tables/stats_bwa.tsv"
        #table="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/stats_assembly.tsv"
    output:
        report(
            directory("{PROJECT}/runs/{run}/tables/bwa"),
            #directory("{PROJECT}/runs/{run}/tables/assembly_"+config["ASSEMBLER"]+"/quast/tbl/"),
            htmlindex="index.html",
            #category="1. Project Information",
            category="5. Read Mapping",
            #subcategory="Assembler: " + config["ASSEMBLER"],
            #labels={"table", "Trimming results"},
        ),
    wrapper:
        "v4.7.2/utils/datavzrd"


rule total_coverage:
    '''
    This rule prepare files for concoct and maxbin (no diff mtx only in maxbin), by taking the sample 
    and the total coverage, not the "totalAverage" which is col 3, now it takes col 4.   
    '''
    input:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth.txt"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth_avg.txt"
    shell:
        "awk -F'\\t' 'NR>1{{print $1\"\\t\"$4}}' {input} > {output}"
        if config["bwa"]["differential_coverage_matrix"].lower() == "f" else
        "awk -F'\\t' 'NR > 1 {{for(x=1;x<=NF;x++) if(x == 1 || (x >= 4 && x % 2 == 0)) printf \"%s\", $x (x == NF || x == (NF-1) ? \"\\n\":\"\\t\")}}'  {input} > {output}"
rule maxbin_coverage:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth.txt"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth_avg_maxbin.txt"
    shell:
        "awk -F '\\t' -v col=\"CONTIGS_SPADESvs_{wildcards.sample}_mapped_against_cross-assembly_sorted.bam\" 'NR==1{{for (i=1; i<=NF; i++) if ($i == col){{c=i; break}}}} NR>1{{print $1\"\\t\"$c}}' {input} > {output}"


if config["BINNING"] == "METABAT" or config["BINNING"] == "DAS":
    rule metabat:
        input:
            depth="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth.txt",
            assembly="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_scaffolds.fasta"
            if config["ANALYSIS"] == "SCAFFOLDS" else "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_contigs.fasta"
        output:
        #"{PROJECT}/runs/{run}/{sample}_data/metabat2/bin/metabat2.log",
        #"{PROJECT}/runs/{run}/{sample}_data/metabat2/bin{number}.fa",
        #expand("{PROJECT}/runs/{run}/{sample}_data/metabat2/bin{{number}}.fa",PROJECT=config["PROJECT"],sample=config["SAMPLES"], run=run),
        #dynamic("{PROJECT}/runs/{run}/{sample}_data/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/bin.{number}.fa")
        #"{PROJECT}/runs/{run}/{sample}_data/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/bin.{number}.fa"
            "{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/metabat.log"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/bin"
        benchmark:
            "{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/metabat.benchmark"
        threads:
            int(config["metabat2"]["threads"])
        shell:
            "metabat2 -o {params} -i {input.assembly} -t {config[metabat2][threads]} "
            "-m {config[metabat2][min_contig]} -a {input.depth} "
            "--minS {config[metabat2][min_score]} --maxP {config[metabat2][maxP]} "
            "--maxEdges {config[metabat2][maxEdge]} -s {config[metabat2][min_bin_size]} "
            "{config[metabat2][extra_params]}  > {output}"
if config["BINNING"] == "MAXBIN" or (config["BINNING"] == "DAS" and config["das"]["maxbin"]["run"]=="T" ):
    rule maxbin:
        """
        Please make sure that your abundance information is provided in the following format:
        (contig header)\t(abundance)
        """
        input:
            depth="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth_avg.txt"
            if config["bwa"]["differential_coverage_matrix"].lower() == "f" else
            "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth_avg_maxbin.txt",
            assembly="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_scaffolds.fasta"
            if config["ANALYSIS"] == "SCAFFOLDS" else "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_contigs.fasta"
        output:
            log="{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/maxbin.log"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/bin"
        benchmark:
            "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/maxbin.benchmark"
        threads:
            int(config["maxbin"]["threads"])
        shell:
            "run_MaxBin.pl -contig {input.assembly} "
            "-abund  {input.depth} -out {params} -thread {config[maxbin][threads]} "
            "-prob_threshold {config[maxbin][prob_threshold]} -markerset {config[maxbin][markerset]} "
            "-min_contig_length {config[maxbin][min_contig_length]}  "
            "{config[maxbin][plotmarker]} {config[maxbin][extra_params]} > {output.log}"

elif config["BINNING"] == "DAS" and config["das"]["maxbin"]["run"]!="T":
   rule skip_maxbin:
        output:
            log="{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/maxbin.log"
        shell:
            "touch {output}"

if config["BINNING"] == "CONCOCT" or ( config["BINNING"] == "DAS"  and config["das"]["concoct"]["run"]=="T"):

    rule activate_concoct:
        output:
            concoct_activation="{PROJECT}/runs/{run}/concoct_activation.log"
        shell:
            "{config[concoct][activation_cmd]} && echo {config[concoct][activation_cmd]} > {output}"
    """
    This rules try to follow the steps recommended by using CONCOCT according to
    the following documentation: https://concoct.readthedocs.io/en/latest/complete_example.html
    """
    rule concoct:
        input:
            depth="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth_avg.txt",
            assembly="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_scaffolds.fasta"
            if config["ANALYSIS"] == "SCAFFOLDS" else "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_contigs.fasta"
            #concoct_activation="{PROJECT}/runs/{run}/concoct_activation.log"
        output:
            #log="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/concoct.log",
            clustering="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/bin_clustering_gt"+config["concoct"]["min_contig_length"]+".csv"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/bin"
        threads:
            int(config["concoct"]["threads"])
        shell:
            "concoct -l {config[concoct][min_contig_length]} -i {config[concoct][max_iteration]} -t {config[concoct][threads]} "
            "--coverage_file {input.depth} --composition_file {input.assembly} {config[concoct][extra_params]}   -b {params}  > /dev/null 2>&1"
    """
    We only need the bins in case that the user is running concoct alone, otherwise
    we only use the clustering file to DAS to create new bins. NOT anymore! now we run 
    checkM for all the methods
    """
#    if config["BINNING"] == "CONCOCT":
    rule extract_concoct_bins:
        input:
            assembly="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_scaffolds.fasta"
            if config["ANALYSIS"] == "SCAFFOLDS" else "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_contigs.fasta",
            clustering="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/bin_clustering_gt"+config["concoct"]["min_contig_length"]+".csv"
        output:
            log="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/concoct.log"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/"
        shell:
            #"/usr/local/bioinf/concoct-0.4/scripts/extract_fasta_bins.py  --output_path  {params} {input.assembly} {input.clustering} > {output.log}"
            #"/export/data/aabdala/utils/CONCOCT/scripts/extract_fasta_bins.py  --output_path  {params} {input.assembly} {input.clustering} > {output.log}"
            #"python /opt/biolinux/anaconda2.2019.07/pkgs/concoct-1.1.0-py27h88e4a8a_0/bin/extract_fasta_bins.py  --output_path  {params} {input.assembly} {input.clustering} > {output.log}"
            #"python /export/data01/tools/concoct_scripts/concoct-0.4/extract_fasta_bins.py  --output_path  {params} {input.assembly} {input.clustering} > {output.log}"
            #"/opt/biolinux/anaconda/73/3.2020.11/envs/concoct/bin/extract_fasta_bins.py --output_path  {params} {input.assembly} {input.clustering} > {output.log}"
            "extract_fasta_bins.py --output_path  {params} {input.assembly} {input.clustering} > {output.log}"
elif config["BINNING"] == "DAS" and config["das"]["concoct"]["run"]!="T":
   rule skip_concoct:
        output:
            log="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/concoct.log",
            log2="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/bin_clustering_gt"+config["concoct"]["min_contig_length"]+".csv"
        shell:
            "touch {output.log} && touch {output.log2}"

if config["BINNING"] == "BINSANITY" or (config["BINNING"] == "DAS" and config["das"]["binsanity"]["run"]=="T"):
    rule log_transform_coverage:
        '''
        We used to work with col 3 "totalAvgDepth" but now we use col4 CONTIGS_SPADESvs_{sample}_mapped_against_cross-assembly_sorted.bam 
        '''
        input:
            depth="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth.txt"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth_avg_log.txt"
        shell:
            "cat {input} | awk 'NR>1 && $4>1{{ if($4 <= 1) a = 0; else  a = log($4)/log(10); printf(\"%s\\t%0.4f\\n\",$1,a)}}' > {output}"
            if config["bwa"]["differential_coverage_matrix"].lower() == "f" else
            "cat {input} | awk -F '\\t' 'NR > 1 {{for(x=1;x<=NF;x++) if(x == 1 || (x >= 4 && x % 2 == 0)) {{if($x <= 1) a = 0; else if(x == 1) a = $x; else  a = log($x)/log(10);  printf \"%s\", a (x == NF || x == (NF-1) ? \"\\n\":\"\\t\")}}}}' > {output}"
            #"cat {input} | awk 'NR>1 {{print $1,log($3+1)/log(10)}}' > {output}"

    rule filter_fasta_by_coverage:
        input:
            depth="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth_avg_log.txt",
            fasta="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_scaffolds.fasta"
            if config["ANALYSIS"] == "SCAFFOLDS" else "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_contigs.fasta"
        output:
            temp("{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/assembly_coverage_gt0.fasta")
        shell:
            "cat {input.depth} | cut -f1 | grep -w -A1 --no-group-separator -F -f - {input.fasta}  > {output}" 

    rule bin_sanity:
        input:
            depth="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth_avg_log.txt",
            assembly="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_scaffolds.fasta"
            if config["ANALYSIS"] == "SCAFFOLDS" else "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_contigs.fasta"
        params:
            contig_directory="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"],
            bin_directory="{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/",
            #final bins will be at  {bin_directory}/BinSanity-Final-bins/final_Bin-xx.fna
            #also they can be found in the form final_Bin-xx_refined-xx.fna
            assembly="{sample}_scaffolds.fasta"
            if config["ANALYSIS"] == "SCAFFOLDS" else "{sample}_contigs.fasta"
        output:
            log="{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanityWf.log"
        threads:
            int(config["binsanity"]["threads"])
        shell:
            "Binsanity-wf -f {params.contig_directory} -l {params.assembly} -c {input.depth} -o {params.bin_directory} --binPrefix  final "
            #"/export/data/aabdala/utils/BInSanity/BinSanity-master/bin/Binsanity-wf -f {params.contig_directory} -l {params.assembly} -c {input.depth} -o {params.bin_directory} --binPrefix  final "
            "-p {config[binsanity][preference]} -x {config[binsanity][min_contig_length]} --threads {config[binsanity][threads]} "
            "{config[binsanity][extra_params]}"
elif config["BINNING"] == "DAS" and config["das"]["binsanity"]["run"]!="T":
   rule skip_bin_sanity:
        output:
            log="{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanityWf.log"
        shell:
            "touch {output}"

if config["BINNING"] == "DAS":
    rule bins_to_table:
        input:
            maxbin="{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/maxbin.log",
            metabat="{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/metabat.log",
            #concoct="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/concoct.log"
            concoct="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/bin_clustering_gt"+config["concoct"]["min_contig_length"]+".csv",
            binsanity="{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanityWf.log" 
        output:
            metabat_out="{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv",
            maxbin_out="{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv",
            concoct_out="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv",
            binsanity_out="{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv"
        params:
            output_dir_metabat="{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"],
            output_dir_maxbin="{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"],
            output_dir_concoct="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"],
            output_dir_binsanity="{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanity-Final-bins/",
            file_ext_metabat="fa",
            file_ext_maxbin="fasta",
            file_ext_concoct="fa",
            file_ext_binsanity="fna"
        script:
            "Scripts/tableBins.py"
    rule das:
        input:#https://stackoverflow.com/questions/51362210/snakemake-optional-input-for-rules
            metabat_bin2t="{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv",
            maxbin_bin2t="{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv",
            concoct_bin2t="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv",
            binsanity_bin2t="{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv",
            assembly="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_scaffolds.fasta"
            if config["ANALYSIS"] == "SCAFFOLDS" else "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_contigs.fasta"
        params:
            das_out_dir="{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/DasOut",
            bs_input=",{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv" if config["das"]["binsanity"]["run"]=="T" else "",
            mx_input=",{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv" if config["das"]["maxbin"]["run"]=="T" else "",
            cc_input=",{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv" if config["das"]["concoct"]["run"]=="T" else "",
            bs_l=",binsanity" if config["das"]["binsanity"]["run"]=="T" else "",
            mx_l=",maxbin"  if config["das"]["maxbin"]["run"]=="T" else "",
            cc_l=",concoct"  if config["das"]["concoct"]["run"]=="T" else ""

        output:
            log="{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/das.log",
            t2bin="{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/DasOut_DASTool_summary.tsv"
        threads:
            int(config["das"]["threads"])
        shell:
            #"DAS_Tool -i {input.metabat_bin2t},{input.maxbin_bin2t},{input.concoct_bin2t},{input.binsanity_bin2t} "
            #"-l metabat,maxbin,concoct,binsanity -c {input.assembly} -t {config[das][threads]} "
            #"--write_bins 1 --db_directory {config[das][db]} --search_engine {config[das][search_engine]} "
            #"--create_plots 1 {config[das][extra_params]} -o {params} > {output.log}"
            "DAS_Tool -i {input.metabat_bin2t}{params.mx_input}{params.cc_input}{params.bs_input} "
            "-l metabat{params.mx_l}{params.cc_l}{params.bs_l} -c {input.assembly} -t {config[das][threads]} "
            "--write_bins  --dbDirectory {config[das][db]} --search_engine {config[das][search_engine]} "
            " {config[das][extra_params]} -o {params.das_out_dir} > {output.log} 2>&1"
else:
    rule skip_das:
        params:
            "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]
        output:
            "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/das.log"
        shell:
            "touch {output}"


rule checkM_metabat2:
   input:
       "{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/metabat.log"
   params:
       bin_folder="{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/",
       out_folder="{PROJECT}/runs/{run}/{sample}_data/binning/checkM_metabat2/",
       bin_ext="fa"
   output:
       out_file="{PROJECT}/runs/{run}/{sample}_data/binning/checkM_metabat2/summary.txt"
   threads:
       int(config["checkM"]["threads"])
   shell:
       "checkm lineage_wf -f {output.out_file} -t  {config[checkM][threads]} -x {params.bin_ext} {config[checkM][extra_params]} {params.bin_folder} {params.out_folder} "

rule checkM_maxbin:
   input:
       "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/maxbin.log"
   params:
       bin_folder="{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/",
       out_folder="{PROJECT}/runs/{run}/{sample}_data/binning/checkM_maxbin/",
       bin_ext="fasta"
   output:
       out_file="{PROJECT}/runs/{run}/{sample}_data/binning/checkM_maxbin/summary.txt"
   threads:
       int(config["checkM"]["threads"])
   shell:
       "checkm lineage_wf -f {output.out_file} -t  {config[checkM][threads]} -x {params.bin_ext} {config[checkM][extra_params]} {params.bin_folder} {params.out_folder} "

rule checkM_concoct:
   input:
       "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/concoct.log"
   params:
       bin_folder="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/",
       out_folder="{PROJECT}/runs/{run}/{sample}_data/binning/checkM_concoct/",
       bin_ext="fa",
       res_folder="/export/lv3/scratch/workshop_2021/S12_Pipelines/{sample}/checkM_concoct/"
   output:
       out_file="{PROJECT}/runs/{run}/{sample}_data/binning/checkM_concoct/summary.txt"
   shell:
       "checkm lineage_wf -f {output.out_file} -t  {config[checkM][threads]} -x {params.bin_ext} {config[checkM][extra_params]} {params.bin_folder} {params.out_folder} "

rule checkM_binsanity:
   input:
       "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanityWf.log"
   params:
       bin_folder="{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanity-Final-bins/",
       out_folder="{PROJECT}/runs/{run}/{sample}_data/binning/checkM_binsanity/",
       bin_ext="fna"
   threads:
       int(config["checkM"]["threads"])         
   output:
       out_file="{PROJECT}/runs/{run}/{sample}_data/binning/checkM_binsanity/summary.txt"
   shell:
       "checkm lineage_wf -f {output.out_file} -t  {config[checkM][threads]} -x {params.bin_ext} {config[checkM][extra_params]} {params.bin_folder} {params.out_folder} "

rule checkM_das:
   input:
       "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/das.log",
       "{PROJECT}/runs/{run}/{sample}_data/binning/checkM_metabat2/summary.txt",
       "{PROJECT}/runs/{run}/{sample}_data/binning/checkM_maxbin/summary.txt" if config["das"]["maxbin"]["run"]=="T" and  config["das"]["maxbin"]["checkm_analysis"]=="T"  else "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/das.log",
       "{PROJECT}/runs/{run}/{sample}_data/binning/checkM_concoct/summary.txt" if config["das"]["concoct"]["run"]=="T" and config["das"]["concoct"]["checkm_analysis"]=="T"  else "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/das.log",
       "{PROJECT}/runs/{run}/{sample}_data/binning/checkM_binsanity/summary.txt" if config["das"]["binsanity"]["run"]=="T" and config["das"]["binsanity"]["checkm_analysis"]=="T" else "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/das.log"
   params:
       bin_folder="{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/DasOut_DASTool_bins/",
       out_folder="{PROJECT}/runs/{run}/{sample}_data/binning/checkM_das/",
       bin_ext="fa"
   output:
       out_file="{PROJECT}/runs/{run}/{sample}_data/binning/checkM_das/summary.txt"
   threads:
       int(config["checkM"]["threads"])
   shell:
      "checkm lineage_wf -f {output.out_file} -t  {config[checkM][threads]} -x {params.bin_ext} {config[checkM][extra_params]} {params.bin_folder} {params.out_folder} "

rule checkM_bins:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/binning/checkM_metabat2/summary.txt"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/checkM_maxbin/summary.txt"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/checkM_concoct/summary.txt"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/checkM_binsanity/summary.txt"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/checkM_das/summary.txt"
    params:
        dir="checkM_metabat2"
        if config["BINNING"] == "METABAT" else
        "checkM_maxbin"
        if config["BINNING"] == "MAXBIN" else
        "checkM_concoct"
        if config["BINNING"] == "CONCOCT" else
        "checkM_binsanity"
        if config["BINNING"] == "BINSANITY" else
        "checkM_das"
    output:
        out_file="{PROJECT}/runs/{run}/{sample}_data/binning/checkM/summary.txt"
    shell:
        "ln -s ../{params.dir}/summary.txt {output}"

rule gtdbtk_metabat2:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/metabat.log"
    params:
        bin_folder="{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/",
        out_folder="{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_metabat2/",
        bin_ext="fa"
    output:
        out_file="{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_metabat2/summary.txt"
    threads:
        int(config["gtdbtk"]["cpus"])
    conda:
        "envs/gtdbtk.yaml"
    shell:
        "gtdbtk classify_wf --genome_dir {params.bin_folder} --out_dir {params.out_folder}  -x {params.bin_ext} --cpus {config[gtdbtk][cpus]} {config[gtdbtk][extra_params]} > {output}"

rule gtdbtk_maxbin:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/maxbin.log"
    params:
        bin_folder="{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/",
        out_folder="{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_maxbin/",
        bin_ext="fasta"
    output:
        out_file="{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_maxbin/summary.txt"
    threads:
        int(config["gtdbtk"]["cpus"])
    conda:
        "envs/gtdbtk.yaml"
    shell:
        "gtdbtk classify_wf --genome_dir  {params.bin_folder} --out_dir {params.out_folder}  -x {params.bin_ext} --cpus {config[gtdbtk][cpus]} {config[gtdbtk][extra_params]}  > {output}"

rule gtdbtk_concoct:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/concoct.log"
    params:
        bin_folder="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/",
        out_folder="{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_concoct/",
        bin_ext="fa"
    output:
        out_file="{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_concoct/summary.txt"
    threads:
        int(config["gtdbtk"]["cpus"])
    conda:
        "envs/gtdbtk.yaml"
    shell:
        "gtdbtk classify_wf --genome_dir  {params.bin_folder} --out_dir {params.out_folder}  -x {params.bin_ext} --cpus {config[gtdbtk][cpus]} {config[gtdbtk][extra_params]}  > {output}"

rule gtdbtk_binsanity:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanityWf.log"
    params:
        bin_folder="{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanity-Final-bins/",
        out_folder="{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_binsanity/",
        bin_ext="fna"
    output:
        out_file="{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_binsanity/summary.txt"
    threads:
        int(config["gtdbtk"]["cpus"])
    conda:
        "envs/gtdbtk.yaml"
    shell:
        "gtdbtk classify_wf --genome_dir  {params.bin_folder} --out_dir {params.out_folder}  -x {params.bin_ext} --cpus {config[gtdbtk][cpus]} {config[gtdbtk][extra_params]}  > {output} "

rule gtdbtk_das:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/das.log",
        "{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_metabat2/summary.txt",
        "{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_maxbin/summary.txt" if config["das"]["maxbin"]["run"]=="T" and config["das"]["maxbin"]["gtdbtk_analysis"]=="T" else "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/das.log",
        "{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_concoct/summary.txt" if config["das"]["concoct"]["run"]=="T" and config["das"]["concoct"]["gtdbtk_analysis"]=="T" else "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/das.log",
        "{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_binsanity/summary.txt" if config["das"]["binsanity"]["run"]=="T" and config["das"]["binsanity"]["gtdbtk_analysis"]=="T"  else "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/das.log"
    params:
        bin_folder="{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/DasOut_DASTool_bins/",
        out_folder="{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_das/",
        bin_ext="fa"
    output:
        out_file="{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_das/summary.txt"
    threads:
        int(config["gtdbtk"]["cpus"])
    conda:
        "envs/gtdbtk.yaml"
    shell:
        "gtdbtk classify_wf --genome_dir  {params.bin_folder} --out_dir {params.out_folder}  -x {params.bin_ext} --cpus {config[gtdbtk][cpus]} {config[gtdbtk][extra_params]} > {output} "

rule gtdbtk_bins:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_metabat2/summary.txt"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_maxbin/summary.txt"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_concoct/summary.txt"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_binsanity/summary.txt"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_das/summary.txt"
    params:
        dir="gtdbtk_metabat2"
        if config["BINNING"] == "METABAT" else
        "gtdbtk_maxbin"
        if config["BINNING"] == "MAXBIN" else
        "gtdbtk_concoct"
        if config["BINNING"] == "CONCOCT" else
        "gtdbtk_binsanity"
        if config["BINNING"] == "BINSANITY" else
        "gtdbtk_das"
    output:
        out_file="{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk/summary.txt"
    shell:
        "ln -s ../{params.dir}/summary.txt {output}"



rule summarize_gtdbtk:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_metabat2/summary.txt"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_maxbin/summary.txt"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_concoct/summary.txt"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_binsanity/summary.txt"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_das/summary.txt"
    params:
        search_path="{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk_\\*/summary.txt"
    output:
        temp("{PROJECT}/runs/{run}/{sample}_data/binning/summary_gtdb.tsv")
    shell:
        "Scripts/summary_gtdb.sh {params.search_path}"

rule summarize_checkM:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/binning/checkM_metabat2/summary.txt"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/checkM_maxbin/summary.txt"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/checkM_concoct/summary.txt"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/checkM_binsanity/summary.txt"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/checkM_das/summary.txt"
    params:
        search_path="{PROJECT}/runs/{run}/{sample}_data/binning/checkM_\\*/summary.txt"
    output:
        temp("{PROJECT}/runs/{run}/{sample}_data/binning/summary_checkM.tsv")
    shell:
        "Scripts/summary_checkM.sh {params.search_path}"

rule merge_checkM_gtdb_results:
    input:
        checkM="{PROJECT}/runs/{run}/{sample}_data/binning/summary_checkM.tsv",
        gtdb="{PROJECT}/runs/{run}/{sample}_data/binning/summary_gtdb.tsv"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/binning/summary.tsv"
    shell:
        "cat {input.gtdb} | "
        "awk -F\"\\t\" 'BEGIN{{OFS=\"\\t\"}} NR==FNR{{if(NR==1){{ header=\"GTDB_\"$3\"\\tGTDB_\"$4\"\\tGTDB_\"$5\"\\tGTDB_\"$6}} else{{bin[$1$2]=$3\"\\t\"$4\"\\t\"$5\"\\t\"$6}};next}}  BEGIN{{OFS=\"\\t\"}} {{if(FNR==1){{print $0,header}} else{{if(bin[$1$2]){{print $0,bin[$1$2]}}else{{print $0,\"Filtered\",\"-\",\"-\",\"-\"}}}}}}' "
        "- {input.checkM} > {output}"

rule bin_cvg_metabat2:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/metabat.log"
    params:
        bin_folder="{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/",
        out_folder="{PROJECT}/runs/{run}/{sample}_data/binning/",
        bin_ext="fa",
        coverage="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth.txt"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/binning/abundance.metabat.tsv"
    shell:
        "Scripts/summary_coverage_metabat.sh {params.bin_folder} {params.bin_ext} {params.coverage} {output}" 

rule bin_cvg_concoct:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/concoct.log"
    params:
        bin_folder="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/",
        out_folder="{PROJECT}/runs/{run}/{sample}_data/binning/",
        bin_ext="fa",
        coverage="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth.txt"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/binning/abundance.concoct.tsv" 
    shell:
        "Scripts/summary_coverage_concoct.sh {params.bin_folder} {params.bin_ext} {params.coverage} {output}" 
       
       
rule bin_cvg_maxbin:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/maxbin.log"
    params:
        bin_folder="{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/",
        out_folder="{PROJECT}/runs/{run}/{sample}_data/binning/",
        bin_ext="fasta",
        coverage="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth.txt"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/binning/abundance.maxbin.tsv"
    shell:
        "Scripts/summary_coverage_maxbin.sh {params.bin_folder} {params.bin_ext} {params.coverage} {output}" 
 
       
rule bin_cvg_binsanity:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanityWf.log"
    params:
        bin_folder="{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanity-Final-bins/",
        out_folder="{PROJECT}/runs/{run}/{sample}_data/binning/",
        bin_ext="fna",
        coverage="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth.txt"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/binning/abundance.binsanity.tsv"
    shell:
        "Scripts/summary_coverage_bs.sh {params.bin_folder} {params.bin_ext} {params.coverage} {output}" 
 

rule bin_cvg_das:
     input:
        "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/das.log",
        "{PROJECT}/runs/{run}/{sample}_data/binning/abundance.metabat.tsv",
        "{PROJECT}/runs/{run}/{sample}_data/binning/abundance.maxbin.tsv",
        "{PROJECT}/runs/{run}/{sample}_data/binning/abundance.concoct.tsv",
        "{PROJECT}/runs/{run}/{sample}_data/binning/abundance.binsanity.tsv"
     params:
        bin_folder="{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/DasOut_DASTool_bins/",
        out_folder="{PROJECT}/runs/{run}/{sample}_data/binning/",
        bin_ext="fa",
        coverage="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth.txt"
     output:
        "{PROJECT}/runs/{run}/{sample}_data/binning/abundance.das.tsv"
     shell:
        "Scripts/summary_coverage_das.sh {params.bin_folder} {params.bin_ext} {params.coverage} {output}" 
 
rule summarize_coverage:
    input:
        abundance="{PROJECT}/runs/{run}/{sample}_data/binning/abundance.metabat.tsv"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/abundance.maxbin.tsv"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/abundance.concoct.tsv"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/abundance.binsanity.tsv"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/abundance.das.tsv",
        summary="{PROJECT}/runs/{run}/{sample}_data/binning/summary.tsv"
    params:
        abundance_folder="{PROJECT}/runs/{run}/{sample}_data/binning/"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/binning/summary_abundance.tsv"
    shell:
        "cat {params.abundance_folder}abundance*.tsv | grep  -v num_contigs | "
        "awk -F \"\\t\" 'FNR==NR{{if(NR>1){{h[$1$2]=$3\"\\t\"$4\"\\t\"$5}};next }} BEGIN{{OFS=\"\\t\"}} "
        "{{if(FNR==1){{print $0,\"num_contigs\",\"total_length\",\"avg_depth\" }}else{{print $0,h[$1$2]}} }}' "
        " - {input.summary} > {output}"
        
rule gc_prc_metabat2:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/metabat.log"
    params:
        bin_folder="{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/",
        out_folder="{PROJECT}/runs/{run}/{sample}_data/binning/",
        bin_ext="fa"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/binning/gc_prc.metabat.tsv"
    shell:
        "Scripts/computeGC.sh {params.bin_folder} {params.bin_ext} {output}"

rule gc_prc_concoct:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/concoct.log"
    params:
        bin_folder="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/",
        out_folder="{PROJECT}/runs/{run}/{sample}_data/binning/",
        bin_ext="fa",
        coverage="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth.txt"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/binning/gc_prc.concoct.tsv"
    shell:
        "Scripts/computeGC.sh {params.bin_folder} {params.bin_ext}  {output}"


rule gc_prc_maxbin:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/maxbin.log"
    params:
        bin_folder="{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/",
        out_folder="{PROJECT}/runs/{run}/{sample}_data/binning/",
        bin_ext="fasta",
        coverage="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth.txt"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/binning/gc_prc.maxbin.tsv"
    shell:
        "Scripts/computeGC.sh {params.bin_folder} {params.bin_ext} {output}"


rule gc_prc_binsanity:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanityWf.log"
    params:
        bin_folder="{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanity-Final-bins/",
        out_folder="{PROJECT}/runs/{run}/{sample}_data/binning/",
        bin_ext="fna",
        coverage="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth.txt"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/binning/gc_prc.binsanity.tsv"
    shell:
        "Scripts/computeGC.sh {params.bin_folder} {params.bin_ext} {output}"
        
rule gc_prc_das:
     input:
        "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/das.log",
        "{PROJECT}/runs/{run}/{sample}_data/binning/gc_prc.metabat.tsv",
        "{PROJECT}/runs/{run}/{sample}_data/binning/gc_prc.maxbin.tsv",
        "{PROJECT}/runs/{run}/{sample}_data/binning/gc_prc.concoct.tsv",
        "{PROJECT}/runs/{run}/{sample}_data/binning/gc_prc.binsanity.tsv"
     params:
        bin_folder="{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/DasOut_DASTool_bins/",
        out_folder="{PROJECT}/runs/{run}/{sample}_data/binning/",
        bin_ext="fa",
        coverage="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth.txt"
     output:
        "{PROJECT}/runs/{run}/{sample}_data/binning/gc_prc.das.tsv"
     shell:
        "Scripts/computeGC.sh {params.bin_folder} {params.bin_ext} {output}"

rule summarize_gc_prc:
    input:
        coverage="{PROJECT}/runs/{run}/{sample}_data/binning/gc_prc.metabat.tsv"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/gc_prc.maxbin.tsv"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/gc_prc.concoct.tsv"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/gc_prc.binsanity.tsv"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/gc_prc.das.tsv",
        summary="{PROJECT}/runs/{run}/{sample}_data/binning/summary_abundance.tsv"
    params:
        prc_folder="{PROJECT}/runs/{run}/{sample}_data/binning/"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/binning/summary_abundance_coverage.tsv"
    shell:
        "cat {params.prc_folder}gc_prc*.tsv |  "
        "awk -F \"\\t\" 'FNR==NR{{if(NR>1){{h[$1$2]=$3}};next }} BEGIN{{OFS=\"\\t\"}} "
        "{{if(FNR==1){{print $0,\"avg_gc\" }}else{{print $0,h[$1$2]}} }}' "
        " - {input.summary} > {output}"

if config["CREATE_UNBINNED"] == "T":
    rule get_unbinned_contigs:
        input:
            bin_table="{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv"
            if config["BINNING"] == "METABAT" else
            "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv"
            if config["BINNING"] == "MAXBIN" else
            "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv"
            if config["BINNING"] == "CONCOCT" else
            "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv"
            if config["BINNING"] == "BINSANITY" else
            "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/DasOut_DASTool_summary.tsv",
            assembly="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_scaffolds.fasta"
            if config["ANALYSIS"] == "SCAFFOLDS" else "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_contigs.fasta"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/unbinned/unbinned_contigs_list.txt"
        shell:
            "cat {input.bin_table} | cut -f1 | grep -v -F -w -f - {input.assembly} | grep \"^>\" | sed 's/^>//' > {output}"
    
    rule create_unbinned_fasta:
        input:
            unbinned_list="{PROJECT}/runs/{run}/{sample}_data/unbinned/unbinned_contigs_list.txt",
            assembly="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_scaffolds.fasta"
            if config["ANALYSIS"] == "SCAFFOLDS" else "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/{sample}_contigs.fasta"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/unbinned/unbinned.fasta"
        shell:
            "seqtk subseq {input.assembly} {input.unbinned_list} > {output}"

else:
    rule skip_create_unbinned:
        output:
            "{PROJECT}/runs/{run}/{sample}_data/unbinned/unbinned.txt"
        shell:
            "echo 'CREATE_UNBINNED == F' > {output}"

rule prokka_bins:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/metabat.log"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/maxbin.log"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/concoct.log"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanityWf.log"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/das.log",
        "{PROJECT}/runs/{run}/{sample}_data/binning/checkM/summary.txt"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
    params:
        output_dir="{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanity-Final-bins/"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/DasOut_DASTool_bins/",
        file_ext= "fa"
        if config["BINNING"] == "METABAT" or config["BINNING"] == "CONCOCT" or config["BINNING"] == "DAS" else
        "fasta"
        if config["BINNING"] == "BINSANITY" else
        "fna"
    benchmark:
        #"{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.benchmark"
        #if config["BINNING"] == "METABAT" else
        #"{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.benchmark"
        "{PROJECT}/runs/{run}/{sample}_data/binning/prokka"+config["BINNING"]+"_bins.benchmark"
    script:
            "Scripts/annotateProkkaBins.py"
rule diamond_bins:
    input:
        #"{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
        #if config["BINNING"] == "METABAT" else
        #"{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
        "{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/diamond.log"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/diamond.log"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/diamond.log"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/diamond.log"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/diamond.log"
    params:
        output_dir="{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanity-Final-bins/"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/DasOut_DASTool_bins/",
        file_ext= "faa"
    benchmark:
        "{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.benchmark"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.benchmark"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.benchmark"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.benchmark"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.benchmark"
    script:
        "Scripts/diamondProkkaBins.py"


rule rename_Final_bins:
    input:
        bin_table="{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/DasOut_DASTool_summary.tsv"
    params:
        bin_folder="{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/" if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/" if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/" if config["BINNING"] == "CONCOCT" else 
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanity-Final-bins/"  if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/DasOut_DASTool_bins/", 
        bin_ext="fa"  if config["BINNING"] == "METABAT" else
        "fasta"  if config["BINNING"] == "MAXBIN" else
        "fa" if config["BINNING"] == "CONCOCT" else
        "fna" if config["BINNING"] == "BINSANITY" else
        "fa",
        smp="{sample}",
        out_dir="{PROJECT}/runs/{run}/{sample}_data/binning/FinalBins/" 
    output:
        "{PROJECT}/runs/{run}/{sample}_data/binning/FinalBins/new_names.txt"
    shell:
        "Scripts/renameFinalBins.sh {params.bin_folder}  {params.bin_ext} {params.smp} {params.out_dir} {output}"  

rule coverage_contigs_final_bins:
    input:
        final_bins="{PROJECT}/runs/{run}/{sample}_data/binning/FinalBins/new_names.txt",
        bwa="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth.txt"
    params:
        bin_folder="{PROJECT}/runs/{run}/{sample}_data/binning/FinalBins/"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/binning/FinalBins/contig_coverage.txt"
    shell:
        "cat  {params.bin_folder}*.fna  | grep \"^>\" | sed 's/>// ; s/ /\t/g' | "
        " awk -F\"\\t\" 'NR==FNR{{h[$2]=$1;next}}BEGIN{{OFS=\"\\t\"; FS=\"\\t\"}}{{if(h[$1]){{print h[$1],$3}}}}' "
        "- {input.bwa} > {output}"

rule summarize_final_bins:
    input:
        summary="{PROJECT}/runs/{run}/{sample}_data/binning/summary_abundance_coverage.tsv",
        new_names="{PROJECT}/runs/{run}/{sample}_data/binning/FinalBins/new_names.txt"
    params:
        prc_folder="{PROJECT}/runs/{run}/{sample}_data/binning/"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/binning/FinalBins.summary.tsv"
    shell:
        "cat {input.new_names} |  "
        "awk -F \"\\t\" 'FNR==NR{{ parts=split($2,n,\".\");name=n[1]; for(i=2;i<parts;i++){{name=name\".\"n[i]}}; "
        "nparts=split($3,nn,\".\");nname=nn[1]; for(i=2;i<nparts;i++){{nname=nname\".\"nn[i]}}; h[$1name]=nname;next }} "
        "BEGIN{{OFS=\"\\t\"}} "
        "{{if(FNR==1){{print \"New_BinID\",$0}}else if(h[$1$2]){{print h[$1$2],$0}} }}' "
        " - {input.summary} > {output}"

rule create_yaml_bins_tbl:
    output:
        "{PROJECT}/runs/{run}/tables/bins.yaml"
    shell:
        "cp resources/datavzrd/bins.yaml {output}"

rule datavzrd_bins:
    input:
        config="{PROJECT}/runs/{run}/tables/bins.yaml",
        table=expand("{PROJECT}/runs/{run}/{sample}_data/binning/FinalBins.summary.tsv", PROJECT=config["PROJECT"],sample=config["SAMPLES"], run=run)
        #table="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/stats_assembly.tsv"
    output:
        report(
            directory("{PROJECT}/runs/{run}/tables/bins/{sample}"),
            htmlindex="index.html",
            #category="1. Project Information",
            category="6. Binning",
            subcategory="{sample}",
            labels={"sample":"{sample}"}, 
            #labels={"table", "Trimming results"},
        ),
    wrapper:
        "v4.7.2/utils/datavzrd"


rule report:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_mapped_against_cross-assembly_sorted.flagstat"
        if config["bwa"]["differential_coverage_matrix"].lower() == "f" else
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/sam-flagstat_cmds.log",
        ##"{PROJECT}/runs/{run}/{sample}_data/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log",
        #"{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/diamond.log"
        #if config["BINNING"] == "METABAT" else
        #"{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/diamond.log"
        #if config["BINNING"] == "MAXBIN" else
        #"{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/diamond.log"
        #if config["BINNING"] == "CONCOCT" else
        #"{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/diamond.log"
        #if config["BINNING"] == "BINSANITY" else
        #"{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/diamond.log",
        #"{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/maxbin.log",
        "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kraken.taxonomy.report"
         if config["TAXONOMY"]["PROFILING"] == "KRAKEN" else
         "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kaiju.taxonomy.out.report"
         if config["TAXONOMY"]["PROFILING"] == "KAIJU" else
         "{PROJECT}/runs/{run}/{sample}_data/taxonomy/all.taxonomy.out"
         if config["TAXONOMY"]["PROFILING"] == "ALL" else
         "{PROJECT}/runs/{run}/{sample}_data/no_tax.txt",
         "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/das.log",
        #"{PROJECT}/runs/{run}/{sample}_data/metabat2/prokka.out"
        #"{PROJECT}/runs/{run}/{sample}_data/metabat2/bin/metabat2.log",
         #"{PROJECT}/runs/{run}/{sample}_data/binning/checkM/summary.txt",
         #"{PROJECT}/runs/{run}/{sample}_data/binning/gtdbtk/summary.txt",
         "{PROJECT}/runs/{run}/{sample}_data/binning/summary.tsv",
         "{PROJECT}/runs/{run}/{sample}_data/binning/summary_abundance.tsv",
         "{PROJECT}/runs/{run}/{sample}_data/unbinned/unbinned.fasta"
         if config["CREATE_UNBINNED"] == "T" else
         "{PROJECT}/runs/{run}/{sample}_data/unbinned/unbinned.txt",
         "{PROJECT}/runs/{run}/{sample}_data/binning/FinalBins.summary.tsv",
         "{PROJECT}/runs/{run}/{sample}_data/binning/FinalBins/contig_coverage.txt",
         "{PROJECT}/runs/{run}/tables/trimmomatic",
         "{PROJECT}/runs/{run}/tables/bwa",
         "{PROJECT}/runs/{run}/tables/bins/{sample}"
    output:
        temp("{PROJECT}/runs/{run}/{sample}_data/report_f.html")
    shell:
        "touch {output}"
#rule tune_report:
#    input:
#        "{PROJECT}/runs/{run}/{sample}_data/report.html"
#    output:
#        "{PROJECT}/runs/{run}/{sample}_data/report_f.html"
#    script:
#        "Scripts/tuneReport.py"
