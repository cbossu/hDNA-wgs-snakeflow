# use SeqPrep2 to remove adaptors and merge paired-end reads with a minimum length of 28 bp and that overlap by > 19bp
# note leeHom is the recommended alternative to SeqPrep2 for this job. Consider using that program instead. There is no conda env, so install the binary from github in projects and call that.
rule seqprep2_trim_merge:
    input:
        unpack(get_fastq)
    output:
        merged="results/bqsr-round-{bqsr_round}/trim-merged/{sample}---{unit}.merged.fastq.gz"
    params:
        adapter_a="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", # same as in config file
        adapter_b="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        minlen="28",
        overlap="19"
    log:
        out="results/bqsr-round-{bqsr_round}/logs/seqprep2_trim_merge/{sample}---{unit}.log",
        err="results/bqsr-round-{bqsr_round}/logs/seqprep2_trim_merge/{sample}---{unit}.err"
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/seqprep2_trim_merge/{sample}---{unit}.bmk"
    shell: """
        /projects/cbossu@colostate.edu/SeqPrep2/SeqPrep2 \
            -f {input.r1} -r {input.r2} \
            -1 /dev/null -2 /dev/null -s {output.merged} \
            -A {params.adapter_a} -B {params.adapter_b} \
            -L {params.minlen} -o {params.overlap} > {log.out} 2> {log.err}
    """
        
        
# filter out low complexity reads with prinseq dust, threshold = 7
rule prinseqpp_filter_dust:
    input:
        fq="results/bqsr-round-{bqsr_round}/trim-merged/{sample}---{unit}.merged.fastq.gz"
    output:
        good="results/bqsr-round-{bqsr_round}/filtered_fq/{sample}---{unit}.filtered.fastq.gz"
    log:
        "results/bqsr-round-{bqsr_round}/logs/prinseq/{sample}---{unit}.log"
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/prinseq/{sample}---{unit}.bmk"
    threads: 2
    conda:
        "prinseqpp"
    shell:
        """
        prinseq++ \
            -fastq {input.fq} \
            -out_good {output.good} -out_gz \
            -out_bad /dev/null \
            -lc_dust=0.07 \
            &> {log}
        """

#align the filtered historical dna to the reference genome using the aln (-l 1024) and samse commands in BWA.
# I split this into two functions using wrappers so I could sort and add read groups and filter by -q 30
#rule align_reads:
#    input:
#        ref="resources/genome.fasta",
#        index=expand("resources/genome.fasta.{ext}", ext=["amb", "ann", "bwt", "pac", "sa"]),
#        reads="results/bqsr-round-{bqsr_round}/filtered_fq/{sample}---{unit}.filtered.fastq"
#    output:
#        sam="results/bqsr-round-{bqsr_round}/mapped/{sample}---{unit}.sam"
#    conda:
#        "bwa"
#    log:
#        "results/bqsr-round-{bqsr_round}/logs/map_aln_samse/{sample}---{unit}.log"
#    benchmark:
#        "results/bqsr-round-{bqsr_round}/benchmarks/map_aln_samse/{sample}---{unit}.bmk"
#    threads: 4
#    shell:
#        """
#        bwa aln -l 1024 -t {threads} {input.ref} {input.reads} > {output.sam}.sai 2>> {log}
#        bwa samse {input.ref} {output.sam}.sai {input.reads} > {output.sam} 2>> {log}
#        rm {output.sam}.sai
#        """

rule align_reads:
    input:
        idx=multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
        fastq="results/bqsr-round-{bqsr_round}/filtered_fq/{sample}---{unit}.filtered.fastq.gz"
    output:
        "results/bqsr-round-{bqsr_round}/sai/{sample}---{unit}.sai",
    conda:
        "bwa"
    params:
        extra="-l 1024",
    log:
        "results/bqsr-round-{bqsr_round}/logs/bwa_aln/{sample}---{unit}.log",
    threads: 8
    wrapper:
        "v1.23.3/bio/bwa/aln"        

rule bwa_samse:
    input:
        fastq="results/bqsr-round-{bqsr_round}/filtered_fq/{sample}---{unit}.filtered.fastq.gz",
        sai="results/bqsr-round-{bqsr_round}/sai/{sample}---{unit}.sai",
        idx=multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        temp("results/bqsr-round-{bqsr_round}/mapped_samse/{sample}---{unit}.sorted.bam"),
    conda:
        "bwa"
    params:
        extra=get_read_group, # r"-r '@RG\tID:{sample}\tSM:{sample}'",  # optional: Extra parameters for bwa.
        sort="samtools",  # optional: Enable sorting. Possible values: 'none', 'samtools' or 'picard'`
        sort_order="coordinate",  # optional: Sort by 'queryname' or 'coordinate'
        sort_extra="-q 30",  # optional: extra arguments for samtools/picard
    threads: 4
    log:
        "results/bqsr-round-{bqsr_round}/logs/bwa_samse/{sample}--{unit}.log",
    wrapper:
        "v1.23.3/bio/bwa/samse"    

# use samtools to convert sam to bam files with minimum mapping quality of 30 (for each unit) and sort.
rule filt_sorted_bam:
    input:
        "results/bqsr-round-{bqsr_round}/mapped_samse/{sample}---{unit}.sorted.bam"
    output:
        "results/bqsr-round-{bqsr_round}/mapped_samse_filt/{sample}---{unit}.sorted.filt2.bam"
    log:
        "results/bqsr-round-{bqsr_round}/logs/samtools/filtered_bams/{sample}---{unit}.log"
    conda:
        "bioinf"
    threads: 10
    shell:
        """
        samtools view -b -q 30 -F 4 {input}  -@ {threads} -o {output} 2> {log}
        """
#-bq 30 bam files had some unmapped sites that weren't 0, to remove unmapped flag -F 4

        
# holden removed this in historical processing workflow, but to remove duplicates, you just 
# need to add --REMOVE_DUPLICATES true to config file, so I kept this rule because I want to first know the percent dups per sample, using
# --until flag, and then remove the duplicates by changing the config file

rule mark_duplicates:
    input:
        get_all_bams_of_common_sample
    output:
        bam="results/bqsr-round-{bqsr_round}/mkdup/{sample}.bam",
        bai="results/bqsr-round-{bqsr_round}/mkdup/{sample}.bai",
        metrics="results/bqsr-round-{bqsr_round}/qc/mkdup/{sample}.metrics.txt",
    log:
        "results/bqsr-round-{bqsr_round}/logs/picard/mkdup/{sample}.log",
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/mark_duplicates/{sample}.bmk"
    conda:
        "picard"
    params:
        extra=config["params"]["picard"]["MarkDuplicates"],
    resources:
        cpus = 4
    wrapper:
        "v1.1.0/bio/picard/markduplicates"

# use samtools to merge all bams from common sample together, remove duplicates, and index the resulting bam files.
#rule rmdup:
#    input:
#        get_all_bams_of_common_sample
#    output:
#        bam="results/bqsr-round-{bqsr_round}/rmdup/{sample}_rmdup.bam",
#        bai="results/bqsr-round-{bqsr_round}/rmdup/{sample}_rmdup.bai"
#    log:
#        "results/bqsr-round-{bqsr_round}/logs/samtools/rmdup/{sample}.log"
#    conda:
#        "bioinf"
#    shell:
#        """
#        samtools merge -f {wildcards.sample}_merged.bam {input} && \
#        samtools rmdup {wildcards.sample}_merged.bam {output.bam} && \
#        samtools index -o {output.bai} {output.bam} && \
#        rm {wildcards.sample}_merged.bam
#        """

# use mapDamage2 to rescale base quality scores in accordance with their probability of being damaged
rule rescale_base_quality_scores:
    input:
        ref="resources/genome.fasta",
        bam="results/bqsr-round-{bqsr_round}/mkdup/{sample}.bam"
    output:
        db="results/bqsr-round-{bqsr_round}/mapdamage2/{sample}",
        log="results/bqsr-round-{bqsr_round}/mapdamage2/{sample}/Runtime_log.txt",
        GtoA3p="results/bqsr-round-{bqsr_round}/mapdamage2/{sample}/3pGtoA_freq.txt",
        CtoT5p="results/bqsr-round-{bqsr_round}/mapdamage2/{sample}/5pCtoT_freq.txt",
        dnacomp="results/bqsr-round-{bqsr_round}/mapdamage2/{sample}/dnacomp.txt",
        frag_misincorp="results/bqsr-round-{bqsr_round}/mapdamage2/{sample}/Fragmisincorporation_plot.pdf",
        len="results/bqsr-round-{bqsr_round}/mapdamage2/{sample}/Length_plot.pdf",
        lg_dist="results/bqsr-round-{bqsr_round}/mapdamage2/{sample}/lgdistribution.txt",
        misincorp="results/bqsr-round-{bqsr_round}/mapdamage2/{sample}/misincorporation.txt",
        rescaled_bam="results/bqsr-round-{bqsr_round}/mapdamage2/rescaled_bams/{sample}.rescaled.bam"
    params:
        extra="--rescale"
    log:
        "results/bqsr-round-{bqsr_round}/logs/mapdamage2/{sample}.log"
    threads: 1
    shell: """
        mapDamage --input {input.bam} --reference {input.ref} \
        --folder {output.db} --rescale-out {output.rescaled_bam} --rescale
        """
#    wrapper:
#        "v6.2.0/bio/mapdamage2"

### need to define absolute path to tmpdir in hpcc_profile config file! Otherwise this snakemake rule failed regardless of shell or wrapper usage (Bayesian estimation never completed)
### NOTE: mkdir rescaled_bam inside mapdamage2 folder first- this might not be necessary once the tmpdir was fixed, but it was done anyway.
