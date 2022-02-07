
"""
Workflow to map sequencing data from 4 runs to the genome and mark duplicates
inputs: 
    sra
outputs:
    mods
run:
module load snakemake/5.26.1-foss-2019b-Python-3.7.4 
snakemake -np --use-envmodules --max-status-checks-per-second 0.01 --snakefile code/snakemake_wg_lrtaps.smk --cluster "qsub -o lrtaps.log -e lrtaps.err  -P ludwig.prjc -q long.qc@@long.hge -cwd -V -S /bin/bash -N taps -pe shmem 4" -j 4
samtools view -bS bc1001_cmb.fwd.trim.minimap2.bam 4kb:246-4260 >bc1001_cmb.fwd.trim.minimap2.4kb.bam
samtools view -bS -s 0.01 bc1001_cmb.fwd.trim.minimap2.4kb.bam >bc1001_cmb.fwd.trim.minimap2.4kb.0.01.bam
samtools view -bS bc1001_cmb.rev.trim.minimap2.bam 4kb:246-4260 >bc1001_cmb.rev.trim.minimap2.4kb.bam
samtools view -bS -s 0.01 bc1001_cmb.rev.trim.minimap2.4kb.bam >bc1001_cmb.rev.trim.minimap2.4kb.0.01.bam
samtools index bc1001_cmb.fwd.trim.minimap2.4kb.0.01.bam;samtools index bc1001_cmb.rev.trim.minimap2.4kb.0.01.bam
"""

from snakemake.utils import min_version
min_version("5.26")


SAMPLES = ["bc1001_cmb"] #["bc1001_reseq"] #["bc1001", "bc1002","bc1003"]
REF = "/gpfs2/well/ludwig/users/cfo155/whole_genome_lrTAPS/resource/mm9_genome.spikeins.fasta"
REF_SIZE = "/gpfs2/well/ludwig/users/cfo155/whole_genome_lrTAPS/resource/mm9_genome.spikeins.fasta.fai"
# SAMPLES = ["bc1001_cmb_realign"]
# REF = "/gpfs2/well/ludwig/users/cfo155/whole_genome_lrTAPS/resource/C57BL_6J_Eve_chromosome.spike_in.fasta"
# REF_SIZE = "/gpfs2/well/ludwig/users/cfo155/whole_genome_lrTAPS/resource/C57BL_6J_Eve_chromosome.spike_in.fasta.fai"
# SAMPLES = ["bc1001_cmb_e14"]
# REF = "/gpfs2/well/ludwig/users/cfo155/whole_genome_lrTAPS/resource/e14_genome.spikeins.fasta"
# REF_SIZE = "/gpfs2/well/ludwig/users/cfo155/whole_genome_lrTAPS/resource/e14_genome.spikeins.fasta.fai"
# SAMPLES = ["bc1001_cmb_mm10"]
# REF = "/gpfs2/well/ludwig/users/cfo155/whole_genome_lrTAPS/resource/mm10_genome.spikeins.fasta"
# REF_SIZE = "/gpfs2/well/ludwig/users/cfo155/whole_genome_lrTAPS/resource/mm10_genome.spikeins.fasta.fai"


print(SAMPLES)

rule all:
    input:
        # expand("trim/{sample}.fwd.trim.fastq.gz", sample = SAMPLES),
        # expand("trim/{sample}.fwd.trim_fastqc.html", sample = SAMPLES),
        # expand("align/{sample}.fwd.trim.minimap2.bam", sample = SAMPLES),
        # expand("align/{sample}.tg.bam", sample = SAMPLES),      
        # expand("align/{sample}.tg.uniflag.bam", sample = SAMPLES),
        # expand("align/{sample}.uniflag.bam", sample = SAMPLES),
        # expand("align/{sample}.uniflag.md.bam", sample = SAMPLES),
        # expand("meth/{sample}.uniflag.md_CpG.bedGraph", sample = SAMPLES),
        # expand("sta/{sample}.all.sta", sample = SAMPLES),
        # expand("align/{sample}.uniflag.gatk_haplotypecaller.vcf.gz", sample = SAMPLES)
        # expand("fastq/{sample}.read_len.txt", sample = SAMPLES),
        # expand("meth/{sample}_gene.gz",  sample = SAMPLES),
        # expand("meth/{sample}_cgi.gz",  sample = SAMPLES),
        expand("plots/{sample}_gene.avg.pdf",sample = SAMPLES )
        # expand("align/{sample}.cov.bedGraph", sample = SAMPLES),
        # expand("align/{sample}.tg.len.txt", sample = SAMPLES),
        # expand("align/{sample}.uniflag.gc_bias_summ.txt", sample = SAMPLES),
        # expand("align/{sample}.chrHMM.cov.sta", sample = SAMPLES),
        # expand("meth/{sample}.chrHMM.meth.sta", sample = SAMPLES),
        # expand("align/{sample}.chrHMM.len.sta", sample = SAMPLES),
        # expand("align/{sample}.rm.sta", sample = SAMPLES),
        # expand("align/{sample}.uniflag.sort.bam", sample = SAMPLES),
        # expand("align/{sample}.ca.cgi.cov.sta", sample = SAMPLES),
        # expand("meth/{sample}.uniflag_CpG.cgi.sta", sample = SAMPLES),
        # expand("meth/{sample}.uniflag_CpG.merge.bedGraph", sample = SAMPLES),
        # expand("meth/{sample}.uniflag.rmdup_CpG.bedGraph", sample = SAMPLES)
        # expand("meth/{sample}.tg.uniflag_CpG.bedGraph", sample = SAMPLES),
        # expand("fastq/{sample}.first20.txt", sample = SAMPLES),        
        # expand("align/{sample}.uniflag.cov_metrics.txt", sample = SAMPLES),        
        # expand("meth/{sample}.uniflag_CpG.context.sta", sample = SAMPLES)
        # expand("align/{sample}.tg.spikein.bam", sample = SAMPLES),


rule trim:
    input:
        "fastq/{sample}.fastq.gz"
    output:
        fwd_read="trim/{sample}.fwd.trim.fastq.gz", 
        rev_read="trim/{sample}.rev.trim.fastq.gz",
        untrim=temp("trim/{sample}.untrim.fastq.gz")
    shell:
        """
        /gpfs2/well/ludwig/users/cfo155/miniconda2/bin/cutadapt -g CCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT {input} --untrimmed-output {output.untrim} |\
        /gpfs2/well/ludwig/users/cfo155/miniconda2/bin/cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTA -o {output.fwd_read} -
        /gpfs2/well/ludwig/users/cfo155/miniconda2/bin/cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGG {output.untrim} --discard-untrimmed |\
        /gpfs2/well/ludwig/users/cfo155/miniconda2/bin/cutadapt -g TACGAGATACATCGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -o {output.rev_read} -
        """


rule fastqc:
    input:
        read="fastq/{sample}.fastq.gz",
        fwd_read="trim/{sample}.fwd.trim.fastq.gz", 
        rev_read="trim/{sample}.rev.trim.fastq.gz"
    output:
        fwd_qc="trim/{sample}.fwd.trim_fastqc.html", 
        rev_qc="trim/{sample}.rev.trim_fastqc.html"
    shell:
        """
        /gpfs2/well/ludwig/users/cfo155/miniconda2/opt/fastqc-0.11.8/fastqc -i {input.read}
        /gpfs2/well/ludwig/users/cfo155/miniconda2/opt/fastqc-0.11.8/fastqc -i {input.fwd_read}
        /gpfs2/well/ludwig/users/cfo155/miniconda2/opt/fastqc-0.11.8/fastqc -i {input.rev_read}
        """

rule align:
    input:
        fwd_read="trim/{sample}.fwd.trim.fastq.gz", 
        rev_read="trim/{sample}.rev.trim.fastq.gz"
    output:
        fwd_bam="align/{sample}.fwd.trim.minimap2.bam", 
        rev_bam="align/{sample}.rev.trim.minimap2.bam"
    params:
        ref=REF,
        prefix="{sample}"
    threads: 3
    shell:
        """
        /gpfs2/well/ludwig/users/cfo155/miniconda2/bin/minimap2 -a -x map-pb {params.ref} {input.fwd_read} -t {threads} --MD | \
        /gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools sort -T {params.prefix} -o {output.fwd_bam}
        /gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools index {output.fwd_bam}
        /gpfs2/well/ludwig/users/cfo155/miniconda2/bin/minimap2 -a -x map-pb {params.ref} {input.rev_read} -t {threads} --MD | \
        /gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools sort -T {params.prefix} -o {output.rev_bam}
        /gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools index {output.rev_bam}
        """

rule bam_merge:
    input:
        fwd_bam="align/{sample}.fwd.trim.minimap2.bam", 
        rev_bam="align/{sample}.rev.trim.minimap2.bam"
    output:
        fwd_bam_flag16="align/{sample}.fwd.flag16.bam",
        rev_bam_flag16="align/{sample}.rev.flag16.bam",
        fwd_bam_flag0="align/{sample}.fwd.flag0.bam", 
        rev_bam_flag0="align/{sample}.rev.flag0.bam",
        fwd_bam="align/{sample}.tg.bam", 
        rev_bam="align/{sample}.ca.bam"
    shell:
        """
        /gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools view -h {input.fwd_bam}|awk '$0~/^@/||$2=="16"'|samtools view -bS - >{output.fwd_bam_flag16}
        /gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools view -h {input.rev_bam}|awk '$0~/^@/||$2=="16"'|samtools view -bS - >{output.rev_bam_flag16}
        /gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools view -h {input.fwd_bam}|awk '$0~/^@/||$2=="0"' |samtools view -bS ->{output.fwd_bam_flag0}
        /gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools view -h {input.rev_bam}|awk '$0~/^@/||$2=="0"' |samtools view -bS ->{output.rev_bam_flag0}
        /gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools merge {output.rev_bam} {output.fwd_bam_flag16} {output.rev_bam_flag0}
        /gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools merge {output.fwd_bam} {output.rev_bam_flag16} {output.fwd_bam_flag0}
        """

rule bam_change_flag:
    input:
        fwd_bam="align/{sample}.tg.bam", 
        rev_bam="align/{sample}.ca.bam"
    output:
        fwd_bam="align/{sample}.tg.uniflag.bam", 
        rev_bam="align/{sample}.ca.uniflag.bam"
    shell:
        """
        /gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools view {input.fwd_bam} |awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{gsub("[0-9]*","0",$2);print $0}}'|cat <(/gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools view -H {input.fwd_bam}) -|/gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools view -bS - >{output.fwd_bam}
        /gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools view {input.rev_bam} |awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{gsub("[0-9]*","16",$2);print $0}}'|cat <(/gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools view -H {input.rev_bam}) -|/gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools view -bS - >{output.rev_bam}
        """

rule merge_flag_changed_bam:
    input:
        fwd_bam="align/{sample}.tg.uniflag.bam", 
        rev_bam="align/{sample}.ca.uniflag.bam"
    output:
        bam="align/{sample}.uniflag.bam"
    shell:
        """
        /gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools merge {output.bam} {input.fwd_bam} {input.rev_bam}
        """

rule remove_duplicates:
    input:
        "align/{sample}.uniflag.bam"
    output:
        nrm="align/{sample}.uniflag.sort.bam",
        rm="align/{sample}.uniflag.rmdup.bam",
        md="align/{sample}.uniflag.md.bam",
        md_matrix="align/{sample}.uniflag.md_metrics.txt"
    envmodules: "picard/2.23.0-Java-11"
    shell:
        """
        /gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools sort {input} -o {output.nrm}
        /gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools rmdup -s {input} {output.rm}
        java -jar /apps/eb/skylake/software/picard/2.23.0-Java-11/picard.jar MarkDuplicates \
            I={input} \
            O={output.md} \
            M={output.md_matrix}
        """


rule call_merge:
    input:
        nrm_bam="align/{sample}.uniflag.sort.bam",
        rm_bam="align/{sample}.uniflag.rmdup.bam",
        md_bam="align/{sample}.uniflag.md.bam"
    output: 
        nrm_mcall="meth/{sample}.uniflag.sort_CpG.bedGraph",
        rm_mcall="meth/{sample}.uniflag.rmdup_CpG.bedGraph",
        md_mcall="meth/{sample}.uniflag.md_CpG.bedGraph",
        nrm_mcall_merge="meth/{sample}.uniflag.sort_CpG.merge.bedGraph",
        rm_mcall_merge="meth/{sample}.uniflag.rmdup_CpG.merge.bedGraph",
        md_mcall_merge="meth/{sample}.uniflag.md_CpG.merge.bedGraph"
    envmodules: 
        "matplotlib/3.1.1-foss-2019b-Python-3.7.4"
    params:
        pars="-q 10 -p 5 -t 2 --CHG --CHH",
        nrm_prefix="meth/{sample}.uniflag.sort",
        rm_prefix="meth/{sample}.uniflag.rmdup",
        md_prefix="meth/{sample}.uniflag.md",
        ref=REF,
        ref_size=REF_SIZE,
        genome_cg="resource/mm9_genome.cg.info",
        spikes_cg="resource/spike_in.cg.info"
    shell:
        """
        /users/ludwig/cfo155/cfo155/cfDNA/022020_cfDNA/test/mbias_spikein/MethylDackel/bin/MethylDackel extract \
            {params.ref} {input.nrm_bam} \
            {params.pars} -o {params.nrm_prefix} 
        /users/ludwig/cfo155/cfo155/cfDNA/022020_cfDNA/test/mbias_spikein/MethylDackel/bin/MethylDackel extract \
            {params.ref} {input.rm_bam} \
            {params.pars} -o {params.rm_prefix} 
        /users/ludwig/cfo155/cfo155/cfDNA/022020_cfDNA/test/mbias_spikein/MethylDackel/bin/MethylDackel extract \
            {params.ref} {input.md_bam} \
            {params.pars} -o {params.md_prefix} 
        bedtools intersect -nonamecheck -a <(cat {params.genome_cg} {params.spikes_cg}|grep -v seqID|awk 'BEGIN{{OFS="\\t"}}{{print $1,$5,$5+2}}'|sort -k1,1 -k2,2n ) -b <(tail -n +2 {output.nrm_mcall}|sort -k1,1 -k2,2n) -wa -wb -sorted |\
            awk '{{mC[$1"_"$2]+=$9;aC[$1"_"$2]+=$8+$9}}END{{for(i in mC)print i"\\t"mC[i]"\\t"aC[i]}}' |\
            sort -k1,1  |awk 'BEGIN{{OFS="\\t"}}{{if($3==0)print $1,$2,$3,"0";else print $1,$2,$3,$2/$3}}' |\
            join -1 1 <(cat {params.genome_cg} {params.spikes_cg}|grep -v seqID|awk 'BEGIN{{OFS="\\t"}}{{print $1"_"$5}}'|sort -k1,1 )  -2 1 - -a1 -t$'\\t'|\
            awk 'BEGIN{{OFS="\\t"}}{{if($2=="")print $1,"NA","NA","NA";else print $1,$2,$3,$4}}' |\
            sed 's/_/\\t/g'|sort -k1,1 -k2,2n |sed 1i"chr\\tstart\\tmC\\taC\\tmeth" >{output.nrm_mcall_merge}
        bedtools intersect -nonamecheck -a <(cat {params.genome_cg} {params.spikes_cg}|grep -v seqID|awk 'BEGIN{{OFS="\\t"}}{{print $1,$5,$5+2}}'|sort -k1,1 -k2,2n ) -b <(tail -n +2 {output.rm_mcall}|sort -k1,1 -k2,2n) -wa -wb -sorted |\
            awk '{{mC[$1"_"$2]+=$9;aC[$1"_"$2]+=$8+$9}}END{{for(i in mC)print i"\\t"mC[i]"\\t"aC[i]}}' |\
            sort -k1,1  |awk 'BEGIN{{OFS="\\t"}}{{if($3==0)print $1,$2,$3,"0";else print $1,$2,$3,$2/$3}}' |\
            join -1 1 <(cat {params.genome_cg} {params.spikes_cg}|grep -v seqID|awk 'BEGIN{{OFS="\\t"}}{{print $1"_"$5}}'|sort -k1,1 )  -2 1 - -a1 -t$'\\t'|\
            awk 'BEGIN{{OFS="\\t"}}{{if($2=="")print $1,"NA","NA","NA";else print $1,$2,$3,$4}}' |\
            sed 's/_/\\t/g'|sort -k1,1 -k2,2n |sed 1i"chr\\tstart\\tmC\\taC\\tmeth" >{output.rm_mcall_merge}
        bedtools intersect -nonamecheck -a <(cat {params.genome_cg} {params.spikes_cg}|grep -v seqID|awk 'BEGIN{{OFS="\\t"}}{{print $1,$5,$5+2}}'|sort -k1,1 -k2,2n ) -b <(tail -n +2 {output.md_mcall}|sort -k1,1 -k2,2n) -wa -wb -sorted |\
            awk '{{mC[$1"_"$2]+=$9;aC[$1"_"$2]+=$8+$9}}END{{for(i in mC)print i"\\t"mC[i]"\\t"aC[i]}}' |\
            sort -k1,1  |awk 'BEGIN{{OFS="\\t"}}{{if($3==0)print $1,$2,$3,"0";else print $1,$2,$3,$2/$3}}' |\
            join -1 1 <(cat {params.genome_cg} {params.spikes_cg}|grep -v seqID|awk 'BEGIN{{OFS="\\t"}}{{print $1"_"$5}}'|sort -k1,1 )  -2 1 - -a1 -t$'\\t'|\
            awk 'BEGIN{{OFS="\\t"}}{{if($2=="")print $1,"NA","NA","NA";else print $1,$2,$3,$4}}' |\
            sed 's/_/\\t/g'|sort -k1,1 -k2,2n |sed 1i"chr\\tstart\\tmC\\taC\\tmeth" >{output.md_mcall_merge}
        """

rule sta_new:
    input:
        raw="fastq/{sample}.fastq.gz",
        nrm_bam="align/{sample}.uniflag.sort.bam",
        rm_bam="align/{sample}.uniflag.md.bam",
        fwd_read="trim/{sample}.fwd.trim.fastq.gz", 
        rev_read="trim/{sample}.rev.trim.fastq.gz",
        nrm_mcall_merge="meth/{sample}.uniflag.sort_CpG.merge.bedGraph",
        rm_mcall_merge="meth/{sample}.uniflag.md_CpG.merge.bedGraph",
        nrm_mcall_merge_ch=expand("meth/{{sample}}.uniflag.sort_{context}.bedGraph", context=['CHG','CHH']),
    output:
        "sta/{sample}.all.sta"
    params:
        pos="resource/spike_in.cg.info",
        genome_cg="resource/mm9_genome.cg.info"
    shell:
        """
        nraw=`zcat {input.raw}|awk 'NR%4==1'|wc -l`
        nfwd_read=`zcat {input.fwd_read}|awk 'NR%4==1'|wc -l`
        nrev_read=`zcat {input.rev_read}|awk 'NR%4==1'|wc -l`
        nq20map_bam=`/gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools view -q20 {input.nrm_bam}|cut -f1|sort -u|wc -l`
        nq20rm_bam=`/gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools view -q20 {input.rm_bam}|cut -f1|sort -u|wc -l`
        nq1map_bam=`/gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools view -q1 {input.nrm_bam}|cut -f1|sort -u|wc -l`
        nq1rm_bam=`/gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools view -q1 {input.rm_bam}|cut -f1|sort -u|wc -l`
        
        kb4_meth=`awk '$7=="CCGG" && $1=="4kb:246-4260"' {params.pos} |awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{print $1"_"$5"\\n"$1"_"$5+1}}' |sort -k1,1 |join -1 1 - -2 1 <(awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{print $1"_"$2,$0}}' {input.nrm_mcall_merge}|sort -k1,1) -t$'\\t'|awk '{{mC+=$4;aC+=$5}}END{{print mC/aC";("mC"/"aC")"}}'`
        lambda_meth=`awk '$7=="CCGG" && $1=="J02459.1:44132-47728"' {params.pos} |awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{print $1"_"$5"\\n"$1"_"$5+1}}' |sort -k1,1 |join -1 1 - -2 1 <(awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{print $1"_"$2,$0}}' {input.nrm_mcall_merge}|sort -k1,1) -t$'\\t'|awk '{{mC+=$4;aC+=$5}}END{{print mC/aC";("mC"/"aC")"}}'`
        kb4_unmeth_cg=`awk '$7!="CCGG" && $1=="4kb:246-4260"' {params.pos} |awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{print $1"_"$5"\\n"$1"_"$5+1}}' |sort -k1,1 |join -1 1 - -2 1 <(awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{print $1"_"$2,$0}}' {input.nrm_mcall_merge}|sort -k1,1) -t$'\\t'|awk '{{mC+=$4;aC+=$5}}END{{print mC/aC";("mC"/"aC")"}}'`
        lambda_unmeth_cg=`awk '$7!="CCGG" && $1=="J02459.1:44132-47728"' {params.pos} |awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{print $1"_"$5"\\n"$1"_"$5+1}}' |sort -k1,1 |join -1 1 - -2 1 <(awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{print $1"_"$2,$0}}' {input.nrm_mcall_merge}|sort -k1,1) -t$'\\t'|awk '{{mC+=$4;aC+=$5}}END{{print mC/aC";("mC"/"aC")"}}'`
        kb4_unmeth_ch=`awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{if($1=="4kb:246-4260")print $0}}' {input.nrm_mcall_merge_ch} |awk '{{mC+=$6;aC+=$5+$6}}END{{print mC/aC";("mC"/"aC")"}}'`
        lambda_unmeth_ch=`awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{if($1=="J02459.1:44132-47728")print $0}}' {input.nrm_mcall_merge_ch} |awk '{{mC+=$6;aC+=$5+$6}}END{{print mC/aC";("mC"/"aC")"}}'`
        
        genome_meth=`awk '$1!="4kb:246-4260"&&$1!="J02459.1:44132-47728"' {params.genome_cg} |awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{print $1"_"$5"\\n"$1"_"$5+1}}'|sort -k1,1 |join -1 1 - -2 1 <(awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{print $1"_"$2,$0}}' {input.rm_mcall_merge}|sort -k1,1) -t$'\\t'|awk '{{mC+=$4;aC+=$5}}END{{print mC/aC";("mC"/"aC")"}}'`

        echo -e "id\\tnraw\\tnfwd_read\\tnrev_read\\tnq20map_bam\\tnq20rm_bam\\tnq1map_bam\\tnq1rm_bam\\tkb4_meth\\tlambda_meth\\tkb4_unmeth_cg\\tlambda_unmeth_cg\\tkb4_unmeth_ch\\tlambda_unmeth_ch\\tgenome_meth" >{output}
        echo -e "{output}\\t${{nraw}}\\t${{nfwd_read}}\\t${{nrev_read}}\\t${{nq20map_bam}}\\t${{nq20rm_bam}}\\t\\t${{nq1map_bam}}\\t${{nq1rm_bam}}\\t${{kb4_meth}}\\t${{lambda_meth}}\\t${{kb4_unmeth_cg}}\\t${{lambda_unmeth_cg}}\\t${{kb4_unmeth_ch}}\\t${{lambda_unmeth_ch}}\\t${{genome_meth}}" >>{output}
        """

rule gatk_HaplotypeCaller:
    input:
        "align/{sample}.uniflag.rmdup.bam"
    output:
        bam="align/{sample}.uniflag.rmdup.rf.bam",
        vcf="align/{sample}.uniflag.gatk_haplotypecaller.vcf.gz"
    params:
        ref=REF,
        prams="--read-filter NotSupplementaryAlignmentReadFilter --pcr-indel-model AGGRESSIVE --native-pair-hmm-threads 1"
    envmodules: "GATK/4.1.7.0-GCCcore-8.3.0-Java-11"
    shell:
        """
        gatk CreateSequenceDictionary --REFERENCE {params.ref}
        java -jar /apps/eb/skylake/software/picard/2.23.0-Java-11/picard.jar AddOrReplaceReadGroups \
            I={input}  O={output.bam} RGID=4 RGLB=lib1 RGPL=PACBIO RGPU=bc1001 RGSM=bc1001_cmb
        samtools index {output.bam}
        gatk HaplotypeCaller --input {output.bam} --output {output.vcf} --reference {params.ref} {params.prams}
        """


rule gc_bias:
    input:
        fwd_bam="align/{sample}.tg.bam", 
        rev_bam="align/{sample}.ca.bam",
        bam="align/{sample}.uniflag.bam"
    output:
        fwd_mat="align/{sample}.tg.gc_bias_metrics.txt", 
        rev_mat="align/{sample}.ca.gc_bias_metrics.txt",
        mat="align/{sample}.uniflag.gc_bias_metrics.txt",
        fwd_pdf="align/{sample}.tg.gc_bias_metrics.pdf", 
        rev_pdf="align/{sample}.ca.gc_bias_metrics.pdf",
        pdf="align/{sample}.uniflag.gc_bias_metrics.pdf",
        fwd_summ="align/{sample}.tg.gc_bias_summ.txt", 
        rev_summ="align/{sample}.ca.gc_bias_summ.txt",
        summ="align/{sample}.uniflag.gc_bias_summ.txt",
    params:
        ref=REF
    envmodules: "picard/2.23.0-Java-11"
    shell:
        """
        java -jar /apps/eb/skylake/software/picard/2.23.0-Java-11/picard.jar CollectGcBiasMetrics \
                I={input.fwd_bam} \
                O={output.fwd_mat}  \
                CHART={output.fwd_pdf} \
                S={output.fwd_summ}  \
                R={params.ref}
        java -jar /apps/eb/skylake/software/picard/2.23.0-Java-11/picard.jar CollectGcBiasMetrics \
                I={input.rev_bam} \
                O={output.rev_mat}  \
                CHART={output.rev_pdf} \
                S={output.rev_summ}  \
                R={params.ref}
        java -jar /apps/eb/skylake/software/picard/2.23.0-Java-11/picard.jar CollectGcBiasMetrics \
                I={input.bam} \
                O={output.mat}  \
                CHART={output.pdf} \
                S={output.summ}  \
                R={params.ref}
        """

rule readlen_sta:
    input:
        "fastq/{sample}.fastq.gz"
    output:
        bin_sta="fastq/{sample}.read_len.txt",
        all_sta="fastq/{sample}.read_all.txt"
    shell:
        """
        zcat {input} |awk 'NR%4==2'|awk '{{print length($0)}}' |tee {output.all_sta}|awk '{{print (int($0/100)*100)}}'|sort|uniq -c |sed 's/^ \+//g;s/ /\\t/g'|awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{print $2,$1}}'|sort -k1,1 |sed 1i"len\t{output.bin_sta}" >{output.bin_sta}
        """


rule fastq_sta:
    input:
        "fastq/{sample}.fastq.gz"
    output:
        "fastq/{sample}.first20.txt"
    shell:
        """
        zcat {input} |awk 'NR%4==2'|awk '{{len=length($0); print substr($0,1,20)";"substr($0,len-20,len)}}'|sort|uniq -c |sed 's/^ \+//g;s/ /\\t/g'|awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{print $2,$1}}'|sort -k2,2n|sed 1i"first20\t{output}" >{output}
        """

rule alignlen_sta:
    input:
        fwd_bam="align/{sample}.tg.bam", 
        rev_bam="align/{sample}.ca.bam",
        bam="align/{sample}.uniflag.bam"
    output:
        fwd_len="align/{sample}.tg.len.txt", 
        rev_len="align/{sample}.ca.len.txt",
        fwd_sta="align/{sample}.tg.sta.txt", 
        rev_sta="align/{sample}.ca.sta.txt",
        lens="align/{sample}.uniflag.len.txt",
        sta="align/{sample}.uniflag.sta.txt"
    shell:
        """
        samtools view {input.fwd_bam}|cut -f6|sed 's/[0-9]*[S,H,I]//g;s/[M,D]/+/g;s/+$//g'|bc|tee {output.fwd_len}| awk '{{print (int($1/10))*10}}'|sort|uniq -c |sed 's/^ \+//g;s/ /\\t/g'|awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{print $2,$1}}'|sort -k1,1 |sed 1i"len\t{output.fwd_len}" >{output.fwd_sta}
        samtools view {input.rev_bam}|cut -f6|sed 's/[0-9]*[S,H,I]//g;s/[M,D]/+/g;s/+$//g'|bc|tee {output.rev_len}| awk '{{print (int($1/10))*10}}'|sort|uniq -c |sed 's/^ \+//g;s/ /\\t/g'|awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{print $2,$1}}'|sort -k1,1 |sed 1i"len\t{output.rev_len}" >{output.rev_sta}
        samtools view {input.bam}|cut -f6|sed 's/[0-9]*[S,H,I]//g;s/[M,D]/+/g;s/+$//g'|bc|tee {output.lens}| awk '{{print (int($1/10))*10}}'|sort|uniq -c |sed 's/^ \+//g;s/ /\\t/g'|awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{print $2,$1}}'|sort -k1,1 |sed 1i"len\t{output.lens}" >{output.sta}
        """

rule alignlen_feature:
    input:
        "align/{sample}.uniflag.bam"
    output:
        lens="align/{sample}.chrHMM.len.txt",
        stas="align/{sample}.chrHMM.len.sta"
    params:
        feature="resource/mESC_cStates_HMM"
    shell:
        """
        paste <(samtools view {input}|cut -f3,4) <(samtools view {input} |cut -f6|sed 's/[0-9]*[S,H,I]//g;s/[M,D]/+/g;s/+$//g'|bc)|\
            grep ^chr|awk 'BEGIN{{OFS="\\t"}}{{print $1,int(($2+$3)/2),int(($2+$3)/2)+1,$3}}'|\
            intersectBed -a - -b {params.feature} -wa -wb| tee {output.lens} |\
            awk '{{len[$8]+=$4;cnt[$8]+=1}}END{{for(i in len)print i"\\t"cnt[i]"\\t"len[i]/cnt[i]}}'|sort -k1,1|\
            sed 1i"feature\tread\tmean_len" >{output.stas}
        """


rule coverage:
    input:
        fwd_bam="align/{sample}.tg.bam", 
        rev_bam="align/{sample}.ca.bam",
        bam="align/{sample}.uniflag.bam"
    output:
        fwd_cov="align/{sample}.tg.cov.bedGraph", 
        rev_cov="align/{sample}.ca.cov.bedGraph",
        cov="align/{sample}.cov.bedGraph"
    params:
        ref_size=REF_SIZE
    shell:
        """
        bedtools coverage -a <( bedtools makewindows -w 1000  -g {params.ref_size} |sort -k1,1 -k2,2n) -b {input.fwd_bam} > {output.fwd_cov}
        bedtools coverage -a <( bedtools makewindows -w 1000  -g {params.ref_size} |sort -k1,1 -k2,2n) -b {input.rev_bam} > {output.rev_cov}
        bedtools coverage -a <( bedtools makewindows -w 1000  -g {params.ref_size} |sort -k1,1 -k2,2n) -b {input.bam} > {output.cov}
        """

rule chrHMM_cov_sta:
    input:
        bam="align/{sample}.uniflag.bam"
    output:
        cov="align/{sample}.chrHMM.cov.bedGraph",
        sta="align/{sample}.chrHMM.cov.sta", 
    params:
        bed="resource/mESC_cStates_HMM"
    shell:
        """
        bedtools coverage -a <( grep ^chr {params.bed} |sort -k1,1 -k2,2n) -b {input.bam}|
            tee {output.cov}|awk '{{depth[$4]+=$10;covp[$4]+=$13;cov[$4]+=$11}}END{{for(i in depth)print i"\t"depth[i]"\t"covp[i]/NR"\t"cov[i]}}'|sort -k1,1 >{output.sta}
        """

rule chrHMM_meth_sta:
    input:
        "meth/{sample}.uniflag_CpG.bedGraph"
    output:
        meth="meth/{sample}.chrHMM.meth.bedGraph",
        sta="meth/{sample}.chrHMM.meth.sta"
    params:
        bed="resource/mESC_cStates_HMM"
    shell:
        """
        bedtools intersect -a <( grep ^chr {params.bed} |sort -k1,1 -k2,2n) -b {input} -wa -wb |\
            tee {output.meth}|awk '{{mC[$4]+=$14;aC[$4]+=$15}}END{{for(i in mC)print i"\\t"mC[i]/aC[i]"\\t"mC[i]"\\t"aC[i]}}' >{output.sta}
        """

rule cgi_cov_sta:
    input:
        fwd_bam="align/{sample}.tg.bam", 
        rev_bam="align/{sample}.ca.bam",
        bam="align/{sample}.uniflag.bam"
    output:
        fwd_cov="align/{sample}.tg.cgi.cov.bedGraph", 
        rev_cov="align/{sample}.ca.cgi.cov.bedGraph",
        fwd_sta="align/{sample}.tg.cgi.cov.sta", 
        rev_sta="align/{sample}.ca.cgi.cov.sta",
        cov="align/{sample}.cgi.cov.bedGraph",
        sta="align/{sample}.cgi.cov.sta", 
    params:
        bed="resource/cpgIslandExt.txt"
    shell:
        """
        bedtools coverage -a <( cut -f2-4 {params.bed} |sort -k1,1 -k2,2n) -b {input.fwd_bam} |tee {output.fwd_cov}|awk '{{depth+=$4;covp+=$7;cov+=$5}}END{{print depth"\\t"covp/NR"\\t"cov}}'|sort -k1,1 >{output.fwd_sta}
        bedtools coverage -a <( cut -f2-4 {params.bed} |sort -k1,1 -k2,2n) -b {input.rev_bam} |tee {output.rev_cov}|awk '{{depth+=$4;covp+=$7;cov+=$5}}END{{print depth"\\t"covp/NR"\\t"cov}}'|sort -k1,1 >{output.rev_sta}
        bedtools coverage -a <( cut -f2-4 {params.bed} |sort -k1,1 -k2,2n) -b {input.bam} |tee {output.cov}|awk '{{depth+=$4;covp+=$7;cov+=$5}}END{{print depth"\\t"covp/NR"\\t"cov}}'|sort -k1,1 >{output.sta}
        """

rule cgi_meth_sta:
    input:
        "meth/{sample}.uniflag_CpG.bedGraph"
    output:
        meth="meth/{sample}.uniflag_CpG.cgi.bedGraph",
        sta="meth/{sample}.uniflag_CpG.cgi.sta"
    params:
        bed="resource/cpgIslandExt.txt"
    shell:
        """
        bedtools intersect -a <( cut -f2-4 {params.bed} |sort -k1,1 -k2,2n) -b {input} -wa -wb |tee {output.meth}|awk '{{mC+=$8;aC+=$9}}END{{print mC/aC"\\t"mC"\\t"aC}}' >{output.sta}
        """


rule align_metrics:
    input:
        bam="align/{sample}.uniflag.bam"
    output:
        mat="align/{sample}.uniflag.cov_metrics.txt",
        pdf="align/{sample}.uniflag.cov_metrics.pdf"
    params:
        ref=REF
    envmodules: "picard/2.23.0-Java-11"
    shell:
        """
        java -jar /apps/eb/skylake/software/picard/2.23.0-Java-11/picard.jar CollectWgsMetricsWithNonZeroCoverage  \
                I={input.bam} \
                O={output.mat} \
                CHART={output.pdf}  \
                R={params.ref}
        """




rule rmdup_sta:
    input:
        nrm="align/{sample}.uniflag.sort.bam",
        rm="align/{sample}.uniflag.rmdup.bam"
    output:
        "align/{sample}.rm.sta"
    shell:
        """
        nrm=`/gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools view {input.nrm} |awk '$3~/chr/' |wc -l`
        rm=`/gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools view {input.rm}|awk '$3~/chr/' |wc -l`
        prm=`echo "scale=3;${{rm}}/${{nrm}}"|bc`
        echo -e "${{nrm}}\\t${{rm}}\\t${{prm}}" >{output}
        """

rule preseq:
    input:
        "align/{sample}.uniflag.sort.bam"
    output:
        c_curve="align/{sample}.c_curve.txt",
        lc_extrap="align/{sample}.lc_extrap.txt"
    shell:
        """
        /users/ludwig/cfo155/miniconda2/envs/preseq/bin/preseq c_curve -B {input} -v > {output.c_curve} 2>&1 
        /users/ludwig/cfo155/miniconda2/envs/preseq/bin/preseq lc_extrap -B -o {output.lc_extrap} {input}
        """

rule context_meth:
    input:
        "meth/{sample}.uniflag_CpG.bedGraph"
    output:
        "meth/{sample}.uniflag_CpG.context.sta"
    params:
        cg="resource/mm9_genome.cg.info"
    shell:
        """
        tail -n +2 {params.cg}|grep -v N|\
        awk 'BEGIN{{OFS="\\t"}}{{print $1,$5+1,toupper($7)"\\n"$1,$5+$2,toupper($7)}}'|\
        awk 'BEGIN{{OFS="\\t"}}{{print $1"_"$2,$3}}'|sort -k1,1 |\
        join -1 1 - -2 1 <(awk 'BEGIN{{OFS="\\t"}}{{print $1"_"$2,$5,$6}}' {input}|sort -k1,1 ) -t$'\\t'|\
        awk '{{mC[$2]+=$3;aC[$2]+=$4}}END{{for(i in mC)print i,mC[i],aC[i],mC[i]/aC[i]}}' >{output}
        """

rule bam_spikeins:
    input:
        fwd_bam="align/{sample}.tg.bam", 
        rev_bam="align/{sample}.ca.bam"
    output:
        fwd_bam="align/{sample}.tg.spikein.bam", 
        rev_bam="align/{sample}.ca.spikein.bam"
    params:
        ref_size=REF_SIZE
    shell:
        """
        /gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools view -bS {input.fwd_bam} -L <(awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{if($1!~/chr/)print $1,"1",$2}}' {params.ref_size})  >{output.fwd_bam}
        /gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools index {output.fwd_bam}
        /gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools view -bS {input.rev_bam} -L <(awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{if($1!~/chr/)print $1,"1",$2}}' {params.ref_size})  >{output.rev_bam}
        /gpfs2/well/ludwig/users/cfo155/tools/samtools-1.11/bin/samtools index {output.rev_bam}
        """

rule call:
    input:
        fwd_bam="align/{sample}.tg.uniflag.bam", 
        rev_bam="align/{sample}.ca.uniflag.bam"
    output: 
        fwd_mcall="meth/{sample}.tg.uniflag_CpG.bedGraph",
        rev_mcall="meth/{sample}.ca.uniflag_CpG.bedGraph"
    envmodules: 
        "matplotlib/3.1.1-foss-2019b-Python-3.7.4"
    params:
        pars="-q 10 -p 5 -t 2",
        fwd_prefix="meth/{sample}.tg.uniflag",
        rev_prefix="meth/{sample}.ca.uniflag",
        ref=REF,
        ref_size=REF_SIZE
    shell:
        """
        /users/ludwig/cfo155/cfo155/cfDNA/022020_cfDNA/test/mbias_spikein/MethylDackel/bin/MethylDackel extract \
            {params.ref} {input.fwd_bam} \
            {params.pars} -o {params.fwd_prefix} 
        /users/ludwig/cfo155/cfo155/cfDNA/022020_cfDNA/test/mbias_spikein/MethylDackel/bin/MethylDackel extract \
            {params.ref} {input.rev_bam} \
            {params.pars} -o {params.rev_prefix}
        """

rule deeptools_compute_matrix_gene: 
    input:
        bw="meth/{sample}.meth.bw"
    output:
        outfilename="meth/{sample}_gene.gz",
        outfilenamemat="meth/{sample}_gene.tab",
        outfilereg="meth/{sample}_gene.bed"
    log:
        "logs/{sample}.bw.log"
    params:
        bw="meth/{sample}.meth.bw",
        region="resource/mm9.refGene.gtf",
        other=" -b 5000 -a 5000 --regionBodyLength 5000 --binSize 100 ",
        addbw="meth/short_read.meth.bw"
    envmodules:
        "deepTools/3.3.1-foss-2018b-Python-3.6.6","plotly.py/4.4.1-foss-2018b-Python-3.6.6"
    shell:
        """
        computeMatrix scale-regions \
            -R {params.region} --transcriptID transcript\
            -S {params.bw} {params.addbw} {params.other} \
            -o {output.outfilename} \
            --outFileNameMatrix {output.outfilenamemat} \
            --outFileSortedRegions {output.outfilereg}
        """


rule deeptools_compute_matrix_cgi: 
    input:
        bw="meth/{sample}.meth.bw"
    output:
        outfilename="meth/{sample}_cgi.gz",
        outfilenamemat="meth/{sample}_cgi.tab",
        outfilereg="meth/{sample}_cgi.bed"
    log:
        "logs/{sample}.bw.log"
    params:
        bw="meth/{sample}.meth.bw",
        region="resource/cpgIslandExt.bed",
        other=" -b 2000 -a 2000 --regionBodyLength 500 --binSize 100 ",
        addbw="meth/short_read.meth.bw"
    envmodules:
        "deepTools/3.3.1-foss-2018b-Python-3.6.6","plotly.py/4.4.1-foss-2018b-Python-3.6.6"
    shell:
        """
        computeMatrix scale-regions \
            -R {params.region} \
            -S {params.bw} {params.addbw} {params.other} \
            -o {output.outfilename} \
            --outFileNameMatrix {output.outfilenamemat} \
            --outFileSortedRegions {output.outfilereg}
        """

rule deeptools_plotheatmap: 
    input:
        expand("meth/{{sample}}_{region}.gz", region=['gene','cgi'])
    output:
        ht=expand("plots/{{sample}}_{region}.heatmap.png", region=['gene','cgi']),
        av=expand("plots/{{sample}}_{region}.avg.pdf", region=['gene','cgi'])
    log:
        "logs/{{sample}}.bw.log"
    params:
        other="--colorMap jet --missingDataColor \"#999999\" --heatmapHeight 15 --whatToShow 'heatmap and colorbar' --yMin 0 --yMax 1 "
    envmodules:
        "deepTools/3.3.1-foss-2018b-Python-3.6.6","plotly.py/4.4.1-foss-2018b-Python-3.6.6"
    shell:
        """
        plotHeatmap -m {input[0]} -out {output.ht[0]} {params.other}
        plotProfile -m {input[0]} -out {output.av[0]} --yMin 0 --yMax 1 
        plotHeatmap -m {input[1]} -out {output.ht[1]} {params.other}
        plotProfile -m {input[1]} -out {output.av[1]} --yMin 0 --yMax 1 
        """
        
