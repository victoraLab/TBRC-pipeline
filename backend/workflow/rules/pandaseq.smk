#run pandaseq
rule pandaseq:
    input:
        expand("{source}/{sample}/{sample}.R1.fastq", source=DEST, sample=SAMPLE),
        expand("{source}/{sample}/{sample}.R2.fastq", source=DEST, sample=SAMPLE)
    output:
        expand("{destiny}/{sample}/{sample}.fasta", sample=SAMPLE, destiny=DEST)
    shell:
        "pandaseq -f {input[0]} -r {input[1]} -w {output}"
