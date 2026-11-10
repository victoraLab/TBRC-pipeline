#run pandaseq
# -T 1: pandaseq's default threading writes the same reads in a different order each run, and fastx_collapser breaks ties by input order, so `head` kept different equally common reads.
rule pandaseq:
    input:
        expand("{source}/{sample}/{sample}.R1.fastq", source=DEST, sample=SAMPLE),
        expand("{source}/{sample}/{sample}.R2.fastq", source=DEST, sample=SAMPLE)
    output:
        expand("{destiny}/{sample}/{sample}.fasta", sample=SAMPLE, destiny=DEST)
    shell:
        "pandaseq -T 1 -f {input[0]:q} -r {input[1]:q} -w {output:q}"
