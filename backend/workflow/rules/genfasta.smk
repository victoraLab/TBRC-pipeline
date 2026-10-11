rule genfasta:
    input:
        fasta=expand("{source}/{sample}/{sample}.fasta", source=DEST, sample=SAMPLE),
        barcode1=expand("data/barcodes/{barset}/BC1.txt", barset=BARC, sample=SAMPLE),
        barcode2=expand("data/barcodes/{barset}/BC2.txt", barset=BARC, sample=SAMPLE)
    params:
        destiny=expand("{destiny}/{sample}/", destiny=DEST, sample=SAMPLE),
        meth=METH,
        head=HEAD[0],
        trimming_seq=TRIM
    output:
        expand("{destiny}/{sample}/split1/{sample}", destiny=DEST, sample=SAMPLE)
    log:
        std=expand("logs/{sample}.genfasta.out.log", sample=SAMPLE),
        error=expand("logs/{sample}.genfasta.err.log", sample=SAMPLE)
    shell:
        """
        workflow/scripts/split_fasta.sh -f {input.fasta} -1 {input.barcode1} -2 {input.barcode2} -o {output} -d '{params.destiny}' -m '{params.meth}' -h '{params.head}' -t '{params.trimming_seq}' > {log.std} 2> {log.error} 
        """

