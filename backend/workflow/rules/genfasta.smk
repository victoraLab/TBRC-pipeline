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
        workflow/scripts/split_fasta.sh -f {input.fasta:q} -1 {input.barcode1:q} -2 {input.barcode2:q} -o {output:q} -d {params.destiny:q} -m {params.meth:q} -h {params.head:q} -t ''{params.trimming_seq:q} > {log.std:q} 2> {log.error:q} 
        """

