# The Illumina front half, shared by every Illumina assembly backend.
#
# prepare -> merge -> demultiplex -> cluster each well. Everything here is
# independent of what a pipeline then does with a well's clusters: tiago_vsearch
# ranks them and takes the top N per well, tiago_sulkseq keeps every cluster and
# counts clones by the wells they appear in. Splitting the two apart means the
# second pipeline reuses this rather than copying it, and a fix to demultiplexing
# reaches both.
#
# What differs from classic_luka, and why:
#
#   prep     trims the Illumina read-through adapter off the RAW mates, so a
#            zero-insert amplicon becomes an ~85 nt pair instead of a pair that
#            runs into adapter and mis-merges into a ~504 nt chimera. This is the
#            source of the 24-32% adapter-bearing records in delivered output.
#   merge    vsearch with a real overlap and length window, emitting FASTQ.
#   demux    mismatch tolerance derived from barcode separation on BOTH passes,
#            rather than --exact on BC1 (which discarded 35.6% of merged reads to
#            single-base errors) and an unsafe default on BC2.
#   wells    clusters at 99% identity so a read's error cloud counts towards its
#            own clone instead of splintering into singletons. Centroids, not
#            synthetic consensus: every delivered sequence is still a real read.
#   rank     deterministic total ordering, and a per-well verdict rather than a
#            silent top-N.

_DEST = DEST[0]
_SAMPLE = SAMPLE
_RUNDIR = "{}/{}".format(_DEST, _SAMPLE)


rule prep_reads_vsearch:
    input:
        r1=READ1_GZ,
        r2=READ2_GZ,
    output:
        r1=temp("{}/{}.prep.R1.fastq.gz".format(_RUNDIR, _SAMPLE)),
        r2=temp("{}/{}.prep.R2.fastq.gz".format(_RUNDIR, _SAMPLE)),
        report="{}/qc/prep/cutadapt.json".format(_RUNDIR),
    threads: 4
    params:
        spacer=param("prep", "spacer_nt", 2),
        a3_r1=param("prep", "adapter3_r1", "AGATCGGAAGAGC"),
        a3_r2=param("prep", "adapter3_r2", "AGATCGGAAGAGC"),
        overlap=param("prep", "adapter_overlap", 8),
        err=param("prep", "adapter_error_rate", 0.15),
        minlen=param("prep", "min_read_len", 60),
    log:
        "logs/{}.prep.log".format(_SAMPLE),
    shell:
        """
        mkdir -p "$(dirname {output.report})"
        cutadapt -j {threads} \
            -u {params.spacer} -U {params.spacer} \
            -a {params.a3_r1} -A {params.a3_r2} \
            -O {params.overlap} -e {params.err} \
            -m {params.minlen} --pair-filter=any \
            --json {output.report} \
            -o {output.r1} -p {output.r2} \
            {input.r1} {input.r2} > {log} 2>&1
        """


rule spacer_probe:
    """Verify the library really has the spacer we are about to cut."""
    input:
        reads=READ1_GZ,
        barcodes=expand("data/barcodes/{barset}/BC1.txt", barset=BARC),
    output:
        report="{}/qc/prep/spacer_probe.tsv".format(_RUNDIR),
    params:
        expect=param("prep", "spacer_nt", 2),
        maxoff=param("prep", "probe_max_offset", 6),
        limit=param("prep", "probe_reads", 200000),
    log:
        "logs/{}.spacer_probe.log".format(_SAMPLE),
    shell:
        """
        mkdir -p "$(dirname {output.report})"
        python3 workflow/scripts/probe_spacer.py \
            --reads {input.reads} --barcodes {input.barcodes} \
            --expect {params.expect} --max-offset {params.maxoff} \
            --limit {params.limit} --report {output.report} > {log} 2>&1
        """


rule merge_vsearch:
    input:
        r1=rules.prep_reads_vsearch.output.r1,
        r2=rules.prep_reads_vsearch.output.r2,
        probe=rules.spacer_probe.output.report,
    output:
        merged=temp("{}/{}.merged.fastq".format(_RUNDIR, _SAMPLE)),
        log="{}/qc/merge/vsearch_merge.log".format(_RUNDIR),
    threads: 4
    params:
        minov=param("merge", "min_overlap", 20),
        maxdiffs=param("merge", "max_diffs", 30),
        maxdiffpct=param("merge", "max_diff_pct", 25),
        minlen=param("merge", "min_len", 300),
    shell:
        """
        mkdir -p "$(dirname {output.log})"
        vsearch --fastq_mergepairs {input.r1} --reverse {input.r2} \
            --fastqout {output.merged} \
            --fastq_allowmergestagger \
            --fastq_minovlen {params.minov} \
            --fastq_maxdiffs {params.maxdiffs} \
            --fastq_maxdiffpct {params.maxdiffpct} \
            --fastq_minmergelen {params.minlen} \
            --threads {threads} \
            --log {output.log} 2>&1
        """


rule demux_vsearch:
    input:
        merged=rules.merge_vsearch.output.merged,
        barcode1=expand("data/barcodes/{barset}/BC1.txt", barset=BARC),
        barcode2=expand("data/barcodes/{barset}/BC2.txt", barset=BARC),
    output:
        sentinel="{}/demux.done".format(_RUNDIR),
        stats="{}/BC.stats.txt".format(_RUNDIR),
    params:
        destiny="{}/".format(_RUNDIR),
        minreads=param("demux", "min_reads_per_well", 10),
        keepunmatched=str(param("demux", "keep_unmatched", False)).lower(),
    log:
        "logs/{}.demux.log".format(_SAMPLE),
    shell:
        """
        bash workflow/scripts/demux_wells.sh \
            -f {input.merged} -1 {input.barcode1} -2 {input.barcode2} \
            -d {params.destiny} -o {output.sentinel} \
            -n {params.minreads} -u {params.keepunmatched} > {log} 2>&1
        """


rule wells_vsearch:
    """Dereplicate and cluster each well. Centroids keep every delivered
    sequence a real observed read; what changes is that its size= now aggregates
    the error cloud instead of counting exact duplicates only."""
    input:
        sentinel=rules.demux_vsearch.output.sentinel,
    output:
        sentinel="{}/wells.done".format(_RUNDIR),
    threads: 8
    params:
        split2="{}/split1/split2/".format(_RUNDIR),
        clustered="{}/split1/split2/clustered/".format(_RUNDIR),
        ident=param("wells", "cluster_id", 0.99),
        iddef=param("wells", "iddef", 2),
        minlen=param("merge", "min_len", 300),
        topn=param("wells", "max_clusters_per_well", TOPN),
        preprocess=str(METH[0]) if METH else "standard",
        light_5p=str(TRIM[0]) if TRIM else "",
        light_3p=param("wells", "light_chain_3p", "CACTGACACTCCCACGAGCTC"),
        light_footprint=param("wells", "light_chain_primer_nt", 23),
    log:
        "logs/{}.wells.log".format(_SAMPLE),
    shell:
        """
        mkdir -p {params.clustered}
        # Light-chain rescue, applied per well AFTER demultiplexing -- the same
        # place classic_luka applies it, and the only place it can go: the 5'
        # handle sits immediately after BC1, so stripping it any earlier would
        # take the barcode with it.
        light_rescue() {{
            f="$1"; five="$2"; three="$3"; foot="$4"
            tmp="$(mktemp)"; tmp2="$(mktemp)"
            cutadapt -g "$five" -e 0.2 -O 15 --discard-untrimmed -o "$tmp" "$f" >/dev/null 2>&1 || return 1
            [ -s "$tmp" ] || {{ rm -f "$tmp" "$tmp2"; return 1; }}
            fastx_trimmer -Q33 -f "$((foot + 1))" -i "$tmp" -o "$tmp2" 2>/dev/null || {{ rm -f "$tmp" "$tmp2"; return 1; }}
            cutadapt -a "$three" -e 0.2 -O 10 --minimum-length 1 -o "$f" "$tmp2" >/dev/null 2>&1
            rc=$?; rm -f "$tmp" "$tmp2"; return $rc
        }}
        export -f light_rescue

        cluster_one() {{
            f="$1"; out="$2"; ident="$3"; iddef="$4"; minlen="$5"; topn="$6"
            prep="${{7:-standard}}"; five="${{8:-}}"; three="${{9:-}}"; foot="${{10:-23}}"
            w="$(basename "$f" .fastq)"
            if [ "$prep" = "luka_light" ]; then
                light_rescue "$f" "$five" "$three" "$foot" || return 0
                [ -s "$f" ] || return 0
            fi
            # No length floor here. The merge stage already applied
            # --fastq_minmergelen, so everything reaching this point cleared it
            # BEFORE trimming. Re-applying it after light-chain rescue -- which
            # removes ~75 nt of handle, primer and tail -- discarded most of the
            # library and left clusters of one or two reads.
            vsearch --fastx_uniques "$f" --sizeout --fastaout "${{out}}${{w}}.uniq.fa" \
                    --minseqlength 32 --threads 1 --quiet || return 0
            [ -s "${{out}}${{w}}.uniq.fa" ] || return 0
            vsearch --cluster_size "${{out}}${{w}}.uniq.fa" --id "$ident" --iddef "$iddef" \
                    --strand plus --qmask none --sizein --sizeout \
                    --centroids "${{out}}${{w}}.cent.fa" --threads 1 --quiet || return 0
            vsearch --sortbysize "${{out}}${{w}}.cent.fa" --topn "$topn" \
                    --sizein --sizeout --output "${{out}}${{w}}.top.fa" --threads 1 --quiet || return 0
            rm -f "${{out}}${{w}}.uniq.fa" "${{out}}${{w}}.cent.fa"
        }}
        export -f cluster_one
        ls {params.split2}*.fastq 2>/dev/null \
          | xargs -P {threads} -I{{}} bash -c \
            'cluster_one "$@"' _ {{}} {params.clustered} {params.ident} {params.iddef} {params.minlen} {params.topn} \
            '{params.preprocess}' '{params.light_5p}' '{params.light_3p}' {params.light_footprint} \
          > {log} 2>&1
        touch {output.sentinel}
        """


