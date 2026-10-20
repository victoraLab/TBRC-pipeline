# tiago_nanopore -- single-cell Ig/TCR amplicons from Oxford Nanopore.
#
# The other two pipelines assemble a short Illumina pair per well. Nanopore
# reads span the whole amplicon on their own, so there is nothing to merge; the
# work is separating wells and turning each well's ~270 noisy reads into one
# accurate sequence.
#
# Shape:
#   demux_nanopore    plate cassette + well/FR1 anchor -> per-read well identity
#   consensus_nanopore  per well: cluster to the dominant clone, align with
#                     abPOA, polish with racon
#   emit_nanopore     write the same output contract every other pipeline emits,
#                     so final.smk, the archive and the run report are unchanged
#
# Two decisions worth stating, both measured on this lab's data rather than
# assumed:
#
# Clone identity comes from clustering each well's own reads, not from
# annotating them. The reference pipeline runs IgBLAST on every exact
# representative -- and the exact-collapse ahead of it is a no-op here, since
# 100% of reads are unique at Q24 over 1.6kb, so that is IgBLAST on every read:
# about 1.04M queries and roughly 11.6 hours per run, against 8 minutes for the
# entire alignment stage. A well holds one cell and clones differ sharply across
# CDR3, so ~270 reads cluster locally in milliseconds. IgBLAST then runs once
# per well on the finished consensus. Measured paired on 887 wells: productive
# 92.0% against 92.2%, and yield rose from 890 wells to 949.
#
# abPOA replaces mafft. Partial order alignment is not more accurate here -- at
# this depth a column majority is already near-optimal and there is no headroom
# for any aligner to recover -- but it is 8.2x faster end to end at parity,
# validated on three libraries across both loci and both library modes.

_DEST = DEST[0]
_SAMPLE = SAMPLE
_RUNDIR = "{}/{}".format(_DEST, _SAMPLE)
_SCRIPTS = "workflow/scripts"
# Prepended to PATH inside each rule. Empty falls back to whatever is already
# on PATH, so a machine with the tools installed globally still works.
_NP_BIN = param("tools", "env_bin", "")
_NP_PATH = ("export PATH='{}':\"$PATH\"; ".format(_NP_BIN) if _NP_BIN else "")

_NP_INPUTS = NANOPORE_INPUTS if NANOPORE_INPUTS else [READ1_GZ]

# Barcode plates and well/FR1 anchors live with the barcode set, alongside the
# BC1/BC2 files the Illumina pipelines use, so a run selects them the same way.
_BC_DIR = "data/barcodes/{}".format(PRIMER_SET) if PRIMER_SET else "data/barcodes"

# Reverse primers the plate cassette follows; per barcode set, see
# pipelines.yaml. Checked here because they are pasted into a shell command, and
# a typo should stop the run at parse time rather than demultiplex nothing.
_NP_REV_PRIMERS = param("demux", "inner_reverse_primers", []) or []
for _p in _NP_REV_PRIMERS:
    if not re.fullmatch(r"[ACGTacgt]{15,40}", str(_p)):
        raise ValueError("tiago_nanopore: inner_reverse_primers entry {!r} is not a 15-40 nt ACGT primer".format(_p))
_NP_REV_ARGS = " ".join("--inner-reverse-primer {}".format(str(p).upper()) for p in _NP_REV_PRIMERS)


rule demux_nanopore:
    input:
        reads=_NP_INPUTS,
    output:
        clean="{}/nanopore/{}.platewell_clean.fasta".format(_RUNDIR, _SAMPLE),
        partial="{}/nanopore/{}.well_only_partial.fasta".format(_RUNDIR, _SAMPLE),
        summary="{}/qc/nanopore/{}.demux_summary.tsv".format(_RUNDIR, _SAMPLE),
    log:
        "logs/{}.nanopore_demux.log".format(_SAMPLE),
    params:
        prefix="{}/nanopore/{}".format(_RUNDIR, _SAMPLE),
        sample=_SAMPLE,
        anchors="{}/BCFR1_anchors.tsv".format(_BC_DIR),
        plates="{}/plate_barcodes.tsv".format(_BC_DIR),
        well_mm=param("demux", "well_barcode_mm", 2),
        plate_mm=param("demux", "plate_barcode_mm", 1),
        primer_mm=param("demux", "primer_mm", 2),
        min_insert=param("demux", "min_insert_length", 250),
        rev_primers=_NP_REV_ARGS,
        extra=lambda w, input: " ".join(
            "--extra-fastq '{}'".format(p) for p in list(input.reads)[1:]
        ),
    shell:
        r"""
        set -euo pipefail
        {_NP_PATH}
        mkdir -p "$(dirname {output.clean})" "$(dirname {output.summary})"
        python3 {_SCRIPTS}/demux_nanopore_wells.py \
            "$(echo {input.reads} | cut -d' ' -f1)" \
            '{params.anchors}' '{params.plates}' '{params.prefix}' \
            --sample '{params.sample}' \
            --well-barcode-mm {params.well_mm} \
            --plate-barcode-mm {params.plate_mm} \
            --primer-mm {params.primer_mm} \
            --min-insert-length {params.min_insert} \
            {params.rev_primers} \
            {params.extra} > {log} 2>&1
        # A library sequenced one plate per file legitimately has no plate
        # cassette, and its reads land in the well-only file instead. Both
        # outputs are declared so neither path is a missing-output failure.
        touch {output.clean} {output.partial}
        cp "{params.prefix}.platewell_summary.tsv" {output.summary} 2>/dev/null || \
            printf 'section\tid\tmetric\tvalue\n' > {output.summary}
        """


rule library_qc_nanopore:
    """What the PCR made: full-length products, short products, primer dimers,
    concatemers, and reads too broken to tell.

    The demultiplexer only keeps what it can use, and its counters do not say
    why a library yields little: a library whose PCR made mostly primer dimers
    reads as a run with too few reads until this is measured. Descriptive only:
    it changes nothing downstream."""
    input:
        reads=_NP_INPUTS,
    output:
        summary="{}/qc/nanopore/{}.library_qc.tsv".format(_RUNDIR, _SAMPLE),
        inserts="{}/qc/nanopore/{}.library_qc_inserts.tsv".format(_RUNDIR, _SAMPLE),
        wells="{}/qc/nanopore/{}.library_qc_wells.tsv".format(_RUNDIR, _SAMPLE),
    log:
        "logs/{}.nanopore_library_qc.log".format(_SAMPLE),
    threads: 4
    params:
        prefix="{}/qc/nanopore/{}".format(_RUNDIR, _SAMPLE),
        anchors="{}/BCFR1_anchors.tsv".format(_BC_DIR),
        plates="{}/plate_barcodes.tsv".format(_BC_DIR),
        well_mm=param("demux", "well_barcode_mm", 2),
        plate_mm=param("demux", "plate_barcode_mm", 1),
        primer_mm=param("demux", "primer_mm", 2),
        min_insert=param("demux", "min_insert_length", 250),
        dimer_max=param("qc", "dimer_max_insert", 100),
        rev_primers=_NP_REV_ARGS,
        reads=lambda w, input: " ".join("'{}'".format(p) for p in input.reads),
    shell:
        r"""
        set -euo pipefail
        {_NP_PATH}
        mkdir -p "$(dirname {output.summary})"
        python3 {_SCRIPTS}/nanopore_library_qc.py {params.reads} \
            --well-anchors '{params.anchors}' --plate-barcodes '{params.plates}' \
            --prefix '{params.prefix}' {params.rev_primers} \
            --primer-mm {params.primer_mm} --well-barcode-mm {params.well_mm} \
            --plate-barcode-mm {params.plate_mm} --min-insert {params.min_insert} \
            --dimer-max {params.dimer_max} --jobs {threads} > {log} 2>&1
        """


rule consensus_nanopore:
    input:
        clean=rules.demux_nanopore.output.clean,
        partial=rules.demux_nanopore.output.partial,
    output:
        fasta="{}/nanopore/{}.consensus_per_well.fasta".format(_RUNDIR, _SAMPLE),
        table="{}/qc/nanopore/{}.consensus_per_well.tsv".format(_RUNDIR, _SAMPLE),
    log:
        "logs/{}.nanopore_consensus.log".format(_SAMPLE),
    threads: 8
    params:
        prefix="{}/nanopore/{}".format(_RUNDIR, _SAMPLE),
        min_reads=param("consensus", "min_reads_per_well", 10),
        rescue_min=param("consensus", "rescue_min_reads", 5),
        min_dominance=param("consensus", "min_dominance", 0.70),
        max_unique=param("consensus", "max_unique", 200),
        cluster_id=param("consensus", "cluster_id", 0.85),
        engine=param("consensus", "engine", "abpoa"),
    shell:
        r"""
        set -euo pipefail
        {_NP_PATH}
        # Whichever demultiplex output actually carries reads is the input:
        # plate-multiplexed libraries fill the platewell file, one-plate-per-file
        # libraries fill the well-only file. Choosing by size rather than by
        # configuration means a mislabelled run still processes correctly.
        # Prefer the plate-resolved reads whenever there are any. Choosing by
        # file size is wrong: a plate-multiplexed library still produces a large
        # well-only file from reads whose plate cassette was not read, and on
        # one heavy-chain library that file (45 MB) is bigger than the plate-resolved one (24 MB).
        # Selecting on size processed 92 wells instead of 949, silently.
        src={input.clean}
        if [ ! -s "$src" ]; then
            src={input.partial}
        fi
        python3 {_SCRIPTS}/nanopore_consensus.py "$src" '{params.prefix}' \
            --cluster-wells --engine '{params.engine}' \
            --min-reads {params.min_reads} --rescue-min-reads {params.rescue_min} \
            --min-dominance {params.min_dominance} --max-unique {params.max_unique} \
            --cluster-id {params.cluster_id} --jobs {threads} > {log} 2>&1
        cp "{params.prefix}.consensus_per_well.tsv" {output.table}
        """


rule annotate_long_nanopore:
    """One IgBLAST query per well, on the finished consensus.

    This is the whole point of the redesign. The reference pipeline annotates
    every read to get clone identity and VDJ coordinates -- about 1.04M queries
    and roughly 11.6 hours. Clone identity now comes from clustering, so the only
    thing still needed from a germline database is where V starts and J ends, and
    that can be read off the consensus itself."""
    input:
        fasta=rules.consensus_nanopore.output.fasta,
    output:
        airr="{}/qc/nanopore/{}.long_consensus.airr.tsv".format(_RUNDIR, _SAMPLE),
    log:
        "logs/{}.nanopore_annotate.log".format(_SAMPLE),
    threads: 8
    params:
        # Relative to the pipeline root, which is snakemake's working directory.
        # This was an environment variable and did not survive the hop into the
        # SLURM job: the path resolved to /igblast and IgBLAST reported a missing
        # germline database rather than a missing variable.
        igroot=param("annotate", "igblast_root", "igblast"),
        # The germline panel follows the RUN, not the pipeline default. A light
        # chain library annotated against the heavy chain database produces
        # wrong V and J calls, and every downstream number -- productive rate,
        # clonality, mutation load -- inherits that silently. The frontend
        # already collects this per run as igblast_panel.
        panel=(str(IGBLAST_PANEL[0]).strip()
               if IGBLAST_PANEL and str(IGBLAST_PANEL[0]).strip().lower() not in ("", "none", "nan")
               else param("annotate", "panel", "igh")),
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.airr})"
        W="{params.igroot}/work/mouse/{params.panel}"
        export IGDATA="{params.igroot}"
        # IgBLAST wants a D database even for loci that have no D gene. Without
        # one it falls back to a default name that does not exist and dies with
        # "No alias or index file found for [mouse_gl_D]". Light chain therefore
        # borrows the heavy chain D database, which is what it is for.
        if [ -f "$W/mouse_{params.panel}_D.nhr" ]; then
            D="-germline_db_D $W/mouse_{params.panel}_D"
        else
            D="-germline_db_D {params.igroot}/work/mouse/igh/mouse_igh_D"
        fi
        igblastn -germline_db_V "$W/mouse_{params.panel}_V"                  -germline_db_J "$W/mouse_{params.panel}_J" $D                  -auxiliary_data "{params.igroot}/optional_file/mouse_gl.aux"                  -organism mouse -ig_seqtype Ig -domain_system imgt                  -query {input.fasta} -outfmt 19 -num_threads {threads}                  -out {output.airr} > {log} 2>&1
        """


rule vdj_consensus_nanopore:
    """The second sequence per well: VDJ only, called more strictly.

    Separate from the long consensus rather than cropped from it, because
    mutation analysis needs a higher bar per base -- a miscalled base here is
    reported as a somatic hypermutation that did not happen."""
    input:
        clean=rules.demux_nanopore.output.clean,
        partial=rules.demux_nanopore.output.partial,
        long_fa=rules.consensus_nanopore.output.fasta,
        airr=rules.annotate_long_nanopore.output.airr,
    output:
        fasta="{}/nanopore/{}.vdj_per_well.fasta".format(_RUNDIR, _SAMPLE),
        table="{}/qc/nanopore/{}.vdj_per_well.tsv".format(_RUNDIR, _SAMPLE),
    log:
        "logs/{}.nanopore_vdj.log".format(_SAMPLE),
    threads: 8
    params:
        prefix="{}/nanopore/{}".format(_RUNDIR, _SAMPLE),
        min_q=param("vdj", "min_mean_q", 12.0),
        min_reads=param("vdj", "min_reads", 10),
        max_unique=param("vdj", "max_unique", 300),
        min_depth=param("vdj", "min_depth", 10),
        min_frac=param("vdj", "min_base_fraction", 0.70),
    shell:
        r"""
        set -euo pipefail
        {_NP_PATH}
        src={input.clean}
        if [ ! -s "$src" ]; then
            src={input.partial}
        fi
        python3 {_SCRIPTS}/nanopore_vdj_consensus.py "$src" {input.long_fa} {input.airr}             '{params.prefix}' --min-mean-q {params.min_q} --min-reads {params.min_reads}             --max-unique {params.max_unique} --min-depth {params.min_depth}             --min-base-fraction {params.min_frac} --jobs {threads} > {log} 2>&1
        cp "{params.prefix}.vdj_per_well.tsv" {output.table}
        """


rule isotype_nanopore:
    """Constant-region tail and IgA locus per well.

    The tail is a presence test on the consensus. The locus is not: wild-type
    against locus-swapped turns on a single discriminating base, which is exactly
    where nanopore error becomes a wrong answer rather than a noisy one. That
    base is counted from the well's own reads mapped back to its consensus, and
    the call requires 20-fold coverage and 90% agreement -- below either, the
    well is reported as MIXED or INSUFFICIENT_COVERAGE rather than guessed."""
    input:
        clean=rules.demux_nanopore.output.clean,
        partial=rules.demux_nanopore.output.partial,
        long_fa=rules.consensus_nanopore.output.fasta,
    output:
        table="{}/qc/nanopore/{}.isotype_locus.tsv".format(_RUNDIR, _SAMPLE),
    log:
        "logs/{}.nanopore_isotype.log".format(_SAMPLE),
    threads: 4
    shell:
        r"""
        set -euo pipefail
        {_NP_PATH}
        src={input.clean}
        if [ ! -s "$src" ]; then src={input.partial}; fi
        python3 {_SCRIPTS}/nanopore_isotype.py "$src" {input.long_fa} {output.table}             --jobs {threads} > {log} 2>&1
        """


rule clonality_nanopore:
    """Clonal grouping across wells, at three identity thresholds.

    Two sequences are the same clone when they share a V gene, a J gene and a
    CDR3 length, and their CDR3s are within the threshold. Reported both within
    a sample and globally, because the same clone appearing in two samples means
    something different from it appearing twice in one."""
    input:
        airr=rules.annotate_long_nanopore.output.airr,
        isotype=rules.isotype_nanopore.output.table,
        cons=rules.consensus_nanopore.output.table,
    output:
        table="{}/qc/nanopore/{}.per_well_final.tsv".format(_RUNDIR, _SAMPLE),
        csv="{}/qc/nanopore/{}.per_well_final.csv".format(_RUNDIR, _SAMPLE),
    log:
        "logs/{}.nanopore_clonality.log".format(_SAMPLE),
    shell:
        r"""
        set -euo pipefail
        Rscript {_SCRIPTS}/nanopore_clonality.R {input.airr} {input.isotype}             {input.cons} {output.table} > {log} 2>&1
        """


rule emit_nanopore:
    input:
        fasta=rules.consensus_nanopore.output.fasta,
        table=rules.consensus_nanopore.output.table,
        vdj=rules.vdj_consensus_nanopore.output.fasta,
        final_table=rules.clonality_nanopore.output.table,
        # Not used by this rule: listed so the library QC exists before final
        # builds the report and the archive.
        library_qc=rules.library_qc_nanopore.output.summary,
    output:
        marker="{}/split1/{}".format(_RUNDIR, _SAMPLE),
        final="{}/split1/split2/collapsed/head/final/final.fasta".format(_RUNDIR),
        collapsed="{}/collapsed.stats.txt".format(_RUNDIR),
    log:
        "logs/{}.nanopore_emit.log".format(_SAMPLE),
    shell:
        r"""
        set -euo pipefail
        # The marker is a sentinel FILE, matching what final.smk consumes from
        # every other pipeline. Creating it as a directory raises
        # ImproperOutputException: snakemake requires directory outputs to be
        # declared with directory(), and this one is not a directory.
        mkdir -p "$(dirname {output.final})" "$(dirname {output.marker})"
        touch {output.marker}
        python3 {_SCRIPTS}/nanopore_emit_fasta.py {input.fasta} {output.final}
        # collapsed.stats.txt is what plate_stats.R reads to draw the plate maps.
        # One consensus per well means the per-well count IS the read support,
        # which is the number worth plotting.
        python3 {_SCRIPTS}/nanopore_emit_stats.py '{input.table}' '{output.collapsed}' > {log} 2>&1
        n=$(grep -c '^>' {output.final} || true)
        echo "tiago_nanopore emitted ${{n}} per-well consensus sequences" >> {log}
        """
