# tiago_sulkseq -- MIBS libraries, where clone size is incidence.
#
# Multiwell Incidence-Based Sequencing sorts ~100 cells into each of 16-32 wells
# and measures a clone by the number of wells it turns up in, not by its read
# depth. That single change makes the other pipelines the wrong shape: they all
# reduce a well to its dominant sequence, which assumes one cell per well.
#
# The front half is rules/illumina_core.smk, unchanged and shared with
# tiago_vsearch. What differs starts at the well:
#
#   rank     keeps every cluster above a read floor instead of the top few,
#            because a hundred-cell well has no dominant clone to find.
#   annotate IgBLAST over all of them; clone assignment needs V, J and junction.
#   clones   groups ACROSS wells within a sample -- the opposite direction from
#            every other clonality step here -- and then applies the paper's two
#            filters: at least 6 reads in some well, present in at least 2 wells.
#            The second filter is the method's error control. A chimera appears
#            once; a clone sampled repeatedly from an expanded population does
#            not.
#
# The clone rule was validated against the 586,821 annotated sequences the MIBS
# repository publishes alongside their own clone calls: of the pairs this rule
# joins, 100.00% are joined by the published calls, and it recovers 99.78% of
# theirs. After the published filters the wells-per-clone distribution matches
# at the median, 90th and 99th percentile, with the same maximum of 32.
#
# Read that with the germline panel in mind. It was measured on sequences
# annotated against the IgH panel, which is what this pipeline used for every
# run until annotate_sulkseq below started following the submission form. V and
# J calls are the buckets clone assignment works inside, so a run submitted with
# a different panel is not the run those numbers describe. The report says which
# panel was submitted and which was applied, and after that change the two are
# the panel this rule annotated with.

_SULK_QC = "{}/qc/sulkseq".format(_RUNDIR)


def _sulk_platemap():
    """The well -> sample map, if the submitter put one beside the FASTQs.

    Without it every well belongs to one sample, which is right for a
    single-sample library and wrong for a real MIBS run holding several mice.
    The clones step says so loudly rather than quietly producing one giant
    sample."""
    for name in ("platemap.tsv", "platemap.csv", "plate_map.tsv", "plate_map.csv",
                 "groups.tsv", "groups.csv", "groups.txt",
                 "experimental_groups.tsv", "experimental_groups.csv"):
        candidate = os.path.join(INPUTPATH, name)
        if os.path.exists(candidate):
            return candidate
    return ""


SULK_PLATEMAP = _sulk_platemap()


# The submitter told the form a design was attached and none of the names above
# turned up beside the FASTQs, so the transfer lost it. Running anyway pools
# every well into one sample and reports a clone count that looks entirely
# normal, which is worse than either working or failing.
if PLATEMAP_DECLARED and not SULK_PLATEMAP:
    raise ValueError(
        "sample.json declares an experimental design (platemap_attached), but no "
        "plate map was found in {}. Expected one of: platemap.tsv, platemap.csv, "
        "plate_map.tsv, plate_map.csv, groups.tsv, groups.csv, groups.txt, "
        "experimental_groups.tsv, experimental_groups.csv.".format(INPUTPATH)
    )


rule rank_sulkseq:
    """Every cluster in every well, with its read count, rather than the top N.

    The read floor is a cost control, not the paper's filter: annotating each
    well's singleton error cloud would multiply the IgBLAST input by an order of
    magnitude to add sequences that cannot reach the 6-read threshold on their
    own. The published filter is applied later, on reads summed per clone per
    well, which is where it belongs."""
    input:
        sentinel=rules.wells_vsearch.output.sentinel,
    output:
        final="{}/split1/split2/collapsed/head/final/final.fasta".format(_RUNDIR),
        collapsed="{}/collapsed.stats.txt".format(_RUNDIR),
        verdicts="{}/qc/assembly/well_verdicts.tsv".format(_RUNDIR),
    params:
        clustered="{}/split1/split2/clustered/".format(_RUNDIR),
        sample=_SAMPLE,
        funnel="{}/qc/demux/funnel.tsv".format(_RUNDIR),
        # TOPN is the form's "Contigs or isoforms to retrieve" -- for this
        # pipeline, how many cells went into each well. pipelines.yaml can pin
        # it, but by default the submitter decides.
        topn=param("wells", "max_clusters_per_well", TOPN),
        mincluster=param("wells", "min_reads_per_cluster", 2),
        mindom=param("wells", "min_dominant_reads", 2),
        ratio=param("wells", "multiclonal_ratio", 1.0),
    log:
        "logs/{}.sulk_rank.log".format(_SAMPLE),
    shell:
        """
        python3 workflow/scripts/rank_wells.py \
            --wells-dir {params.clustered} --sample {params.sample} \
            --final-fasta {output.final} --collapsed-stats {output.collapsed} \
            --verdicts {output.verdicts} --funnel {params.funnel} \
            --topn {params.topn} --min-cluster-reads {params.mincluster} \
            --min-dominant-reads {params.mindom} \
            --multiclonal-ratio {params.ratio} > {log} 2>&1
        """


rule annotate_sulkseq:
    """IgBLAST over every kept cluster. Clone assignment needs V, J and junction,
    so this cannot wait for the optional annotation step in final.script.sh."""
    input:
        fasta=rules.rank_sulkseq.output.final,
    output:
        airr="{}/{}.airr.tsv".format(_SULK_QC, _SAMPLE),
    params:
        outdir=_SULK_QC,
        sample=_SAMPLE,
        species=IGBLAST_SPECIES[0] if IGBLAST_SPECIES else "mouse",
        # The germline panel follows the RUN, not the pipeline default. yaml
        # defines annotate.panel, and param() returns a key that is present over
        # any fallback, so putting the form's value in the fallback slot meant a
        # light-chain MIBS library was annotated against IgH germlines and every
        # downstream number inherited that silently. Same precedence, and the
        # same reason, as tiago_nanopore.smk.
        panel=(str(IGBLAST_PANEL[0]).strip()
               if IGBLAST_PANEL and str(IGBLAST_PANEL[0]).strip().lower() not in ("", "none", "nan")
               else param("annotate", "panel", "igh")),
        organism=IGBLAST_ORGANISM[0] if IGBLAST_ORGANISM else "mouse",
        aux=IGBLAST_AUX[0] if IGBLAST_AUX else "",
        data=IGBLAST_DATA[0] if IGBLAST_DATA else "",
    log:
        "logs/{}.sulk_annotate.log".format(_SAMPLE),
    shell:
        """
        mkdir -p {params.outdir}
        IGBLAST_DATA='{params.data}' bash workflow/scripts/run_igblast.sh \
            {input.fasta} {params.outdir} {params.sample} \
            '{params.species}' '{params.panel}' '' '' '' '{params.aux}' '{params.organism}' \
            > {log} 2>&1
        mv {params.outdir}/{params.sample}.igblast.tsv {output.airr}
        """


rule clones_sulkseq:
    """Clones across wells, then incidence, then the paper's two filters."""
    input:
        airr=rules.annotate_sulkseq.output.airr,
    output:
        clones="{}/{}.clones.tsv".format(_SULK_QC, _SAMPLE),
        by_well="{}/{}.clone_by_well.tsv".format(_SULK_QC, _SAMPLE),
        public="{}/{}.public_clonotypes.tsv".format(_SULK_QC, _SAMPLE),
        sequences="{}/{}.sequences.tsv".format(_SULK_QC, _SAMPLE),
        summary="{}/{}.summary.tsv".format(_SULK_QC, _SAMPLE),
        stats="{}/{}.sample_stats.tsv".format(_SULK_QC, _SAMPLE),
        resolved="{}/{}.platemap_resolved.tsv".format(_SULK_QC, _SAMPLE),
    params:
        outdir=_SULK_QC,
        sample=_SAMPLE,
        platemap=SULK_PLATEMAP,
        maxdist=param("clones", "junction_max_distance", 4),
        minreads=param("clones", "min_reads_per_well", 6),
        minwells=param("clones", "min_wells_per_clone", 2),
        productive="--productive-only" if param("clones", "productive_only", False) else "",
        grouping=param("clones", "gene_grouping", "strict"),
        # Passed for the first time here, alongside the flag and the summary row
        # the clones script gains in the same change. 'sample' is what the
        # script did unconditionally before any of this, so a run that does not
        # ask for anything else calls exactly the clones it called yesterday.
        scope=param("clones", "clone_scope", "sample"),
    log:
        "logs/{}.sulk_clones.log".format(_SAMPLE),
    shell:
        """
        python3 workflow/scripts/sulkseq_clones.py \
            --airr {input.airr} --sample {params.sample} --outdir {params.outdir} \
            --platemap '{params.platemap}' \
            --max-distance {params.maxdist} --min-reads {params.minreads} \
            --min-wells {params.minwells} --gene-grouping {params.grouping} \
            --clone-scope {params.scope} \
            {params.productive} > {log} 2>&1
        """


rule emit_sulkseq:
    """The sentinel final.smk consumes, once the clone tables exist."""
    input:
        clones=rules.clones_sulkseq.output.clones,
        final=rules.rank_sulkseq.output.final,
    output:
        marker="{}/split1/{}".format(_RUNDIR, _SAMPLE),
    shell:
        """
        touch {output.marker}
        """
