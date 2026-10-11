# tiago_vsearch -- take the top N clusters per well.
#
# The front half is rules/illumina_core.smk, shared with tiago_sulkseq. What is
# specific to this pipeline is the last step: a deterministic total ordering over
# each well's clusters and an explicit per-well verdict, replacing the silent
# top-N cut of `fastx_collapser | head -n`.

rule rank_and_emit:
    input:
        sentinel=rules.wells_vsearch.output.sentinel,
    output:
        # The sentinel final.smk consumes, at the path it already expects.
        marker="{}/split1/{}".format(_RUNDIR, _SAMPLE),
        final="{}/split1/split2/collapsed/head/final/final.fasta".format(_RUNDIR),
        collapsed="{}/collapsed.stats.txt".format(_RUNDIR),
        verdicts="{}/qc/assembly/well_verdicts.tsv".format(_RUNDIR),
    params:
        clustered="{}/split1/split2/clustered/".format(_RUNDIR),
        sample=_SAMPLE,
        funnel="{}/qc/demux/funnel.tsv".format(_RUNDIR),
        topn=TOPN,
        mindom=param("wells", "min_dominant_reads", 3),
        droplow="--drop-low-depth" if param("wells", "drop_low_depth", False) else "",
        ratio=param("wells", "multiclonal_ratio", 0.2),
    log:
        "logs/{}.rank.log".format(_SAMPLE),
    shell:
        """
        python3 workflow/scripts/rank_wells.py \
            --wells-dir {params.clustered} --sample {params.sample} \
            --final-fasta {output.final} --collapsed-stats {output.collapsed} \
            --verdicts {output.verdicts} --funnel {params.funnel} \
            --topn {params.topn} --min-dominant-reads {params.mindom} \
            --multiclonal-ratio {params.ratio} {params.droplow} > {log} 2>&1
        touch {output.marker}
        """
