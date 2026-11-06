# Pipeline diagrams

These diagrams are drawn from the code: the Snakefile's `PIPELINE_RULES`, the rule files in `backend/workflow/rules/` and the scripts they call. White boxes are Snakemake rules, pale boxes are steps inside one rule, dashed shapes are files, and amber diamonds and hexagons are the settings that change the path. Every assembly pipeline ends in the same file, `split1/split2/collapsed/head/final/final.fasta`, and the final stage takes it from there.

## How a run flows

A submission from the app or from `./tbrc` runs the same launch scripts, and every pipeline ends in the same assembled FASTA, so everything after assembly is shared.

```mermaid
%% Kept identical to the "How a Run Flows" diagram in README.md; change both together.
flowchart TD
    subgraph submit["Submission"]
        app["Shiny app<br/>writes sample.json"]
        cli["./tbrc run<br/>same sample.json, its own checks"]
        pre["preflight.bash<br/>names, source folder, FASTQ pairing and integrity"]
        launch["goServer.bash: FASTQ in a folder on the NAS,<br/>checks free space on the DTN first<br/>goUpload.bash: FASTQ uploaded in the browser"]
        app -->|"server mode"| pre
        app -->|"upload"| launch
        cli --> launch
        pre --> launch
    end

    nas[("Synology NAS<br/>raw FASTQ")]
    scratch[("HPC scratch<br/>data/runs/USER/RUN/")]
    nas -->|"goServer: rsync, run on the DTN"| scratch
    launch -->|"sample.json, and uploaded FASTQ"| scratch

    launch -->|"ssh to the login node, detached"| iso["run_isolated.sh USER-RUN<br/>its own .parallel/ working directory"]
    iso --> sm["snakemake.sh and the Snakefile<br/>read sample.json, submit SLURM jobs"]
    scratch --> sm
    sm --> choice{"pipeline"}

    choice -->|"classic_luka"| classic["trim_reads, pandaseq, genfasta"]
    choice -->|"tiago_vsearch or tiago_sulkseq"| core["illumina_core<br/>adapter trim, vsearch merge,<br/>demultiplex, cluster each well"]
    core -->|"tiago_vsearch"| vs["rank_and_emit<br/>top N clusters per well"]
    core -->|"tiago_sulkseq"| sulk["rank, annotate, clones, emit<br/>clone size = wells it appears in"]
    choice -->|"tiago_nanopore"| np["demultiplex, abPOA consensus per well,<br/>annotate, VDJ, isotype, clonality, emit"]

    contract[/"final.fasta, the output every pipeline produces<br/>split1/split2/collapsed/head/final/final.fasta"/]
    classic --> contract
    vs --> contract
    sulk --> contract
    np --> contract
    contract --> final

    subgraph final["rule final: final.script.sh"]
        direction TB
        f1["FASTA QC, trim the 5' handle"]
        f2["plate QC"]
        f3["receptor and species check<br/>detect_receptor.py"]
        f4["IgBLAST<br/>per locus when loci need different germline sets"]
        f5["cut the primer footprint per locus,<br/>annotate again"]
        f6["well verdicts per locus"]
        f7["clonality, plate clonotype QC"]
        f8["run report, HTML and PDF"]
        f9["archive: zip or tar.gz, light or full scope"]
        f1 --> f2 --> f3 --> f4 --> f5 --> f6 --> f7 --> f8 --> f9
    end

    f9 -->|"rsync push, or the NAS pulls"| results[("NAS results<br/>results/USER/ or results/USER/GROUP/")]
    rerun["rerun_final.sh USER-RUN"] -.->|"rule refinal, no reassembly"| final
```

## classic_luka

The original method: fastx_trimmer and pandaseq make one merged read per pair, and `split_fasta.sh` demultiplexes, collapses identical reads and keeps the top sequences of each well. On luka_light runs the handle, the light-chain primer and the 3' tail are cut before collapsing, so the final stage does not trim them again.

```mermaid
flowchart TD
    classDef rule fill:#ffffff,stroke:#5f7965,stroke-width:1.5px,color:#294034
    classDef step fill:#f7faf6,stroke:#b8cbb6,color:#35523c
    classDef file fill:#f3f7f1,stroke:#8fad8d,stroke-dasharray:4 3,color:#35523c
    classDef param fill:#fbf3e4,stroke:#c9a46a,color:#5c4420
    classDef contract fill:#e3edfa,stroke:#1c4f8f,stroke-width:1.5px,color:#1c4f8f

    fq[/"R1 and R2 FASTQ"/]:::file
    bc[/"BC1.txt and BC2.txt<br/>of the run's primer_set"/]:::file
    trim["trim_reads<br/>fastx_trimmer -f 3"]:::rule
    panda["pandaseq<br/>pandaseq"]:::rule
    gen["genfasta<br/>split_fasta.sh"]:::rule
    contract[/"final.fasta → Final stage"/]:::contract

    fq --> trim
    trim -->|"SAMPLE.R1.fastq, SAMPLE.R2.fastq<br/>first 2 nt of each mate cut"| panda
    panda -->|"SAMPLE.fasta<br/>one merged read per pair"| gen
    bc --> gen
    gen --> sf

    style sf fill:#f7faf6,stroke:#cdddcf,color:#294034
    subgraph sf["inside genfasta"]
        b1["split on BC1, 5' end<br/>fastx_barcode_splitter<br/>exact match"]:::step
        b2["split on BC2, 3' end<br/>1 mismatch, or exact<br/>if the BC2s are close"]:::step
        meth{"luka_light_protocol"}:::param
        l1["cutadapt: trim_seq handle<br/>reads without it dropped"]:::step
        l2["fastx_trimmer: next 23 nt,<br/>the light-chain primer"]:::step
        l3["cutadapt: 3' primer tail"]:::step
        col["fastx_collapser<br/>identical reads per well"]:::step
        top["head -n sequence_number<br/>top N records per well"]:::step
        b1 -->|"wells under 4 kB dropped"| b2
        b2 -->|"wells under 4 kB dropped"| meth
        meth -->|"luka_light"| l1
        l1 --> l2 --> l3 --> col
        meth -->|"standard"| col
        col --> top
    end

    stats[/"BC.stats.txt, collapsed.stats.txt<br/>read by plate QC"/]:::file
    b2 -.-> stats
    col -.-> stats
    top -->|"wells concatenated,<br/>headers WELL_rank-count"| contract
```

## tiago_vsearch

Adapter trimming, vsearch merging, error-tolerant demultiplexing and 99% clustering in each well, after which the clusters are ranked and the top N of each well delivered as whole, untrimmed reads. Clustering sees a copy of each read with the primer footprint hidden, so a degenerate primer cannot split one clone into several clusters.

```mermaid
flowchart TD
    classDef rule fill:#ffffff,stroke:#5f7965,stroke-width:1.5px,color:#294034
    classDef step fill:#f7faf6,stroke:#b8cbb6,color:#35523c
    classDef file fill:#f3f7f1,stroke:#8fad8d,stroke-dasharray:4 3,color:#35523c
    classDef param fill:#fbf3e4,stroke:#c9a46a,color:#5c4420
    classDef contract fill:#e3edfa,stroke:#1c4f8f,stroke-width:1.5px,color:#1c4f8f

    fq[/"R1 and R2 FASTQ"/]:::file
    bc[/"BC1.txt and BC2.txt<br/>of the run's primer_set"/]:::file
    prep["prep_reads_vsearch<br/>cutadapt"]:::rule
    probe["spacer_probe<br/>probe_spacer.py"]:::rule
    merge["merge_vsearch<br/>vsearch --fastq_mergepairs"]:::rule
    demux["demux_vsearch<br/>demux_wells.sh"]:::rule
    wells["wells_vsearch<br/>vsearch"]:::rule
    rank["rank_and_emit<br/>rank_wells.py"]:::rule
    contract[/"final.fasta → Final stage"/]:::contract

    fq --> prep
    fq -->|"R1"| probe
    bc --> probe
    prep -->|"SAMPLE.prep.R1/R2.fastq.gz<br/>2 nt spacer and<br/>Illumina adapter cut"| merge
    probe -->|"qc/prep/spacer_probe.tsv<br/>stops the run if the spacer is off"| merge
    merge -->|"SAMPLE.merged.fastq<br/>overlap ≥ 9 nt<br/>length ≥ 300"| demux
    bc --> demux
    demux -->|"split1/split2/WELL.fastq<br/>wells with ≥ 2 reads"| wells
    demux -.-> dq[/"BC.stats.txt<br/>qc/demux/funnel.tsv"/]:::file
    wells --> wv

    style wv fill:#f7faf6,stroke:#cdddcf,color:#294034
    subgraph wv["wells_vsearch, per well"]
        pp{"luka_light_protocol"}:::param
        lr["copy: handle cut, reads<br/>without it dropped, then<br/>23 nt primer and 3' tail cut"]:::step
        mk["copy: handle and the next<br/>cluster_mask_nt bases cut"]:::step
        un["vsearch --fastx_uniques"]:::step
        cl["vsearch --cluster_size<br/>--id 0.99, centroids"]:::step
        tp["vsearch --sortbysize<br/>top N clusters"]:::step
        rs["restore_centroid_reads.py<br/>most common whole read"]:::step
        pp -->|"luka_light"| lr
        pp -->|"standard"| mk
        lr --> un
        mk --> un
        un --> cl --> tp --> rs
    end

    mask{{"cluster_mask_nt<br/>23, or 0 for the T set"}}:::param
    mask -.-> mk
    rs -->|"split1/split2/clustered/WELL.top.fa"| rank
    rp{{"N = sequence_number / 2<br/>min_dominant_reads 3<br/>multiclonal_ratio 0.2"}}:::param
    rp -.-> rank
    rank -.-> rq[/"collapsed.stats.txt<br/>qc/assembly/well_verdicts.tsv"/]:::file
    rank -->|"top N per well, untrimmed"| contract
```

## tiago_sulkseq

For MIBS libraries, where a well holds about a hundred cells and a clone's size is the number of wells it appears in. It shares the front half with tiago_vsearch, keeps every cluster above a read floor, annotates them itself and groups clones across wells with the `mibs_*` settings from the form.

```mermaid
flowchart TD
    classDef rule fill:#ffffff,stroke:#5f7965,stroke-width:1.5px,color:#294034
    classDef step fill:#f7faf6,stroke:#b8cbb6,color:#35523c
    classDef file fill:#f3f7f1,stroke:#8fad8d,stroke-dasharray:4 3,color:#35523c
    classDef param fill:#fbf3e4,stroke:#c9a46a,color:#5c4420
    classDef contract fill:#e3edfa,stroke:#1c4f8f,stroke-width:1.5px,color:#1c4f8f

    fq[/"R1 and R2 FASTQ"/]:::file
    bc[/"BC1.txt and BC2.txt<br/>of the run's primer_set"/]:::file

    style core fill:#f7faf6,stroke:#cdddcf,color:#294034
    subgraph core["illumina_core.smk"]
        prep["prep_reads_vsearch<br/>cutadapt"]:::rule
        probe["spacer_probe<br/>probe_spacer.py"]:::rule
        merge["merge_vsearch<br/>vsearch --fastq_mergepairs"]:::rule
        demux["demux_vsearch<br/>demux_wells.sh"]:::rule
        wells["wells_vsearch<br/>vsearch --cluster_size"]:::rule
        prep -->|"prep.R1/R2.fastq.gz"| merge
        probe -->|"spacer_probe.tsv"| merge
        merge -->|"merged.fastq<br/>overlap ≥ 20 nt, ≤ 5 diffs"| demux
        demux -->|"WELL.fastq<br/>wells with ≥ 10 reads"| wells
    end

    fq --> core
    bc --> core

    rank["rank_sulkseq<br/>rank_wells.py"]:::rule
    ann["annotate_sulkseq<br/>run_igblast.sh"]:::rule
    clo["clones_sulkseq<br/>sulkseq_clones.py"]:::rule
    emit["emit_sulkseq<br/>touch"]:::rule
    contract[/"final.fasta → Final stage"/]:::contract

    wells -->|"clustered/WELL.top.fa<br/>whole reads on standard,<br/>trimmed on luka_light"| rank
    m1{{"mibs_min_reads_per_cluster<br/>default 2; at most N kept,<br/>N = sequence_number / 2"}}:::param
    m1 -.-> rank
    rank -->|"final.fasta, every<br/>cluster above the floor"| ann
    pan{{"igblast_panel from the form,<br/>igh when blank"}}:::param
    pan -.-> ann
    ann -->|"qc/sulkseq/SAMPLE.airr.tsv"| clo
    m2{{"mibs_* fields, defaults:<br/>junction_max_distance 4<br/>min_reads_per_well 6<br/>min_wells_per_clone 2<br/>gene_grouping strict<br/>clone_scope sample<br/>productive_only false"}}:::param
    m2 -.-> clo
    pm[/"plate map beside the FASTQ<br/>declared but missing:<br/>the run stops"/]:::file
    pm -.-> clo
    clo -->|"clones.tsv, clone_by_well.tsv,<br/>public_clonotypes.tsv"| emit
    rank -->|"final.fasta"| emit
    emit -->|"marker split1/SAMPLE"| contract
```

## tiago_nanopore

Single-end long reads span the whole amplicon, so there is nothing to merge: reads are assigned to wells by their plate and well barcodes, and each well's dominant clone becomes one abPOA consensus polished with racon. Its tools come from a separate conda environment put first on the PATH.

```mermaid
flowchart TD
    classDef rule fill:#ffffff,stroke:#5f7965,stroke-width:1.5px,color:#294034
    classDef step fill:#f7faf6,stroke:#b8cbb6,color:#35523c
    classDef file fill:#f3f7f1,stroke:#8fad8d,stroke-dasharray:4 3,color:#35523c
    classDef param fill:#fbf3e4,stroke:#c9a46a,color:#5c4420
    classDef contract fill:#e3edfa,stroke:#1c4f8f,stroke-width:1.5px,color:#1c4f8f

    fq[/"nanopore FASTQ<br/>one or more files"/]:::file
    bc[/"BCFR1_anchors.tsv, plate_barcodes.tsv<br/>of the run's primer_set"/]:::file
    rp{{"inner_reverse_primers<br/>per primer set"}}:::param
    demux["demux_nanopore<br/>demux_nanopore_wells.py"]:::rule
    lqc["library_qc_nanopore<br/>nanopore_library_qc.py"]:::rule
    pick{"plate-resolved<br/>reads present?"}:::param
    reads[/"reads per well"/]:::file
    cons["consensus_nanopore<br/>nanopore_consensus.py"]:::rule
    ann["annotate_long_nanopore<br/>igblastn"]:::rule
    vdj["vdj_consensus_nanopore<br/>nanopore_vdj_consensus.py"]:::rule
    iso["isotype_nanopore<br/>nanopore_isotype.py"]:::rule
    clo["clonality_nanopore<br/>nanopore_clonality.R"]:::rule
    emit["emit_nanopore<br/>nanopore_emit_fasta.py"]:::rule
    contract[/"final.fasta → Final stage"/]:::contract

    fq --> demux
    bc --> demux
    rp -.-> demux
    fq --> lqc
    demux -->|"platewell_clean.fasta<br/>well_only_partial.fasta"| pick
    pick -->|"yes: platewell_clean<br/>no: well_only_partial"| reads
    reads --> cons
    cons --> cw

    style cw fill:#f7faf6,stroke:#cdddcf,color:#294034
    subgraph cw["per well"]
        c1["mean Q ≥ 15, 250-3000 nt"]:::step
        c2["vsearch at 0.85 identity,<br/>keep the dominant clone"]:::step
        c3["tiers: ≥ 10 reads,<br/>rescue ≥ 5, flag<br/>dominance under 0.70"]:::step
        c4["abPOA draft<br/>at most 200 unique reads"]:::step
        c5["racon polish<br/>reads mapped by minimap2"]:::step
        c1 --> c2 --> c3 --> c4 --> c5
    end

    lc[/"consensus_per_well.fasta<br/>consensus_per_well.tsv"/]:::file
    c5 --> lc
    lc --> ann
    lc -->|"with the well's reads"| vdj
    lc -->|"with the well's reads"| iso
    ann -->|"long_consensus.airr.tsv"| vdj
    ann --> clo
    lc --> clo
    iso -->|"isotype_locus.tsv"| clo
    vdj -->|"vdj_per_well.fasta"| emit
    clo -->|"per_well_final.tsv"| emit
    lc -->|"headers rewritten"| emit
    lqc -->|"library_qc.tsv"| emit
    emit -->|"final.fasta and<br/>collapsed.stats.txt"| contract
    pan{{"igblast_panel from the form,<br/>igh when blank"}}:::param
    pan -.-> ann
```

## Final stage

`final.script.sh` takes the assembled FASTA from any pipeline through trimming, QC, annotation, clonality, the run report and delivery, and the diamonds are the settings that change its path. `rerun_final.sh` reruns the same stage on a finished run from its kept `SAMPLE.final.fasta`, without assembling again.

```mermaid
flowchart TD
    classDef rule fill:#ffffff,stroke:#5f7965,stroke-width:1.5px,color:#294034
    classDef step fill:#f7faf6,stroke:#b8cbb6,color:#35523c
    classDef file fill:#f3f7f1,stroke:#8fad8d,stroke-dasharray:4 3,color:#35523c
    classDef param fill:#fbf3e4,stroke:#c9a46a,color:#5c4420
    classDef contract fill:#e3edfa,stroke:#1c4f8f,stroke-width:1.5px,color:#1c4f8f

    contract[/"final.fasta → Final stage"/]:::contract
    fin["rule final<br/>final.script.sh"]:::rule
    rerun["rerun_final.sh USER-RUN<br/>previous outputs copied<br/>to previous-DATE/"]:::rule
    refinal["rule refinal<br/>final.script.sh"]:::rule
    ff[/"SAMPLE.final.fasta, always archived<br/>untrimmed unless luka_light trimmed it"/]:::file

    contract --> fin
    fin -->|"final.fasta moved here"| ff
    rerun --> refinal
    refinal -.->|"reused, no reassembly"| ff

    ff --> ht{"handle trim?"}:::param
    defs{{"PRETRIMMED: luka_light,<br/>not tiago_vsearch<br/>HANDLE_REMOVED: tiago_nanopore"}}:::param
    defs -.-> ht
    ht -->|"trim_primers on,<br/>handle still on"| cut["cutadapt -g trim_seq<br/>reads without it dropped<br/>luka_light: 3' tail cut too"]:::step
    ht -->|"trim_primers off, or<br/>PRETRIMMED, HANDLE_REMOVED"| cp["copied unchanged"]:::step
    cut -->|"none kept: run fails"| imgt
    cp --> imgt[/"SAMPLE.imgt.ready.fasta"/]:::file

    imgt --> fqc["FASTA QC<br/>fasta_qc.R"]:::step
    fqc --> pqc["plate QC: plate_stats.R, plate_depth_qc.R<br/>skipped when their stats files are gone"]:::step
    pqc --> ig{"run_igblast"}:::param

    ig -->|"true"| rc{"receptor_check"}:::param
    rc -->|"auto"| dra["detect_receptor.py<br/>corrects panel and species"]:::step
    rc -->|"warn, or auto with<br/>pinned V/J databases"| drw["detect_receptor.py<br/>reports only"]:::step
    rc -->|"off"| rig
    dra --> plan{"per-locus plan?"}:::param
    plan -->|"yes: loci need<br/>different germline sets"| apl["annotate_per_locus.py<br/>each locus on its own set"]:::step
    plan -->|"no"| rig["run_igblast.sh<br/>igblastn"]:::step
    drw --> rig
    apl --> ex{"exit status"}:::param
    rig --> ex
    ex -->|"64: not set up here"| skip["IgBLAST skipped,<br/>run continues"]:::step
    ex -->|"other non-zero"| fail["run fails"]:::step
    ex -->|"0"| vid["median V identity<br/>flagged under 85%"]:::step

    vid --> fp{"footprint cut?"}:::param
    fp -->|"primer_footprint set,<br/>trim_primers on,<br/>not PRETRIMMED"| fpc["trim_primer_footprint.py<br/>IGK=23, by each V call"]:::step
    fpc -->|"records cut"| again["annotate again,<br/>same route"]:::step
    fpc -->|"none cut"| wv
    fp -->|"otherwise"| wv
    again --> wv["per-locus well verdicts<br/>well_verdicts_by_locus.py"]:::step

    wv --> rcl{"run_clonality"}:::param
    skip --> rcl
    ig -->|"false"| rcl
    rcl -->|"true, with a table"| clon["clonality<br/>run_clonality.R"]:::step
    rcl -->|"false"| pcq
    clon --> pcq["plate clonotype QC<br/>plate_clonotype_qc.R"]:::step
    pcq --> rpt["run report, HTML and PDF<br/>run_report.R"]:::step

    rpt --> sc{"archive_scope"}:::param
    sc -->|"light, or unset<br/>with keep_raw false"| al["archive without<br/>per-plate PNGs"]:::step
    sc -->|"full, or unset<br/>with keep_raw true"| af["archive with<br/>per-plate PNGs"]:::step
    al --> zip[/"SAMPLE.zip or SAMPLE.tar.gz"/]:::file
    af --> zip

    zip --> cln["cleanup if keep_raw is false:<br/>split1/, SAMPLE.fasta,<br/>BC and collapsed stats"]:::step
    cln --> tm{"transfer_mode"}:::param
    tm -->|"push"| push["rsync to the storage server,<br/>run_group subfolder if set"]:::step
    tm -->|"pull, or refinal<br/>without --push"| pull["archive stays on the cluster"]:::step
    push --> mail["notification email"]:::step
    pull --> mail
```
