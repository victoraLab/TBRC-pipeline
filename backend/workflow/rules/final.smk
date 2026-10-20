FINAL_PARAMS = dict(
    base_path=expand("{destiny}/{sample}", destiny=DEST, sample=SAMPLE),
    base_name=SAMPLE,
    trimming_seq=TRIM,
    should_trim=TRIM_PRIMERS,
    preprocess_method=METH,
    server_user=SAMP,
    server_storage=SERVER,
    server_results_path=list(samples_df['server_results_path']),
    server_ssh_port=SERVER_SSH_PORT,
    server_transfer_mode=SERVER_TRANSFER_MODE,
    run_igblast=RUN_IGBLAST,
    run_clonality=RUN_CLONALITY,
    archive_format=ARCHIVE_FORMAT,
    igblast_species=IGBLAST_SPECIES,
    igblast_panel=IGBLAST_PANEL,
    igblast_bin=IGBLAST_BIN,
    igblast_organism=IGBLAST_ORGANISM,
    igblast_db_v=IGBLAST_DB_V,
    igblast_db_d=IGBLAST_DB_D,
    igblast_db_j=IGBLAST_DB_J,
    igblast_data=IGBLAST_DATA,
    igblast_aux=IGBLAST_AUX,
    notification_email=EMAIL,
    should_keep=KEEP,
    archive_scope=ARCHIVE_SCOPE,
    receptor_check=RECEPTOR_CHECK,
    run_group=RUN_GROUP,
    # 'none', not an empty string: see resolve_archive_scope in the Snakefile.
    primer_footprint=param("final", "primer_footprint", "") or "none",
    light_tail=param("wells", "light_chain_3p", "CACTGACACTCCCACGAGCTC"),
    # Carried through purely so the run report can name them.
    pipeline_name=PIPELINE,
    primer_set=PRIMER_SET,
)

FINAL_SHELL = """
        workflow/scripts/final.script.sh -i {input} -o {params.archive} -s '{params.base_name}' -d '{params.base_path}' -t '{params.trimming_seq}' -y '{params.should_trim}' -L '{params.preprocess_method}' -u '{params.server_user}'  -e '{params.notification_email}' -k '{params.should_keep}' -p '{params.server_results_path}' -j '{params.server_storage}' -r '{params.server_ssh_port}' -m '{params.server_transfer_mode}' -g '{params.run_igblast}' -q '{params.run_clonality}' -f '{params.archive_format}' -z '{params.igblast_species}' -l '{params.igblast_panel}' -b '{params.igblast_bin}' -c '{params.igblast_organism}' -v '{params.igblast_db_v}' -w '{params.igblast_db_d}' -x '{params.igblast_db_j}' -I '{params.igblast_data}' -a '{params.igblast_aux}' -A '{params.archive_scope}' -P '{params.pipeline_name}' -B '{params.primer_set}' -R '{params.receptor_check}' -G '{params.run_group}' -F '{params.primer_footprint}' -T '{params.light_tail}' > {log.std} 2> {log.error}
        """

FINAL_ARCHIVE = expand("{destiny}/{sample}/{sample}.{archive_ext}", destiny=DEST, sample=SAMPLE, archive_ext=ARCHIVE_EXT)

rule final:
    input:
        expand("{destiny}/{sample}/split1/{sample}", destiny=DEST, sample=SAMPLE)
    output:
        FINAL_ARCHIVE
    log:
        std=expand("logs/{sample}.final.out.log", sample=SAMPLE),
        error=expand("logs/{sample}.final.err.log", sample=SAMPLE)
    params:
        archive=FINAL_ARCHIVE,
        **FINAL_PARAMS
    shell:
        FINAL_SHELL


# The same final stage again on a finished run, from the assembled FASTA alone.
# `final` depends on the per-well tree its own cleanup deletes, so asking
# Snakemake to rerun it rebuilt the assembly from FASTQ -- and on the vsearch
# path failed outright (handoff 2026-09-17, section 2). Everything final.script.sh
# needs after assembly survives the cleanup: <sample>.final.fasta, which it
# already reuses, and the QC tables and figures of the first pass. Run through
# backend/rerun_final.sh, which keeps a copy of the previous outputs first. The
# archive is pushed to the storage server only when that script is asked to
# (config refinal_push=1); otherwise it stays on the cluster.
rule refinal:
    input:
        expand("{destiny}/{sample}/{sample}.final.fasta", destiny=DEST, sample=SAMPLE)
    output:
        touch(expand("{destiny}/{sample}/refinal.done", destiny=DEST, sample=SAMPLE))
    log:
        std=expand("logs/{sample}.refinal.out.log", sample=SAMPLE),
        error=expand("logs/{sample}.refinal.err.log", sample=SAMPLE)
    params:
        archive=FINAL_ARCHIVE,
        **dict(FINAL_PARAMS, server_transfer_mode=(
            SERVER_TRANSFER_MODE if str(config.get("refinal_push", "")) == "1" else ["pull"]))
    shell:
        FINAL_SHELL
