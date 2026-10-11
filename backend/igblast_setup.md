# IGBlast And Clonality Setup Notes

If you want the optional Shiny `Run IGBlast annotation` and `Run clonality analysis from IGBlast` checkboxes to work, install and configure these pieces on the backend side.

The UI now lets you choose:

- species: `human` or `mouse`
- receptor scope:
  - `IgH only`
  - `Ig light only`
  - `All immunoglobulins`
  - `TCR alpha only`
  - `TCR beta only`
  - `All TCR`
  - `All TCR and immunoglobulins`

So the backend should be organized around a fixed database build layout, not manual per-run V/D/J paths.

## Required in the Snakemake environment

1. `igblastn`
2. `makeblastdb`
3. `edit_imgt_file.pl` from the standalone IgBlast distribution
4. the `internal_data/` directory from the standalone IgBlast distribution
5. R packages:
   - `clonality`
   - `ggplate`
   - `ggplot2`
   - `dplyr`
   - `readr`

Example Conda-side R package install:

```bash
conda activate snakemake
Rscript -e "install.packages(c('ggplate','ggplot2','dplyr','readr'), repos='https://cloud.r-project.org')"
```

If your `clonality_0.10.tar.gz` tarball is on the server, install it into the same environment with:

```bash
conda activate snakemake
R CMD INSTALL /path/to/clonality_0.10.tar.gz
```

You can also install it directly from the Victora Lab GitHub repository:

```bash
conda activate snakemake
Rscript -e "install.packages('remotes', repos='https://cloud.r-project.org')"
Rscript -e "remotes::install_github('victoraLab/clonality')"
```

## Required IGBlast reference setup

1. Put the raw IMGT germline FASTA sets under:
   - [`backend/igblast/db/human`](backend/igblast/db/human)
   - [`backend/igblast/db/mouse`](backend/igblast/db/mouse)
2. Put `edit_imgt_file.pl` under:
   - [`backend/igblast/bin/edit_imgt_file.pl`](backend/igblast/bin)
   or make it available in `PATH`
3. Put the NCBI IgBlast `internal_data/` directory under:
   - [`backend/igblast/internal_data`](backend/igblast)
   or export `IGBLAST_DATA` to that directory on the server
4. Put optional auxiliary files under:
   - [`backend/igblast/refs/human.gl.aux`](backend/igblast/refs)
   - [`backend/igblast/refs/mouse.gl.aux`](backend/igblast/refs)
   or point config at the NCBI `optional_file/` directory
5. Build the combined panel databases with:

```bash
cd /path/to/TBRC/backend/igblast
bash setup_igblast_databases.sh
```

You can also build one subset only:

```bash
bash setup_igblast_databases.sh human ig_all
bash setup_igblast_databases.sh mouse tcr_all
```

The build script now:
- concatenates the relevant raw IMGT FASTAs for the selected panel
- runs `edit_imgt_file.pl` on the combined FASTA
- then runs `makeblastdb` on the edited FASTA

Recommended validation flow on the server:

```bash
cd /path/to/TBRC/backend/igblast
./check_igblast_setup.sh mouse ig_all
./smoke_test_igblast.sh /path/to/sample.imgt.ready.fasta mouse ig_all
```

Use a small real FASTA first. If that works, the pipeline-side IGBlast step should work too.

What `ready` should mean before you redeploy:

1. `./check_igblast_setup.sh mouse ig_all` exits `0`
2. `./check_igblast_setup.sh human ig_all` exits `0`
3. `./check_igblast_setup.sh mouse tcr_all` exits `0`
4. `./smoke_test_igblast.sh /path/to/sample.imgt.ready.fasta mouse ig_all` writes a non-empty TSV
5. `igblastn` in the `snakemake` environment no longer complains about missing `internal_data`

Notes:
- `edit_imgt_file.pl` may live either in [`backend/igblast/bin`](backend/igblast/bin) or in the active `PATH`
- `internal_data` may live either in [`backend/igblast/internal_data`](backend/igblast/internal_data) or wherever `IGBLAST_DATA` points
- species `.gl.aux` files are recommended and used automatically when present
- `cluster.igblast_auxiliary_data` may be either a single aux file or the whole NCBI `optional_file/` directory
- if no aux override is set, the backend tries [`backend/igblast/refs`](backend/igblast/refs) first and then the `optional_file/` directory next to `IGBLAST_DATA`

6. The built databases will appear under:
   - [`backend/igblast/work`](backend/igblast/work)
7. Fill these fields in [`frontend/server_config.json`](frontend/server_config.json):
   - `cluster.igblast_bin`
   - optionally `cluster.igblast_data` if your `internal_data` directory lives outside [`backend/igblast/internal_data`](backend/igblast/internal_data)
   - optionally `cluster.igblast_auxiliary_data` if you want to override the species-specific aux resolution with either one file or an `optional_file/` directory

The pipeline no longer needs manual per-run `V/D/J` database paths in config if the panel builds under [`backend/igblast/work`](backend/igblast/work) are present. Species and receptor scope come from the Shiny form, and the backend resolves the matching built databases automatically.

Minimal command shape:

```bash
igblastn \
  -germline_db_V /path/to/V_db \
  -germline_db_J /path/to/J_db \
  -auxiliary_data /path/to/optional_file \
  -query sample.imgt.ready.fasta \
  -outfmt 19 \
  -out sample.igblast.tsv
```

The D database is optional and is only passed for scopes that include D segments, such as `IgH`, `TRB`, `All TCR`, or `All TCR and immunoglobulins`.

Raw-reference naming convention expected by the build script:

```bash
backend/igblast/db/human/IGHV.fasta
backend/igblast/db/human/IGHD.fasta
backend/igblast/db/human/IGHJ.fasta
backend/igblast/db/human/IGKV.fasta
backend/igblast/db/human/IGKJ.fasta
backend/igblast/db/human/IGLV.fasta
backend/igblast/db/human/IGLJ.fasta
backend/igblast/db/human/TRAV.fasta
backend/igblast/db/human/TRAJ.fasta
backend/igblast/db/human/TRBV.fasta
backend/igblast/db/human/TRBD.fasta
backend/igblast/db/human/TRBJ.fasta
```

## What the pipeline now does

- Produces `*.imgt.ready.fasta`
- Runs IGBlast optionally and writes `igblast/<sample>.igblast.tsv`
- Runs clonality optionally from the IGBlast AIRR-style TSV and writes:
  - `clonality/<sample>.clonality.tsv`
  - `clonality/<sample>.clonality.summary.tsv`
- Resolves IGBlast databases from the chosen species and receptor scope
- Generates plate read-depth heatmaps from FASTA header well IDs under `qc/plates/`
- Skips redundant `pre_trim` contig-length plots when primer trimming was not requested

## What still needs your server-specific decisions

- final auxiliary data files for human and mouse
- whether you want to expose additional narrow scopes like `TRG only` or `TRD only` later
- whether the clonality package tarball will live in the repo, home directory, or a shared software path

## Fast server checklist

```bash
conda activate snakemake
which igblastn
which makeblastdb
perl backend/igblast/bin/edit_imgt_file.pl 2>/dev/null || true

cd /path/to/TBRC/backend/igblast
./setup_igblast_databases.sh mouse ig_all
./setup_igblast_databases.sh human ig_all
./setup_igblast_databases.sh mouse tcr_all
./check_igblast_setup.sh mouse ig_all
./check_igblast_setup.sh human ig_all
./check_igblast_setup.sh mouse tcr_all
./smoke_test_igblast.sh /path/to/sample.imgt.ready.fasta mouse ig_all
```
