# Define Shiny UI
shinyUI(fluidPage(theme = shinytheme("flatly"),
  # Custom CSS
  tags$style(HTML("
    :root {
      --sage-050: #f3f7f1;
      --sage-100: #edf4ea;
      --sage-200: #dde8db;
      --sage-300: #cdddcf;
      --sage-400: #8fad8d;
      --sage-500: #7f9e7c;
      --sage-700: #35523c;
      --ink-700: #294034;
      --paper: rgba(255, 255, 255, 0.78);
    }

    body {
      background:
        radial-gradient(circle at top left, rgba(143, 173, 141, 0.18), transparent 28%),
        linear-gradient(180deg, #f3f7f1 0%, #e7efe5 100%);
      color: var(--ink-700);
    }

    .well {
      background-color: transparent;
      border: none;
      box-shadow: none;
      padding: 0;
    }

    .main-panel {
      padding-left: 22px;
    }

    h1, h2, h3, h4, .control-label {
      color: var(--sage-700);
    }

    .hero {
      position: relative;
      overflow: hidden;
      background:
        linear-gradient(0deg, rgba(255,255,255,0.56), rgba(255,255,255,0.56)),
        url('bcr.png') right center / contain no-repeat,
        linear-gradient(180deg, #e7f3fb 0%, #f3f7f1 100%);
      border: 1px solid var(--sage-300);
      border-radius: 28px;
      padding: 30px 34px 28px 34px;
      margin-bottom: 22px;
      box-shadow: 0 20px 45px rgba(71, 99, 78, 0.09);
      min-height: 340px;
    }

    .hero-grid {
      position: relative;
      z-index: 2;
      display: flex;
      align-items: flex-end;
      min-height: 240px;
    }

    .hero-title {
      font-size: 32px;
      font-weight: 700;
      letter-spacing: 0.015em;
      margin: 0 0 10px 0;
      max-width: 760px;
    }

    .hero-subtitle {
      color: #4d6954;
      font-size: 15px;
      margin: 0;
      max-width: 720px;
      line-height: 1.55;
    }

    .hero-badge {
      display: none;
    }

    .hero-copy {
      display: inline-block;
      max-width: 620px;
      background: rgba(255, 255, 255, 0.78);
      border: 1px solid rgba(205, 221, 207, 0.72);
      border-radius: 24px;
      padding: 18px 22px;
      backdrop-filter: blur(4px);
      box-shadow: 0 10px 24px rgba(71, 99, 78, 0.08);
    }

    .panel-card {
      background: var(--paper);
      border: 1px solid var(--sage-300);
      border-radius: 20px;
      padding: 18px 18px 14px 18px;
      margin-bottom: 16px;
      box-shadow: 0 14px 28px rgba(71, 99, 78, 0.08);
      backdrop-filter: blur(4px);
    }

    .card-title {
      font-size: 16px;
      font-weight: 700;
      margin: 0 0 6px 0;
      color: var(--sage-700);
    }

    .card-note {
      font-size: 13px;
      color: #5f7965;
      margin-bottom: 14px;
      line-height: 1.45;
    }

    .metric-card {
      background: linear-gradient(135deg, rgba(143,173,141,0.18), rgba(255,255,255,0.92));
      border: 1px solid var(--sage-300);
      border-radius: 18px;
      padding: 16px 18px;
      margin-bottom: 16px;
    }

    .metric-label {
      text-transform: uppercase;
      letter-spacing: 0.08em;
      font-size: 11px;
      color: #5f7965;
      margin-bottom: 6px;
    }

    .metric-note {
      font-size: 11px;
      color: #5f7965;
      margin-top: 6px;
      line-height: 1.35;
    }

    .metric-value {
      font-size: 22px;
      font-weight: 700;
      color: var(--sage-700);
    }

    /* Readiness checklist above the Run button. Muted when satisfied, so the
       eye lands on what is still outstanding rather than on a wall of ticks. */
    .rdy-box {
      border: 1px solid #DDE2E8; border-radius: 8px; padding: 10px 12px;
      margin-bottom: 12px; background: #FBFCFD;
    }
    .rdy-item { display: flex; align-items: baseline; gap: 8px; padding: 3px 0; font-size: 13.5px; }
    .rdy-mark { width: 14px; flex: 0 0 14px; text-align: center; font-weight: 700; }
    .rdy-done .rdy-mark { color: #1B7A4B; }
    .rdy-done .rdy-label { color: #79818E; }
    .rdy-todo .rdy-mark { color: #A85B10; }
    .rdy-todo .rdy-label { color: #16191F; font-weight: 600; }
    .rdy-todo-text { color: #A85B10; font-size: 12.5px; }
    .rdy-blocked {
      background: #EDF1F4; color: #79818E; border: 1px solid #DDE2E8;
      border-radius: 6px; padding: 11px; font-weight: 600; cursor: not-allowed;
    }
    .rdy-hint { color: #79818E; font-size: 12px; text-align: center; margin-top: 6px; }
    .preview-box {
      background: rgba(248, 251, 246, 0.9);
      border: 1px solid var(--sage-300);
      border-radius: 16px;
      padding: 12px 14px;
      min-height: 140px;
      white-space: pre-wrap;
    }

    .btn-default,
    .bootstrap-select > .dropdown-toggle {
      background-color: var(--sage-100);
      border-color: #c8d8c7;
      color: var(--sage-700);
      border-radius: 10px;
    }

    .bootstrap-select .dropdown-menu {
      border-radius: 12px;
      border-color: #c8d8c7;
    }

    .form-control {
      border-radius: 12px;
      border-color: var(--sage-300);
      box-shadow: none;
    }

    .form-control:focus {
      border-color: var(--sage-400);
      box-shadow: 0 0 0 3px rgba(143, 173, 141, 0.18);
    }

    #run {
      background-color: var(--sage-400);
      border: none;
      color: white;
      padding: 15px 32px;
      text-align: center;
      display: inline-block;
      font-size: 16px;
      margin: 4px 2px;
      cursor: pointer; /* Add a mouse pointer on hover */
      border-radius: 4px; /* Rounded corners */
    }

    #run:hover {
      background-color: var(--sage-500);
    }

    .progress {
      background-color: var(--sage-200);
    }

    .irs-bar,
    .irs-bar-edge,
    .irs-single {
      background: var(--sage-400);
      border-color: var(--sage-400);
    }

    table.dataTable tbody tr {
      background-color: rgba(255, 255, 255, 0.6);
    }

    .status-box {
      margin-top: 14px;
      background: rgba(248, 251, 246, 0.92);
      border: 1px solid var(--sage-300);
      border-radius: 14px;
      padding: 12px 14px;
      color: #45604c;
      white-space: pre-wrap;
      font-size: 13px;
      line-height: 1.5;
    }

    .progress-layout {
      display: grid;
      grid-template-columns: minmax(280px, 0.9fr) minmax(320px, 1.1fr);
      gap: 16px;
      align-items: stretch;
    }

    .progress-stack {
      display: flex;
      flex-direction: column;
      gap: 16px;
      height: 100%;
    }

    @media (max-width: 991px) {
      .main-panel {
        padding-left: 15px;
        margin-top: 12px;
      }

      .hero-title {
        font-size: 24px;
      }

      .progress-layout {
        grid-template-columns: 1fr;
      }

      .hero {
        min-height: 250px;
        padding: 20px;
        background-position: center center, center center, center center;
      }

      .hero-grid {
        min-height: 190px;
      }

    }
  ")),

  useShinyalert(),
  tags$div(
    class = "hero",
    tags$div(
      class = "hero-grid",
      tags$div(
        class = "hero-copy",
        tags$div(
          class = "hero-title",
          "TBRCa Pipeline - T and B cell Receptor Clonality analysis"
        ),
        tags$p(
          class = "hero-subtitle",
          "Plant-based BCR/TCR receptor assembly, QC and annotation"
        )
      )
    )
  ),

  sidebarLayout(
    sidebarPanel(
      tags$div(
        class = "panel-card",
        tags$div(class = "card-title", "Submit Run"),
        tags$div(class = "card-note", "Set the owner, source, and destination for this sequencing run."),
        selectizeInput("user_id", "User name", choices = NULL, selected = NULL, options = list(create = TRUE, placeholder = "Choose an existing user or type a new one")),
        textInput("email", "User email", value = NULL, placeholder = "name@lab.org"),
        app_radio_buttons("type_input", "Method of input", choices = c("Upload", "Server"), selected = "Server", status = "success", outline = TRUE, animation = "pulse"),
        app_radio_buttons("server_name", "Server to send/get results", choices = c("empty"), status = "success", outline = TRUE, animation = "pulse"),
        textInput("input_folder",  "Folder in server", value = "", placeholder = "plate_001"),
        fileInput("input_files", "Upload paired FASTQ files (R1 and R2)", accept = c(".fastq.gz", ".fastq", ".gz"), multiple = TRUE)
      ),
      tags$div(
        class = "panel-card",
        tags$div(class = "card-title", "Pipeline Options"),
        tags$div(class = "card-note", "Choose barcode presets and processing behavior."),
        app_picker_input("primer", "Barcode set", choices = list("Loading..." = "placeholder"), options = list(`live-search` = TRUE)),
        app_picker_input("sequencing_kit", "Sequencing run",
                         choices = list("Loading..." = "placeholder")),
        tags$div(
          class = "card-note",
          style = "margin-top:-8px;",
          paste(
            "Recorded with the submission and not used by the pipeline. How many",
            "barcode plates the library used is counted from the delivered wells",
            "once the run finishes, so there is nothing to fill in here."
          )
        ),
        app_radio_buttons("keep_intermediate", "Keep intermediate files on the compute host", choices = c("FALSE", "TRUE"), selected = "FALSE", status = "success", outline = TRUE, animation = "pulse"),
        tags$div(class = "card-note", "Intermediates are the per-well FASTA tree and trimmed FASTQ left in the run directory. They are never put in the archive."),
        app_radio_buttons(
          "archive_scope",
          "Results archive contents",
          choices = c(
            "Light (recommended)" = "light",
            "Full" = "full"
          ),
          selected = "light",
          status = "success",
          outline = TRUE,
          animation = "pulse"
        ),
        tags$div(class = "card-note", "Light ships the imgt.ready FASTA, every QC/IGBlast/clonality table and the multi-page QC PDFs. Full adds the untrimmed assembled FASTA and one PNG per plate."),
        app_radio_buttons("trim_primers", "Trim primers", choices = c(TRUE, FALSE), selected = TRUE, status = "success", outline = TRUE, animation = "pulse"),
        uiOutput("run_local_control"),
        # On by default. Annotation is why most people submit a run, and leaving
        # it off was the path of least resistance: 12 of the 30 runs still on
        # disk were delivered with no annotation at all, and a delivered archive
        # with no igblast/ folder reads as a fault rather than a choice.
        checkboxInput("run_igblast", "Run IGBlast annotation", TRUE),
        checkboxInput("run_clonality", "Run clonality analysis from IGBlast", TRUE),
        tags$div(
          class = "card-note",
          style = "margin-top:-8px;",
          paste(
            "Both are on by default. Clonality reads the IgBLAST table, so it needs",
            "the box above ticked. Untick them only if you want assembled sequences",
            "alone: without annotation there is no V or J usage, no CDR3 lengths,",
            "no mutation load and no clones."
          )
        ),
        app_radio_buttons(
          "archive_format",
          "Final results archive format",
          choices = c(
            "ZIP (Recommended)" = "zip",
            "tar.gz" = "tar.gz"
          ),
          selected = "zip",
          status = "success",
          outline = TRUE,
          animation = "pulse"
        ),
        app_radio_buttons(
          "igblast_species",
          "IGBlast species",
          choices = c(
            "Human" = "human",
            "Mouse" = "mouse"
          ),
          selected = "mouse",
          status = "success",
          outline = TRUE,
          animation = "pulse"
        ),
        app_picker_input(
          "igblast_panel",
          "IGBlast receptor scope",
          choices = igblast_panel_choices,
          selected = "ig_all",
          options = list(`live-search` = TRUE)
        ),
        tags$div(class = "card-note", style = "margin-bottom:6px;", textOutput("chainHint", inline = TRUE)),
        app_radio_buttons(
          "receptor_check",
          "If the sequences disagree with that scope",
          choices = c(
            "Correct it, and say so in the report" = "auto",
            "Report it, but annotate as chosen" = "warn",
            "Do not check" = "off"
          ),
          selected = "auto",
          status = "success",
          outline = TRUE,
          animation = "pulse"
        ),
        tags$div(
          class = "card-note",
          style = "margin-top:-8px;",
          paste(
            "Before annotating, a sample of the assembled sequences is aligned against both",
            "the immunoglobulin and the T-cell receptor germline databases. Annotating a",
            "T-cell library against immunoglobulin germlines does not fail: it produces an",
            "empty table and a clonality run that finds no clones, and both report success."
          )
        ),
        textInput("trim_sequence", "Trim sequence", value = "GGGAATTCGAGGTGCAGCTGCAGGAGTCTGG"),
        tags$div(
          class = "card-note",
          style = "margin-top:-8px;",
          "Filled in from the barcode set. Change it only for a library that uses a different 5' handle."
        ),
        tags$div(class = "card-note", "The 5' handle that sits between the plate barcode and the V region. The default is the B-cell handle (EcoRI site + IgH FR1 consensus) and matches barcode sets B, B2, B cell TA, LukaP and Carla_Barcodes. A T-cell library uses a different handle: leaving this default with barcode set T discards the whole library, because trimming keeps only reads where this sequence is found."),
        sliderInput("n_read", "Contigs or isoforms to retrieve", min = 1, max = 500, step = 1, value = 3),
        tags$div(class = "card-note", style = "margin-top:-8px;",
                 textOutput("nReadHint", inline = TRUE)),
        app_radio_buttons(
          "method",
          "Core assembly pipeline",
          choices = c(
            "vsearch cluster assembly (tiago_vsearch)" = "tiago_vsearch",
            "Classic clonality assembly (classic_luka)" = "classic_luka",
            "Nanopore long-read assembly (tiago_nanopore)" = "tiago_nanopore",
            "MIBS / sulkseq multiwell incidence (tiago_sulkseq)" = "tiago_sulkseq"
          ),
          selected = "tiago_vsearch",
          status = "success",
          outline = TRUE,
          animation = "pulse"
        ),
        tags$div(
          class = "card-note",
          style = "margin-top:-6px;",
          paste(
            "tiago_vsearch is the default. Benchmarked on seven datasets it removes",
            "the Illumina adapter contamination that made up 24-32% of heavy-chain",
            "and 3-19% of TCR records, and gives each delivered sequence 2-12x more",
            "read support. It reports fewer wells, but 97% of the ones it drops held",
            "a single read. IgBlast puts 12-15 more percentage points of heavy-chain",
            "sequences in-frame and productive, and 18 more for light chain.",
            "Choose tiago_nanopore for Oxford Nanopore runs: one FASTQ instead of",
            "a read pair, wells resolved by plate cassette and FR1 anchor, and one",
            "consensus per well. It is not an option for Illumina data."
          )
        ),
        # MIBS clone calling, exposed only for the pipeline that reads it. These
        # numbers decide what counts as a clone and how many wells it has to
        # reach, which for an incidence assay is the experiment rather than the
        # method; until now they lived only in backend/config/pipelines.yaml on
        # the cluster, so in practice nobody ever changed them.
        #
        # conditionalPanel hides with CSS rather than unbinding, so every input
        # below still exists and still holds its default on a classic_luka
        # submission. server.R depends on that: it writes these columns for
        # every pipeline and subscripts the frame by name, and a column that
        # vanished with the panel would be an error rather than a blank.
        #
        # The experimental design upload lives here too. It was in the barcode
        # card, offered to all four pipelines, and no rule outside MIBS reads it.
        conditionalPanel(
          condition = "input.method == 'tiago_sulkseq'",
          tags$div(class = "card-title", "MIBS options"),
          fileInput("platemap_file", "Experimental groups (plate map)",
                    accept = c(".tsv", ".txt", ".csv")),
          downloadLink("platemapTemplate", "Download a template"),
          tags$div(
            class = "card-note",
            style = "margin-top:6px;",
            paste(
              "Two columns, like a barcode file. Either one row per plate",
              "(P01<tab>GF_M1) with duplicates where plates share a group, or one",
              "row per well (P01A01<tab>GF_M1) when a plate is split. A third",
              "column, if you add one, is the broader experimental group."
            )
          ),
          uiOutput("platemapPreview"),
          tags$div(
            class = "card-note",
            style = paste0("border-top:1px solid rgba(0,0,0,.08); margin-top:12px; ",
                           "padding-top:10px; font-size:16px; font-weight:700; ",
                           "color:var(--sage-700);"),
            "Clone calling"
          ),
          numericInput("mibs_min_reads_per_well", "Minimum reads for a clone in a well",
                       value = 6, min = 1, max = 10000, step = 1),
          numericInput("mibs_min_wells_per_clone", "Minimum wells a clone must appear in",
                       value = 2, min = 1, max = 96, step = 1),
          tags$div(
            class = "card-note",
            style = "margin-top:-8px;",
            paste(
              "The published MIBS filters: at least 6 reads for a clone in some well,",
              "and present in at least 2 wells of the same sample. The second is the",
              "error control -- a PCR artefact or a chimera appears once. Lowering",
              "either one admits artefacts as clones."
            )
          ),
          numericInput("mibs_junction_max_distance", "Junction Hamming distance for one clone",
                       value = 4, min = 0, max = 30, step = 1),
          tags$div(
            class = "card-note",
            style = "margin-top:-8px;",
            paste(
              "Two sequences join one clone when they share a V gene, share a J gene,",
              "and their junctions differ by no more than this many nucleotides. The",
              "published rule is 4."
            )
          ),
          numericInput("mibs_min_reads_per_cluster", "Minimum reads to keep a sequence in a well",
                       value = 2, min = 1, max = 1000, step = 1),
          tags$div(
            class = "card-note",
            style = "margin-top:-8px;",
            paste(
              "A sequence seen in fewer reads than this is dropped before annotation.",
              "Raise it to about 10 for a deeply sequenced library; leave it low",
              "otherwise, because the filter above works on reads summed per clone per",
              "well and this one drops sequences before they can contribute to that sum."
            )
          ),
          app_radio_buttons(
            "mibs_gene_grouping",
            "How gene calls decide clone membership",
            choices = c(
              "Strict (a shared gene call is required)" = "strict",
              "Transitive (a chain of ambiguous calls may link)" = "transitive"
            ),
            selected = "strict",
            status = "success",
            outline = TRUE,
            animation = "pulse"
          ),
          tags$div(
            class = "card-note",
            style = "margin-top:-8px;",
            paste(
              "Strict cannot invent a clone. Transitive splits less when V calls are",
              "noisy, and on de novo assembled sequences they are."
            )
          ),
          app_radio_buttons(
            "mibs_clone_scope",
            "Where clones are assigned",
            choices = c(
              "Within each sample" = "sample",
              "Across the whole run, then attributed to a sample (as published)" = "global"
            ),
            selected = "sample",
            status = "success",
            outline = TRUE,
            animation = "pulse"
          ),
          tags$div(
            class = "card-note",
            style = "margin-top:-8px;",
            paste(
              "The published pipeline clustered every sample together and only then",
              "attributed each clone to a sample. Clustering inside each sample instead",
              "is stricter and gives finer clones. This changes the answer, not the",
              "presentation. It also changes what the step costs: clone assignment",
              "compares every pair of junctions inside a V, J and length bucket, and",
              "pooling forty samples makes each bucket forty times larger, so a run",
              "that finishes comfortably per sample can outlast its wall clock",
              "globally."
            )
          ),
          checkboxInput("mibs_productive_only", "Use productive rearrangements only", FALSE),
          tags$div(
            class = "card-note",
            style = "margin-top:-8px;",
            paste(
              "Off by default. A de novo assembled sequence can be called unproductive",
              "for a frameshift the assembly introduced, and dropping it loses a well",
              "the clone was really in. Turn it on only for a library clean enough that",
              "an out-of-frame junction means what it says."
            )
          )
        ),
        # Light-chain rescue is a barcode-splitting mode read by genfasta.smk
        # and illumina_core.smk. No nanopore rule reads METH, and a wrong answer
        # there is silent rather than a failed run, so the question is withdrawn
        # instead of being left to be answered wrongly.
        conditionalPanel(
          condition = "input.method != 'tiago_nanopore'",
          app_radio_buttons(
            "luka_light",
            "Barcode preprocessing mode",
            choices = c(
              "Standard preprocessing" = "standard",
              "Light-chain rescue preprocessing" = "luka_light"
            ),
            selected = "standard",
            status = "success",
            outline = TRUE,
            animation = "pulse"
          ),
          tags$div(class = "card-note", "Light-chain rescue removes the 5' handle, the gene-specific forward-primer bases and the reverse-primer tail during barcode splitting, so the primer trim above is not applied a second time at the end of the run.")
        ),
        # Rendered server-side so the button can state its own readiness. See
        # output$run_control: when something is missing it is replaced by a
        # disabled control naming how many things, above a checklist saying
        # which. A button that is always pressable teaches people to press it in
        # order to find out what is wrong.
        uiOutput("run_control")
      ),
      tags$div(
        class = "panel-card",
        tags$div(class = "card-title", "Upload New Barcodes"),
        tags$div(class = "card-note", "Add a BC1 and BC2 pair and sync it to the backend barcode store."),
        textInput("bc_name", "Barcode name"),
        textInput("folder_name", "Folder name"),
        fileInput("BC_files", "Upload BC1 and BC2 files", multiple = TRUE),
        app_action_button("bcupload", "Upload", style = "material-flat", color = "success", size = "sm", block = TRUE)
      )
    ),

    mainPanel(class = "main-panel",
      tabsetPanel(
        id = "main_tabs",
        tabPanel(
          "Submit a run",
          fluidRow(
            column(
              3,
              tags$div(
                class = "metric-card",
                tags$div(class = "metric-label", "This submission"),
                tags$div(class = "metric-value", textOutput("datasetCountDisplay", inline = TRUE)),
                tags$div(class = "metric-note", textOutput("datasetCountNote", inline = TRUE))
              )
            ),
            column(
              3,
              tags$div(
                class = "metric-card",
                tags$div(class = "metric-label", "Estimated Manual Time Saved"),
                tags$div(class = "metric-value", textOutput("savedTimeDisplay", inline = TRUE)),
                tags$div(class = "metric-note", "Estimate only \u2014 assumptions in server_config.json")
              )
            ),
            column(
              3,
              tags$div(
                class = "metric-card",
                tags$div(class = "metric-label", "Processing Estimator"),
                tags$div(class = "metric-value", textOutput("runtimeEstimateDisplay", inline = TRUE))
              )
            ),
            column(
              3,
              tags$div(
                class = "metric-card",
                tags$div(class = "metric-label", "Estimated Cost, All Libraries"),
                tags$div(class = "metric-value", textOutput("libraryCostDisplay", inline = TRUE)),
                tags$div(class = "metric-note", textOutput("libraryCostNote", inline = TRUE)),
                tags$div(class = "metric-note", textOutput("libraryCostPerPlate", inline = TRUE)),
                tags$div(class = "metric-note", textOutput("libraryCostAllRuns", inline = TRUE))
              )
            )
          ),
          fluidRow(
            column(
              12,
              tags$div(
                class = "progress-layout",
                tags$div(
                  class = "panel-card",
                  tags$div(class = "card-title", "Run Progress"),
                  tags$div(class = "card-note", "Submission progress in the app. Detailed live Snakemake rule tracking is not wired yet."),
                  pgPaneUI(
                    pane_id = "thispg",
                    titles = c("Validate", "Stage", "Dispatch", "Workflow", "Export"),
                    pg_ids = c("validate", "stage", "dispatch", "workflow", "export"),
                    title_main = c("Run Progress"),
                    opened = TRUE,
                    top = "0%",
                    right = "0%"
                  )
                ),
                tags$div(
                  class = "progress-stack",
                  tags$div(
                    class = "panel-card",
                    tags$div(class = "card-title", "Status Detail"),
                    tags$div(class = "card-note", "Narrative status updates from the app while the workflow is submitted and packaged."),
                    tags$div(class = "status-box", textOutput("progressDetail"))
                  ),
                  tags$div(
                    class = "panel-card",
                    tags$div(class = "card-title", "Submission Preview"),
                    tags$div(class = "card-note", "This is the metadata that will be written to the run folder and passed to the pipeline."),
                    DT::dataTableOutput("submission_form")
                  )
                )
              )
            )
          ),
          fluidRow(
            column(
              5,
              tags$div(
                class = "panel-card",
                tags$div(class = "card-title", "Barcode Library"),
                tags$div(class = "card-note", "Review the preset barcode folders currently available to the app."),
                app_picker_input("folder_dropdown", "Select a folder", choices = NULL, options = list(`live-search` = TRUE)),
                textOutput("confirmation")
              )
            ),
            column(
              7,
              tags$div(
                class = "panel-card",
                tags$div(class = "card-title", "Barcode Preview"),
                fluidRow(
                  column(6, tags$div(class = "preview-box", verbatimTextOutput("barcode_display1"))),
                  column(6, tags$div(class = "preview-box", verbatimTextOutput("barcode_display2")))
                )
              )
            )
          )
        ),
        tabPanel(
          "Run log",
          tags$div(style = "padding-top:14px;", runlogUI("runlog"))
        ),
        tabPanel(
          "What changed",
          tags$div(style = "padding-top:14px;", changelogUI("changelog"))
        )
      )
    )
  )
))
