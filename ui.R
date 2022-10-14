#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(shinythemes)
library(shinyalert)
library(spsComps)

shinyUI(fluidPage(theme = shinytheme("cerulean"),
                  useShinyalert(),
                  titlePanel("tbrc pipeline - T and B Receptor Clonality"),

    # Generate form to run pipeline
    sidebarLayout(
        sidebarPanel(
            tags$div(title="The username becomes the name of the folder containing results.",
            textInput   ("useid", "User name:", value = NULL, width = NULL, placeholder = NULL)
            ),
            
            tags$div(title="You will receive an email when the run finished.",
            textInput   ("email", "User email:", value = NULL, width = NULL, placeholder = NULL)
            ),
            
            tags$div(title="If your fastq files are bigger than 500Mb, use the server option.",
            radioButtons("typei", "Method of input", choices = c("Upload", "Server"), selected = "Server")
            ),
            
           tags$div(title="Name of the folder containing fastq files in the server.",
                    p("Your folder name inside the server[Mucida/Victora]/pipelines/tbrc/rawdata/"),
           textInput   ("inpuf",  "Folder in server", value = "")
           ),
            
            # tags$div(title="Deprecated. It self recognizes it.",
            # textInput   ("basename",  "Fastq basename (Deprecated. It self recognizes it.)", value = ""),
            # ),
            # 
            tags$div(title="Which server the data should be sent.",
            radioButtons("serve", "Server to send/get results", choices = c("storage2", "storage1")),
            ),
            
            tags$div(title="Only for upload option. Don't use if fastq is bigger than 500Mb.",
            fileInput   ("files", "Upload fastq files. Don't use if fastq is bigger than 500Mb.", accept = c(".gz", ".fastq"), multiple = T),
            ),
            
            tags$div(title="The cell option are barcode sets. Check if you are using the correct barcode set.",
            radioButtons("prime", "Primer Sets:", choices = c("T cell"  = "T",
                                                              "B cell"  = "B", 
                                                              "B cell2" = "B2", 
                                                              "Upload"  = "New"), selected = "B2")
            ),
            
            tags$div(title="If you are using a different barcode set. Upload it here.",
            fileInput   ("barco", "Upload barcode files:", accept = c(".txt"), multiple = TRUE)
            ),
            
            tags$div(title="If you wish to trim primer sequences used in the library construction.",
            radioButtons("trimp", "Trim Primers:", choices = c(TRUE, FALSE), selected = TRUE),
            ),
            
            tags$div(title="How many contigs you want to get per barcode or well.",
            sliderInput ("nread", "Sequence number to retrieve:", min = 1, max = 500, step = 1, value = 3)
            ),
            
            tags$div(title="The pipeline you want to use for contig assembly.",
            radioButtons("metho", "Method:", choices = c("classic_luka", "trust (Under Construction)", "flash (Under Construction)"))
            ),
            
            tags$div(title="Run the pipeline",
            actionButton("run", "Run")
            ),
        ),

        # Show a plot of the generated distribution
        mainPanel(
            h2("Pipeline submission form"),
            h4("This information will be send to the server. Don't overthink it too much."),
            DT::dataTableOutput("sub.form"),
            
            pgPaneUI(
                pane_id = "thispg",
                titles = c("Creating folders", "Running pipeline", "Finished"),
                pg_ids = c("cre", "run", "fin"),
                title_main = c("Run Progress"),
                opened = TRUE,
                top = "0%",
                right = "0%"
            ),
            
            h2("Check if your barcodes correspond to these presets:"),
            tabsetPanel(
                tabPanel(title = "B cells",
            h4("Barcode 1"),             
            pre(includeText("bc/B/BC1.txt")),
            
            br(),
            h4("Barcode 2"),
            pre(includeText("bc/B/BC2.txt")),
            ),
            
            tabPanel(title = "B cells2",
                     h4("Barcode 1"),
                     pre(includeText("bc/B2/BC1.txt")),
                     br(),
                     h4("Barcode 2"),  
                     pre(includeText("bc/B2/BC2.txt")),
            ),
            
            tabPanel(title = "T cells",
                     h4("Barcode 1"),
                     pre(includeText("bc/T/BC1.txt")),
                     br(),
                     h4("Barcode 2"),
                     pre(includeText("bc/T/BC2.txt")),
            )
            
            
            ),
            
            h3("Recent updates:"),
            h5("Added mouse hover tips."),
            h5("Fastq basename now is self detected."),
            
           
            
        )
    )
))
