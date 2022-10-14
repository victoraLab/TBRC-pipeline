
library(shiny)
library(DT)

server <- function(input, output, session) {

  #Create a submission form with selected options for the pipeline
  form  <- reactive({
    
    fname <- gsub("*.*\\/","", input$inpuf)
    fname <- paste0(input$useid,"/", fname)
    if(grepl("uploading", fname)){
      fname <- paste0(input$useid, "/", gsub("\\.gz","", input$files$name[1]))
    }
        df <- data.frame(row.names = "Submission",
            user_id         = input$useid,
            user_email      = input$email,
            input_method    = input$typei,
            server_name     = input$serve,
            primer_set      = input$prime,
            trim_primers    = input$trimp,
            sequence_number = input$nread*2,
            pipeline        = input$metho,
            source          = fname,
            dest            = paste0("results/", fname)
            #basename        = input$basename
            )
        df <- as.data.frame(t(df))
        
        return(df)
    })
  
  #If the user uploads any fastq file, automatically changes the option from Server to Upload.
  observe({
    if(!is.null(input$files$datapath)){
      updateRadioButtons (session, inputId = "typei",
                         choices = c("Upload", "Server"),
                         selected = "Upload")}
    
  })
  
  #When in upload mode, inutilizes the text input for the folder name with data in the server.
  observe({
    if(input$typei == "Upload"){
      updateTextInput (session, inputId = "inpuf",
                       value = "Not used when uploading...")}
    
    if(input$typei == "Server"){
      updateTextInput (session, inputId = "inpuf",
                       value = "")}
    
  })
  
  #If barcodes were uploaded, automatically changed the option for primer set selection.
  observe({
    if(!is.null(input$barco$datapath)){
      updateRadioButtons (session, inputId = "prime",
                         # choices = c("Upload", "Server"),
                          selected = "New")}
    
  })


  #Render form table
  output$sub.form <- renderDataTable(form())
  
  
  #When clicking Run button, run this code.
  
  observeEvent(input$run, {
    #Gives a warning if Upload mode is selected but no fastq files were uploaded.
    if(input$typei == "Upload" & is.null(input$files)){
    shinyalert(title = "Oops!", text = "No fastq files detected! Upload or change input method.", type = "error")}

    #Gives a warning if there are arguments missing in the form
    if(any(form()$Submission == "")){
      shinyalert(title = "Oops!", text = "There are arguments missing. Please complete the form.", type = "error")}
    
    #If Upload mode is on, and fastq files were uploaded, runs the pipeline.
    if(input$typei == "Upload" & !is.null(input$files)){
      
      path <- paste0("runs/", input$useid, "/", gsub("\\.gz","", input$files$name[1] ))
      dir.create(path, showWarnings = T, recursive = T)
      write.table(t(form()), file = paste0(path,"/", "sample.tsv"), quote = F, row.names = F, sep = "\t")
      #$1 = user name
      #$2 = 
      system(paste("./scripts/goUpload.bash", 
                   input$useid, 
                   input$files$datapath[1], input$files$datapath[2], 
                   input$files$name[1], input$files$name[2]))
    }
    
    #If Server mode is on, and all fields were filled, runs the pipeline.
    if(input$typei == "Server" & !any(form()$Submission == "")){
      
      pgPaneUpdate("thispg", "cre", 25)
      pgPaneUpdate("thispg", "cre", 50)
      pgPaneUpdate("thispg", "cre", 100)
      #Creates a run folder within the user id Folder with the folder name
      fname <- gsub("*.*\\/","", input$inpuf) #folder on server name (Run name)
      path <- paste0("runs/", input$useid, "/", fname) #full path to the local run
      dir.create(path, showWarnings = T, recursive = T) #create the folder
      #Write the submission form on the run folder
      write.table(t(form()), file = paste0(path, "/", "sample.tsv"), quote = F, row.names = F, sep = "\t") 
      
      pgPaneUpdate("thispg", "run", 50)
      #script + 1userid + 2server name + 3full path in folder in server name + 4folder name
      system(paste("./scripts/goServ.bash", input$useid, input$serve, input$inpuf, fname)) #Run the snakemake script on the cluster
      pgPaneUpdate("thispg", "run", 100)
      pgPaneUpdate("thispg", "fin", 100)
    }
    
    
    

    
  })
  
  
  
  }
