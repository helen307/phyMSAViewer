library(shiny)
library(ape)
library(msa)
library(ggtree)
library(ggmsa)
library(seqinr)

ui <- fluidPage(

    # App title ----
    titlePanel("Visualize your phylogenetic trees with MSA"),

    fluidRow(


        column(2,
               h3("Description"),
               helpText("phyMSAViewer is a helpful tool that enables users to have a general overview on multiple selected sequences. It also provides multiple sequence alignment on the side and enables users to zoom on a particular range of the multiple sequence alignment to have a better view of the amino acids.",
               )),

        column(1,
               textInput("text", h3("Enter Uniprot ID"))),

        actionButton("goButton", "Go!"),

        column(1,
               h3("Note"),
               helpText("Please enter the Uniprot entries like:  ",
                        "AC=P19838 OR AC=Q00653 OR AC=Q01201. ",
                        "Please wait for 20 seconds for the ",
                        "program to run.")),

        sliderInput("range",
                    label = "Range of interest (MSA):",
                    min = 0, max = 1600, value = c(0, 100))
    ),

    # Sidebar layout with input and output definitions ----
    sidebarLayout(

        # Sidebar panel for inputs ----
        sidebarPanel(

            # Input: Slider for the number of bins ----


        ),

        # Main panel for displaying outputs ----
        mainPanel(
            tabsetPanel(type = "tabs",
                        tabPanel("PHY + MSA (full)", plotOutput("phyPlot")),
                        tabPanel("MSA (zoomed)", plotOutput("msaPlot"))
            )
        )
    )
)

server <- function(input, output) {


    shinyjs::disable("goButton")
    Sys.sleep(20)
    shinyjs::enable("goButton")

    output$phyPlot <- renderPlot({
        uniID <- input$text
        mybank <- choosebank(bank = "swissprot")
        seq1 <- query("relSeq", uniID)
        # get sequence
        seq2 <- getSequence(seq1)

        write.fasta(sequences = seq2,
                    names = getName(seq1),
                    nbchar = 80, file.out = "seqs.fasta")
        # read sequence from the fasta file
        mySeqs <- readAAStringSet("seqs.fasta")   # from package Biostrings

        # perform multiple sequence alignment
        to_align <- msa(mySeqs)

        # Build tree
        my_align <- msaConvert(to_align, type="seqinr::alignment")

        # write into fasta
        write.fasta(as.list(my_align$seq),my_align$nam,file.out="msa.fasta")

        d <- dist.alignment(my_align, "identity")
        myTree <- nj(d) # neighbor-joining

        # use ggtree to plot
        ggtree::msaplot(p=ggtree(myTree) + geom_tiplab(hjust=1,vjust=-1),
                        fasta="msa.fasta")
    })


    output$msaPlot <- renderPlot({
        ggmsa("relseqs.fasta", start = input$range[1], end = input$range[2])
    })
}

# Run the application
shinyApp(ui = ui, server = server)
