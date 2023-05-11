suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinyjs))
suppressPackageStartupMessages(library(shinyAce))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(shinycssloaders))
suppressPackageStartupMessages(library(shinydashboard))
suppressPackageStartupMessages(library(shinyalert))

options(shiny.error = function() {
  stop("An error has occurred")
})
tablestyle  = "height:500px; overflow-y: scroll;"
filterstyle = "padding: 15px;font-size=17px;"

sidebar <- dashboardSidebar(
    disable = T,
    sidebarMenu(
        id="nav",
#        menuItem("Summary", icon = icon("chart-pie"), tabName = "summary"),
        menuItem("Dashboard", icon = icon("table"), tabName = "dashboard"),
        menuItem("circRNA reconstruction details", icon = icon("circle-notch"), tabName = "circRNA")
    )
)

body <- dashboardBody(
            id = "page",
            tags$head(
                tags$link(rel = "stylesheet", type = "text/css", href = "css/app.css"),
                tags$script(src="js/app.js")
            ),
            tabItems(

#------------------------------------------------------------------------------
# DASHBOARD TAB 
#------------------------------------------------------------------------------

  
                  tabItem(tabName = "dashboard",

                    fluidRow( # COUNT BOXES
                      valueBoxOutput("step_first", width="12"),
                      valueBoxOutput("count_circ", width="6"),
                      valueBoxOutput("count_mirna", width="3"),
                      valueBoxOutput("count_gene", width="3")
                    ),

                    fluidRow( # FILTER TABLE CIRCRNA-MIRNA-GENE
                        box(title = "circRNA", status = "primary", height = "595" ,
                            solidHeader = T, width="6", 
                            column(width=12, withSpinner(DT::dataTableOutput("circ_tbl")), 
                                   style = tablestyle)),
                        box(title = "miRNA", status = "primary", height = "595", 
                            solidHeader = T, width="3",
                            column(width=12, withSpinner(DT::dataTableOutput("mirna_tbl")), 
                                   style = tablestyle)),
                        box(title = "Gene", status = "primary", height = "595", 
                            solidHeader = T, width="3",
                            column(width=12, withSpinner(DT::dataTableOutput("gene_tbl")), 
                                   style = tablestyle))
                    ),
                    
                    
                    fluidRow(
                        valueBoxOutput("step_second", width="12"),
                        box(title="Filter by one or more circRNA, miRNA and Gene", status="primary",
                            subtitle = "Select circRNA, miRNA and/or genes to compose circRNA-miRNA-gene inte",
                            height="auto",
                            solidHeader = T, width="12",
                            fluidRow(
                            box(title = "circRNA", status = "primary", height = "auto" ,
                                solidHeader = F, width="3", 
                                column(width=12, withSpinner(DT::dataTableOutput("circ_tbl2")), 
                                       style = tablestyle)),
                            box(title = "miRNA", status = "primary", height = "auto", 
                                solidHeader = F, width="3",
                                column(width=12, withSpinner(DT::dataTableOutput("mirna_tbl2")), 
                                       style = tablestyle)),
                            box(title = "Gene", status = "primary", height = "auto", 
                                solidHeader = F, width="3",
                                column(width=12, withSpinner(DT::dataTableOutput("gene_tbl2")), 
                                       style = tablestyle)),
                            box(title = "Active filters", status = "primary", height = "auto", 
                                solidHeader = T, width="3",
                                fluidRow(column(width=12, 
                                      helpText("Filtered circRNA"),
                                       htmlOutput("filters_circ"))),
                                fluidRow(column(width=12, 
                                       helpText("Filtered miRNA"),
                                       htmlOutput("filters_mirna")
                                       )),
                                fluidRow(column(width=12, 
                                       helpText("Filtered Gene"),htmlOutput("filters_gene"))))
                            )
                        )
                    ),
                    fluidRow( # FILTERS
                        box(title="Filter by circRNA, miRNA and Gene features",  status = "primary", height = "auto",
                            solidHeader = T, width="12",
                            #strong("circRNAs/miRNAs/Genes", style=filterstyle),
                            #column(width=12, htmlOutput("filters"),
                            #       helpText("Data from AT&T (1961) The World's Telephones.")),
                            #br(),br(),br(),br(),br(),
                            column(width=2, offset=0, uiOutput("flt_type"),
                                   helpText("Options: 'validated', 'predicted' or
                                            'disease.drug'.")),
                            column(width=3, offset=0, uiOutput("flt_database"),
                                   helpText("Databases used to retrieve miRNA
                                            and gene interaction in
                                            multiMiR.")),
                            #column(width=2,offset=0,uiOutput("flt_support_type")),
                            column(width=2, offset=0, uiOutput("flt_site_type"),
                                   helpText("miRNA seed regions.")),
                            column(width=3, offset=0, uiOutput("flt_category_H"),
                                   helpText("Gene set categories")), 
                            column(width=2,offset=0,uiOutput("flt_logFC_max"),
                                   helpText("Minimum absolute log2-fold-change.")),
                            br(),br(),br(),br(),br(),br(),br(),br(),br(),
                            #olumn(width=1,offset=0,uiOutput("flt_logFC_min")),
                          
                            #column(width=1,offset=1,uiOutput("flt_pval_min")),
                            column(width=2,offset=0,uiOutput("flt_pval_max"),
                                   helpText("Only genes with p-values under the
                                            cutoff are listed.")),
                            column(width=2, offset=0, uiOutput("flt_chr"),
                                   helpText("filtering by circRNA chromosome.")),
                            column(width=2, offset=0, uiOutput("flt_start_pos"),
                                   helpText("filtering by circRNA starting position.")),
                            column(width=2, offset=0, uiOutput("flt_ending_pos"),
                                   helpText("filtering by circRNA ending position.")),
#                            br(),
#                            column(width=2, offset=2, actionButton("flt_clear", "Reset filters")),
                        )
                    ),
                    fluidRow( # PLOTS
                        box(title="Distribution of miRNA target site types", 
                            status="primary", solidHeader=T, width ="4",
                            column(width=12,textOutput("error_pie"),
                                   plotOutput("plot_sitetype", height = 500))),
                        box(title="Gene set distribution", status="primary",
                            solidHeader=T,width ="8",column(width=12,
                                                            textOutput("error_barplot"),
                                                            plotOutput("plot_hallmark", height = 500)))),
                    
                    fluidRow( # ALL TABLE
                        valueBoxOutput("step_third", width="12"),
                        box(title = "circRNA-miRNA-Gene interaction", status = "primary", height = "auto" , solidHeader = T, width="12",
                            #column(width = 1, offset=0, uiOutput("clip")),
                            column(width = 1, offset=0, downloadButton("downloadCSV", "Download CSV")),
                            column(width = 1, offset=1, downloadButton("downloadXLSX", "Download XLSX")),
                            br(), br(),
                            column(width=12, withSpinner(DT::dataTableOutput("all_tbl"))))
                    )
                )
))

#-------------------------------------------------------------------------------
#  UI generation
#-------------------------------------------------------------------------------

shinyUI(
    dashboardPage(
        dashboardHeader(title = "EasyCircR"),
        sidebar,
        body,
        title = "EasyCircR"
        #tags$style(type="text/css",
        #           ".shiny-output-error { visibility: hidden; }",
        #           ".shiny-output-error:before { visibility: hidden; }")
    )
)

