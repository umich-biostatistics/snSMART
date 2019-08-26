# install.packages("shiny")
# install.packages("shinythemes")
# install.packages('rsconnect')
# library(rsconnect)
# rsconnect::setAccountInfo(name='kidwell',
#                           token='5E5F29B01B30E95B7A3A4B08A97A098E',
#                           secret='WdSG61dQRyIBuwGuFHpnNPweY/cIL5oP1NKbj+Qe')
library(truncdist)
library(rmutil)
library(rootSolve)
library(dplyr)
library(truncnorm)
library(cubature)
library(VGAM)
library(EnvStats)
library(shinythemes)
library(shinyBS)

source("functions.R")

options(encoding = 'UTF-8')

shinyUI(

  fluidPage(theme = shinytheme("cosmo"),
            tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "customize.css")),
            
    fluidRow(

      column(width=12,style="background-color:#E0E0E0;",
             align="center",
             h1("snSMART Sample Size App", 
                tags$p("\t"))
      )
    ),
    hr(),
    
    fluidRow(
      
      # column allocation for widgets
      column(width=2,style='margin-top:-15px;',h2("Description:")),
      # column(width=2,h3("Description:")),
      suppressWarnings(column(width=10,
               tags$ul(
                 tags$li(h4("What is an snSMART"),h4("The small n Sequential, Multiple Assignment, Randomized Trial (snSMART) is a two-stage clinical trial design for rare diseases.  
                            In Stage 1, all patients are randomized to one of the multiple treatments. 
                            In Stage 2, patients who respond to their initial treatment receive the same treatment again, 
                            while those who fail to respond are re-randomized to one of the remaining treatments.  
                            The creators of a Bayesian joint stage model for estimating the response rates of
                            each individual treatment in a three-arm snSMART demonstrated efficiency gains for
                            a given sample size relative to other existing frequentist approaches."),style=""),
                 tags$li(h4("The goal for the snSMART sample size App"), h4("Our App is to conduct Bayesian sample size calculation for an snSMART designed to distinguish the best treatment from
                                                          the second-best treatment using the Bayesian joint stage model."),style="")
             )
      ))
    ),
    hr(),
    
    fluidRow(
      
      # column allocation for widgets
      column(width=2,style='margin-top:-15px;',h2("Design Map:")),
      # column(width=2,h3("Description:")),
      suppressWarnings(column(width=10,
                              tags$ul(
                                tags$li(h4("Figure of the study design of an snSMART"),
                                                                   tags$div(
                                                     h4("Patients are randomized (R) to one of the treatment arms, A, B or C equally (1:1:1) 
                                                        and followed for six months. The responders continue the same treatment for another six months, 
                                                        while the non-responders are re-randomized to one of the remaining treatments for an additional six months.")
                                                     ),
                                                     tags$div(
                                                       tags$img(src = "snSMART_design.png", width = "800px")
                                                     ),style="")
                                                   )
                                                   ))
                                ),
    hr(),
    
    fluidRow(
      id = "form_design",
      # column allocation for widgets
      column(width=2,style='margin-top:-15px;',h2("Design Settings:")),
      # column(width=2,h3("Description:")),
      suppressWarnings(column(width=10,
                              tags$ul(
                                tags$li(h4("Set Design Parameters"),
                                                                   tags$div(
                                                                     fluidRow(
                                                                       tags$style(type = "text/css",
                                                                                  "label {font-size: 17px; font-weight: normal}"
                                                                       ),
                                                                       column(width=4,
                                                                              h4("Treatment Expectation Setting"),
                                                                              tags$div(
                                                                                numericInput("piA",label = HTML("&pi;<sub>A</sub>: <p>The response rate (ranges from 0.01 to 0.99) for treatment A.</p>"),0.50,
                                                                                             min = 0.01, max = 0.99, step = 0.01,width = NULL),
                                                                                
                                                                                numericInput("piB",label = HTML("&pi;<sub>B</sub>: <p>The response rate (ranges from 0.01 to 0.99) for treatment B.</p>"),0.25,
                                                                                             min = 0.01, max = 0.99, step = 0.01,width = NULL),
                                                                                
                                                                                numericInput("piC",label = HTML("&pi;<sub>C</sub>: <p>The response rate (ranges from 0.01 to 0.99) for treatment C.</p>"),0.25,
                                                                                             min = 0.01, max = 0.99, step = 0.01,width = NULL)
                                                                                
                                                                              )
                                                                       ),
                                                                       column(width=4,
                                                                              h4("Linkage Parameter Setting"),
                                                                              tags$div(
                                                                                
                                                                                numericInput("beta1",label = HTML("&beta;<sub>1</sub>: <p>The linkage parameter (ranges from 1.00 to 1/largest response rate) for first stage responders. (A smaller value leads to more conservative sample size calculation because two stages are less correlated.)</p>"),1.5,
                                                                                             min = 1.00, max = 100, step = 0.01,width = NULL),
                                                                                
                                                                                numericInput("beta0",label = HTML("&beta;<sub>0</sub>: <p>The linkage parameter (ranges from 0.01 to 0.99) for first stage non-responders. A larger value leads to a more conservative sample size calculation because two stages are less correlated.)</p>"),0.5,
                                                                                             min = 0.01, max = 0.99, step = 0.01,width = NULL)
                                                                              )
                                                                       ),
                                                                       column(width=4,
                                                                              h4("Coverage and Power Setting"),
                                                                              tags$div(
                                                                                
                                                                                numericInput("coverage",label = HTML("1-&alpha;</sub>: <p>The coverage rate (ranges from 0.01 to 0.99) for the posterior difference of top two treatments.</p>"),0.9,
                                                                                             min = 0.01, max = 0.99, step = 0.01,width = NULL),
                                                                                
                                                                                numericInput("power",label = HTML("1-&xi;: <p>The probability (ranges from 0.01 to 0.99) for identify the best treatment.</p>"),0.8,
                                                                                             min = 0.01, max = 0.99, step = 0.01,width = NULL)
                                                                              )
                                                                       )
                                                                     )
                                                                   )
                                                                   ,style="")
                                )
      ))
  ),
  hr(),
  
  fluidRow(
    id = "form_prior",
    # column allocation for widgets
    column(width=2,style='margin-top:-7px;',h2("Prior Settings:")),
    # column(width=2,h3("Description:")),
    suppressWarnings(column(width=10,
                            tags$ul(
                              tags$li(h4("Set Priors for Treatment Response Rates"),
                                                                 tags$div(
                                                                   fluidRow(
                                                                     tags$style(type = "text/css",
                                                                                "label {font-size: 17px; font-weight: normal}"
                                                                     ),
                                                                     column(width=6,
                                                                            h4("Treatment Prior mean Settings"),
                                                                            tags$div(
                                                                              
                                                                              numericInput("muA",label = HTML("&mu;<sub>A</sub>: <p>The prior mean (ranges from 0.01 to 0.99) for treatment A.</p>"),0.50,
                                                                                           min = 0.01, max = 0.99, step = 0.01,width = NULL),
                                                                              
                                                                              numericInput("muB",label = HTML("&mu;<sub>B</sub>: <p>The prior mean (ranges from 0.01 to 0.99) for treatment B.</p>"),0.25,
                                                                                           min = 0.01, max = 0.99, step = 0.01,width = NULL),
                                                                              
                                                                              numericInput("muC",label = HTML("&mu;<sub>C</sub>: <p>The prior mean (ranges from 0.01 to 0.99) for treatment C.</p>"),0.25,
                                                                                           min = 0.01, max = 0.99, step = 0.01,width = NULL)
                                                                            )
                                                                     ),
                                                                     column(width=6,
                                                                            h4("Treatment Prior Sample Size Settings*"),
                                                                            tags$div(
                                                                              
                                                                              numericInput("nA",label = HTML("n<sub>A</sub>: <p>The prior sample size (larger than 0) for treatment A.</p>"),2,
                                                                                           min = 1, max = 999, step = 1,width = NULL),
                                                                              
                                                                              numericInput("nB",label = HTML("n<sub>B</sub>: <p>The prior sample size (larger than 0) for treatment B.</p>"),2,
                                                                                           min = 1, max = 999, step = 1,width = NULL),
                                                                              
                                                                              numericInput("nC",label = HTML("n<sub>C</sub>: <p>The prior sample size (larger than 0) for treatment C.</p>"),2,
                                                                                           min = 1, max = 999, step = 1,width = NULL),
                                                                              h5("*The larger value of prior sample size, the higher degree of belief for prior means.")
                                                                            )
                                                                     )
                                                                   )
                                                                 )
                                                                 ,style="")
                            )
    ))
  ),
  hr(),
    
    fluidRow(
      
      # column allocation for widgets
      column(width=3,offset = 1,
             tags$div(
               br(),
               actionButton("goButton", "Go!")
             )
      ),
      column(width=3,offset = 1,
             tags$div(
               br(),
               actionButton("refreshButton", "Refresh")
             )
      ),
      column(width=3,offset = 1,
             tags$div(
               br(),
               downloadButton("downloadButton", "Download Results")
             )
      )
      
    ),
    hr(),
    
    fluidRow(
      
      # column allocation for widgets
      column(width=12,style='margin-top:-7px;',  tags$div(style='margin-top:0px;',
                                                             htmlOutput("myResult")))


    ),
    hr(),
    
        
    fluidRow(
      
      # column allocation for widgets
      column(width=2,style='margin-top:-15px;',h2("Reference:")),
      # column(width=2,h3("Description:")),
      suppressWarnings(column(width=10,
                              tags$ul(
                                tags$li(h4("Tamura, R. N., Krischer, J. P., Pagnoux, C., Micheletti, R., Grayson, P. C., Chen, Y.-F.,
and Merkel, P. A. (2016). A small n sequential multiple assignment randomized trial
design for use in rare disease research. Contemporary Clinical Trials 46, 48-51."),style=""),
                                tags$li(h4("Wei, B., Braun, T. M., Tamura, R. N., and Kidwell, K. M. (2018). A Bayesian analysis of small n Sequential Multiple Assignment Randomized Trials (snSMARTs). Statistics in
                                                                                                                      Medicine 37, 3723-3732."),style="")
                                                   )
                                                   ))
                                ),
    hr()
  )
)
