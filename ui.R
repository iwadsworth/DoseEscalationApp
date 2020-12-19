library(shiny)

# Define UI
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Remifentanil study"
              ),
  
  sidebarPanel(
    h3("To input new patient data check the box below, input the new data, then press 'Update':"),
    br(),
    checkboxInput(inputId = "inputdata", 
                  label = strong("To input new patient data, click here."),
                  value = FALSE),
    br(),
    numericInput("dadd", "What dose was given to the new patient? (mcg/kg/min)","",value=""),
    numericInput("radd", "Did the patient have a successful outcome? (1=yes, 0=no)","",value=""),
    textInput("date", "Date of this dosing (yyyy.mm.dd)",value=""),
    textInput("time", "Time of this dosing (24 hour: hh.mm)",value="")  , 
    
    #conditionalPanel(condition = "input.inputdata == true"     
    #                ),
    
    submitButton("Update")
                ),
  
  
  mainPanel(
    
    tabsetPanel(
      tabPanel(title = "Observed data", 
               h3("Summary of the dataset"),
               dataTableOutput('view'), # Show a summary of the dataset
               #div(style="display:inline-block",plotOutput("bar")),        #div() to make plots side by side
               #div(style="display:inline-block",plotOutput("respplot"))
               plotOutput("bar"),
               plotOutput("respplot")
      ),
      tabPanel(title = "Dose recommendation and posterior summaries", 
                  p(h4(uiOutput('TEXT1'),style = "color:red")),
               br(),
               h5(uiOutput('TEXT2')),
               h5(uiOutput('TEXT3')),
               h5(uiOutput('TEXT4')),
               br(),
               plotOutput("DRPlot"),
               br(),
               br(),
               br(),
               br(),
               br(),
               plotOutput("fDRPlot"),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               h5(uiOutput('TEXT5'))
               ),
      tabPanel(title= "Final output",             
                  h4(uiOutput('FINAL')),
                  h5(uiOutput('FINAL1')),
                  h5(uiOutput('FINAL2')),
                  h5(uiOutput('FINAL3')),
                  h5(uiOutput('FINAL4')),
                  br(),
                  h5(uiOutput('FINAL5')),
                  h5(uiOutput('FINAL6')),
                  h5(uiOutput('FINAL7')),
                  h5(uiOutput('FINAL8'))
               )
              )   
    )
))
