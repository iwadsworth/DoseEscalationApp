
library(shiny)
library(ggplot2)
source("functions.R")
source("Startfunctions.R")

#############################################
# Define server logic required
shinyServer(function(input, output) {
  

  observe({
    if (input$inputdata == 0){                      ### maybe 'inputdata' instead of 'save'
      table1 <- read.csv("datanew.csv",header=TRUE)
      output$view = renderDataTable({data.frame(ID=table1[1],Dose=table1[2],Response=table1[3])},
                    options = list(iDisplayLength = 10))
    }else{}
  })
 
  observe({
    if (input$inputdata == 0){                      ### maybe 'inputdata' instead of 'save'
      Soutput1 <- Sfun1()        
      Soutput2 <- Sfun2(Soutput1)           
      output$DRPlot <- renderPlot(function(){
        print(Soutput2[14])
      }  , height = 500, width = 900)
      
      output$fDRPlot <- renderPlot(function(){
        print(Soutput2[17])
      }  , height = 500, width = 900)
    
      if(length(Soutput1$subjectID)==0){
        output$bar <- renderPlot(function(){}  , height = 400, width = 920, bg="transparent")
        output$respplot <- renderPlot(function(){}  , height = 400, width = 920, bg="transparent")
      }else{
        output$bar <- renderPlot(function(){
          print(Soutput2[15])
        }  , height = 400, width = 920, bg="transparent") ## bg= for transparent background
        
        output$respplot <- renderPlot(function(){
          plot(Soutput1$subjectID,Soutput1$Dosenew,pch=Soutput1$Responsenew*19,ylab="Dose given",xlab="Patient number",main="Plot of current doses given",ylim=c(0.1,0.3), xaxt = "n")
          axis(side = 1, at = seq(from = 1, to = max(Soutput1$subjectID), by = 1)) 
          abline(h=0.1, col="gray");abline(h=0.15, col="gray");abline(h=0.2, col="gray");abline(h=0.25, col="gray");abline(h=0.3, col="gray");abline(h=0.11, col="gray90")
          abline(h=0.12, col="gray90");abline(h=0.13, col="gray90");abline(h=0.14, col="gray90");abline(h=0.16, col="gray90");abline(h=0.17, col="gray90")
          abline(h=0.18, col="gray90");abline(h=0.19, col="gray90");abline(h=0.21, col="gray90");abline(h=0.22, col="gray90");abline(h=0.23, col="gray90")
          abline(h=0.24, col="gray90");abline(h=0.26, col="gray90");abline(h=0.27, col="gray90");abline(h=0.28, col="gray90");abline(h=0.29, col="gray90")
          legend("topleft",c("No","Yes"),pch=c(0,19), title = "Successful outcome")
        }  , height = 400, width = 920)
      }

      output$TEXT1 = renderUI({Soutput2[1]})
      output$TEXT2 = renderUI({Soutput2[2]})
      output$TEXT3 = renderUI({Soutput2[3]})
      output$TEXT4 = renderUI({Soutput2[4]})
      output$FINAL = renderUI({Soutput2[5]})
      output$FINAL1 = renderUI({Soutput2[6]})
      output$FINAL2 = renderUI({Soutput2[7]})
      output$FINAL3 = renderUI({Soutput2[8]})
      output$FINAL4 = renderUI({Soutput2[9]})
      output$FINAL5 = renderUI({Soutput2[10]})
      output$FINAL6 = renderUI({Soutput2[11]})
      output$FINAL7 = renderUI({Soutput2[12]})
      output$FINAL8 = renderUI({Soutput2[13]})
      output$TEXT5 = renderUI({Soutput2[16]})
    }else{
         
    isolate({
      
      output1 <- reactive({
        fun1(input$dadd, input$radd, input$date, input$time)  
      })
      
      output2 <- reactive({
        fun2(output1())
      })
      
      
      output$DRPlot <- reactivePlot(function(){
        print(output2()[14])
      }  , height = 500, width = 900)
      
      output$fDRPlot <- renderPlot(function(){
        print(output2()[17])
      }  , height = 500, width = 900)
      
      output$bar <- reactivePlot(function(){
        print(output2()[15])
      }  , height = 400, width = 920, bg="transparent") ## bg= for transparent background
      
      output$respplot <- reactivePlot(function(){
        plot(output1()$subjectID,output1()$Dosenew,pch=output1()$Responsenew*19,ylab="Dose given",xlab="Patient number",main="Plot of current doses given",ylim=c(0.1,0.3), xaxt = "n")
axis(side = 1, at = seq(from = 1, to = max(output1()$subjectID), by = 1)) 
abline(h=0.1, col="gray");abline(h=0.15, col="gray");abline(h=0.2, col="gray");abline(h=0.25, col="gray");abline(h=0.3, col="gray");abline(h=0.11, col="gray90")
abline(h=0.12, col="gray90");abline(h=0.13, col="gray90");abline(h=0.14, col="gray90");abline(h=0.16, col="gray90");abline(h=0.17, col="gray90")
abline(h=0.18, col="gray90");abline(h=0.19, col="gray90");abline(h=0.21, col="gray90");abline(h=0.22, col="gray90");abline(h=0.23, col="gray90")
abline(h=0.24, col="gray90");abline(h=0.26, col="gray90");abline(h=0.27, col="gray90");abline(h=0.28, col="gray90");abline(h=0.29, col="gray90")
        legend("topleft",c("No","Yes"),pch=c(0,19), title = "Successful outcome")
      }  , height = 400, width = 920)
      
      output$TEXT1 = renderUI({output2()[1]})
      output$TEXT2 = renderUI({output2()[2]})
      output$TEXT3 = renderUI({output2()[3]})
      output$TEXT4 = renderUI({output2()[4]})
      output$FINAL = renderUI({output2()[5]})
      output$FINAL1 = renderUI({output2()[6]})
      output$FINAL2 = renderUI({output2()[7]})
      output$FINAL3 = renderUI({output2()[8]})
      output$FINAL4 = renderUI({output2()[9]})
      output$FINAL5 = renderUI({output2()[10]})
      output$FINAL6 = renderUI({output2()[11]})
      output$FINAL7 = renderUI({output2()[12]})
      output$FINAL8 = renderUI({output2()[13]})
      output$TEXT5 = renderUI({output2()[16]})
      
      
      output$view = renderDataTable({data.frame(ID=output1()[29],Dose=output1()[30],Response=output1()[31])},options = list(iDisplayLength = 10))
      
          
    })
  }
  })
   
  
})




