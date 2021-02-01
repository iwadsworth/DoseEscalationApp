# DoseEscalationApp

An R Shiny app to run the adaptive Bayesian logistic regression dose escalation approach for the study found in Arnold et al. (under review).

The files 'ui.R' and 'server.R' contain the main Shiny app UI and server code. The files 'functions.R' and 'Startfunctions.R' contain functions called by the Shiny app files, used to run the dose escalation approach, read-in data and output data. The two csv files, 'datanew.csv' and 'dataprev.csv' are used by the app to read-in existing data and save newly added data - To begin with, these two files will be empty (other than headings) with new rows added wih each subject input.

Opening the 'ui.R' or 'server.R' files in RStudio and using the 'RunApp' buttun or calling the function 'runApp()' in R with the working directory set to the location of the app files will run the app. A step-by-step guide of the running of the app can be found in file 'Help document.pdf'.

Reference: Arnold et al. (under review)
