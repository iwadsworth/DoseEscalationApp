fun1 <- function(dose,resp,DATE,TIME){
 
  
  ## load dataset
  data = read.csv("datanew.csv",header=TRUE)
  write.csv(data,"dataprev.csv",row.names=F)
  
  datanew <- data
  
  if(length(data[,1])==0){
    datanew[length(data[,1])+1,1] <- 1
  }else{
    datanew[length(data[,1])+1,1] <- datanew[length(data[,1]),1] + 1
  }
  
  datanew[length(data[,1])+1,2] <- dose
  datanew[length(data[,1])+1,3] <- resp
  
  write.csv(datanew,"datanew.csv",row.names=F)
  write.csv(datanew,paste("data",DATE,"-",TIME,".csv",sep=''),row.names=F)
  
  ## define dosing set
  doseSet = seq(0.1, 0.3, by=0.01)
  logdoseSet = log(doseSet)
  doseSum = vector(mode="numeric", length(doseSet))
  TD = 0.8
  cstop = 1.3
  
  ##########################***************************########################################
  # Below holds the pseudo data prior information
  ##########################***************************########################################
  ## w1 contains number of patients, doses contains the dose given, outcome contains #successes
  w1 <- c(3,3,rep(1,times=length(datanew[,1])))
  doses <- c(0.13,0.18,datanew[,2])
  outcome <- c(1.5,2.4,datanew[,3])
  
  ##########################***************************########################################
  
  index <- NULL
  for(i in 1:length(w1)){
    index[i] = which(mapply(function(x, y) {isTRUE(all.equal(x, y))}, doseSet, doses[i]))
  }
  
  
  ## Find modal estimate of dose-response model parameters
  FIT <- glm(outcome/w1~log(doses),family=binomial(link="logit"),weights=w1)
  SF <- summary(FIT)
  theta1 <- coef(SF)[1,1]
  theta2 <- coef(SF)[2,1]
  se1 <- coef(SF)[1,2]
  se2 <- coef(SF)[2,2]
  logdstar <- (log(TD/(1-TD)) - theta1)/theta2
  dstar = exp(logdstar)
  
  t1CIlow95  = theta1-qnorm(0.975)*se1
  t1CIupp95 = theta1+qnorm(0.975)*se1
  t2CIlow95  = theta2-qnorm(0.975)*se2
  t2CIupp95 = theta2+qnorm(0.975)*se2
  
  ndose <- vector(mode="numeric", length = length(doseSet))
  for(j in 1:length(index)){
    ndose[index[j]] = ndose[index[j]] + w1[j]
  }
  
  pj <- 1.0/(1.0 + exp(-(theta1 + theta2*logdoseSet))) 
  
  ## Calculate CI for ED80
  asympVar1 = sum(ndose*pj*(1-pj)*((logdoseSet - logdstar)^2))
  sum1 = sum(ndose*pj*(1-pj))
  gbar = sum(ndose*pj*(1-pj)*logdoseSet)/sum(ndose*pj*(1-pj))
  sum2 = sum(ndose*pj*(1-pj)*((logdoseSet - gbar)^2))
  asympVar = asympVar1/((theta2^2)*sum1*sum2)
  
  if(theta2 <= 0 & length(w1) < 62){
    CIlow95  = 1
    CIupp95 = 2.3
  }else{
    CIlow95  = dstar*exp(-qnorm(0.975)*sqrt(asympVar))
    CIupp95 = dstar*exp(qnorm(0.975)*sqrt(asympVar))
  }
  
  CIlow90  = dstar*exp(-qnorm(0.95)*sqrt(asympVar))
  CIupp90 = dstar*exp(qnorm(0.95)*sqrt(asympVar))
  
  CIlow50  = dstar*exp(-qnorm(0.75)*sqrt(asympVar))
  CIupp50 = dstar*exp(qnorm(0.75)*sqrt(asympVar))
  
  gain <- abs(pj - TD)
  d <- which.min(gain)
  RecDose <- doseSet[d]
  
  
  #############################################
  #For the plots
  
  ## Calculate CI for ED100pi
  dpi <- seq(0.1,0.3,by=0.001)
  logdpi <- log(dpi)
  pi1 <- exp(theta2*logdpi + theta1)/(1+exp(theta2*logdpi + theta1))
  
  asympVarpi1 <- NULL
  
  sum1 = sum(ndose*pj*(1-pj))
  gbar = sum(ndose*pj*(1-pj)*logdoseSet)/sum(ndose*pj*(1-pj))
  sum2 = sum(ndose*pj*(1-pj)*((logdoseSet - gbar)^2))
  
  for(k in 1:length(pi1)){
    asympVarpi1[k] = sum(ndose*pj*(1-pj)*((logdoseSet - logdpi[k])^2))
  }
  asympVarpi = asympVarpi1/((theta2^2)*sum1*sum2) #vector containing the asymptotic variance for each pi
  
  #95% CI
  CIlowpi95  = dpi*exp(-qnorm(0.975)*sqrt(asympVarpi))  #CI for dose
  CIupppi95 = dpi*exp(qnorm(0.975)*sqrt(asympVarpi))
  PMElowpi95 <- 1.0/(1.0 + exp(-(theta1 + theta2*log(CIlowpi95)))) #Corresponding prob of response
  PMEupppi95 <- 1.0/(1.0 + exp(-(theta1 + theta2*log(CIupppi95))))
  
  #50% CI
  CIlowpi50  = dpi*exp(-qnorm(0.75)*sqrt(asympVarpi))
  CIupppi50 = dpi*exp(qnorm(0.75)*sqrt(asympVarpi))
  PMElowpi50 <- 1.0/(1.0 + exp(-(theta1 + theta2*log(CIlowpi50)))) 
  PMEupppi50 <- 1.0/(1.0 + exp(-(theta1 + theta2*log(CIupppi50))))
  
  
  #If failure occurs
  
  dosesgiven <- as.numeric(names(table(datanew[,2])))  #vector of the dose level which have been given
  num_doses <- table(datanew[,2])[] #vector containing the number of times each dose has been given
  
  num_out <- vector(mode="numeric",length(dosesgiven)) #vector to contain the number of responses each dose level has had
  for(i in 1:length(dosesgiven)){
    num_out[i] <- sum((datanew[,3])[which(mapply(function(x, y) {isTRUE(all.equal(x, y))}, datanew[,2], dosesgiven[i]))])
  }
  prob_outcome <- num_out/num_doses  #probability of a response for each dose level given
  
  #############################################
  # For the final frequentist analysis
  
  ## Find modal estimate of dose-response model parameters
  fFIT <- glm(datanew[,3]~log(datanew[,2]),family=binomial(link="logit"))
  fSF <- summary(fFIT)
  
  if(fFIT$converged==FALSE){
    stop.fFit.Fail <- 1}else{stop.fFit.Fail <- 0
    }
  
  if( length(coef(fSF)[,1])==2 ){
    ftheta1 <- coef(fSF)[1,1]
    ftheta2 <- coef(fSF)[2,1]
    fse1 <- coef(fSF)[1,2]
    fse2 <- coef(fSF)[2,2]
    
    Ft1CIlow95  = ftheta1-qnorm(0.975)*fse1
    Ft1CIupp95 = ftheta1+qnorm(0.975)*fse1
    Ft2CIlow95  = ftheta2-qnorm(0.975)*fse2
    Ft2CIupp95 = ftheta2+qnorm(0.975)*fse2
    
    flogdstar <- (log(TD/(1-TD)) - ftheta1)/ftheta2
    fdstar = exp(flogdstar)
    
    ni <- num_doses
    ti <- num_out
    di <- dosesgiven
    expthet <- exp(ftheta1+ftheta2*log(di))
    negexpthet <- exp(-ftheta1-ftheta2*log(di))
    
    #Second derivative with respect to theta1
    d2ldtheta1 <- sum(ti*(expthet/((expthet+1)^2)) + 
                        (ni-ti)*(negexpthet/((negexpthet+1)^2)))
    
    #Second derivative with respect to theta2
    d2ldtheta2 <- sum(ti*(((log(di)^2)*expthet)/((expthet+1)^2)) +
                        (ni-ti)*(((log(di)^2)*negexpthet)/((negexpthet+1)^2)))
    
    #Second derivative with respect to theta1 and theta2
    d2ldtheta12 <- sum(ti*((log(di)*expthet)/((expthet+1)^2)) +
                         (ni-ti)*((log(di)*negexpthet)/((negexpthet+1)^2)))
    #Information matrix
    Io <- matrix(data=c(d2ldtheta1,d2ldtheta12,d2ldtheta12,d2ldtheta2),nrow=2,ncol=2)
    
    #Derivatives of the logED100pi function
    dgtheta1 <- -1/ftheta2                                ###*exp((log(4)-ftheta1)/ftheta2) for ED100pi
    dgtheta2 <- -(log(0.8/(1-0.8))-ftheta1)/(ftheta2^2)   ###*exp((log(4)-ftheta1)/ftheta2)
    
    del_gtheta <- c(dgtheta1,dgtheta2)
    vargtheta <- t(del_gtheta) %*% solve(Io) %*% (del_gtheta)
    segtheta <- sqrt(vargtheta)
    
    ## Calculate CI for ED80  
    FlogCIlow95  = flogdstar - qnorm(0.975)*segtheta
    FlogCIupp95 = flogdstar + qnorm(0.975)*segtheta
    FCIlow95 = exp(FlogCIlow95)
    FCIupp95 = exp(FlogCIupp95)
    
    checkFCIlow95  = fdstar*exp(-qnorm(0.975)*segtheta)
    checkFCIupp95 = fdstar*exp(qnorm(0.975)*segtheta)
    
    #For final frequentist plot
    
    fpj <- 1.0/(1.0 + exp(-(ftheta1 + ftheta2*logdoseSet))) #probs of response for the strict set of doses
    
    
    fpi1 <- exp(ftheta2*logdpi + ftheta1)/(1+exp(ftheta2*logdpi + ftheta1)) #probs of response for a range of possible doses
    
    #Derivatives of the logED100pi function
    fdgtheta1 <- -1/ftheta2                
    fdgtheta2 <- -(log(fpi1/(1-fpi1))-ftheta1)/(ftheta2^2)
    
    fdel_gtheta <- matrix(ncol=2,nrow=length(fpi1))
    for(i in 1:length(fpi1)){
      fdel_gtheta[i,1] <- fdgtheta1
      fdel_gtheta[i,2] <- fdgtheta2[i]  
    }
    
    fvargtheta <- NULL
    for(i in 1:length(fpi1)){
      fvargtheta[i] <- t(fdel_gtheta[i,]) %*% solve(Io) %*% (fdel_gtheta[i,])
    }
    
    fsegtheta <- sqrt(fvargtheta) #asyptotic SE for a range of possible ED100pi's
    
    
    ## Calculate 95% CI for ED100fpi  
    fFlogdCIlow95  = logdpi - qnorm(0.975)*fsegtheta
    fFlogdCIupp95 = logdpi + qnorm(0.975)*fsegtheta
    fFdCIlow95 = exp(fFlogdCIlow95)
    fFdCIupp95 = exp(fFlogdCIupp95)
    fFCIlow95 <- 1.0/(1.0 + exp(-(ftheta1 + ftheta2*log(fFdCIlow95)))) #Corresponding prob of response
    fFCIupp95 <- 1.0/(1.0 + exp(-(ftheta1 + ftheta2*log(fFdCIupp95))))
    
    
    #50% CI
    fFlogdCIlow50  = logdpi - qnorm(0.75)*fsegtheta
    fFlogdCIupp50 = logdpi + qnorm(0.75)*fsegtheta
    fFdCIlow50 = exp(fFlogdCIlow50)
    fFdCIupp50 = exp(fFlogdCIupp50)
    fFCIlow50 <- 1.0/(1.0 + exp(-(ftheta1 + ftheta2*log(fFdCIlow50)))) #Corresponding prob of response
    fFCIupp50 <- 1.0/(1.0 + exp(-(ftheta1 + ftheta2*log(fFdCIupp50))))
    
    
  }else{
    ftheta1 <- NULL
    ftheta2 <- NULL
    fse1 <- NULL
    fse2 <- NULL
    Ft1CIlow95  = NULL
    Ft1CIupp95 = NULL
    Ft2CIlow95  = NULL
    Ft2CIupp95 = NULL
    flogdstar <- NULL
    fdstar = NULL
    ni <- NULL
    ti <- NULL
    di <- NULL
    expthet <- NULL
    negexpthet <- NULL
    d2ldtheta1 <- NULL
    d2ldtheta2 <- NULL
    d2ldtheta12 <- NULL
    Io <- NULL
    dgtheta1 <- NULL
    dgtheta2 <- NULL
    del_gtheta <- NULL
    vargtheta <- NULL
    segtheta <- NULL
    FlogCIlow95  = NULL
    FlogCIupp95 = NULL
    FCIlow95 = NULL
    FCIupp95 = NULL
    checkFCIlow95  = NULL
    checkFCIupp95 = NULL
    fpj <- NULL
    fpi1 <- NULL
    fdgtheta1 <- NULL                
    fdgtheta2 <- NULL
    fdel_gtheta <- NULL
    fvargtheta <- NULL
    fsegtheta <- NULL
    fFlogdCIlow95  = NULL
    fFlogdCIupp95 = NULL
    fFdCIlow95 = NULL
    fFdCIupp95 = NULL
    fFCIlow95 <- NULL
    fFCIupp95 <- NULL
    fFlogdCIlow50  = NULL
    fFlogdCIupp50 = NULL
    fFdCIlow50 = NULL
    fFdCIupp50 = NULL
    fFCIlow50 <- NULL
    fFCIupp50 <- NULL
  }

  
  
  
  #############################################
  #Statements for the displayed text
  
  if(CIupp95/CIlow95 <= cstop){
    stop.Accuracy.Reached <- 1}else{stop.Accuracy.Reached <- 0
    }
  
  if(theta2 <= 0){
    stop.PME.Fail <- 1}else{stop.PME.Fail <- 0
    }
  
  if((sum(w1)-6)>=60){
    stop.Max.Subj <- 1}else{stop.Max.Subj <- 0
    }

  if(FIT$converged==FALSE){
    stop.Fit.Fail <- 1}else{stop.Fit.Fail <- 0
    }
  

  
  subjectID <- datanew[,1]
  Dosenew <- datanew[,2]
  Responsenew <- datanew[,3]
  
  
  
result1 <- list(stop.Accuracy.Reached=stop.Accuracy.Reached,  #1
                stop.PME.Fail=stop.PME.Fail,                  #2
                stop.Max.Subj=stop.Max.Subj,                  #3
                dosesgiven=dosesgiven,                        #4
                RecDose=RecDose,                              #5
                doses=doses,                                  #6
                dstar=dstar,                                  #7
                CIlow90=CIlow90,                              #8
                CIupp90=CIupp90,                              #9
                CIlow50=CIlow50,                              #10
                
                CIupp50=CIupp50,                              #11
                fdstar=fdstar,                                #12
                FCIlow95=FCIlow95,                            #13
                FCIupp95=FCIupp95,                            #14
                ftheta1=ftheta1,                              #15
                ftheta2=ftheta2,                              #16
                Ft2CIlow95=Ft2CIlow95,                        #17
                Ft2CIupp95=Ft2CIupp95,                        #18
                Ft1CIlow95=Ft1CIlow95,                        #19
                Ft1CIupp95=Ft1CIupp95,                        #20
                
                CIlow95=CIlow95,                              #21
                CIupp95=CIupp95,                              #22
                theta2=theta2,                                #23
                theta1=theta1,                                #24
                t2CIlow95=t2CIlow95,                          #25
                t2CIupp95=t2CIupp95,                          #26
                t1CIlow95=t1CIlow95,                          #27
                t1CIupp95=t1CIupp95,                          #28
                subjectID=subjectID,                          #29
                Dosenew=Dosenew,                              #30
                
                Responsenew=Responsenew,                      #31
                prob_outcome=prob_outcome,                    #32
                doseSet=doseSet,                              #33
                pj=pj,                                        #34
                PMElowpi95=PMElowpi95,                        #35
                PMEupppi95=PMEupppi95,                        #36
                PMElowpi50=PMElowpi50,                        #37
                PMEupppi50=PMEupppi50,                        #38
                dpi=dpi,                                      #39
                stop.Fit.Fail=stop.Fit.Fail,                  #40
                
                stop.fFit.Fail=stop.fFit.Fail,                #41
                fpj=fpj,                                      #42
                fFCIlow95=fFCIlow95,                          #43
                fFCIupp95=fFCIupp95,                          #44
                fFCIlow50=fFCIlow50,                          #45
                fFCIupp50=fFCIupp50)                          #46
return(result1)
}

#output1 <- fun1(0.2,1,12,12)


fun2 <- function(output1){
  if(output1[1]==0 & output1[2]==0 & output1[3]==0 & output1[40]==0){
    text1 <- paste("The current recommended dose for the next patient is:",output1[5],"mcg/kg/min.",sep=' ')
  }else if(output1[1]==0 & output1[2]==1 & output1[3]==0 & output1[40]==0 ){
    text1 <- paste("Repeat the previous dose of:",(output1$doses)[length(output1$doses)],"mcg/kg/min.",sep=' ')
  }else if(output1[1]==0 & output1[2]==0  & output1[3]==0 & output1[40]==1){
    text1 <- paste("Repeat the previous dose of:",(output1$doses)[length(output1$doses)],"mcg/kg/min.",sep=' ')
  }else if(output1[1]==0  & output1[2]==0 & output1[3]==1 & output1[40]==0){
    text1 <- paste("Stop recruitment, the maximum number of patients has been reached.")
  }else if(output1[1]==1 & output1[2]==0 & output1[3]==0 & output1[40]==0){
    text1 <- paste("Stop recruitment, the required level of accuracy has been reached.")
  }else if(output1[1]==1 & output1[2]==1 & output1[3]==0 & output1[40]==0){
    text1 <- paste("Stop recruitment, the required level of accuracy has been reached. The slope of
                   the dose-response model is negative.")
  }else if(output1[1]==0 & output1[2]==1 & output1[3]==1 & output1[40]==0){
    text1 <- paste("Stop recruitment, the maximum number of patients has been reached. The slope of
                   the dose-response model is negative.")
  }else if(output1[1]==1 & output1[2]==0 & output1[3]==1 & output1[40]==0){
    text1 <- paste("Stop recruitment, the maximum number of patients has been reached and the
                   required level of accuracy has been reached.")
  }else if(output1[1]==1 & output1[2]==1 & output1[3]==1 & output1[40]==0){
    text1 <- paste("Stop recruitment, the maximum number of patients and the required level of accuracy
                   has been reached. The slope of the dose-response model is negative.")
  }else if(output1[1]==1 & output1[2]==0 & output1[3]==0 & output1[40]==1){
    text1 <- paste("Stop recruitment, the required level of accuracy has been reached. Model failed
                   to converge.")
  }else if(output1[1]==0 & output1[2]==0 & output1[3]==1 & output1[40]==1){
    text1 <- paste("Stop recruitment, the maximum number of patients has been reached. Model failed
                   to converge.")
  }else if(output1[1]==1 & output1[2]==0 & output1[3]==1 & output1[40]==1){
    text1 <- paste("Stop recruitment, the maximum number of patients and the required level of
                   accuracy has been reached. Model failed to converge.")
  }
  
  if(output1[1]==0 & output1[2]==1 & output1[3]==0 & output1[40]==0){
    text2 <- ""
    text3 <- ""
    text4 <- ""
  }else{
    text2 <- paste("Given prior opinion and the data from the last",length(output1$doses)-2,"patients, the current most likely value of the ED80 for remifentanil is",
                   round(output1$dstar,digits=3),"mcg/kg/min. Our uncertainty about the true value of the ED80 is summarised by the following credibility intervals:")
    text3 <- paste("* 95% Credibility Interval for the ED80: (",round(output1$CIlow95,digits=3),", ",round(output1$CIupp95,digits=3),"),",sep='')
    text4 <- paste("* 50% Credibility Interval for the ED80: (",round(output1$CIlow50,digits=3),", ",round(output1$CIupp50,digits=3),").",sep='')
  }
  
  if(output1[2]==1 & output1[3]==0 & output1[40]==0){
    text5 <- "Warning: The slope of the dose-response relationship is negative. As a result the previous
              dose should be repeated and the graph shown above is a plot of the doses given to
              patients vs. the ratio of number of responses on each dose to the number of times
              each dose was given."
  }else if(output1[2]==0 & output1[3]==0 & output1[40]==1){
    text5 <- "Warning: The dose-response model failed to converge. As a result the previous
              dose should be repeated and the graph shown above is a plot of the doses given to
              patients vs. the ratio of number of responses on each dose to the number of times
              each dose was given."
  }else if(output1[2]==1 & output1[3]==0 & output1[40]==1){
    text5 <- "Warning: The slope of the dose-response relationship is negative and the dose-response model
              failed to converge.  As a result the previous dose should be repeated and the graph
              shown above is a plot of the doses given to patients vs. the ratio of number of responses
              on each dose to the number of times each dose was given."
  }else if(output1[2]==1 & output1[3]==1 & output1[40]==0){
    text5 <- "Warning: The slope of the dose-response relationship is negative."
  }else if(output1[2]==0 & output1[3]==1 & output1[40]==1){
    text5 <- "Warning: The dose-response model failed to converge."
  }else if(output1[2]==1 & output1[3]==1 & output1[40]==1){
    text5 <- "Warning: The slope of the dose-response relationship is negative and the dose-response model
              failed to converge."
  }else{
    text5 <- ""
  }
  
  if(output1[1]==1 | output1[3]==1){
    final <- "Final analysis:"
    final1 <- paste("The primary frequentist analysis for the final ",length(output1$doses)-2," patients, gives the maximum likelihood estimates of the ED80, slope and intercept parameters to be:",sep='')
    final2 <- paste("* ED80 = ",round(output1$fdstar,digits=3)," mcg/kg/min, with 95% Confidence Interval: (",round(output1$FCIlow95,digits=3),", ",round(output1$FCIupp95,digits=3),").",sep='')
    final3 <- paste("* slope = ",round(output1$ftheta2,digits=3)," with 95% Confidence Interval: (",round(output1$Ft2CIlow95,digits=3),", ",round(output1$Ft2CIupp95,digits=3),").",sep='')
    final4 <- paste("* intercept = ",round(output1$ftheta1,digits=3)," with 95% Confidence Interval: (",round(output1$Ft1CIlow95,digits=3),", ",round(output1$Ft1CIupp95,digits=3),").",sep='')
    final5 <- paste("The secondary Bayesian analysis for the final ",length(output1$doses)-2," patients, gives the current most likely value of the ED80 for remifentanil, slope and intercept parameters to be:",sep='')
    final6 <- paste("* ED80 = ",round(output1$dstar,digits=3)," mcg/kg/min, with 95% Confidence Interval: (",round(output1$CIlow95,digits=3),", ",round(output1$CIupp95,digits=3),").",sep='')
    final7 <- paste("* slope = ",round(output1$theta2,digits=3)," with 95% Confidence Interval: (",round(output1$t2CIlow95,digits=3),", ",round(output1$t2CIupp95,digits=3),").",sep='')
    final8 <- paste("* intercept = ",round(output1$theta1,digits=3)," with 95% Confidence Interval: (",round(output1$t1CIlow95,digits=3),", ",round(output1$t1CIupp95,digits=3),").",sep='')
  }else{
    final <- ""
    final1 <- ""
    final2 <- ""
    final3 <- ""
    final4 <- ""
    final5 <- ""
    final6 <- ""
    final7 <- ""
    final8 <- ""
  }
    
  if((output1[2]==1 | output1[40]==1) & output1[3]==0){
    
    data.df <- data.frame(dosesgiven=output1$dosesgiven,prob_outcome=output1$prob_outcome)
    a <- ggplot(data.df, aes(x = dosesgiven, y = prob_outcome)) + geom_line() + xlab("Dose (mcg/kg/min)") + ylab("Probability of successful outcome") + ggtitle("Plot of predicted probability of successful outcome") + xlim(0.1, 0.3) + ylim(0, 1)
    
  }else{
    
    data.df <- data.frame(doseSet=output1$doseSet,pj=output1$pj)
    data.df2 <- data.frame(PMElowpi95=output1$PMElowpi95,PMEupppi95=output1$PMEupppi95,PMElowpi50=output1$PMElowpi50,PMEupppi50=output1$PMEupppi50,dpi=output1$dpi)
    data.df3 <- data.frame(RecDose=output1$RecDose)
    a <- ggplot(data.df, aes(x = doseSet, y = pj,colour="Probability of successful outcome for a given dose, p(d)")) + 
         geom_line() + 
         geom_line(data=data.df2, aes(x = dpi, y = PMElowpi95, colour= 'The 95% Credibility Interval for p(d)')) +  
         geom_line(data=data.df2, aes(x = dpi, y = PMEupppi95, colour= 'The 95% Credibility Interval for p(d)')) + 
         geom_line(data=data.df2, aes(x = dpi, y = PMElowpi50, colour= 'The 50% Credibility Interval for p(d)')) +
         geom_line(data=data.df2, aes(x = dpi, y = PMEupppi50, colour= 'The 50% Credibility Interval for p(d)')) + 
         geom_vline(data=data.df3, aes(xintercept=RecDose, linetype = "Recommended dose")) +
         xlab("Dose (mcg/kg/min)") + 
         ylab("Probability of successful outcome") + 
         ggtitle("Dose response curve") + 
         xlim(0.1, 0.3) + ylim(0, 1) +
         scale_colour_manual("", values=c("Probability of successful outcome for a given dose, p(d)"="black",
				"The 95% Credibility Interval for p(d)"="red","The 50% Credibility Interval for p(d)"= "green")) +
      scale_linetype_manual("",values=c("Recommended dose"="longdash")) +
      scale_x_continuous(breaks = seq(0.1, 0.3, by = 0.02),limits=c(0.1,0.3)) +
  	 scale_y_continuous(breaks = round(seq(0, 1, by = 0.1),1)) 
  }
  
  if(output1[1]==1 | output1[3]==1){
    data.df <- data.frame(doseSet=output1$doseSet,fpj=output1$fpj)
    data.df2 <- data.frame(fFCIlow95=output1$fFCIlow95,fFCIupp95=output1$fFCIupp95,fFCIlow50=output1$fFCIlow50,fFCIupp50=output1$fFCIupp50,dpi=output1$dpi)
    c <- ggplot(data.df, aes(x = doseSet, y = fpj,colour="Probability of successful outcome for a given dose, p(d)")) + 
      geom_line() + 
      geom_line(data=data.df2, aes(x = dpi, y = fFCIlow95, colour= 'The 95% Credibility Interval for p(d)')) +  
      geom_line(data=data.df2, aes(x = dpi, y = fFCIupp95, colour= 'The 95% Credibility Interval for p(d)')) + 
      geom_line(data=data.df2, aes(x = dpi, y = fFCIlow50, colour= 'The 50% Credibility Interval for p(d)')) +
      geom_line(data=data.df2, aes(x = dpi, y = fFCIupp50, colour= 'The 50% Credibility Interval for p(d)')) + 
      xlab("Dose (mcg/kg/min)") + 
      ylab("Probability of successful outcome") + 
      ggtitle("Frequentist Dose response curve") + 
      xlim(0.1, 0.3) + ylim(0, 1) +
      scale_colour_manual("", values=c("Probability of successful outcome for a given dose, p(d)"="black",
                                       "The 95% Credibility Interval for p(d)"="red","The 50% Credibility Interval for p(d)"= "green")) +
      scale_x_continuous(breaks = seq(0.1, 0.3, by = 0.02),limits=c(0.1,0.3)) +
      scale_y_continuous(breaks = round(seq(0, 1, by = 0.1),1)) 
  }else{
    c <- NULL
  }
  
  ### For the bar chart
  table2 <- cbind(output1$Responsenew[order(output1$Dosenew)],output1$Dosenew[order(output1$Dosenew)])
  counts <- table(table2[,1],table2[,2])
  
  b <- qplot(factor(table2[,2]), geom="bar", fill=factor(table2[,1]))+ xlab("Dose (mcg/kg/min)") + 
	ylab("Number of patients given dose") +
       	theme(text = element_text(size=15)) + 
	scale_fill_discrete(name="Successful \noutcome",breaks=c("0", "1"),labels=c("No", "Yes")) +
	ggtitle("Number of patients treated on each dose") + 
	theme(
    	panel.grid.minor = element_blank(), 
    	panel.grid.major = element_blank(),
   	panel.background = element_blank(),
    	plot.background = element_blank()
   	)					### + theme makes background transparent
  


  
  result2 <- list(text1=text1,      #1
                  text2=text2,      #2
                  text3=text3,      #3
                  text4=text4,      #4
                  final=final,      #5
                  final1=final1,    #6
                  final2=final2,    #7
                  final3=final3,    #8
                  final4=final4,    #9
                  final5=final5,    #10
                  final6=final6,    #11
                  final7=final7,    #12
                  final8=final8,    #13
                  a=a,              #14
                  b=b,              #15
                  text5=text5,      #16
                  c=c)              #17
  return(result2)
}
