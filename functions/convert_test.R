convert.test <- function(conversion.model, src, tgt) {
  #Tests a LASSO model generated using convert.train()
  #conversion.model is a named list of glmnet cvfit models from convert.train()
  #src is an MxL dataframe with replicates on rows and probes/genes on columns
  #tgt is an MxN dataframe with replicates on rows and probes/genes on columns
  require(glmnet)
  
  #Get names of each g/st in tgt to test
  tgt.names <- colnames(tgt)
  N <- length(tgt.names)
  
  #Initialize returned error
  rep.names <- row.names(tgt)
  M <- length(rep.names)
  
  predicted <- convert.predict(conversion.model, src, tgt)
    
  #Determine the Pearson Correlation between the log tgt and predicted
  pearson <- rep(NA,M)
  for (m in 1:M) {
    #Handle NA and -Inf from log
    pred <- log(predicted[m,])
    targ <- log(tgt[m,])
    compat <- pred!=-Inf&!is.na(pred)&targ!=-Inf&!is.na(targ)
    
    #Pearson Correlation calculation    
    pearson[m] <- cor.test(log(predicted[m,compat]),log(tgt[m,compat]),alternative = "t",method="pearson")[["estimate"]]
  }
  
  #Return RMSE
  return(list("predicted" = predicted, "pearson" = pearson))
}

convert.test.long <- function(conversion.model, src, tgt) {
  #Tests a LASSO model generated using convert.train()
  #conversion.model is a named list of glmnet cvfit models from convert.train()
  #src is an MxL dataframe with replicates on rows and probes/genes on columns
  #tgt is an MxN dataframe with replicates on rows and probes/genes on columns
  require(glmnet)
  
  #Get names of each g/st in tgt to test
  tgt.names <- colnames(tgt)
  N <- length(tgt.names)
  
  #Initialize returned error
  rep.names <- row.names(tgt)
  M <- length(rep.names)
  
  
  predicted <- convert.predict(conversion.model, src, tgt)
  
  #Validation phase
  E <- predicted - tgt #Error from fit
  SE <- E^2 #Squared Error
  SSE <- colSums(SE) #Sum of the squared errors by g/st
  MSE <- SSE/M #Mean squared error by g/st
  RMSE <- sqrt(MSE) #Length N vector of RMSE values for each probe/gene
  
  #Determine RMSE Between Technical Replicates
  RMSE.intra <- matrix(data=NA,nrow=M,ncol=N,dimnames=list(rep.names, tgt.names))
  for (m in 1:M) {
    current.rep <- tgt[m,] #Get mth replicate
    current.tgt <- tgt[-m,] #Remove mth replicate from tgt
    E.intra  <- current.tgt - current.rep
    SE.intra <- E.intra^2
    SSE.intra <- colSums(SE.intra)
    MSE.intra <- SSE.intra/(M-1)
    RMSE.intra[m,] <- sqrt(MSE.intra) #Populate RMSE matrix
  }
  
  #Determine Wilcoxon Rank Sum p-value between RMSEs
  wRS <- wilcox.test(RMSE,c(RMSE.intra))[[3]]
  
  #Determine the Pearson Correlation between the tgt and predicted/intra
  pearson <- matrix(data = NA,nrow = M,ncol = M+1,dimnames = list(rep.names,c(rep.names, 'Pred')))
  for (m in 1:M) {
    row <- tgt[m,]
    for (p in 1:(M+1)) {
      if (p == M+1) {
        col <- predicted[m,]
      }
      else {
        col <- tgt[p,]
      }
      pearson[m,p] <- cor.test(log(row),log(col),alternative = "two.sided",method="pearson")[['estimate']]
    }
  }
  
  #Return RMSE
  return(list("predicted" = predicted, "RMSE" = RMSE, "RMSE.intra" = RMSE.intra,
              "wilcoxon" = wRS, "pearson" = pearson))
}
