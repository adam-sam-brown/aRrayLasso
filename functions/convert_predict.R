convert.predict <- function(conversion.model, src, tgt) {
  #Tests a LASSO model generated using convert.train()
  #conversion.model is a named list of glmnet cvfit models from convert.train()
  #src is an MxL dataframe with replicates on rows and probes/genes on columns
  #tgt is an MxN dataframe with replicates on rows and probes/genes on columns
  #tgt is used to format the returned data matrix of predicted values (can be empty, but must have correct dim and dimnames)
  #subset is an optional character containing probes or genes from tgt that should be tested, default is NA
  require(glmnet)
  
  #Get names of each g/st in tgt to predict
  tgt.names <- colnames(tgt)
  N <- length(tgt.names)
  
  #Get names and number of replicates to predict
  rep.names <- row.names(src)
  M <- length(rep.names)
  
  #Initialize prediction matrix
  predicted <- matrix(data=NA,nrow=M,ncol=N,dimnames=list(rep.names, tgt.names))
  
  #Loop through each g/st in tgt, predict values for all replicates based on cvfit model conversion.model
  for (n in 1:N) {
    tgt.name.current <- tgt.names[n]
    model.current <- conversion.model[[tgt.name.current]]
    if (is.na(model.current)) {
      predicted.current <- rep(0,M)
    }
    else {
      predicted.current <- predict(model.current, src, s="lambda.min") #Predict using the nth model with lambda = lambda.min  
    }
    predicted[,n] <- predicted.current
  }
  
  #Return predicted values
  return(predicted)
}