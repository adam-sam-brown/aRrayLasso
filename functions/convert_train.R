convert.train <- function(src, tgt, folds=10) {
  #Trains a LASSO model (GLM with alpha = 0) for each gene/sequence tag in tgt
  #src is an MxL matrix with replicates on rows and probes/genes on columns
  #tgt is an MxN matrix with replicates on rows and probes/genes on columns
  #M, the number of replicates must be the same for both src and tgt, L/N, the number of probes/genes may vary
  #folds is an integer number corresponding to the number of folds for cross-validation, default is 10
  require(glmnet)
  
  #Get names
  tgt.names <- colnames(tgt)
  N <- length(tgt.names)
  
  #Initialize returned datastructure
  conversion.model <- list()
  #skip <- 0
  
  #Loop through each g/st in tgt, generate cvfit model for each based on src
  for (n in 1:N) {
    tgt.name.current <- tgt.names[n]
    if (all(tgt[,n] == 0)) {
      model.current <- NA
      #skip <- skip + 1
    }
    else {
      model.current <- cv.glmnet(src,as.vector(tgt[,n]),nfolds=folds)  
    }
    conversion.model[[tgt.name.current]] <- model.current
  }
  
  #Return completed cv models
  #print(skip)
  return(conversion.model)
}