convert.eSet <- function(input) {
  #Converts the typical expression set output of RMA to glmnet compatible matrices
  #Input is an object of the class expressionSet
  data <- exprs(input)
  samples <- sampleNames(input)
  features <- featureNames(input)
  compatibility.matrix <- matrix(data=data,nrow=length(samples),ncol=length(features),dimnames=list(samples,features))
  
  return(compatibility.matrix)
}

convert.GEO <- function(GSE, platform) {
  #Converts a GEO accession to an expression set
  #GSE is a GEO excession number
  #Platform is any platform in the GSE
  eSet.list <- getGEO(GSE, GSEMatrix=TRUE)
  plat.string <- paste(GSE,"-",platform,"_series_matrix.txt.gz", sep='')
  eSet <- eSet.list[[plat.string]]
  
  return(eSet)
}