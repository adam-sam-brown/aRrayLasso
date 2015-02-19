# aRrayLasso
aRrayLasso in R
Adam Brown
19 February 2015

Robust conversion between microarray platforms is needed to leverage the wide variety of microarray expression studies that have been conducted to date. Currently available conversion methods rely heavily on manufacturer annotations, which are often incomplete, or on direct alignment of probes from different platforms, which are computationally intensive. Here, we describe aRrayLasso, which uses the Lasso-penalized generalized linear model to model the relationships between individual probes in different probe sets. We have implemented aRrayLasso in a set of five open-source R functions that allow the user to acquire data from public sources such as GEO, train a set of Lasso models on that data, and directly map one microarray platform to another. aRrayLasso significantly predicts expression levels with higher fidelity than technical replicates of the same RNA pool, demonstrating its utility in the integration of data sets from different platforms.

aRrayLasso is composed of three R functions: 
1) convert.train, which trains an aRrayLasso model based on a source and target dataset containing experiments performed with different microarrays on the same samples
2) convert.predict, which predicts expression levels in a target microarray platform given an aRrayLasso model and a source experiment
3) convert.test, which predicts expression levels in a target platform from a source platform and calculates Pearson's product-moment correlation coefficients between the prediction and given target data.

aRrayLasso also comes with two functions to convert existing datasets into a format usable by aRrayLasso:
1. convert.eSet, which converts an existing ExpressionSet object to an aRrayLasso compatible matrix
2. convert.GEO, which takes in a GEO Accession and the desired platform and returns an aRrayLasso compatible matrix
