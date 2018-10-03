rm(list=ls())
#### Library to install if not present
suppressMessages(library(Hmisc))
#### Function to compute the weighted average and the uncertainty measures
source('weighted_ensemble_method.R')
#### Function to compute the distribution of local identifiability
source('uncertainty.r')
############################################# Data explanation
#################### path to data and results
wd.name = getwd()
#args = commandArgs(trailingOnly=TRUE)
#### If run from terminal
#ModelName = args[1]
#### if run from R studio then uncomment one of the following:
#ModelName = 'bacteria_ts_bermuda2011'
ModelName = 'bacteria_ts_hawaii'
data.data = paste('data_from_', ModelName, '.RData', sep = '')
### file to save
results.data = paste('results_from_', ModelName, '.RData', sep = '')
################################# Compute the weighted average and save the outputfile
weighted.average(results.data, data.data, methods.pesi = 'R2')
################################# Now look at the confidence intervals and the coefficient of variation
dist.of.coef.var =practical.identifiability(results.data, plt = TRUE, save.plt = F)
nomi.interactions = unlist(lapply(1:length(dist.of.coef.var), function(x,X)
  names(X[[x]][1]), dist.of.coef.var))
################################# Now plot the identifiability of coefficients
### Distribution:
dist.of.coef.var = (do.call(cbind, dist.of.coef.var))
boxplot(dist.of.coef.var, ylim = c(0,3.), ylab = 'distribution of coefficient of variation', 
        names = nomi.interactions)

colnames(dist.of.coef.var) = nomi.interactions
rownames(dist.of.coef.var) = NULL
### Temporal evolution
dist.of.coef.var[is.na(dist.of.coef.var)] = 0
colori = c('red', 'orange', 'blue', 'green')
plot(dist.of.coef.var[,1], type = 'l', lwd = 2, col = colori[1],
     xlab = 'Time', ylab = 'coefficient of variation', ylim = c(0,max(dist.of.coef.var)),
     main = 'Uncertainty in time')
for(k in 2:ncol(dist.of.coef.var)){
  lines(dist.of.coef.var[,k], lwd = '2', col = colori[k])
}
