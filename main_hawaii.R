rm(list=ls(all=TRUE)) 
suppressMessages(library(Matrix))
suppressMessages(library(quantreg))
suppressMessages(library(parallel))
suppressMessages(library(compiler))
suppressMessages(library(lars))
suppressMessages(library(elasticnet))
suppressMessages(library(caret))
options(warn=-1)
#####################################################################################
source('Auxiliar.r')
source('elastic_net_fit.r')
source('LOOCV.r')
source('KernelFunctions.r')
source('OutOfSample.r')
source('TrainingError.r')
#####################################################################################
############## Test different prior distributions and look at the inferred parameters
#####################################################################################
ShowPlot = FALSE
lags = TRUE
ModelName = 'bacteria_ts_hawaii'
FileName = paste(ModelName, '.txt', sep = '')
#### get directory name
wd.name = getwd()
data.name = paste(wd.name, '/outputFiles/data_from_', ModelName, '.RData', sep = '')
###################################
logspace <- function(d1, d2, n) exp(log(10)*seq(d1, d2, length.out=n)) 
std_err <- function(x) sd(x)/sqrt(length(x))
############# Choose the kernel
Kernel.Options = c('Exponential.Kernel', 'Epanechnikov.Kernel', 'TriCubic.Kernel')
Regression.Kernel = Kernel.Options[3]
############# Parameters for cross validation
if(Regression.Kernel == 'Exponential.Kernel'){
  lambda = logspace(-3,0,15)                       
  tht = seq(from = 0., to = 6, length = 20)         
}else{
  lambda = logspace(-3,0,15)                       
  tht = seq(from = 0.1, to = 3, length = 20)     
}
parameters_on_grid = expand.grid(tht, lambda)     
### Read Time series
d = as.matrix(read.table(FileName, header= T))
######################
original.Embedding = c('pbact', 'ebact')
original.TargetList = original.Embedding
d = d[, original.Embedding]
#### Here you take combinations of lags
x.lag = 2; y.lag = 2; z.lag = 1
sp.lag.selection = c(x.lag, y.lag)
lagged.time.series = make.lagged.ts(d, sp.lag.selection)
add.temperature = FALSE
if(isTRUE(add.temperature)){
  s = as.matrix(read.table(FileName, header= T))
  d = cbind(lagged.time.series$time.series, s[(max(x.lag,y.lag)+1):nrow(s), 'temp'])
  cat('Temperature is also a predictor\n')
}else{
d = lagged.time.series$time.series
}
original.col = lagged.time.series$original.variables
if(lags == TRUE){ var.sel = original.col; }else{ var.sel = colnames(d)}
##### Names and embedding in the laged dataset
if(lags == TRUE){ colnames(d) = Embedding =  TargetList = LETTERS[1:ncol(d)]}else{
  Embedding =  TargetList = original.Embedding
}
##### length of training and test set
length.testing = 2
length.training = nrow(d) - length.testing
#### Preserve training for the interactions
ts.train.preserved = d[1:length.training, var.sel]
std.ts.train = Standardizza(ts.train.preserved)
#### Preserve testing for the test (you want your algorithm to learn the real structure of the model)
ts.test.preserved = d[(length.training + 1):nrow(d), var.sel]
#### Training set:
d.training = Standardizza(d[1:length.training, ])
#### You now need to standardize the test set using mean and sd of the training set
d.testing = Standardizza.test(ts.test.preserved,ts.train.preserved)
############## Prepare for parallel computing
Lavoratori = detectCores() - 1
cl <- makeCluster(Lavoratori, type = "FORK")
####
iterations = 40
jac.coefficients = list()
training.error.rho = test.error.rho = naive.pred.rmse = 
training.error.rmse = test.error.rmse = 
Explained.Variance.test = Explained.Variance.train =
rep(0,iterations)
RegressionType = 'ELNET_fit'
alpha = 0.6
cat('Species in the system:', original.Embedding, '\n')
##################################### Main loop #########################################
for(i in 1:iterations){
 	cat('alpha:', alpha, '\n')
  BestModel = BestModelLOOCV(cl, d.training, TargetList, Embedding, parameters_on_grid, RegressionType,alpha)
  BestCoefficients = BestModel$BestCoefficients
	BestParameters = BestModel$BestParameters
	### Forecast
	out.of.samp.ELNET = out_of_sample_sequence(cl, BestCoefficients, 
	                                           BestParameters$BestTH,
        	                                   BestParameters$BestLM, 
                	                           d.training, length.testing)
	prd = out.of.samp.ELNET$out_of_samp
	###############################
	###### Now take the coefficients of the unlagged variables
	jacobiano.inference = take.coeff(BestCoefficients, var.sel, original.Embedding)
	###### Now check the in-sample error of the unlagged variables
	TrainErr = ComputeTrainingError(d.training, BestCoefficients, var.sel)
	Reconstruction = ReconstructionOfTrainingSet(d.training, BestCoefficients)
	Reconstruction = Reconstruction[,var.sel]
  ##############################################################
	###### Now check the quality of the out-of-sample forecast
	###### of only the unlaged variables
	prd = prd[,var.sel]
	rmse.test = compute.rmse.test(d.testing, prd)
	rmse.test.naive = naive.forecast(std.ts.train[nrow(d.training),],d.testing)
	R2.training = as.numeric(postResample(pred = Reconstruction, 
	                                      obs = d.training[1:nrow(d.training)-1,var.sel])['Rsquared'])
  R2.testing = as.numeric(postResample(pred = prd, obs = d.testing)['Rsquared'])
	jac.coefficients[[i]] = jacobiano.inference
	training.error.rmse[i] = TrainErr$rmse
	test.error.rmse[i] = rmse.test
	naive.pred.rmse[i] = rmse.test.naive
	Explained.Variance.train[i] = R2.training
	Explained.Variance.test[i] = R2.testing
	alpha = alpha + 0.01
}
num.species = ncol(std.ts.train)
#### Now save the variables to pass to analyze
save(jac.coefficients, naive.pred.rmse, 
     Explained.Variance.train, Explained.Variance.test, 
     training.error.rmse, test.error.rmse, num.species, ts.train.preserved, file = data.name)
cat('... Done. Results saved in /outputFiles/Data\n')
###############################
if(ShowPlot == TRUE){
  source('PlotFunctions.r')
  Plot_out(d.testing, prd)
  Plot_AllSpeciesTraining(std.ts.train[1:nrow(d.training)-1,], Reconstruction)
  Plot_SingleSpeciesTraining(std.ts.train[1:nrow(d.training)-1,], 
                             Reconstruction, 1)
}
stopCluster(cl)
