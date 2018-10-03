

weighted.average <- function(file.name, data.name, methods.pesi = 'R2'){
	########################################################################################################################
	#### This function compute the uncertainty and the confidence interval to the coefficient of the S-map
	#### Input: 
	#### file.name = name of the file where to save the results
	#### data.name = name of the file where the ensemble of models is saved
	#### methods.pesi = one of the following: 'R2', 'uniform', 'all.models'
	#### Output:
	#### average Jacobian, standard error and coefficient of variation of its coefficients
	########################################################################################################################
	load(data.name)
	epsilon.round.rmse = 2
	### subset those models with error at maximum 10% bigger than the minimum
	idx.rmse.test = which(round(test.error.rmse,epsilon.round.rmse) < 
	                        1.1*min(round(test.error.rmse,epsilon.round.rmse))) 
	### subset those models with explained variance at minimum 98% of the maximum explained variance
	idx.r2.train = which(round(Explained.Variance.train,epsilon.round.rmse) > 
	                        0.98*max(round(Explained.Variance.train[idx.rmse.test],epsilon.round.rmse)))
	### Find the intersection between the best models in the training and test set
	idx.subset = intersect(idx.r2.train, idx.rmse.test)
	cat('Number of models in the ensemble:', length(idx.subset), '\n')
	####################################################################################################
	species.names = colnames(jac.coefficients[[1]][[1]])
	idx.to.use.for.subset = idx.subset
	#############################################
	##### Coefficients of the ensemble method with their standard error
	##### ensemble.coefficients will be the time series of Jacobian from the ensemble method
	##### ensemble.coefficients.se is the time series of their standard errors
	##### ensemble.coefficients.cv is the the series of coefficient of variations
	ensemble.coefficients = ensemble.coefficients.se =  ensemble.coefficients.cv = list()
	##### The weights are like probability of a model and must sum up to one
	##### Choose which method to use to weight the average
	if(methods.pesi == 'R2'){
	    pesi = Explained.Variance.train[idx.to.use.for.subset]/sum(Explained.Variance.train[idx.to.use.for.subset])
	} else if(methods.pesi == 'uniform'){
	    pesi = rep(1./length(idx.to.use.for.subset), length(idx.to.use.for.subset))
	} else if(methods.pesi == 'all.models'){
	    idx.to.use.for.subset = 1:length(test.error.rmse)
	    pesi = Explained.Variance.train[idx.to.use.for.subset]/sum(Explained.Variance.train[idx.to.use.for.subset])
	}
	### Compute ensemble coefficients, the standard errors and the coefficient of variations
	for(k in 1:length(jac.coefficients[[1]])){
	  ensemble.coefficients[[k]] = matrix(NA,num.species,num.species)
	  ensemble.coefficients.se[[k]] = matrix(NA,num.species,num.species)
	  ensemble.coefficients.cv[[k]] = matrix(NA,num.species,num.species)    
	  for(i in 1:num.species){
	    for(j in 1:num.species){
	      boosted.coeff = unlist(lapply(idx.to.use.for.subset, function(x,X,l,m,n) X[[x]][[n]][l,m], 
	                                    jac.coefficients,i,j,k))
	  #### Weigthed mean over models
	      ensemble.coefficients[[k]][i,j] = wtd.mean(boosted.coeff, pesi)
		#### Weighted standard error
	      ensemble.coefficients.se[[k]][i,j] = 
	      1.96*sqrt(wtd.var(boosted.coeff, weights = pesi, na.rm = T,normwt=T))/sqrt(length(idx.to.use.for.subset))
		#### Weighted coefficient of variation
		    ensemble.coefficients.cv[[k]][i,j] = 
	        sqrt(wtd.var(boosted.coeff, weights = pesi, na.rm = T,normwt=T))/abs(wtd.mean(boosted.coeff, weights = pesi))
	    }
	  }
	}
	save(species.names,
       ensemble.coefficients, 
	     ensemble.coefficients.se,
	     ensemble.coefficients.cv, file = file.name)

}


