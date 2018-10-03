practical.identifiability <- function(file.name, plt = TRUE, save.plt = FALSE){
  ##################################################################################################################################
  #### Function to plot the Jacobian coefficient with their confidence interval and to extract the coefficient of variation
  #### Input:
  #### file.name = the name of the data file save from the weighted.averate() function
  #### plt = TRUE if want to plot the jacobian coefficients
  #### Output:
  #### coefficient of variation in time
  ##################################################################################################################################
  load(file.name)
  if(file.exists(file.name)){
    load(file.name)
    coeff.var.dist = list()
    tmp = 1
    num.species = ncol(ensemble.coefficients[[1]])
    if(plt == TRUE){
      if(save.plt == TRUE){
        nome.del.file=paste(ModelName, '_interactions.pdf', sep = '')
        pdf(nome.del.file,width=6,height=6,paper='special')
      }
	     par(mfrow=c(num.species,num.species),  oma = c( 0, 0, 2, 0 ) )
    }
    for(l in 1:num.species){
      for(m in 1:num.species){
        x = l
        y = m
        coeff = confidence = coef.var = rep(0,length(ensemble.coefficients))
        for(i in 1:length(ensemble.coefficients)){
          coeff[i] = ensemble.coefficients[[i]][x,y]
          confidence[i] = ensemble.coefficients.se[[i]][x,y]
          coef.var[i] = ensemble.coefficients.cv[[i]][x,y]
        }
        element.name = paste('d(', species.names[x], ')/d(', 
              species.names[y], ')', sep = '')
        coeff.var.dist[[tmp]] = coef.var
        names(coeff.var.dist[[tmp]]) = element.name
        tmp = tmp + 1
        if(plt == TRUE){
          plot(coeff, type = 'l', lwd = 1, col = 'red', xlab = 'Time', ylab = 'interactions', 
               main = element.name)
                  polygon(c(1:length(coeff),length(coeff):1),
                  c(coeff+(confidence) ,rev(coeff-(confidence))),col='skyblue',border=NA)
 
        }
      }
    }
    if(plt == TRUE){
      title("Jacobian coefficients with confidence intervals", outer = TRUE, cex = 5)
      if(save.plt == TRUE){
          dev.off()
      }
      par(mfrow=c(1,1))
      }
  } else{
	cat('Problem in reading the file\n')
   }
  return(coeff.var.dist)
}
