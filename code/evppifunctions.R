
library(voi)

evppivar_npreg <- function (outputs, inputs, pars, method, verbose, ...) {
  fitted <- voi:::fitted_npreg_fn(method)(y = outputs, inputs = inputs, 
                                          pars = pars, verbose = verbose, ...)
  list(data.frame(evppi = var(fitted)),fitted)
}

evppivar <- function (outputs, inputs, pars = NULL, method = NULL, nsim = NULL, 
                      verbose = TRUE, ...) {
  inputs <- voi:::check_inputs(inputs, iname = deparse(substitute(inputs)))
  voi:::check_outputs_vector(outputs, inputs)
  if (is.list(pars)) {
    return(evppivar_list(outputs = outputs, inputs = inputs, 
                         pars = pars, method = method, nsim = nsim, verbose = verbose, 
                         ...))
  }
  pars <- voi:::check_pars(pars, inputs)
  opts <- list(...)
  if (is.null(method)) 
    method <- voi:::default_evppi_method(pars)
  if (is.null(nsim)) 
    nsim <- nrow(inputs)
  outputs <- outputs[1:nsim]
  inputs <- inputs[1:nsim, , drop = FALSE]
  if (method %in% voi:::npreg_methods) {
    rese <- evppivar_npreg(outputs = outputs, inputs = inputs, 
                           pars = pars, method = method, verbose = verbose, 
                           ...)
  }
  else stop("Other methods not implemented yet")
  res <- cbind(pars = paste(pars, collapse = ","), rese[[1]])
  res
}

evppifit <- function (outputs, inputs, pars = NULL, method = NULL, nsim = NULL, 
                      verbose = TRUE, ...) {
  inputs <- voi:::check_inputs(inputs, iname = deparse(substitute(inputs)))
  voi:::check_outputs_vector(outputs, inputs)
  if (is.list(pars)) {
    return(evppivar_list(outputs = outputs, inputs = inputs, 
                         pars = pars, method = method, nsim = nsim, verbose = verbose, 
                         ...))
  }
  pars <- voi:::check_pars(pars, inputs)
  opts <- list(...)
  if (is.null(method)) 
    method <- voi:::default_evppi_method(pars)
  if (is.null(nsim)) 
    nsim <- nrow(inputs)
  outputs <- outputs[1:nsim]
  inputs <- inputs[1:nsim, , drop = FALSE]
  if (method %in% voi:::npreg_methods) {
    rese <- evppivar_npreg(outputs = outputs, inputs = inputs, 
                           pars = pars, method = method, verbose = verbose, 
                           ...)
  }
  else stop("Other methods not implemented yet")
  res <- rese[[2]]
  res
}


get_voi = function(outcomes, parameterdf, groupings, scenlabs){
  ydim1 = ncol(outcomes[[1]])
  matchpar = ifelse(class(parameterdf)=='list',T,F)
  voi <- matrix(NA,length(outcomes)*ydim1, length(groupings))
  for(i in 1:length(outcomes)){
    thispar = if(matchpar) parameterdf[[names(outcomes)[i]]] else parameterdf
    for(k in 1:ydim1){
      y = outcomes[[i]][,k]
      vary <- var(y) 
      if(vary>0){
        for(j in 1:length(groupings)){
          sourcesj <- groupings[[j]]
          fittedvalues <- evppifit(y,thispar,pars=colnames(thispar)[sourcesj],method='earth')
          voi[(i-1)*ydim1+k,j] <- (vary - mean((y - fittedvalues) ^ 2)) / vary * 100
        }
      }
    }
  }
  voitab = reshape2::melt(voi)
  pl = ggplot(voitab) +
    geom_tile(aes(x=Var2,y=Var1,fill=value),colour='black',linewidth=0.05) +
    scale_fill_gradient(low = "azure1", high = "darkslategray", na.value = 'white') +
    scale_y_continuous(breaks=1:nrow(voi),labels=rep(scenlabs,length(outcomes)), expand=c(0,0),
                       sec.axis = sec_axis(~ ., breaks = ydim1*c(1:length(outcomes))-3, 
                                           labels = TeX(names(outcomes)),name='Outcomes')) +
    theme_bw(base_size=15) +
    theme(axis.text.y.right = element_text(angle=-90,hjust=.2),
          axis.text.x = element_text(angle=-45,hjust=0,vjust=0.5),
          axis.ticks = element_blank()) +
    scale_x_continuous(breaks=1:ncol(voi),labels=names(groupings), expand=c(0,0)) +
    labs(x='Parameters',y='Scenarios',fill='VOI')
  for(i in 1:length(outcomes))
    pl = pl + annotate('rect',xmin=.5,xmax=ncol(voi)+.5,ymin=ydim1*(i-1)+.5,ymax=ydim1*i+.5,fill=NA,colour='black',linewidth=1)
  return(list(pl, voitab))
}

