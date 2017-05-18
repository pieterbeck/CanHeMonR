#' @title  Fitting Simple Polynomial Models
#' @description Fits a simple polynomial of the form y ~ scale*(x^rate)+interc, using nls,
#' designed to be fit to a time series of a canopy health index that increases(?) in case of anomaly
#' @param data Data frame with two variables: im.dates, the independent variable, and DN.c, the dependent one (e.g. a VI or reflectance)
#' @param scale. The starting value for the scale parameter in the nls estimation of the polynomial. It determines how high the values rise.
#' Default is 2.8.
#' @param rate. The starting value for the rate parameter in the nls estimation of the polynomial.
#' **rate.** determines how abrupt the polynomial increases.
#' **rate.** = 1 generates a linear increase. Default is 1.
#' @param interc. The starting value for the intercept parameter in the nls estimation of the polynomial. Default is 1.
#' @param plott Logical. Should results be plotted. Default is FALSE.
#' @param min.n.time. The minimum nr of unique time observations for the model to run
#' @return An nls model
#' @details The current formulation of the model, combined with a scaling of the tim ethat makes it start at 0,
#' forces the model through VI==interc at the first time observation
#' the rate determines how abrupt the increase/decrease is. Rate=1 generates a linear increase
#' Note that when the rate is 0, the model fails
#  Scale determines how high the values go
#' @note called by decline_trees_RS_analysis.r
#' @export
fit_simple_polynomial <- function(data, scale.=2.8, rate.=1, interc.=0,plott=F,min.n.time.=3){
  time. <- data$im.dates
  VI. <- data$DN.c


  #nls functions
  #http://www.r-bloggers.com/a-better-nls/
  #require(minpack.lm)

  #convert the time. variables (from date) to numeric, if needed
  if (!is.numeric(time.)){
    time. <- as.numeric(time.)
  }

  outp <- NULL
  #if you have enough temporal replicates
  if (length(unique(time.)) >= min.n.time.){
    #browser()
    if (plott){
      plot(time.,VI.,type="p")
    }

    #in the formula, exp(1)*(time.-min(time.))/(max(time.)-min(time.)),
    #scales the time variable to the interval [0,e]
    nls.mod <- try(minpack.lm::nlsLM( formula = VI. ~ scale*(exp(1)*(time.-min(time.))/(max(time.)-min(time.)))^rate+interc,
                          start=list(scale=scale.,rate=rate.,interc=interc.) ),
                   silent=T)

    if(class(nls.mod) != 'try-error'){
      #if the model converged
      outp <- nls.mod
      #the current formulation of the model, combined with a scaling of the time
      #that makes it start at 0, forces the model through VI==interc at the first time observation
      #the rate determines how abrupt the increase/decrease is
      #rate=1 is a linear increase
      #Note that when the rate is 0, the model fails

      #scale determines how high the values go

    }#if the model converged
  }#if there were enough observations

  return(outp)
  #confint(nls.mod) should be possible according to http://www.r-bloggers.com/a-better-nls/
  #but gave error
}

