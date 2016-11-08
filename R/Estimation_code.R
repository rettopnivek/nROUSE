#-----------------------------------------#
# R code for estimating nROUSE parameters #
#-----------------------------------------#

#' The log-likelihood function for the nROUSE model
#'
#' Calculates the sum of the log-likelihoods for the nROUSE 
#' model given a set of data. This function can be used in 
#' an optimization routine, like \code{optim} to obtain 
#' estimates of the nROUSE model parameters.
#'
#' @param par a vector of log-transformed estimates for 
#'   the subset of nROUSE parameters.
#' @param dat a matrix (or data frame) with a set of named 
#'   columns:
#' \describe{
#'   \item{TarDur}{The duration of the target flash in ms.}
#'   \item{MaskDur}{The duration of the mask in ms.}
#'   \item{PrimeDur}{The duration of the prime in ms.}
#'   \item{Type}{The input for the prime. Negative values 
#'     denote foil primes, positive values indicate 
#'     target primes.}
#'   \item{Y}{The number of correct responses per condition.}
#'   \item{N}{The total number of items per condition.}
#' }
#' @param mapping an index for which parameters of the nROUSE 
#'   model to estimate.
#' \enumerate{
#'   \item The feedback scalar.
#'   \item The noise multiplier.
#'   \item The constant leak current.
#'   \item The synaptic depletion rate.
#'   \item The replenishment rate.
#'   \item The inhibition constant.
#'   \item The noise multiplier.
#'   \item The temporal attention parameter.
#'   \item The visual integration rate.
#'   \item The orthographic integration rate.
#'   \item The semantic integration rate.
#' }
#' @param estimate a logical value. If false, returns 
#'   the predicted proportion correct per condition 
#'   instead of the sum of the log-likelihoods.
#' @return Either the sum of the log-likelihoods, or if 
#' \code{estimate} is false, the vector of predicted 
#' proportions for correct responses.
#' @examples
#' # Load in example data set
#' data('priming_ex')
#' # Select a single subject
#' d = priming_ex[ priming_ex$Subject == 1, ]
#' # Specify a set of log-transformed starting values
#' sv = log( c( N = .0302, I = .9844, Ta = 1 ) )
#' # Estimate the parameters using maximum likelihood
#' mle = optim( sv, nROUSE_logLik, dat = d, 
#'   control = list( fnscale = -1, maxit = 10000 ) )
#' # Print the parameter estimates
#' round( exp( mle$par ), 3 )
#' # Compare to default values:
#' # N  = 0.030
#' # I  = 0.9844
#' # Ta = 1.0
#' @export

nROUSE_logLik = function( par, dat, mapping = c(2,6,8), 
                          estimate = T ) {
  
  # Restrict parameters to be positive
  par = exp( par );
  
  # Define default parameters
  param = c(
    Fe = 0.25, 
    N = 0.0302, 
    L = 0.15, 
    D = 0.324, 
    R = 0.022, 
    I = 0.9844, 
    Th = 0.15, 
    Ta = 1, 
    SV = 0.0294, 
    SO = 0.0609, 
    SS = 0.015
  )
  
  # Define function to simulate nROUSE model for given set of conditions
  f = function(x) {
    presentations = c(x['PrimeDur'],x['TarDur'],x['MaskDur'],500);
    primeType = c(0,0);
    if (x['Type'] < 0) primeType[2] = abs( x['Type'] );
    if (x['Type'] > 0) primeType[1] = x['Type'];
    param[ mapping ] = par;
    out = simulate_nROUSE(presentations,primeType,param);
    return( out$Latencies[3] )
  }
  
  theta = apply( dat, 1, f )
  
  if (estimate) {
    logLik = sum( dbinom( dat[,'Y'], dat[,'N'], theta, log = T ) )
    if (is.na(logLik)) logLik = -Inf
    
    return( logLik )
  } else return( theta )
  
}
