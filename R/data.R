#' Response time and choice data from a priming task.
#'
#' A dataset containing response times and choice data for
#' 25 subjects who completed a word identification task with
#' immediate word priming.
#'
#' @format A data frame with 100 rows and 8 variables:
#' \describe{
#'   \item{PrimeDur}{The duration in ms for the prime word.}
#'   \item{TarDur}{The duration in ms for the flash of 
#'     the target word}
#'   \item{MaskDur}{The duration of the backward mask for the 
#'     target word.}
#'   \item{Type}{The type of prime, where -2 indicates a 
#'     double prime of the foil word, and 2 indicates a double 
#'     prime of the target word.}
#'   \item{Subject}{The subject index.}
#'   \item{P}{The proportion correct per subject and condition.}
#'   \item{Y}{The frequency correct per subject and condition.}
#'   \item{N}{The total number of trials per subject and condition.}
#'   }
"priming_ex"
