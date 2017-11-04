#------------------------------------#
# Description for the nROUSE package #
#------------------------------------#

# Useful tools for package development
# library(devtools)
# library(roxygen2)

#' Functions for the nROUSE model
#'
#' @docType package
#' @name nROUSE-package
#' @aliases nROUSE
#' @author Kevin W Potter, \email{kevin.w.potter@@gmail.com}
#' 
#' @useDynLib nROUSE, .registration = TRUE
#' @importFrom Rcpp evalCpp sourceCpp
#' 
#' @description
#'
#' The \pkg{nROUSE} package provides a set of functions 
#' to simulate and estimate the nROUSE model (Huber & 
#' O'Reilly, 2003; Rieth & Huber, 2017). The nROUSE model 
#' is a neural network model that accounts for accuracy 
#' performance in experimental paradigms that involve 
#' perceptual identification with repetition priming. The 
#' model posits that manipulations of prime type and duration 
#' influence accuracy performance via the non-linear 
#' interaction of lingering perceptual activation and neural 
#' habituation.
#' 
#' @references
#' 
#' Huber, D. E., & O'Reilly, R. C. (2003). Persistence and 
#'   accommodation in short-term priming and other perceptual 
#'   paradigms: Temporal segregation through synaptic depression. 
#'   Cognitive Science, 27, 403-430.
#'   
#' Rieth, C. A., & Huber, D. E. (2017). Comparing different 
#'   kinds of words and word-word relations to test an habituation 
#'   model of priming. Cognitive Psychology, 95, 79-104.
#'   DOI: https://doi.org/10.1016/j.cogpsych.2017.04.002
#'
NULL
