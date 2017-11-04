// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]


/*
Purpose:
C++ code for fast simulation of the nROUSE model (Huber and O'Reilly, 
2003; Reith & Huber, 2017).

// Printing to R console
Rcpp::Rcout << "Debugging example" << std::endl;

Index
Lookup - 01:  subplus
Lookup - 02:  update_nROUSE
Lookup - 03:  simulate_nROUSE
*/

// Lookup - 01
// An internal function to calculate when membrane potentials 
// are above a firing threshold.
arma::mat subplus( arma::mat A, double Th ) {
  
  arma::mat B( A.n_rows, A.n_cols, arma::fill::ones );
  B = A - B*Th;
  
  B.elem( find( B <= 0.0 ) ).zeros();
  
  return( B );
}

// Lookup - 02
// An internal function that updates the membrane potential, amplitude, 
// and output for a specified layer in the nROUSE network.
void update_nROUSE( arma::mat& new_mem, arma::mat& new_amp, arma::mat& old_out, 
                    arma::mat inp, arma::mat old_mem, arma::mat old_amp, 
                    double Th, double L, double R, double D, double I, 
                    Rcpp::NumericVector S, int level ) {
  
  // Output is above threshold activation multiplied by available 
  // synaptic resources
  double sums;
  
  old_out = subplus( old_mem, Th ) % old_amp;
  
  // Initialize matrix for inhibition
  arma::mat inhibit(1, old_out.n_cols, arma::fill::zeros);
  if ( level == 1 ) {
    // Visual inhibition
    
    // Prime, target, and mask mutually inhibit each other 
    // (centrally presented)
    sums = old_out(0,0) + old_out(0,1) + old_out(0,2);
    for (int i = 0; i < 3; i++) inhibit(0,i) = sums;
    // Choice words only self inhibit (unique screen location)
    for (int i = 3; i < 5; i++) inhibit(0,i) = old_out(0,i);
    
  } else {
    // For orthography and semantics, everything inhibits everything
    
    sums = accu( old_out );
    int Tot = old_out.n_cols;
    for ( int i = 0; i < Tot; i++ ) inhibit(0,i) = sums;
    
  }
  
  // Define matrix filled with ones
  arma::mat ones_const( old_mem.n_rows, old_mem.n_cols, arma::fill::ones );
  
  // Update membrane potential and synaptic resources
  new_mem = old_mem + S(level-1)*( inp % ( ones_const - old_mem ) - 
    L*old_mem - 
    I*(inhibit % old_mem) );
  new_amp = old_amp + S(level-1)*( R*(ones_const - old_amp) - D*old_out );
  
  // Keep things within bounds
  new_mem.elem( find( new_mem < 0.0 ) ).zeros();
  new_mem.elem( find( new_mem > 1.0 ) ).ones();
  new_amp.elem( find( new_amp < 0.0 ) ).zeros();
  new_amp.elem( find( new_amp > 1.0 ) ).ones();
  
}

//' The nROUSE model
//'
//' Function to simulate the predictions of the nROUSE model (Huber 
//' and O'Reilly, 2003; Reith & Huber, 2017) for perceptual 
//' identification latencies and forced-choice accuracy.
//'
//' @param presentations A vector giving the duration (in ms) of 
//'   the prime, the target, the mask, and the choice alternatives.
//' @param primeInput A vector of two elements used to set the type 
//'   of prime. For instance, [2,0] indicates a double prime for 
//'   targets, while [0,1] indicates a single prime for foils.
//' @param param A vector giving the values for the parameters.
//'  \describe{
//'    \item{\code{param[1]}}{Semantic to orthographic feedback scalar}
//'    \item{\code{param[2]}}{Noise constant}
//'    \item{\code{param[3]}}{Constant leak current}
//'    \item{\code{param[4]}}{Synaptic depletion rate}
//'    \item{\code{param[5]}}{Recovery rate}
//'    \item{\code{param[6]}}{Inhibition constant}
//'    \item{\code{param[7]}}{Activation threshold}
//'    \item{\code{param[8]}}{Temporal attention}
//'    \item{\code{param[9]}}{Integration time constant (Visual)}
//'    \item{\code{param[10]}}{Integration time constant (Orthographic)}
//'    \item{\code{param[10]}}{Integration time constant (Semantic)}
//'     }
//'
//' @section References:
//' Huber, D. E., & O'Reilly, R. C. (2003). Persistence and 
//'   accommodation in short-term priming and other perceptual 
//'   paradigms: Temporal segregation through synaptic depression. 
//'   Cognitive Science, 27(3), 403-430.
//'   
//' Rieth, C. A., & Huber, D. E. (2017). Comparing different 
//'   kinds of words and word-word relations to test an habituation 
//'   model of priming. Cognitive Psychology, 95, 79-104.
//'   DOI: https://doi.org/10.1016/j.cogpsych.2017.04.002
//'
//' @examples
//' # Define duration (in ms) for prime, target, mask, and choices
//' presentations = c( 17, 50, 450, 500 )
//' # Set a double prime for targets
//' primeInput = c(2,0)
//' # Simulate model with default parameters
//' sim = simulate_nROUSE( presentations, primeInput )
//'
//' @export
// [[Rcpp::export]]
Rcpp::List simulate_nROUSE( Rcpp::NumericVector presentations, 
                            Rcpp::NumericVector primeInput, 
                            Rcpp::NumericVector param = 
                              Rcpp::NumericVector::create(
                                .25, .0302, .15, .324, .022, 
                                .9844, .15, 1.0, .0294, .0609, 
                                .015 ) ) {
  // Duration of prime
  int PrimeDur = presentations(0);
  // Duration of target flash
  int TarDur = presentations(1);
  // Duration of target mask
  int MaskDur = presentations(2);
  // Duration of choice options
  int ChoiceDur = presentations(3);
  
  double Fe = param(0); // Semantic to orthographic feedback scalar
  double N = param(1);  // Noise constant
  double L = param(2);  // Constant leak current
  double D = param(3);  // Synaptic depletion rate
  double R = param(4);  // Recovery rate
  double I = param(5);  // Inhibition constant
  double Th = param(6); // Activation threshold
  double Ta = param(7); // Temporal attention
  Rcpp::NumericVector S(3);
  S(0) = param(8);      // Integration time constants at each level
  S(1) = param(9);
  S(2) = param(10);
  
  // Set these as variables to be stored
  arma::mat Latency( 1, 2, arma::fill::zeros );
  double acc;
  
  int VPR = 0;  // Visual prime
  int VTR = 1;  // Visual target
  int VMK = 2;  // Visual mask
  int VTRC = 3; // Visual target choice
  int VFLC = 4; // Visual foil choice
  
  // Assignments for orthographic and semantic layers
  int TARG = 0; // Target
  int FOIL = 1; // Foil
  
  // Weight matrices
  arma::mat VisOrth( 5, 2, arma::fill::zeros );
  VisOrth( 0, 0 ) = primeInput(0);
  VisOrth( 0, 1 ) = primeInput(1);
  VisOrth( 1, 0 ) = 1.0;
  VisOrth( 3, 0 ) = 1.0;
  VisOrth( 4, 1 ) = 1.0;
  
  // Identity matrix for weights
  arma::mat OrthSem(2,2,arma::fill::eye);
  // Same identity matrix for feedback
  arma::mat SemOrth(2,2,arma::fill::eye);
  
  // Determine time when choices are presented
  int SOA = PrimeDur + TarDur + MaskDur;
  // Duration of entire trial
  int TrialDur = SOA + ChoiceDur;
  
  // Create matrix to track membrane potential and activation
  Rcpp::NumericMatrix firing_rate( TrialDur, 2 );
  
  // Initialize neural variables
  
  // Visual layer
  arma::mat mem_vis( 1, 5, arma::fill::zeros );
  arma::mat amp_vis( 1, 5, arma::fill::ones );
  arma::mat out_vis( 1, 5, arma::fill::zeros );
  
  // Orthographic layer
  arma::mat mem_orth( 1, 2, arma::fill::zeros );
  arma::mat amp_orth( 1, 2, arma::fill::ones );
  arma::mat out_orth( 1, 2, arma::fill::zeros );
  
  // Semantic-lexical layer
  arma::mat mem_sem( 1, 2, arma::fill::zeros );
  arma::mat amp_sem( 1, 2, arma::fill::ones );
  arma::mat out_sem( 1, 2, arma::fill::zeros );
  
  // Needed to check for peak output
  arma::mat old_sem( 1, 2, arma::fill::zeros );
  
  // Initialize input
  arma::mat inp_vis( 1, 5);
  
  // Update layers every ms
  for ( int t = 1; t <= TrialDur; t++ ) {
    
    // Present prime
    if ( t == 1 ) {
      inp_vis.zeros();
      inp_vis(0,VPR) = 1.0;
    }
    // Present target
    if ( t == PrimeDur + 1 ) {
      inp_vis.zeros();
      inp_vis(0,VTR) = Ta;
    }
    // Present mask
    if ( t == PrimeDur + TarDur + 1 ) {
      inp_vis.zeros();
      inp_vis(0,VMK) = 1.0;
    }
    // Present choices
    if ( t == SOA+1 ) {
      inp_vis.zeros();
      inp_vis(0,VTRC) = 1.0;
      inp_vis(0,VFLC) = 1.0;
    }
    
    // Update visual layer
    arma::mat new_mem_vis( 1, 5, arma::fill::zeros );
    arma::mat new_amp_vis( 1, 5, arma::fill::zeros );
    update_nROUSE( new_mem_vis, new_amp_vis, out_vis, 
                   inp_vis, mem_vis, amp_vis, 
                   Th, L, R, D, I, S, 1 );
    
    // Determine input from visual layer
    arma::mat inp_orth = out_vis*VisOrth;
    inp_orth = inp_orth + Fe*( out_sem*SemOrth );
    
    // Update orthographic layer
    arma::mat new_mem_orth( 1, 5, arma::fill::zeros );
    arma::mat new_amp_orth( 1, 5, arma::fill::zeros );
    update_nROUSE( new_mem_orth, new_amp_orth, out_orth, 
                   inp_orth, mem_orth, amp_orth, 
                   Th, L, R, D, I, S, 2 );
    
    // Determine input to lexical-semantic layer
    arma::mat inp_sem = out_orth*OrthSem;
    
    // Update semantic-lexical layer
    arma::mat new_mem_sem( 1, 5, arma::fill::zeros );
    arma::mat new_amp_sem( 1, 5, arma::fill::zeros );
    update_nROUSE( new_mem_sem, new_amp_sem, out_sem, 
                   inp_sem, mem_sem, amp_sem, 
                   Th, L, R, D, I, S, 3 );
    
    // Perceptual decision process
    
    // The +50 gives things a chance to get going before peak 
    // activation is checked
    if ( t > SOA + 50 ) {
      for ( int tf = 0; tf < 2; tf++ ) {
        // Step through index for target and foil
        if ( ( out_sem(tf) < old_sem(tf) ) && 
             ( Latency(0,tf) == 0 ) ) {
          // Check for peak activation
          Latency(0,tf) = t - SOA;
        }
      }
      old_sem = out_sem;
    }
    
    // Swap old variables with new variables
    mem_vis = new_mem_vis;
    amp_vis = new_amp_vis;
    mem_orth = new_mem_orth;
    amp_orth = new_amp_orth;
    mem_sem = new_mem_sem;
    amp_sem = new_amp_sem;
    
    firing_rate(t-1,0) = out_sem(0,0);
    firing_rate(t-1,1) = out_sem(0,1);
  }
  
  // Calculate accuracy
  // Average difference between target and foil latency
  double mean_diff = Latency(0,FOIL) - Latency(0,TARG);
  // Variance of difference between target and foil latency
  double var_diff = accu( exp( N * Latency.row(0) ) );
  
  acc = 1.0 - R::pnorm( 0.0, mean_diff, pow(var_diff,0.5), 1, 0 );
  
  if ( ( Latency(0,TARG) == 0 ) && 
       ( Latency(0,FOIL) > 0 ) ) {
    // Target never launched
    acc = 0.0;
  }
  if ( ( Latency(0,TARG) > 0 ) && 
       ( Latency(0,FOIL) == 0 ) ) {
    // Foil never launched
    acc = 1.0;
  }
  if ( ( Latency(0,TARG) == 0 ) && 
       ( Latency(0,FOIL) == 0 ) ) {
    // Neither launched
    acc = 0.5;
  }
  
  Rcpp::NumericVector out(3);
  out(0) = Latency(0,0);
  out(1) = Latency(0,1);
  out(2) = acc;
  out.names() = Rcpp::CharacterVector::create(
    "Target", "Foil", "Accuracy");
  
  colnames(firing_rate) = Rcpp::CharacterVector::create("Target", "Foil");
  
  // Create a list of output
  return Rcpp::List::create( 
    Rcpp::Named("Latencies", out ),
    Rcpp::Named("Activation", firing_rate ) );
}