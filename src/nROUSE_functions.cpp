#include "RcppArmadillo.h"

// [[Rcpp::export]]
arma::mat subplus( arma::mat A, double Th ) {
  
  arma::mat B( A.n_rows, A.n_cols, arma::fill::ones );
  B = A - B*Th;
  
  B.elem( find( B <= 0.0 ) ).zeros();
  
  return( B );
}

// [[Rcpp::export]]
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
    for ( int i = 0; i < old_out.n_cols; i++ ) inhibit(0,i) = sums;
    
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

// [[Rcpp::export]]
Rcpp::List simulate_nROUSE( Rcpp::NumericVector durations, 
                            Rcpp::NumericVector presentations, 
                            Rcpp::NumericVector param, 
                            Rcpp::NumericVector primeInput ) {
  
  // Duration of target flash
  int TarDur = presentations(0);
  // Duration of target mask
  int MaskDur = presentations(1);
  // Duration of choice options
  int ChoiceDur = presentations(2);
  
  double Fe = param(0); // Semantic to orthographic feedback scalar
  double N = param(1);  // Noise constant
  double L = param(2);  // Constant leak current
  double D = param(3);  // Synaptic depletion rate
  double R = param(4);  // Recovery rate
  double I = param(5);  // Inhibition constant
  double Th = param(6); // Activation threshold
  Rcpp::NumericVector S(3);
  S(0) = param(7);      // Integration time constants at each level
  S(1) = param(8);
  S(2) = param(9);
  
  // Set these as variables to be stored
  arma::mat Latency( durations.size(), 2, arma::fill::zeros );
  Rcpp::NumericVector acc( durations.size() );
  
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
  
  
  arma::mat out_sem( 1, 5, arma::fill::zeros );
  for ( int pd = 0; pd < durations.size(); pd++ ) {
    
    // Extract prime durations
    int PrimeDur = durations(pd);
    // Determine time when choices are presented
    int SOA = PrimeDur + TarDur + MaskDur;
    // Duration of entire trial
    int TrialDur = SOA + ChoiceDur;
    
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
        inp_vis(0,VTR) = 1.0;
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
               ( Latency(pd,tf) == 0 ) ) {
            // Check for peak activation
            Latency(pd,tf) = t - SOA;
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
      
    }
    
    // Calculate accuracy
    // Average difference between target and foil latency
    double mean_diff = Latency(pd,FOIL) - Latency(pd,TARG);
    // Variance of difference between target and foil latency
    double var_diff = accu( exp( N * Latency.row(pd) ) );
    
    acc(pd) = 1.0 - R::pnorm( 0.0, mean_diff, pow(var_diff,0.5), 1, 0 );
    
    if ( ( Latency(pd,TARG) == 0 ) && 
         ( Latency(pd,FOIL) > 0 ) ) {
      // Target never launched
      acc(pd) = 0.0;
    }
    if ( ( Latency(pd,TARG) > 0 ) && 
         ( Latency(pd,FOIL) == 0 ) ) {
      // Foil never launched
      acc(pd) = 1.0;
    }
    if ( ( Latency(pd,TARG) == 0 ) && 
         ( Latency(pd,FOIL) == 0 ) ) {
      // Neither launched
      acc(pd) = 0.5;
    }
    
  }
  
  // Create a list of output
  return Rcpp::List::create( 
    Rcpp::Named("Latency", Latency ),
    Rcpp::Named("Accuracy", acc ) );
}