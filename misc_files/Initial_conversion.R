#------------------------------------#
# R code for implementing the nROUSE #
# model (Huber and O'Reilly, 2003)   #
#------------------------------------#

# Assignments for visual layer

VPR = 1;  # Visual prime
VTR = 2;  # Visual target
VMK = 3;  # Visual mask
VTRC = 4; # Visual target choice
VFLC = 5; # Visual foil choice

# Assignments for orthographic and semantic layers

TARG = 1; # Target
FOIL = 2; # Foil

# Weight matrices

VisOrth = rbind(
  c( 0, 0 ), # From VPR
  c( 1, 0 ), # From VTR
  c( 0, 0 ), # From VMK
  c( 1, 0 ), # From VTRC
  c( 0, 1 )  # From VFLC
)
OrthSem = diag(2); # Identity matrix for weights
SemOrth = diag(2); # Same identity matrix for feedback

# Experimental procedures

TarDur = 50; # Duration of target flash (ms)
MaskDur = 500 - TarDur; # Duration of target mask
ChoiceDur = 500; # Duration of choice display
durations = c( 17,50,150,400,2000 ) # Prime durations

# Parameters

Fe = .25;           # Semantic to orthographic feedback scalar
N = 0.0302;         # Noise constant
L = .15;            # Constant leak current
D = .324;           # Synaptic depletion rate
R = .022;           # Recovery rate
I = 0.9844;         # Inhibition constant
Th = .15;           # Activation threshold
S = c( S1 = 0.0294, 
       S2 = 0.0609, 
       S3 = .015 ); # Integration time constants at each level

my_subplus = function( x ) {
  
  x[ x <= 0 ] = 0;
  return( x )
  
}

nROUSE_update = function( old_mem, old_amp, inp, level ) {
  # Purpose:
  # Forthcoming
  # Arguments:
  # old_mem -
  # old_amp -
  # inp     -
  # level   -
  # Returns:
  # Forthcoming
  
  # Output is above threshold activation multiplied by available 
  # synaptic resources
  old_out = my_subplus( old_mem - Th ) * old_amp;
  
  inhibit = numeric(ncol(old_out))
  if ( level == 1 ) { # Visual inhibition
    
    # Prime, target, and mask mutually inhibit each other 
    # (centrally presented)
    inhibit[1:3] = rep( sum( old_out[1,1:3] ), 3 );
    # Choice words only self inhibit (unique screen location)
    inhibit[4:5] = old_out[1,4:5];
  } else {
    # For orthography and semantics, everything inhibits everything
    inhibit = rep( sum(old_out), ncol(old_out) );
  }
  
  # Update membrane potential and synaptic resources
  new_mem = old_mem + ( S[level]*(  ( inp*(1 - old_mem) )  -  
                                      (L*old_mem)  -  
                                      (I*inhibit*old_mem) ) );
  new_amp = old_amp + ( S[level]*(  (  R*(1 - old_amp) )  -  
                                      (D*old_out) ) );
  new_mem[new_mem<0] = 0; # Keep things within bounds
  new_mem[new_mem>1] = 1;
  new_amp[new_amp<0] = 0;
  new_amp[new_amp>1] = 1;
  
  return( list( new_mem, new_amp, old_out ) )
}


nROUSE_simulate = function() {
  
  Latency = matrix( 0, length( durations ), 2 )
  acc = numeric( length( durations ) )
  
  for (pd in 1:length( durations ) ) {
    
    # Extract prime duration
    PrimeDur = durations[pd];
    # Determine time when choices are presented
    SOA = PrimeDur + TarDur + MaskDur;
    # Total duration of the trial
    TrialDur = SOA + ChoiceDur;
    
    # Initialize neural variables
    mem_vis  = matrix( 0, 1, 5 );
    amp_vis  = matrix( 1, 1, 5 );
    out_vis  = matrix( 0, 1, 5 );
    mem_orth = matrix( 0, 1, 2 );
    amp_orth = matrix( 1, 1, 2 );
    out_orth = matrix( 0, 1, 2 );
    mem_sem  = matrix( 0, 1, 2 );
    amp_sem  = matrix( 1, 1, 2 );
    out_sem  = matrix( 0, 1, 2 );
    old_sem  = matrix( 0, 1, 2 ); # Needed to check for peak output
    
    for ( t in 1:TrialDur ) {
      
      # Determine input into visual layer
      
      # Present prime
      if ( t == 1 ) {
        inp_vis = matrix( 0, 1, 5 );
        inp_vis[1,VPR] = 1;
      }
      # Present target
      if ( t == PrimeDur + 1 ) {
        inp_vis = matrix( 0, 1, 5 );
        inp_vis[1,VTR] = 1;
      }
      # Present mask
      if ( t == PrimeDur + TarDur + 1 ) {
        inp_vis = matrix( 0, 1, 5 );
        inp_vis[1,VMK] = 1;
      }
      # Present choices
      if ( t == SOA+1 ) {
        inp_vis = matrix( 0, 1, 5 );
        inp_vis[1,VTRC] = 1;
        inp_vis[1,VFLC] = 1;
      }
      
      # Update visual layer
      out = nROUSE_update( mem_vis, amp_vis, inp_vis, 1 )
      new_mem_vis = out[[1]];
      new_amp_vis = out[[2]];
      out_vis = out[[3]];
      
      # Determine input from visual layer
      inp_orth = out_vis %*% VisOrth;
      inp_orth = inp_orth + Fe*( out_sem %*% SemOrth );
      
      # Update orthographic layer
      out = nROUSE_update( mem_orth, amp_orth, inp_orth, 2 )
      new_mem_orth = out[[1]];
      new_amp_orth = out[[2]];
      out_orth = out[[3]];
      
      # Determine input to lexical-semantic layer
      inp_sem = out_orth %*% OrthSem;
      
      # Update semantic layer
      out = nROUSE_update( mem_sem, amp_sem, inp_sem, 3 )
      new_mem_sem = out[[1]];
      new_amp_sem = out[[2]];
      out_sem = out[[3]];
      
      # Perceptual decision process
      
      # The +50 gives things a chance to get going before peak 
      # activation is checked
      if ( t > SOA + 50 ) {
        for ( tf in 1:2 ) {
          # Step through index for target and foil
          if ( out_sem[tf] < old_sem[tf] & Latency[pd,tf] == 0 ) {
            # Check for peak activation
            Latency[pd,tf] = t - SOA;
          }
        }
        old_sem = out_sem;
      }
      
      # Swap old variables with new variables
      mem_vis = new_mem_vis;
      amp_vis = new_amp_vis;
      mem_orth = new_mem_orth;
      amp_orth = new_amp_orth;
      mem_sem = new_mem_sem;
      amp_sem = new_amp_sem;
      
    }
    
    # Calculate accuracy
    
    # Average difference between target and foil latency
    mean_diff = Latency[pd,FOIL] - Latency[pd,TARG];
    # Variance of difference between target and foil latency
    var_diff = sum( exp(N * Latency[pd,]) );
    
    acc[pd] = 1 - pnorm( 0, mean_diff, sqrt( var_diff ) );
    
    if ( Latency[pd,TARG] == 0 & Latency[pd,FOIL] > 0 ) {
      # Target never launched
      acc[pd] = 0;
    }
    if ( Latency[pd,TARG] > 0 & Latency[pd,FOIL] == 0 ) { 
      # Foil never launched
      acc[pd] = 1;
    }
    if ( Latency[pd,TARG] == 0 & Latency[pd,FOIL] == 0 ) {
      # Neither launched
      acc[pd]=.5;
    }
    
  }
  
  return( list( Latency = Latency, Acc = acc ) )
}

### Test code ###

startTime = Sys.time()

target_lat = matrix( NA, length(durations), 2 )
foil_lat = matrix( NA, length(durations), 2 )
acc = matrix( NA, length(durations), 2 )

for (cd in 1:2) {
  if ( cd == 1 ) { # Target primed
    # Set to 2 because there are two visual copies onscreen
    VisOrth[VPR,] = c(2,0)
  }
  if ( cd == 2 ) { # Foil primed
    VisOrth[VPR,] = c(0,2)
  }
  tst = nROUSE_simulate()
  target_lat[,cd] = tst$Latency[,TARG]
  foil_lat[,cd] = tst$Latency[,FOIL]
  acc[,cd] = tst$Acc
}

R_res = list( target_lat, foil_lat, acc );

endTime = Sys.time() - startTime
print( endTime )

print( target_lat )
print( foil_lat )
print( round( acc, 4 ) )