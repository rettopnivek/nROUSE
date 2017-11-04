#----------------------------------------------#
# R code for plotting/demonstrations of nROUSE #
#----------------------------------------------#

# Index
# Lookup - 01:  nROUSE_predictions
# Lookup - 02:  nROUSE_demonstration

# Lookup - 01
#' Predictions from the nROUSE model
#'
#' Generates a plot of the predicted accuracy over different 
#' prime durations and types for the nROUSE model (output can 
#' be generated as well).
#'
#' @param primeDurations a vector of prime durations between 
#'   17 and 2000 ms.
#' @param primeTypes a character vector with up to 4 types of prime: 
#'   'Target', 'Foil', 'Neither', or 'Both'.
#' @param design a vector with the duration (in ms) of the target flash,
#'   the target mask, and the choice options.
#' @param opt a logical vector indicating if the predicted accuracies 
#'   latencies should be returned, and if a plot should be generated.
#' @param prm an optional named vector for the parameters of the 
#'   nROUSE model, where...
#'   \itemize{
#'     \item Fe = the feedback scalar.
#'     \item N = The noise multiplier.
#'     \item L = The constant leak current.
#'     \item D = The synaptic depletion rate.
#'     \item R = The replenishment rate.
#'     \item I = The inhibition constant.
#'     \item Th = The noise multiplier.
#'     \item Ta = The temporal attention parameter.
#'     \item SV = The visual integration rate.
#'     \item SO = The orthographic integration rate.
#'     \item SS = The semantic integration rate.
#'   }
#' 
#' 
#' @return A list consisting of a matrix of predicted accuracies by 
#'   prime duration and type, and an array with the predicted perceptual 
#'   identification latencies for the target and foil by prime duration 
#'   and prime type.
#' 
#' @examples
#' # Default
#' nROUSE_predictions()
#' 
#' # Specific conditions
#' 
#' # Prime durations
#' pd = c(50,400)
#' # Prime types
#' pt = c('Target','Foil')
#' nROUSE_predictions( primeDurations = pd, primeTypes = pt )
#' 
#' # Different timing for target stimulus flash and subsequent mask
#' nROUSE_predictions( design = c( 84, 416, 500 ) )
#' 
#' # Suppress plotting and return predictions
#' pred = nROUSE_predictions( opt = c( T, F ) )
#' 
#' # Different parameter values (Must be named)
#' prm = c( I = .5, Ta = .8, N = .025 )
#' nROUSE_predictions( prm = prm )
#' 
#' @export
nROUSE_predictions = function( primeDurations = c( 17,50,150,400,2000 ),
                    primeTypes = c('Target','Foil','Neither','Both'),
                    design = c( 50, 450, 500 ),
                    opt = c( F, T ),
                    prm = NULL ) {
  
  # Default parameters (Rieth & Huber, 2017)
  param = c( Fe = .25, N = .0302, L = .15, 
             D = .324, R = .022, I = .9844, 
             Th = .15, Ta = 1.0, SV = .0294, 
             SO = .0609, SS = .015 )
  
  # Check if parameter values should be changed
  if ( length( prm ) > 0 ) {
    
    selPar = names(prm)
    
    for (i in 1:length(selPar)) {
      chk = names(param) == selPar[i]
      if ( sum(chk) == 1 ) 
        param[ chk ] = prm[i]
    }
    
  }
  
  # Sort prime duratiosn
  primeDurations = sort( primeDurations );
  
  # Check if prime types are admissable
  chk = c('Target','Foil','Neither','Both')
  if ( sum(primeTypes %in% chk) < length(primeTypes) ) 
    stop('Check character vector for prime types.')
  
  # Matrices to store output
  Accuracy = matrix( NA, length( primeDurations ),
                     length( primeTypes ) )
  colnames( Accuracy ) = primeTypes;
  rownames( Accuracy ) = as.character( primeDurations )
  
  Latency = array( NA, dim=c( length( primeDurations ),
                              length( primeTypes ),
                              2 ) )
  
  # Simulate nROUSE
  for ( pt in 1:length( primeTypes ) ) {
    for ( pd in 1:length( primeDurations ) ) {
      pres = c( primeDurations[pd], design )
      
      if (primeTypes[pt] == 'Target') prmT = c(2,0)
      if (primeTypes[pt] == 'Foil') prmT = c(0,2)
      if (primeTypes[pt] == 'Neither') prmT = c(0,0)
      if (primeTypes[pt] == 'Both') prmT = c(1,1)
      
      sim = simulate_nROUSE( pres, prmT, param )
      
      Accuracy[pd,pt] = sim$Latencies[3]
      Latency[pd,pt,] = sim$Latencies[1:2]
      
    }
  }
  
  # Plot results
  if (opt[2]) {
    
    clrs = c( Target = 'orange', 
              Foil = 'red', 
              Neither = 'blue', 
              Both = 'purple' )
    
    xl = c( 10, max( primeDurations ) )
    
    plot( log(xl), c(0,1), type = 'n', xaxt = 'n', 
          yaxt = 'n', xlab = ' ', ylab = ' ', 
          bty = 'l' )
    xa = c(10,100,1000,5000)
    axis( 1, log( xa[ xa <= max(xl) ] ), 
          xa[ xa <= max(xl) ], cex.axis = 1.5 )
    mtext('Prime duration (log ms)',side=1,cex=1.5,
          line = 2.25)
    axis(2,seq(0,1,.25),cex.axis=1.5)
    mtext('Accuracy',side=2,cex=1.5,
          line = 2.25)
    abline( h = .5, lty = 2, lwd = 2 )
    mtext('nROUSE predictions', side = 3, line = 1, 
          cex = 1.5 )
    
    whchClr = c()
    for ( pt in 1:length( primeTypes ) ) {
      
      clr = clrs[ names(clrs) %in% primeTypes[pt] ]
      whchClr = c( whchClr, clr )
      
      lines( log( primeDurations ), Accuracy[,pt], col = clr, 
             lwd = 2 )
      points( log( primeDurations ), Accuracy[,pt], bg = clr, 
             pch = 21, cex = 1.2 )
      
    }
    legend('bottomleft',primeTypes, fill = whchClr, bty = 'n', 
           cex = 1.5 )
    
  }
  
  
  if (opt[1]) return( list( Accuracy = Accuracy, Latency = Latency ) )
}
                   
# Lookup - 02
#' Demonstration of the nROUSE model.
#'
#' Generates a plot of the activation at the highest level of the 
#' neural network (i.e. semantic-lexical) for the target and foil 
#' representations under different prime types and an inputted 
#' prime duration.
#'
#' @param primeTypes a character vector with up to 4 types of prime: 
#'   'Target', 'Foil', 'Neither', or 'Both'.
#' @param design a vector with the duration (in ms) of the target flash,
#'   the target mask, and the choice options.
#' @param prm an optional named vector for the parameters of the 
#'   nROUSE model, where...
#'   \itemize{
#'     \item Fe = the feedback scalar.
#'     \item N = The noise multiplier.
#'     \item L = The constant leak current.
#'     \item D = The synaptic depletion rate.
#'     \item R = The replenishment rate.
#'     \item I = The inhibition constant.
#'     \item Th = The noise multiplier.
#'     \item Ta = The temporal attention parameter.
#'     \item SV = The visual integration rate.
#'     \item SO = The orthographic integration rate.
#'     \item SS = The semantic integration rate.
#'   }
#' 
#' @details The user will receive a prompt to input a prime duration. 
#'   The function will then generate a plot showing the time course 
#'   of activation for the target/foil representations for the 
#'   set of prime types. The predicted accuracies and the difference 
#'   distribution from which these accuracies are derived are shown 
#'   on the right.
#' 
#' @examples
#' \dontrun{
#' # Default
#' nROUSE_demonstration()
#' 
#' # Fewer prime types
#' pt = c( 'Target', 'Foil' )
#' nROUSE_demonstration(primeTypes = pt)
#' 
#' # Different timing for target stimulus flash and subsequent mask
#' nROUSE_demonstration(design=c(84,416,500))
#' 
#' # Different parameter values (Must be named)
#' prm = c( I = .5, Ta = .8, N = .025 )
#' nROUSE_demonstration(prm=prm)
#' }
#' 
#' @export

nROUSE_demonstration = function(primeTypes = c('Target','Foil',
                                               'Neither','Both'),
                                design = c( 50, 450, 500 ), 
                                prm = NULL ) {
  
  # Default parameters
  param = c( Fe = .25, N = .0302, L = .15, 
             D = .324, R = .022, I = .9844, 
             Th = .15, Ta = 1.0, SV = .0294, 
             SO = .0609, SS = .015 )
  
  # Check if parameter values should be changed
  if ( length( prm ) > 0 ) {
    
    selPar = names(prm)
    
    for (i in 1:length(selPar)) {
      chk = names(param) == selPar[i]
      if ( sum(chk) == 1 ) 
        param[ chk ] = prm[i]
    }
    
  }
  
  cat("Prime duration? (17 - 2000)","\n")
  cat("Type number (in ms) and hit 'enter'","\n")
  
  pd = scan(n = 1)
  td = design[1]
  md = design[2]
  cd = design[3]
  
  presentations = c( pd, design )
  
  pltStr = c()
  npt = length( primeTypes )
  pltStr = matrix( rep( 1:length( primeTypes ), each = 8 ),
                   length( primeTypes ), 8, byrow = T )
  pltStr = pltStr + 0:(length(primeTypes)-1)
  pltStr = cbind( pltStr, pltStr[,1:2]+1 )

  layout( pltStr )
  
  for (pt in 1:length(primeTypes)) {
    
    if (primeTypes[pt] == 'Target') prmT = c(2,0)
    if (primeTypes[pt] == 'Foil') prmT = c(0,2)
    if (primeTypes[pt] == 'Neither') prmT = c(0,0)
    if (primeTypes[pt] == 'Both') prmT = c(1,1)
    
    sim = simulate_nROUSE( presentations,
                           prmT, param );
    
    par( mar=c( 3, 3, 2, .5 ) )
    if ( pt == nrow(pltStr) ) par( mar=c(4,3,1,.5) )
    
    yl = max( apply( sim$Activation, 2, max ) )
    yl = c(0,ceiling(yl*10)/10)
    xl = nrow(sim$Activation)
    xl = c(0,xl)
    plot( xl, yl, type = 'n', bty = 'n', 
          xlab = ' ', ylab = ' ',
          xaxt = 'n', yaxt = 'n' )
    abline(h=0,lwd=2)
    mtext( primeTypes[pt], side = 2, cex = 1.25, line = -.75 )
    if (pt==nrow(pltStr)) axis(1,round(seq(xl[1],xl[2],length=5)),
                               tick = F, line = -1, cex.axis = 1.4 )
    if (pt == 1) {
      
      xp = grconvertX( 1, from = "npc", to = "user" )*.1
      
      if ( nrow(pltStr) == 4 ) 
        yp = grconvertY( 1, from = "npc", to = "user" )*1.25
      if ( nrow(pltStr) == 3 ) 
        yp = grconvertY( 1, from = "npc", to = "user" )*1.2
      if ( nrow(pltStr) == 2 ) 
        yp = grconvertY( 1, from = "npc", to = "user" )*1.1
      
      legend( xp, yp, c('Prime','Target','Mask','Choices'),
              fill = c('blue','orange','purple','grey'),
              bty = 'n', cex= 1.4, horiz = T, xpd = T )
      
      if ( nrow(pltStr) == 4 ) 
        yp = grconvertY( 1, from = "npc", to = "user" )*.15
      if ( nrow(pltStr) == 3 ) 
        yp = grconvertY( 1, from = "npc", to = "user" )*.1
      if ( nrow(pltStr) == 2 ) 
        yp = grconvertY( 1, from = "npc", to = "user" )*.05
      
      legend( xp, 0 - yp, 
              c('Target activation', 'Foil activation'),
              fill = c('orange','red'),
              bty = 'n', cex= 1.4, horiz = T, xpd = T )
    }
    
    for (i in 1:3) segments( cumsum(presentations[1:i])+.5, 0, 
                             cumsum(presentations[1:i])+.5, yl[2], 
                             col = 'grey80' )
    
    lines( 1:nrow(sim$Activation),
           sim$Activation[,1], lwd = 2, col = 'orange' )
    lines( 1:nrow(sim$Activation),
           sim$Activation[,2], lwd = 2, col = 'red'  )
    
    pck = numeric(2)
    act = sim$Activation[ (pd+td+md+1):(pd+td+md+cd), ]
    for (i in 1:2) {
      pck[i] = which( act[,i] == max( act[,i] ) )
      segments( pck[i] + (pd+td+md), 0,
                pck[i] + (pd+td+md), act[pck[i],i], col = 'black',
                lwd = 2 )
    }
    
    segments( 1, yl[2], pd, yl[2],
              col = 'blue', lwd = 3 )
    segments( pd+1, yl[2], 
              pd+td, yl[2],
              col = 'orange', lwd = 3 )
    segments( pd+td+1, yl[2], 
              pd+td+md, yl[2],
              col = 'purple', lwd = 3 )
    segments( pd+td+md+1, yl[2], 
              pd+td+md+cd, yl[2],
              col = 'grey', lwd = 3 )
    
    mu = diff(sim$Latencies[2:1])
    sigma = sqrt( sum( exp( param['N']*sim$Latencies[1:2] ) ) )
    
    val = seq( -4*sigma, 4*sigma, length = 1000 ) + mu
    dn = dnorm( val, mu, sigma )
    
    par( mar=c(1,.5,.5,3) )
    plot( dn, val, type = 'n', xlab = ' ', ylab = ' ',
          yaxt = 'n', xaxt = 'n', lwd = 2, bty = 'n' )
    
    sel = val < 0
    polygon( c(0, dn[sel], 0 ),
             c( min(val[sel]), val[sel], max(val[sel]) ),
             col = 'grey', border = 'NA' )
    lines( dn, val, lwd = 2 )
    abline( h = 0, lwd = 2, lty = 2 )
    legend( 'bottom', 
            paste( round( sim$Latencies[3], 2 )*100, '%', sep = '' ),
            bty = 'n', cex = 1.5 )
    
  }
  mtext('Activation',side=2,outer=T,cex=1.5,line=-1.75)
  mtext('Duration (ms)',side=1,outer=T,cex=1.5,line=-1.5)
  mtext('Predicted accuracy',side=4,outer=T,cex=1.5,line=-1.5)
  
}

