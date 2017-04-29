# nROUSE
A R package for simulating and fitting the nROUSE model (Huber and O'Reilly, 2003).

## Getting started

### Prerequisites

- The program R ( version >= 3.0.0 )
- A C++ compiler (e.g., via [Rtools](https://cran.r-project.org/bin/windows/Rtools/))
- The R packages Rcpp and RcppArmadillo (for the underlying C++ code and matrix algebra).

### Installation

To easily install the 'nROUSE' package, you'll need the 'devtools' package:  
```
install.packages("devtools")
library(devtools)
```

The 'nROUSE' package can then be installed via the following command:  
```
install_github("rettopnivek/nROUSE")
```

## Using the package

The package can be loaded via:
```
library(nROUSE)
```

A figure with a visual demonstration of how the nROUSE model generates predictions via simulated neural activity can be produced with the command:
```
nROUSE_demonstration()
```

The predicted accuracy over varying prime durations and types can be produced via:
```
nROUSE_predictions()
```

The nROUSE model can be simulated via `simulate_nROUSE`, and the log-likelihood for a set of data can be computed via the function `nROUSE_logLik`. This latter function can be passed into a optimization routine (e.g., `optim`) to carry out maximum likelihood estimation.

Additional details and examples for all four functions can be obtained via:
```
help("function")
```
substituting the appropriate function name for "function".

## Authors

Kevin Potter

## License

MIT
