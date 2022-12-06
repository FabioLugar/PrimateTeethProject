# lower limits for models A and B
PCMParamLowerLimit.BM <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$X0)) {
    o$X0[1] <- 0.0
    o$X0[2] <- 0.0
  }
  
  if(is.Global(o$Sigma_x)) {
    o$Sigma_x[1, 1] <- o$Sigma_x[2, 2] <- 0.0001
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- -.0
    }
  } else {
    for(r in seq_len(R)) {
      o$Sigma_x[1, 1, r] <- o$Sigma_x[2, 2, r] <- 0.0001
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[1, 2, r] <- -2.0
      }
    }
  }
  o
}

# upper limits for models A and B
PCMParamUpperLimit.BM <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$X0)) {
    o$X0[1] <- 5
    o$X0[2] <- 5
  }
  
  if(is.Global(o$Sigma_x)) {
    o$Sigma_x[1, 1] <- o$Sigma_x[2, 2] <- 2.0
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- 2.0
    }
  } else {
    for(r in seq_len(R)) {
      o$Sigma_x[1, 1, r] <- o$Sigma_x[2, 2, r] <- 2.0
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[1, 2, r] <- 2.0
      }
    }
  }
  o
}

# lower limits for models C, ..., F.
PCMParamLowerLimit.OU <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$X0)) {
    o$X0[1] <- 0.0
    o$X0[2] <- 0.0
  }
  
  if(is.Global(o$H)) {
    o$H[1, 1] <- o$H[2, 2] <- 0.001
    if(!is.Diagonal(o$H)) {
      o$H[1, 2] <- -10.0
    }
  } else {
    for(r in seq_len(R)) {
      o$H[1, 1, r] <- o$H[2, 2, r] <- 0.001
      if(!is.Diagonal(o$H)) {
        o$H[1, 2, r] <- -10.0
      }
    }
  }
  
  if(is.Global(o$Theta)) {
    o$Theta[1] <- 0.0
    o$Theta[2] <- 0.0
  } else {
    for(r in seq_len(R)) {
      o$Theta[1, r] <- 0.0
      o$Theta[2, r] <- 0.0
    }
  }
  if(is.Global(o$Sigma_x)) {
    o$Sigma_x[1, 1] <- o$Sigma_x[2, 2] <- 0.0001
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- -2.0
    }
  } else {
    for(r in seq_len(R)) {
      o$Sigma_x[1, 1, r] <- o$Sigma_x[2, 2, r] <- 0.0001
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[1, 2, r] <- -2.0
      }
    }
  }
  o
}

# upper limits for models C, ..., F.
PCMParamUpperLimit.OU <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$X0)) {
    o$X0[1] <- 5
    o$X0[2] <- 5
  } 
  
  if(is.Global(o$H)) {
    o$H[1, 1] <- o$H[2, 2] <- 10.0
    if(!is.Diagonal(o$H)) {
      o$H[1, 2] <- 10.0
    }
  } else {
    for(r in seq_len(R)) {
      o$H[1, 1, r] <- o$H[2, 2, r] <- 10.0
      if(!is.Diagonal(o$H)) {
        o$H[1, 2, r] <- 10.0
      }
    }
  }
  
  if(is.Global(o$Theta)) {
    o$Theta[1] <- 5
    o$Theta[2] <- 5
  } else {
    for(r in seq_len(R)) {
      o$Theta[1, r] <- 5
      o$Theta[2, r] <- 5
    }
  }
  if(is.Global(o$Sigma_x)) {
    o$Sigma_x[1, 1] <- o$Sigma_x[2, 2] <- 2.0
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- 2.0
    }
  } else {
    for(r in seq_len(R)) {
      o$Sigma_x[1, 1, r] <- o$Sigma_x[2, 2, r] <- 2.0
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[1, 2, r] <- 2.0
      }
    }
  }
  o
}