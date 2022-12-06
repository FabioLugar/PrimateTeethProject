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
PCMParamLowerLimit.OU <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$X0)) {
    o$X0[1] <- 0.0
    o$X0[2] <- 0.0
  }
  
  if(is.Global(o$H)) {
    o$H[1, 1] <- o$H[2, 2] <- 0
    if(!is.Diagonal(o$H)) {
      o$H[1, 2] <- -2.0
    }
  } else {
    for(r in seq_len(R)) {
      o$H[1, 1, r] <- o$H[2, 2, r] <- 0
      if(!is.Diagonal(o$H)) {
        o$H[1, 2, r] <- -2.0
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
PCMParamUpperLimit.OU <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$X0)) {
    o$X0[1] <- 5
    o$X0[2] <- 5
  } 
  
  if(is.Global(o$H)) {
    o$H[1, 1] <- o$H[2, 2] <- 2.0
    if(!is.Diagonal(o$H)) {
      o$H[1, 2] <- 2.0
    }
  } else {
    for(r in seq_len(R)) {
      o$H[1, 1, r] <- o$H[2, 2, r] <- 2.0
      if(!is.Diagonal(o$H)) {
        o$H[1, 2, r] <- 2.0
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

PCMParamLowerLimit.BMkappa <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$X0)) {
    o$X0[1] <- 0.0
    o$X0[2] <- 0.0
  }
  
  if(is.ScalarParameter(o$kappa)) {
    o$kappa[1] <- 0.0001
  }
  o
}
PCMParamUpperLimit.BMkappa <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$X0)) {
    o$X0[1] <- 5
    o$X0[2] <- 5
  }
  
  if(is.ScalarParameter(o$kappa)) {
    o$kappa[1] <- 10
  }
  o
}
PCMParamLowerLimit.OUkappa <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$X0)) {
    o$X0[1] <- 0.0
    o$X0[2] <- 0.0
  }
  
  if(is.Global(o$H)) {
    o$H[1, 1] <- o$H[2, 2] <- 0
    if(!is.Diagonal(o$H)) {
      o$H[1, 2] <- -2.0
    }
  } else {
    for(r in seq_len(R)) {
      o$H[1, 1, r] <- o$H[2, 2, r] <- 0
      if(!is.Diagonal(o$H)) {
        o$H[1, 2, r] <- -2.0
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

  if(is.ScalarParameter(o$kappa)) {
    o$kappa[1] <- 0.0001
  }
  o
}
PCMParamUpperLimit.OUkappa <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$X0)) {
    o$X0[1] <- 5
    o$X0[2] <- 5
  } 
  
  if(is.Global(o$H)) {
    o$H[1, 1] <- o$H[2, 2] <- 2.0
    if(!is.Diagonal(o$H)) {
      o$H[1, 2] <- 2.0
    }
  } else {
    for(r in seq_len(R)) {
      o$H[1, 1, r] <- o$H[2, 2, r] <- 2.0
      if(!is.Diagonal(o$H)) {
        o$H[1, 2, r] <- 2.0
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
  
  if(is.ScalarParameter(o$kappa)) {
    o$kappa[1] <- 10
  }
  o
}