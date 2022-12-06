# lower limits for models A and B
PCMParamLowerLimit.BM <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$X0)) {
    o$X0[] <- -5
  }
  
  if(is.Global(o$Sigma_x)) {
    diag(o$Sigma_x) <- 0.0001
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[upper.tri(o$Sigma_x)] <- -2.0
    }
  } else {
    for(r in seq_len(R)) {
      diag(o$Sigma_x[, , r]) <-  0.0001
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[,,r][upper.tri(o$Sigma_x[,,r])] <- -2.0
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
    o$X0[] <- 5
  }
  
  if(is.Global(o$Sigma_x)) {
    diag(o$Sigma_x) <- 2.0
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[upper.tri(o$Sigma_x)] <- 2.0
    }
  } else {
    for(r in seq_len(R)) {
      diag(o$Sigma_x[,,r]) <- 2.0
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[,,r][upper.tri(o$Sigma_x[,,r])] <- 2.0
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
    o$X0[] <- -5
  }
  
  if(is.Global(o$H)) {
    diag(o$H) <- 0
    if(!is.Diagonal(o$H)) {
      o$H[upper.tri(o$H)] <- -2.0
    }
  } else {
    for(r in seq_len(R)) {
      diag(o$H[,,r]) <- 0
      if(!is.Diagonal(o$H)) {
        o$H[,,r][upper.tri(o$H[,,r])] <- -2.0
      }
    }
  }
  
  if(is.Global(o$Theta)) {
    o$Theta[] <- -5
  } else {
    for(r in seq_len(R)) {
      o$Theta[, r] <- -5
    }
  }
  if(is.Global(o$Sigma_x)) {
    diag(o$Sigma_x) <- 0.0001
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[upper.tri(o$Sigma_x)] <- -2.0
    }
  } else {
    for(r in seq_len(R)) {
      diag(o$Sigma_x[, , r]) <-  0.0001
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[,,r][upper.tri(o$Sigma_x[,,r])] <- -2.0
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
    o$X0[] <- 5
  } 
  
  if(is.Global(o$H)) {
    diag(o$H) <- 2.0
    if(!is.Diagonal(o$H)) {
      o$H[upper.tri(o$H)] <- 2.0
    }
  } else {
    for(r in seq_len(R)) {
      diag(o$H[,,r]) <- 2.0
      if(!is.Diagonal(o$H)) {
        o$H[,,r][upper.tri(o$H[,,r])] <- 2.0
      }
    }
  }
  
  
  if(is.Global(o$Theta)) {
    o$Theta[] <- 5
  } else {
    for(r in seq_len(R)) {
      o$Theta[, r] <- 5
    }
  }
  if(is.Global(o$Sigma_x)) {
    diag(o$Sigma_x) <- 2.0
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[upper.tri(o$Sigma_x)] <- 2.0
    }
  } else {
    for(r in seq_len(R)) {
      diag(o$Sigma_x[,,r]) <- 2.0
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[,,r][upper.tri(o$Sigma_x[,,r])] <- 2.0
      }
    }
  }
  o
}

# lower limits for models A and B
PCMParamLowerLimit.BMkappa <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$X0)) {
    o$X0[] <- -5
  }
  
  if(is.ScalarParameter(o$kappa)) {
    o$kappa[1] <- 0.0001
  }
  o
}

# upper limits for models A and B
PCMParamUpperLimit.BMkappa <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$X0)) {
    o$X0[] <- 5
  }
  
  if(is.ScalarParameter(o$kappa)) {
    o$kappa[1] <- 5
  }
  o
}

# lower limits for models C, ..., F.
PCMParamLowerLimit.OUkappa <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$X0)) {
    o$X0[] <- -5
  }
  
  if(is.Global(o$H)) {
    diag(o$H) <- 0
    if(!is.Diagonal(o$H)) {
      o$H[upper.tri(o$H)] <- -2.0
    }
  } else {
    for(r in seq_len(R)) {
      diag(o$H[,,r]) <- 0
      if(!is.Diagonal(o$H)) {
        o$H[,,r][upper.tri(o$H[,,r])] <- -2.0
      }
    }
  }
  
  if(is.Global(o$Theta)) {
    o$Theta[] <- -5
  } else {
    for(r in seq_len(R)) {
      o$Theta[, r] <- -5
    }
  }
  
  if(is.ScalarParameter(o$kappa)) {
    o$kappa[1] <- 0.0001
  }
  o
}

# upper limits for models C, ..., F.
PCMParamUpperLimit.OUkappa <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$X0)) {
    o$X0[] <- 5
  } 
  
  if(is.Global(o$H)) {
    diag(o$H) <- 2.0
    if(!is.Diagonal(o$H)) {
      o$H[upper.tri(o$H)] <- 2.0
    }
  } else {
    for(r in seq_len(R)) {
      diag(o$H[,,r]) <- 2.0
      if(!is.Diagonal(o$H)) {
        o$H[,,r][upper.tri(o$H[,,r])] <- 2.0
      }
    }
  }
  
  if(is.Global(o$Theta)) {
    o$Theta[] <- 5
  } else {
    for(r in seq_len(R)) {
      o$Theta[, r] <- 5
    }
  }
  
  if(is.ScalarParameter(o$kappa)) {
    o$kappa[1] <- 5
  }
  o
}