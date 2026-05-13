# Hamiltonian Monte Carlo

leapfrog_hmc <- function(epsilon,
                         L,
                         niters,
                         U,
                         grad_U,
                         x_init,
                         M,
                         ...){
  
  d <- length(x_init)
  
  M_inv <- chol2inv(chol(M))
  
  x <- matrix(NA, niters, d)
  
  x[1, ] <- x_init
  
  chol_M <- chol(M)
  
  p <- matrix(rnorm(niters * d), niters, d) %*% chol_M
  
  acc <- 0

  
  for(iter in 2:niters){
    
    x_prop <- x[iter - 1, ]
    
    p_prop <- p[iter, ]
    
    # Half momentum step
    
    g <- grad_U(x_prop, ...)
    
    p_prop <- p_prop - 0.5 * epsilon * g

    # Leapfrog steps
    
    for(l in 1:L){
      
      x_prop <- x_prop + epsilon *as.vector(M_inv %*% p_prop)
      
      g <- grad_U(x_prop, ...)
      
      if(l != L){
        
        p_prop <- p_prop - epsilon * g
      }
    }
 
    # Final half momentum step
    
    p_prop <- p_prop - 0.5 * epsilon * g
   
    # Hamiltonians
    
    H_curr <- U(x[iter - 1, ], ...) + 
      0.5 * t(p[iter, ]) %*% (M_inv %*% p[iter, ])
    
    H_prop <- U(x_prop, ... ) + 
      0.5 * t(p_prop) %*% (M_inv %*% p_prop)
  
    # Accept/reject
    
    if(log(runif(1)) < -(H_prop - H_curr)){
      
      x[iter, ] <- x_prop
      acc <- acc + 1
      
    } else {
      
      x[iter, ] <- x[iter - 1, ]
    }
  }
  
  list(
    samples = x,
    acc = acc / niters
  )
}