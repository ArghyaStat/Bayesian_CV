# This file contains the importance sampling proposals of different methods
# "Gelfand", "Bronze", "Zanella" with classical methods
# "TIS and PSIS" on above methods


log_sum_exp <- function(x){
  
  m <- max(x)
  
  m + log(sum(exp(x - m)))
}


# Gelfand (original posterior)


log_q_g <- function(theta,
                    X,
                    y,
                    log_like,
                    log_prior){
  
  ll <- log_like(theta, X, y)
  
  log_prior(theta) + sum(ll)
}


grad_log_q_g <- function(theta,
                         X,
                         y,
                         grad_log_like,
                         grad_log_prior){
  
  grad_ll <- grad_log_like(theta, X, y)
  
  grad_lp <- grad_log_prior(theta)
  
  grad_lp + colSums(grad_ll)
}




# Bronze proposal

log_q_b <- function(theta,
                    X,
                    y,
                    alpha,
                    log_like,
                    log_prior){
  
  ll <- log_like(theta, X, y)
  
  log_prior(theta) + alpha * sum(ll)
}


grad_log_q_b <- function(theta,
                         X,
                         y,
                         alpha,
                         grad_log_like,
                         grad_log_prior){
  
  grad_ll <- grad_log_like(theta, X, y)
  
  grad_lp <- grad_log_prior(theta)
  
  grad_lp + alpha * colSums(grad_ll)
}




# Zanella proposal

log_q_z <- function(theta,
                    X,
                    y,
                    log_like,
                    log_prior){
  
  ll <- log_like(theta, X, y)
  
  log_prior(theta) + sum(ll) + log_sum_exp(-ll)
}


grad_log_q_z <- function(theta,
                         X,
                         y,
                         grad_log_like,
                         grad_log_prior,
                         log_like){
  
  ll <- log_like(theta, X, y)
  
  grad_ll <- grad_log_like(theta, X, y)
  
  grad_lp <- grad_log_prior(theta)
  
  log_w <- -ll
  
  log_w <- log_w - log_sum_exp(log_w)
  
  w <- exp(log_w)
  
  grad_mix <- colSums(grad_ll * w)
  
  grad_lp + colSums(grad_ll) - grad_mix
}


# 
# 
# U_g <- function(theta,
#                 X,
#                 y,
#                 log_like,
#                 log_prior){
#   
#   -log_q_g(
#     theta,
#     X,
#     y,
#     log_like,
#     log_prior
#   )
# }
# 
# 
# grad_g <- function(theta,
#                    X,
#                    y,
#                    grad_log_like,
#                    grad_log_prior){
#   
#   -grad_log_q_g(
#     theta,
#     X,
#     y,
#     grad_log_like,
#     grad_log_prior
#   )
# }
# 
# 
# U_b <- function(theta,
#                 X,
#                 y,
#                 alpha,
#                 log_like,
#                 log_prior){
#   
#   -log_q_b(
#     theta,
#     X,
#     y,
#     alpha,
#     log_like,
#     log_prior
#   )
# }
# 
# 
# grad_b <- function(theta,
#                    X,
#                    y,
#                    alpha,
#                    grad_log_like,
#                    grad_log_prior){
#   
#   -grad_log_q_b(
#     theta,
#     X,
#     y,
#     alpha,
#     grad_log_like,
#     grad_log_prior
#   )
# }
# 
# 
# U_z <- function(theta,
#                 X,
#                 y,
#                 log_like,
#                 log_prior){
#   
#   -log_q_z(
#     theta,
#     X,
#     y,
#     log_like,
#     log_prior
#   )
# }
# 
# 
# grad_z <- function(theta,
#                    X,
#                    y,
#                    grad_log_like,
#                    grad_log_prior,
#                    log_like){
#   
#   -grad_log_q_z(
#     theta,
#     X,
#     y,
#     grad_log_like,
#     grad_log_prior,
#     log_like
#   )
# }
# 
