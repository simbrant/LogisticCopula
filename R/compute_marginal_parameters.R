coefs_to_pars <- function(y, x, xtype, beta_0, betas) {

  # Compute xbar and s2 for each covariate.
  means_x <- apply(x, 2, mean)
  vars_x <- apply(x, 2, var)

  pars_pi_y <- function(pi_y) {

    # Compute mu and sigma parameters
    sigma <- rep(NA, length(xtype))
    mu_0 <- rep(NA, length(xtype))
    mu_1 <- rep(NA, length(xtype))

    if (any(xtype == "c_a")) {

      sigma[which(xtype == "c_a")] <- pmax(sapply(
        which(xtype == "c_a"), function(j) {
          q_j <- betas[j]**2 * pi_y * (1 - pi_y)
          (sqrt(1 + 4 * q_j * vars_x[j]) - 1) / (2 * q_j)
        }), 1e-8)

      mu_0[which(xtype == "c_a")] <- (
        (
          means_x[which(xtype == "c_a")] -
            (
              betas[which(xtype == "c_a")] *
                sigma[which(xtype == "c_a")] *
                pi_y
            )
        )
      )

      mu_1[which(xtype == "c_a")] <- (
        means_x[which(xtype == "c_a")] + (
          betas[which(xtype == "c_a")] *
            (1 - pi_y) * sigma[which(xtype == "c_a")]
        )
      )
    }

    nu_0 <- rep(NA, length(xtype))
    nu_1 <- rep(NA, length(xtype))
    if (any(xtype == "c_p")) {
      nu_0[which(xtype == "c_p")] <- sapply(
        which(xtype == "c_p"), 
        function(j) {
          roots <- polyroot(
            c((1 - pi_y) * betas[j], - (1 + means_x[j]*betas[j]), means_x[j])
          )
          
          roots <- Re(roots[abs(Im(roots)) < 1e-10 & Re(roots) > 0])
          roots <- roots[roots > betas[j]]
          if (length(roots) > 1) {
            
            nu_0_cands <- roots
            nu_1_cands <- (
              nu_0_cands - betas[j]
            )
            marglik <- sapply(
              seq(length(roots)), 
              function(k) {
                sum(log((1 - pi_y) * dexp(x[, j], nu_0_cands[k])
                        + pi_y * dexp(x[, j], nu_1_cands[k])))
              })
            roots[which.max(marglik)]
          } else {
            roots
          }
        }
      )
      nu_1[which(xtype == "c_p")] <- (
        nu_0[which(xtype == "c_p")] - betas[which(xtype == "c_p")]
      )
    }
    
    lambda_0 <- rep(NA, length(xtype))
    lambda_1 <- rep(NA, length(xtype))
    if (any(xtype == "d_c")) {
      lambda_0[which(xtype == "d_c")] <- (
        means_x[which(xtype == "d_c")] / 
          (1 + pi_y * (exp(betas[which(xtype == "d_c")]) - 1))
      )
      lambda_1[which(xtype == "d_c")] <- (
        exp(betas[which(xtype == "d_c")]) * means_x[which(xtype == "d_c")] / 
          (1 + pi_y * (exp(betas[which(xtype == "d_c")]) - 1))
      )
    }
    
    # Compute discrete parameters
    rho_0 <- rep(NA, length(xtype))
    rho_1 <- rep(NA, length(xtype))
    if (any(xtype == "d_b")) {
      
      rho_0[which(xtype == "d_b")] <- sapply(
        which(xtype == "d_b"), 
        function(j) {
          roots <- polyroot(
            c(means_x[j], -(1 + (1 - exp(betas[j])) * (means_x[j] - pi_y)),
              (1 - exp(betas[j])) * (1 - pi_y))
          )
          
          roots <- Re(roots[abs(Im(roots)) < 1e-10 & Re(roots) > 0])
          if(!any(roots > 0 & roots < 1))
          {
            if(which.min(c(abs(roots-1),abs(roots))) < 3)
            {
              roots <- 1e-15
            }
            else
            {
              roots <- 1-1e-15
            }
          }
          else 
          {
            roots <- roots[roots > 0 & roots < 1]
            if (length(roots) > 1) {
              
              rho_0_cands <- roots
              rho_1_cands <- (
                (exp(betas[j]) * rho_0_cands) /
                  (1 - (1 - exp(betas[j])) * rho_0_cands)
              )
              marglik <- sapply(
                1:length(roots), 
                function(k) {
                  sum(log((1 - pi_y) * dbinom(x[, j], 1, rho_0_cands[k])
                          + pi_y * dbinom(x[, j], 1, rho_1_cands[k])))
                })
              roots[which.max(marglik)]
            } else {
              roots
            }
          }
        }
      )
      rho_0[which(xtype == "d_b")] <- pmin(pmax(rho_0[which(xtype == "d_b")],1e-15),1-1e-15)
      
      rho_1[which(xtype == "d_b")] <- pmin(pmax((
        exp(betas[which(xtype == "d_b")]) * rho_0[which(xtype == "d_b")]
      ) / (
        1 - (1 - exp(betas[which(xtype == "d_b")])) * 
          rho_0[which(xtype == "d_b")]
      ),1e-15),1-1e-15)
    }
    list(pi_y = pi_y, mu_1 = mu_1, mu_0 = mu_0, sigma = sigma, nu_1 = nu_1,
         nu_0 = nu_0, lambda_1 = lambda_1, lambda_0 = lambda_0, rho_1 = rho_1,
         rho_0 = rho_0)
  }

  beta_0_pars <- function(pi_y) {
    # beta_0 as a function of all the parameter to the marginal distribution
    # of Y (pi_y), given all of the other parameters.

    pars <- pars_pi_y(pi_y)

    beta_0 <- log(pi_y / (1 - pi_y))

    for (j in which(xtype == "c_a")) {
      beta_0 <- beta_0 + 0.5 * (
        (pars$mu_0[j]**2 - pars$mu_1[j]**2) / pars$sigma[j]
      )
    }

    for (j in which(xtype == "c_p")) {
      beta_0 <- beta_0 + log(pars$nu_1[j] / pars$nu_0[j])
    }
    for (j in which(xtype == "d_c")) {
      beta_0 <- beta_0 + pars$lambda_0[j] - pars$lambda_1[j]
    }
    
    for (j in which(xtype == "d_b")) {
      beta_0 <- beta_0 + log(
        (1 + exp(qlogis(pars$rho_0[j]))) / (1 + exp(qlogis(pars$rho_1[j])))
      )
    }

    beta_0
  }
  
  opt <- optim(
    par = qlogis(mean(y)),
    function(t) sapply(
      t, function(s) abs(beta_0_pars(min(max(plogis(s),1e-15),1-1e-15)) - beta_0)
    ),
    method = "Brent", lower = -15, upper = 15
  )
  
  pi_y <- plogis(opt$par)
  pars <- pars_pi_y(pi_y)
  pars$pi_y <- pi_y
  pars
}

pars_to_linear_pred <- function(x, xtype, pars) {
  lp <- qlogis(pars$pi_y)
  
  for (j in which(xtype == "c_a")) {
    lp <- (lp + dnorm(x[, j], pars$mu_1[j], sqrt(pars$sigma[j]), TRUE) -
             dnorm(x[, j], pars$mu_0[j], sqrt(pars$sigma[j]), TRUE))
  }

  for (j in which(xtype == "c_p")) {
    lp <- (lp + dexp(x[, j], pars$nu_1[j], TRUE) -
             dnorm(x[, j], pars$nu_0[j], TRUE))
  }

  for (j in which(xtype == "d_p")) {
    lp <- (lp + dpois(x[, j], pars$lambda_1[j], TRUE) -
             dpois(x[, j], pars$lambda_0[j], TRUE))
  }

  for (j in which(xtype == "d_b")) {
    lp <- (lp + dbinom(x[, j], 1, pars$rho_1[j], TRUE) -
             dbinom(x[, j], 1, pars$rho_0[j], TRUE))
  }

  lp
}

parameters_to_intercept <- function(xtype, pars){

  if(any(xtype=="c_a")){
    cont_a <- sum(
      sapply(
        which(xtype=="c_a"),
        function(j) ((pars$mu_0[j]**2 -
                       pars$mu_1[j]**2) / pars$sigma[j])/2
      )
    )
  } else {
    cont_a <- 0
  }
  
  if(any(xtype=="c_p")){
    cont_p <- sum(
      sapply(
        which(xtype=="c_p"),
        function(j) log(pars$nu_1[j]) - log(pars$nu_0[j])
      )
    )
  } else {
    cont_p <- 0
  }

  if(any(xtype=="d_b")){
    disc_b <- sum(
      sapply(
        which(xtype=="d_b"),
        function(j) log(1 + exp(qlogis(pars$rho_0[j])))
        - log(1 + exp(qlogis(pars$rho_1[j])))
      )
    )
  } else {
    disc_b <- 0
  }

  if(any(xtype=="d_c")){
    disc_c <- sum(
      sapply(
        which(xtype=="d_c"),
        function(j) pars$lambda_0[j] - pars$lambda_1[j]
      )
    )
  } else {
    disc_c <- 0
  }

  qlogis(pars$pi_y) + cont_a + cont_p + disc_b + disc_c

}
