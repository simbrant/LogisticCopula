
comp_likelihood_model <- function(y, m_obj) {
  eta <- m_obj$beta_0 + m_obj$beta_x_insample + m_obj$copula_eff_insample
  prob <- plogis(eta)
  sum(dbinom(y, 1, prob, log = TRUE))
}

extract_var_type <- function(xtype) {
  stringr::str_split(xtype, "_")[[1]][1]
}

intercept_adjustment <- function(effect, l_p, y){
  
  opt <- nlm(
    function(adj) {
      ind <- which((plogis(l_p + effect + adj) > 0)&(plogis(l_p + effect + adj) < 1))
      if(length(ind) == 0)
      {
        1e15
      }
      else
      {
        -sum(
          dbinom(
            y[ind], 1, plogis(l_p + effect + adj)[ind], T
          ))
      }
      },
    p = 0)
  
  opt$estimate
}

choose_best_effect <- function(y, l_p, effects, adjust_intercept = TRUE) {
  if (adjust_intercept) {
    beta_0_adjustments <- try(apply(effects, 2, intercept_adjustment, l_p = l_p,
                                y = y),silent=TRUE)
    if(length(beta_0_adjustments) == 1)
    { 
      beta_0_adjustments <- rep(0, ncol(effects))
    }
  }
  else {
    beta_0_adjustments <- rep(0, ncol(effects))
  }
  
  scores <-  sapply(seq(ncol(effects)),
                      function(j){
                        sum(y*(beta_0_adjustments[j]+ effects[, j])
                            +log((1+exp(l_p))/(1+exp(l_p + beta_0_adjustments[j]+ effects[, j]))))
                      })
  which.max(scores)
}

transformed_delta_vec <- function(m_obj) {
  t_delta <- c()
  for (y_ind in seq(2)) {
    for (t in seq(length(m_obj$trees))) {
      if (length(m_obj$trees[[t]]$edges) > 0) {
        for (e_j in seq(length(m_obj$trees[[t]]$edges))) {
          t_delta <- c(
            t_delta, 
            inv_par_transform(
              m_obj$trees[[t]]$pair_copulas[[e_j]][[y_ind]]$par,
              m_obj$trees[[t]]$pair_copulas[[e_j]][[y_ind]]$family)
          )
        }
      }
    }
  }

  t_delta
}

inv_par_transform <- function(delta, family) {
  if (family == "gaussian") {
    tau <- (2 / pi) * asin(min(max(delta,-1+1e-8),1-1e-8))
    t <- log((tau + 1) / (1 - tau))
  } else if (family == "gumbel") {
    tau <- min(max(1 - 1 / (delta),1e-8),0.94-1e-8)
    t <- qlogis(tau  / 0.94)
  } else if (family == "clayton") {
    tau <- min(max((delta - 1e-10) / (delta - 1e-10 + 2),-1+1e-8),0.93-1e-8)
    t <- qlogis(tau  / 0.93)
  } else {
    stop(paste0("Transformation for family ", family, " not implemented"))
  }
  t
}

transform_t_to_delta_vec <- function(t_delta, m_obj) {
  delta <- c()
  cnt <- 1
  for (y_ind in seq(2)) {
    for (t in seq(length(m_obj$trees))) {
      if (length(m_obj$trees[[t]]$edges) > 0) {
        for (e_j in seq(length(m_obj$trees[[t]]$edges))) {
          delta <- c(
            delta,
            par_transform(t_delta[cnt], 
                          m_obj$trees[[t]]$pair_copulas[[e_j]][[y_ind]]$family
            )
          )
          cnt <- cnt + 1
        }
      }
    }
  }
  delta
}

par_transform <- function(t, family) {
  if (family == "gaussian") {
    tau <- min(max(2 * plogis(t) - 1, -1 + 1e-8), 1 - 1e-8)
    delta <- sin((pi / 2) * tau)
  } else if (family == "gumbel") {
    tau <- min(max(plogis(t) * 0.94, 1e-8), 0.94-1e-8)
    delta <- 1 / (1 - tau)
  } else if (family == "clayton") {
    tau <- min(max(plogis(t) * 0.93, 1e-8), 0.93-1e-8)
    delta <- 1e-10 + 2 * tau / (1 - tau)
  } else {
    stop(paste0("Transformation for family ", family, " not implemented"))
  }
  delta
}

set_deltas <- function(m_obj, deltas, transform = TRUE) {
  cnt <- 1
  for (y_ind in seq(2)) {
    for (t in seq(length(m_obj$trees))) {
      if (length(m_obj$trees[[t]]$edges) > 0 ){
        for (e_j in seq(length(m_obj$trees[[t]]$edges))) {
          if (transform) {
            m_obj$trees[[t]]$pair_copulas[[e_j]][[y_ind]]$parameters <- matrix(
              par_transform(
                deltas[cnt], m_obj$trees[[t]]$pair_copulas[[e_j]][[y_ind]]$family
              )
            )
          } else {
            m_obj$trees[[t]]$pair_copulas[[e_j]][[y_ind]]$parameters <- (
              matrix(deltas[cnt])
            )
          }
          cnt <- cnt + 1 
        }
      }
    }
  }
  m_obj
}

set_model <- function(y, x, m_obj, beta = NULL, delta = NULL) {
  
  if (is.null(beta) & is.null(delta)) {
    stop(
      "Must provide at least one of beta or delta values to function set_model"
    )
  }
  if (!is.null(beta)) {
    pars <- coefs_to_pars(y, x, m_obj$xtype, beta[1], beta[-1])
    m_obj$parameters <- pars
    m_obj$beta_0 <- beta[1]
    m_obj$beta_vec <- beta
    m_obj$beta_x_insample <- x %*% beta[-1]
  }
  if (!is.null(delta)) {
    m_obj <- set_deltas(m_obj, delta, transform = FALSE)
  }

  m_obj <- set_transformed_vars(x, m_obj$which_include, m_obj)
  m_obj$copula_eff_insample <- compute_g(m_obj)
  m_obj

}