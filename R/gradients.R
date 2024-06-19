derivative_transformation <- function(t, family) {
  if (family == "gaussian") {
    tau <- 2 * plogis(t) - 1
    d_delta_d_tau <- cos((pi / 2) * tau) * (pi / 2)
    d_tau_d_t <- 2 * plogis(t) * (1 - plogis(t))
  } else if (family == "gumbel") {
    tau <- plogis(t) * 0.94
    d_delta_d_tau <- 1 / (1 - tau)**2
    d_tau_d_t <- plogis(t) * (1 - plogis(t)) / 0.94
  } else if (family == "clayton") {
    tau <- plogis(t) * 0.93
    d_delta_d_tau <- 2 / (1 - tau)**2
    d_tau_d_t <- plogis(t) * (1 - plogis(t)) / 0.93
  } else {
    stop(paste0("Transformation for family ", family, " not implemented"))
  }
  
  d_delta_d_tau * d_tau_d_t
  
}

derivative_for_edge_e_j <- function(j_e_j, t, y_cond, fit) {
  
  e_j <- fit$trees[[t]]$edges[[j_e_j]]
  u_e_j <- get_u_upair_for_edge(e_j, fit$u, y_cond)
  # Compute the derivative at e_j
  derivative <- bicop_2_logcopdens_der(
    u = u_e_j, fit$trees[[t]]$pair_copulas[[j_e_j]][[y_cond + 1]], "par"
  )
  # Find the rest
  t_max <- length(fit$trees)
  if (t + 1 <= t_max) {
    for (s in seq(t+1, t_max)) {
      for (j_e_k in seq(1, length(fit$trees[[s]]$edges))) {
        # For each edge e_k in tree s, check if either node in e_k is 
        # e_j, or a parent of e_j
        e_k <- fit$trees[[s]]$edges[[j_e_k]]
        
        # Check whether to differentiate the log copula density of
        # e_k with respect to the first argument
        dir <- check_left_right(e_j, e_k)
        if (dir != "neither") {
          u_e_k <- get_u_upair_for_edge(e_k, fit$u, y_cond)
        }
        
        if (dir %in% c("left", "both")) {
          a <- bicop_2_logcopdens_der(
            u = u_e_k, 
            bicop_obj = fit$trees[[s]]$pair_copulas[[j_e_k]][[y_cond + 1]], 
            deriv = "u1" 
          )
          b <- diff_h_function_recursive(
            e_k, e_j, fit, j_e_k, s, t, y_cond, 1
          )
          contribution <- a * b
          derivative <- derivative + contribution
        }
        
        union_e_k_left <- union(e_k[[1]][[1]], e_k[[1]][[2]])
        
        
        
        if (dir %in% c("right", "both")) {
          contribution <- bicop_2_logcopdens_der(
            u = u_e_k, 
            bicop_obj = fit$trees[[s]]$pair_copulas[[j_e_k]][[y_cond + 1]], 
            deriv = "u2"
          ) * diff_h_function_recursive(
            e_k, e_j, fit, j_e_k, s, t, y_cond, 2
          )
          derivative <- derivative + contribution
        }
      }
    }
  }
  derivative
}

check_left_right <- function(e_j, e_k) {
  # Check whether e_j is an element of either the first or second
  # element of e_k, or both.
  
  union_e_j <- union(
    union(e_j[[1]][[1]], e_j[[1]][[2]]), 
    union(e_j[[2]][[1]], e_j[[2]][[2]])
  )
  union_e_k_left <- union(e_k[[1]][[1]], e_k[[1]][[2]])
  union_e_k_right <- union(e_k[[2]][[1]], e_k[[2]][[2]])
  
  if (all(union_e_j %in% union_e_k_left) & 
      all(union_e_j %in% union_e_k_right)) {
    "both"
  } else if (all(union_e_j %in% union_e_k_left)) {
    "left"
  } else if (all(union_e_j %in% union_e_k_right)) {
    "right"
  } else {
    "neither"
  }
  
}

diff_h_function_recursive <- function(e_k, e_j, fit, j_e_k, s, t, y_cond, 
                                      direction=1) {
  
  j_e_l <- which_edge_is_cond_node(
    e_k[[direction]], fit$trees[[s - 1]]$edges
  )
  e_l <- fit$trees[[s - 1]]$edges[[j_e_l]]
  u_e_l <- get_u_upair_for_edge(e_l, fit$u, y_cond)
  
  if (equal_nodes(e_k[[direction]], condition_node(e_l[[1]], e_l[[2]]))) {
    cond_var <- 2
  } else if (equal_nodes(e_k[[direction]], 
                         condition_node(e_l[[2]], e_l[[1]]))) {
    cond_var <- 1
  }
  
  if (equal_edges(e_j, e_l)) {
    derivative <- bicop_2_hbicop_der(
      u_e_l, 
      fit$trees[[s - 1]]$pair_copulas[[j_e_l]][[y_cond + 1]], 
      cond_var = cond_var, deriv = "par"
    )
    
  } else {
    if (check_left_right(e_j, e_l) == "left") {
      
      a <- bicop_2_hbicop_der(
        u_e_l, fit$trees[[s - 1]]$pair_copulas[[j_e_l]][[y_cond + 1]],
        cond_var, deriv = "u1"
      )
      
      b <- diff_h_function_recursive(
        e_l, e_j, fit, j_e_l, s-1, t, y_cond, 1
      )
      derivative <- a * b
    } else if (check_left_right(e_j, e_l) == "right") {
      derivative <- bicop_2_hbicop_der(
        u_e_l, fit$trees[[s - 1]]$pair_copulas[[j_e_l]][[y_cond + 1]], 
        direction, deriv = "u2"
      ) * diff_h_function_recursive(
        e_l, e_j, fit, j_e_l, s-1, t, y_cond, 2
      )
    } else if (check_left_right(e_j, e_l) == "both") {
      a <- bicop_2_hbicop_der(
        u_e_l, fit$trees[[s - 1]]$pair_copulas[[j_e_l]][[y_cond + 1]],
        cond_var, deriv = "u1"
      ) * diff_h_function_recursive(
        e_l, e_j, fit, j_e_l, s-1, t, y_cond, 1
      )
      b <- bicop_2_hbicop_der(
        u_e_l, fit$trees[[s - 1]]$pair_copulas[[j_e_l]][[y_cond + 1]], 
        direction, deriv = "u2"
      ) * diff_h_function_recursive(
        e_l, e_j, fit, j_e_l, s-1, t, y_cond, 2
      )
      derivative <- a + b
    }
  }
  
  derivative
}

which_edge_is_cond_node <- function(cond_node, edges) {
  # This is a helper function intended to find the index of the
  # edge in the tree on the previous level of the vine that corresponds to one
  # of the nodes of an edge
  which(
    sapply(
      edges, 
      function(e_l) {
        (equal_nodes(condition_node(e_l[[1]], e_l[[2]]), cond_node) |
           equal_nodes(condition_node(e_l[[2]], e_l[[1]]), cond_node))
      }
    )
  )
}

partial_pi_y_partial_beta_k <- function(k, pars, beta, xtype, xbar, S2) {
  
  denom_a <- 1 / (pars$pi_y * (1 - pars$pi_y))
  if (any(xtype == "c_a")) {

    denom_b <- sum(
      sapply(which(xtype == "c_a"), 
             function(j) {
               B_j <- (
                 ((1 - 2*pars$pi_y) / (pars$pi_y * (1 - pars$pi_y))) * 
                   (pars$sigma[j] - S2[j]) / 
                   (1 + 2 * beta[j]**2 * pars$pi_y * 
                      (1 - pars$pi_y) * pars$sigma[j])
               )
               (beta[j]**2 / 2) * (2 * pars$sigma[j] - B_j * 
                                     (1 - 2 * pars$pi_y))
             }
      )
    )
  } else {
    denom_b <- 0
  }
  
  if (any(xtype == "d_b")) {
    denom_c <- sum(
      sapply(
        which(xtype == "d_b"), 
        function(j) {
          num_2 <- (1 - exp(beta[j])) * (pars$rho_0[j] - pars$rho_0[j]**2)
          den_2 <- 1 + (1 - exp(beta[j])) * (
            xbar[j] - pars$pi_y - 2 * pars$rho_0[j] * (1 - pars$pi_y)
          )
          D_jk_0 <- num_2 / den_2
          D_jk_1 <- (exp(beta[j]) /
                       (1 - (1 - exp(beta[j]))*pars$rho_0[j])**2) * D_jk_0
          D_jk_0 / (1 - pars$rho_0[j]) - D_jk_1 / (1 - pars$rho_1[j]) 
        }
      )
    )
  } else {
    denom_c <- 0
  }
  
  if (any(xtype == "d_c")) {
    denom_d <- sapply(
      which(xtype == "d_c"), 
      function(j) {
        F_jk_0 <- (pars$lambda_0[j] * (1 - exp(beta[j])) / 
                     (1 - pars$pi_y * (1 - exp(beta[j]))))
        F_jk_1 <- (pars$lambda_1[j] * (1 - exp(beta[j])) / 
                     (1 - pars$pi_y * (1 - exp(beta[j]))))
        F_jk_0 - F_jk_1
      }
    )
  } else {
    denom_d <- 0
  }
  
  denom <- denom_a + denom_b + denom_c + denom_d
  
  if (k == 0) {
    nom <- 1    
  } else if (xtype[k] == "c_a") {
    
    A_k <- (
      (2 / beta[k]) * (pars$sigma[k] - S2[k]) / 
        (1 + 2 * beta[k]**2 * pars$pi_y * (1 - pars$pi_y) * pars$sigma[k])
    )
    
    nom <- (
      xbar[k] + beta[k] * (1 - 2 * pars$pi_y) * pars$sigma[k] + 
        0.5 * beta[k]**2 * (1 - 2 * pars$pi_y) * A_k
    )
    
  } else if (xtype[k] == "d_b") {
    num_1 <- exp(beta[k]) * (
      pars$rho_0[k] * (xbar[k] - pars$pi_y) - 
        pars$rho_0[k]**2 * (1 - pars$pi_y)
    )
    den_1 <- 1 + (1 - exp(beta[k])) * (
      xbar[k] - pars$pi_y - 2 * pars$rho_0[k] * (1 - pars$pi_y)
    )
    C_kk_0 <- num_1 / den_1
    C_kk_1 <- (exp(beta[k]) / (1 - (1 - exp(beta[k])) * pars$rho_0[k])**2) * (
      C_kk_0 + pars$rho_0[k] - pars$rho_0[k]**2
    )
    nom <- C_kk_1 / (1 - pars$rho_1[k]) - C_kk_0 / (1 - pars$rho_0[k])
  } else if (xtype[k] == "d_c") {
    E_kk_1 <- pars$lambda_1[k]*((1 - pars$pi_y) / 
                                  (1 - pars$pi_y * (1 - exp(beta[k]))))
    E_kk_0 <- (-pars$lambda_0[k] * exp(beta[k]) * pars$pi_y / 
                 (1 - pars$pi_y * (1 - exp(beta[k]))))
    nom <- E_kk_1 - E_kk_0
  }
  nom / denom
}


partial_sigma_j_wrt_beta_k <- function(j, k, pars, beta, xtype, xbar, S2) {

  Q_j <- beta[j]**2 * pars$pi_y * (1 - pars$pi_y)
  d_sigma_jd_Q_j <- (1 / Q_j) * (
    S2[j] / (2 * Q_j * pars$sigma[j] + 1) - pars$sigma[j]
  )

  d_Q_d_beta_k <- (
    beta[j]**2 * (1 - 2 * pars$pi_y) * 
      partial_pi_y_partial_beta_k(
        k, pars, beta, xtype, xbar, S2
      )
  )

  if (j == k) {
    d_Q_d_beta_k <- d_Q_d_beta_k + 2 * beta[j] * pars$pi_y * (1 - pars$pi_y)
  }

  d_sigma_jd_Q_j * d_Q_d_beta_k

}

partial_mu_j_y_wrt_beta_k <- function(j, k, y_cond, pars, beta, xtype, 
                                      xbar, S2) {
  
  if(y_cond == 0) {
    
    der <- (
      -beta[j] * partial_sigma_j_wrt_beta_k(
        j, k, pars, beta, xtype, xbar, S2
      ) * pars$pi_y - 
        beta[j] * pars$sigma[j] * partial_pi_y_partial_beta_k(
          k, pars, beta, xtype, xbar, S2
        )
    )

    if (j == k) {
      der <- der - pars$sigma[j] * pars$pi_y
    }

  } else {
    der <- (
      beta[j] * partial_sigma_j_wrt_beta_k(j, k, pars, beta, xtype, xbar, S2) * 
        (1 - pars$pi_y) - 
        beta[j] * pars$sigma[j] * partial_pi_y_partial_beta_k(
          k, pars, beta, xtype, xbar, S2
        )
    )
    if (j == k) {
      der <- der + pars$sigma[j] * (1 - pars$pi_y)
    }
  }
  der
}

partial_rho_0j_wrt_beta_k <- function(j, k, pars, beta, xtype, xbar, S2) {
  
  partial_pi_y_wrt_beta_k <- partial_pi_y_partial_beta_k(
    k, pars, beta, xtype, xbar, S2
  )
  
  if (j == k) {
    num_1 <- exp(beta[j]) * (
      pars$rho_0[j] * (xbar[j] - pars$pi_y) - 
        pars$rho_0[j]**2 * (1 - pars$pi_y)
    )
    den_1 <- 1 + (1 - exp(beta[j])) * (
      xbar[j] - pars$pi_y - 2 * pars$rho_0[j] * (1 - pars$pi_y)
    )
    C_jk_0 <- num_1 / den_1
  } else {
    C_jk_0 <- 0
  }
  num_2 <- (1 - exp(beta[j])) * (pars$rho_0[j] - pars$rho_0[j]**2)
  den_2 <- 1 + (1 - exp(beta[j])) * (
    xbar[j] - pars$pi_y - 2 * pars$rho_0[j] * (1 - pars$pi_y)
  )
  D_jk_0 <- num_2 / den_2
  
  C_jk_0 + D_jk_0 * partial_pi_y_wrt_beta_k
}

partial_rho_1j_wrt_beta_k <- function(j, k, pars, beta, xtype, xbar, S2) {
  
  partial_pi_y_wrt_beta_k <- partial_pi_y_partial_beta_k(
    k, pars, beta, xtype, xbar, S2
  )
  
  if (j == k) {
    num_1 <- exp(beta[j]) * (
      pars$rho_0[j] * (xbar[j] - pars$pi_y) - 
        pars$rho_0[j]**2 * (1 - pars$pi_y)
    )
    den_1 <- 1 + (1 - exp(beta[j])) * (
      xbar[j] - pars$pi_y - 2 * pars$rho_0[j] * (1 - pars$pi_y)
    )
    C_jk_0 <- num_1 / den_1
    C_jk_1 <- (exp(beta[j]) / (1 - (1 - exp(beta[j])) * pars$rho_0[j])**2) * (
      C_jk_0 + pars$rho_0[j] - pars$rho_0[j]**2
    )
  } else {
    C_jk_1 <- 0
  }
  
  num_2 <- (1 - exp(beta[j])) * (pars$rho_0[j] - pars$rho_0[j]**2)
  den_2 <- 1 + (1 - exp(beta[j])) * (
    xbar[j] - pars$pi_y - 2 * pars$rho_0[j] * (1 - pars$pi_y)
  )
  D_jk_0 <- num_2 / den_2
  D_jk_1 <- (exp(beta[j]) / (1 - (1 - exp(beta[j]))*pars$rho_0[j])**2) * D_jk_0
  
  C_jk_1 + D_jk_1 * partial_pi_y_wrt_beta_k
}

partial_lambda_0j_wrt_beta_k <- function(j, k, pars, beta, xtype, xbar, S2) {
  
  partial_pi_y_wrt_beta_k <- partial_pi_y_partial_beta_k(
    k, pars, beta, xtype, xbar, S2
  )
  
  if (j == k) {
    E_jk_0 <- (-pars$lambda_0[j] * exp(beta[j]) * pars$pi_y / 
                 (1 - pars$pi_y * (1 - exp(beta[j]))))
  } else {
    E_jk_0 <- 0
  }
  F_jk_0 <- (pars$lambda_0[j] * (1 - exp(beta[j])) / 
               (1 - pars$pi_y * (1 - exp(beta[j]))))
  
  E_jk_0 + F_jk_0 * partial_pi_y_wrt_beta_k
  
}

partial_lambda_1j_wrt_beta_k <- function(j, k, pars, beta, xtype, xbar, S2) {
  
  partial_pi_y_wrt_beta_k <- partial_pi_y_partial_beta_k(
    k, pars, beta, xtype, xbar, S2
  )
  
  if (j == k) {
    E_jk_1 <- pars$lambda_1[j]*((1 - pars$pi_y) / 
                                  (1 - pars$pi_y * (1 - exp(beta[j]))))
  } else {
    E_jk_1 <- 0
  }
  F_jk_1 <- (pars$lambda_1[j] * (1 - exp(beta[j])) / 
               (1 - pars$pi_y * (1 - exp(beta[j]))))
  
  E_jk_1 + F_jk_1 * partial_pi_y_wrt_beta_k
  
}

partial_F_wrt_beta <- function(j, k, y_cond, pars, beta, xtype, x) {
  
  if (xtype[j] == "c_a") {
    dmu <- partial_mu_j_y_wrt_beta_k(j, k, y_cond, pars, beta, xtype, apply(x, 2, mean), apply(x, 2, var))
    dsigma <- partial_sigma_j_wrt_beta_k(j, k, pars, beta, xtype, apply(x, 2, mean), apply(x, 2, var))
    if (y_cond == 1) {
      
      -dnorm((x[, j] - pars$mu_1[j]) / sqrt(pars$sigma[j]), 0, 1) * (
        dmu * pars$sigma[j] ** (- 1 / 2) + 
          0.5 * (x[, j] - pars$mu_1[j]) * (
            dsigma * pars$sigma[j] ** (- 3 / 2)
          )
      )
    } else {
      -dnorm((x[, j] - pars$mu_0[j]) / sqrt(pars$sigma[j]), 0, 1) * (
        dmu * pars$sigma[j] ** (- 1 / 2) + 
          0.5 * (x[, j] - pars$mu_0[j]) * (
            dsigma * pars$sigma[j] ** (- 3 / 2)
          )
      )
    }
  } else {
    stop("Margin derivatives w.r.t. betas only implemented for gaussian margins")
  }
}

beta_gradient <- function(y, x, fit, beta_0, betas, update = TRUE) {
  
  if (update) {
    #
    # Update model using the provided beta0 and beta
    #
    fit <- set_model(y, x, fit, beta = c(beta_0, betas))
  }
  eta <- fit$beta_0 + fit$beta_x_insample + fit$copula_eff_insample
  p <- plogis(eta)
  ystar <- (y - p)
  
  sapply(seq(0, ncol(x)), function(k) {
    if (!is.null(fit$trees[[1]]$edges)) {
      copder1 <- derivative_of_log_vine_density_wrt_beta(k, fit, 1, x, betas)
      copder0 <- derivative_of_log_vine_density_wrt_beta(k, fit, 0, x, betas)
    } else {
      copder1 <- 0
      copder0 <- 0
    }
    
    if (k == 0) {
      sum(ystar * (1 + copder1 - copder0))
    } else {
      sum(ystar * (x[, k] + copder1 - copder0))
    }
  })
  
}

derivative_of_h_function_wrt_beta <- function(k, node, fit, y_cond, x, t_prev,
                                              betas) {
  
  edge_ind <- which_edge(node, fit$trees[[t_prev - 1]]$edges)
  edge <- fit$trees[[t_prev - 1]]$edges[[edge_ind]]
  
  cond_var <- which(
    c(equal_nodes(node, condition_node(edge[[2]], edge[[1]])), 
      equal_nodes(node, condition_node(edge[[1]], edge[[2]])))
  )
  
  u_pair <- get_u_upair_for_edge(
    edge, fit$u, y_cond
  )
  
  
  if (t_prev - 1 == 1) {
    
    a <- predict(
      fit$trees[[t_prev - 1]]$pair_copulas[[edge_ind]][[y_cond + 1]],
      newdata = weave_transformed(u_pair$u1, u_pair$u2)
    ) * partial_F_wrt_beta(
      edge[[c(2, 1)[cond_var]]]$vertice,
      k, y_cond, fit$parameters, betas, fit$xtype, x
    )
    
  } else {
    
    a <- predict(
      fit$trees[[t_prev - 1]]$pair_copulas[[edge_ind]][[y_cond + 1]],
      newdata = weave_transformed(u_pair$u1, u_pair$u2)
    ) * derivative_of_h_function_wrt_beta(
      k, edge[[c(2, 1)[cond_var]]], fit, y_cond, x, t_prev - 1, betas
    )
  }
  
  if (t_prev - 1 == 1) {
    
    b <- bicop_2_hbicop_der(
      u_pair, fit$trees[[t_prev - 1]]$pair_copulas[[edge_ind]][[y_cond + 1]],
      cond_var = cond_var, deriv = "u2"
    ) * partial_F_wrt_beta(
      edge[[c(1, 2)[cond_var]]]$vertice,
      k, y_cond, fit$parameters, betas, fit$xtype, x
    )
    
  } else {
    b <- bicop_2_hbicop_der(
      u_pair, fit$trees[[t_prev - 1]]$pair_copulas[[edge_ind]][[y_cond + 1]],
      cond_var = cond_var
    ) * derivative_of_h_function_wrt_beta(
      k, edge[[c(1, 2)[cond_var]]], fit, y_cond, x, t_prev - 1, betas
    )
  }
  return(a + b)
}

derivative_of_log_vine_density_wrt_beta <- function(k, fit, y_cond, x, betas) {

  deriv <- 0
  for (t in seq(length(fit$trees))) {
    for (e_i in seq(length(fit$trees[[t]]$edges))) {
      e <- fit$trees[[t]]$edges[[e_i]]
      u_pair <- get_u_upair_for_edge(e, fit$u, y_cond)
      if (t > 1) {
        a <- bicop_2_logcopdens_der(
          u_pair, fit$trees[[t]]$pair_copulas[[e_i]][[y_cond + 1]], deriv = "u1"
        ) * derivative_of_h_function_wrt_beta(k, e[[1]], fit, y_cond, x, t, 
                                              betas)
        
        b <- bicop_2_logcopdens_der(
          u_pair, fit$trees[[t]]$pair_copulas[[e_i]][[y_cond + 1]], deriv = "u2"
        ) * derivative_of_h_function_wrt_beta(k, e[[2]], fit, y_cond, x, t, 
                                              betas)
        deriv <- deriv + a + b
      } else {
        a <- bicop_2_logcopdens_der(
          u_pair, fit$trees[[t]]$pair_copulas[[e_i]][[y_cond + 1]],
          deriv = "u1"
        ) * partial_F_wrt_beta(
          e[[1]]$vertice, k, y_cond, fit$parameters, betas, fit$xtype, x
        )
        
        b <- bicop_2_logcopdens_der(
          u_pair, fit$trees[[t]]$pair_copulas[[e_i]][[y_cond + 1]], 
          deriv = "u2"
        ) * partial_F_wrt_beta(
          e[[2]]$vertice, k, y_cond, fit$parameters, betas, fit$xtype, x
        )
        deriv <- deriv + a + b
      }
    }    
  }
  deriv
}

copula_gradient <- function(y, x, fit, delta, transformation = TRUE,
                            update = TRUE) {
  
  if (transformation) {
    t_delta <- delta
    delta <- transform_t_to_delta_vec(t_delta, fit)
  }

  if (update) {
    fit <- set_model(y, x, fit, delta = delta)
  }

  eta <- fit$beta_0 + fit$beta_x_insample + fit$copula_eff_insample
  p <- plogis(eta)
  ystar <- y - p
  grad <- c()
  cnt <- 1
  for (t in seq(length(fit$trees))) {
    for (edge_index in seq(length(fit$trees[[t]]$edges))) {
      if (!transformation) {
        grad <- c(
          grad, -sum(derivative_for_edge_e_j(edge_index, t, 0, fit) * ystar)
        )
      } else {
        grad <- c(
          grad, -sum(
            derivative_for_edge_e_j(edge_index, t, 0, fit) * ystar
          ) * derivative_transformation(
            t_delta[cnt],
            family = fit$trees[[t]]$pair_copulas[[edge_index]][[1]]$family
          ) 
        )
      }
      cnt <- cnt + 1
    }
  }
  
  for (t in seq(length(fit$trees))) {
    for (edge_index in seq(length(fit$trees[[t]]$edges))) {
      if (!transformation) {
        grad <- c(
          grad, sum(derivative_for_edge_e_j(edge_index, t, 1, fit) * ystar)
          )
      } else {
        grad <- c(
          grad, sum(
            derivative_for_edge_e_j(edge_index, t, 1, fit) * ystar
          ) * derivative_transformation(
            t_delta[cnt], 
            family = fit$trees[[t]]$pair_copulas[[edge_index]][[2]]$family
          )
        )
      }
      cnt <- cnt + 1
    }
  }

  grad

}


full_gradient <- function(y, x, fit, beta_0, betas, t_delta) {
  
  delta <- transform_t_to_delta_vec(t_delta, fit)
  fit <- set_model(y, x, fit, c(beta_0, betas), delta)
  eta <- fit$beta_0 + fit$beta_x_insample + fit$copula_eff_insample
  p <- plogis(eta)
  ystar <- (y - p)
  
  c(beta_gradient(y, x, fit, beta_0, betas, update = FALSE),
    copula_gradient(y, x, fit, t_delta, transformation = TRUE, update = FALSE))

}


