
transform_margins_to_hash <- function(x, xtype, pars, which_include) {
  
  u_1 <- new.env(hash = TRUE)
  u_0 <- new.env(hash = TRUE)
  if(is.null(which_include)) {
    which_include = seq(length(xtype))
  }
  
  for (j in which_include) {
    
    if (xtype[j] == "c_a") {
      
      assign(cond_node_key(CondNode(j)),
             pnorm(x[, j], pars$mu_0[j], sqrt(pars$sigma[j])), envir = u_0)
      assign(cond_node_key(CondNode(j)),
             pnorm(x[, j], pars$mu_1[j], sqrt(pars$sigma[j])), envir = u_1)
    } else if (xtype[j] == "c_p") {
      
      assign(cond_node_key(CondNode(j)),
             pexp(x[, j], pars$nu_0[j]), envir = u_0)
      assign(cond_node_key(CondNode(j)),
             pexp(x[, j], pars$nu_1[j]), envir = u_1)
    } else if (xtype[j] == "d_b") {
      
      assign(cond_node_key(CondNode(j)),
             cbind(pbinom(x[, j], 1,
                          pars$rho_0[j]),
                   pbinom(x[, j] - 1, 1,
                          pars$rho_0[j])),
             envir = u_0)
      assign(cond_node_key(CondNode(j)),
             cbind(pbinom(x[, j], 1,
                          pars$rho_1[j]),
                   pbinom(x[, j] - 1, 1,
                          pars$rho_1[j])),
             envir = u_1)
    } else if (xtype[j] == "d_c") {
      
      assign(cond_node_key(CondNode(j)),
             cbind(ppois(x[, j], 1, pars$lambda_0[j]),
                   ppois(x[, j] - 1, 1, pars$lambda_0[j])),
             envir = u_0)
      assign(cond_node_key(CondNode(j)),
             cbind(ppois(x[, j], 1, pars$lambda_1[j]),
                   ppois(x[, j] - 1, 1, pars$lambda_1[j])),
             envir = u_1)
    }
  }
  
  list(u_0 = u_0, u_1 = u_1)
  
}


set_transformed_vars <- function(x, which_include, fit) {
  
  #
  # Computes the values of all the conditional U's for the margins and all
  # parts of the copula.
  #
  
  # First compute the u's for the margins
  u_new <- transform_margins_to_hash(
    x, fit$xtype, fit$parameters, which_include
  )

  # Then, for every edge
  for (t in seq(length(fit$trees))) {
    if (!is.null(fit$trees[[t]]$edges)) {
      for (e in seq(length(fit$trees[[t]]$edges))) {
        u_new <- compute_u_for_edge(
          u_new, fit$trees[[t]]$edges[[e]],
          fit$trees[[t]]$pair_copulas[[e]], fit$xtype
        )
      }
    }
  }

  # Update the 'fit' object
  fit$u <- u_new

  # Return fit 'object'
  fit
}

compute_u_for_edge <- function(u_hash, edge, copula_pair, xtype) {
  
  u0 <- weave_transformed(
    get(cond_node_key(edge[[1]]), envir = u_hash$u_0),
    get(cond_node_key(edge[[2]]), envir = u_hash$u_0)
  )
  
  u1 <- weave_transformed(
    get(cond_node_key(edge[[1]]), envir = u_hash$u_1),
    get(cond_node_key(edge[[2]]), envir = u_hash$u_1)
  )
  
  assign(cond_node_key(condition_node(edge[[1]], edge[[2]])),
         bicop_2_hbicop(u0, copula_pair[[1]],
                        return_u_minus = if(xtype[edge[[1]]$vertice] 
                                            %in% c("d_b", "d_c")){
                          TRUE
                        } else {
                          FALSE
                        }), envir = u_hash$u_0)
  
  assign(cond_node_key(condition_node(edge[[2]], edge[[1]])),
         bicop_2_hbicop(u0, copula_pair[[1]], cond_var = 1,
                        return_u_minus = if(xtype[edge[[2]]$vertice] 
                                            
                                            %in% c("d_b", "d_c")){
                          TRUE
                        } else {
                          FALSE
                        }), envir = u_hash$u_0)
  
  assign(cond_node_key(condition_node(edge[[1]], edge[[2]])),
         bicop_2_hbicop(u1, copula_pair[[2]],
                        return_u_minus = if(xtype[edge[[1]]$vertice] %in% 
                                            c("d_b", "d_c")){
                          TRUE
                        } else {
                          FALSE
                        }),
         envir = u_hash$u_1)
  
  assign(cond_node_key(condition_node(edge[[2]], edge[[1]])),
         bicop_2_hbicop(u1, copula_pair[[2]], cond_var = 1,
                        return_u_minus = if(xtype[edge[[2]]$vertice] 
                                            %in% c("d_b", "d_c")){
                          TRUE
                        } else {
                          FALSE
                        }), envir = u_hash$u_1)
  
  u_hash
  
}


weave_transformed <- function(u1, u2) {
  
  # This function weaves two transformed variables (so that the resulting
  # variable has the form [u+, u-], if any of the two margins are discrete,
  # u- only contains columns from discrete variables),
  # where u_j will be a 2xn matrix if the conditioned variable
  # in the j-th transformed variable is discrete, and a numeric vector
  # of length n if continuous.
  
  u <- cbind(u1, u2)
  
  if (is.null(ncol(u1)) || ncol(u1) == 1) {
    u
  } else if (is.null(ncol(u2)) || ncol(u2) == 1) {
    u[, c(1, 3, 2)]
  } else {
    u[, c(1, 3, 2, 4)]
  }
}

bicop_2_hbicop <- function(u, bicop_obj, cond_var=2, return_u_minus=F) {
  
  # Function that lets you compute h-functions, without having to
  # specify the u_1^- when computing C(u_1 | u_2), when u_1 is
  # discrete. This is because u_1^- is redundant in that case, but rvinecopulibs
  # hbicop() function demands that it is provided. In addition, the specifying
  # return_u_minus = T, will make the function output the n x 2 matrix
  # [C(u_2 | u_1), C(u_2^- | u_1)] if cond_var = 1, or
  # [C(u_1 | u_2), C(u_1^- | u_2)] if cond_var = 2.
  
  if (!bicop_obj$var_types[c(2, 1)[cond_var]] == "d") {
    # In this case (that the conditioned variable is continuous), the
    # h-function in rvinecopulib can be called directly, as it then behaves as
    # expected
    
    if (return_u_minus) {
      cbind(rvinecopulib::hbicop(u, cond_var = cond_var,
                                 family = bicop_obj$family,
                                 rotation = bicop_obj$rotation,
                                 parameters = bicop_obj$parameters,
                                 var_types = bicop_obj$var_types),
            rvinecopulib::hbicop(u, cond_var = cond_var,
                                 family = bicop_obj$family,
                                 rotation = bicop_obj$rotation,
                                 parameters = bicop_obj$parameters,
                                 var_types = bicop_obj$var_types))
    } else {
      rvinecopulib::hbicop(u, cond_var = cond_var, family = bicop_obj$family,
                           rotation = bicop_obj$rotation,
                           parameters = bicop_obj$parameters,
                           var_types = bicop_obj$var_types)
    }
  } else {
    
    # This is more complicated. There are four_cases.
    
    u_columns_1 <- switch(1 + 2 * (cond_var - 1) + 1 * (ncol(u) == 4),
                          c(1, 2, 1, 3),
                          c(1, 2, 3, 4),
                          c(1, 2, 3, 2),
                          c(1, 2, 3, 4))
    
    if (return_u_minus) {
      
      u_columns_2 <- switch(1 + 2 * (cond_var - 1) + 1 * (ncol(u) == 4),
                            c(1, 3, 1, 3),
                            c(1, 4, 3, 4),
                            c(3, 2, 3, 2),
                            c(3, 2, 3, 4))
      
      cbind(rvinecopulib::hbicop(u[, u_columns_1], cond_var = cond_var,
                                 family = bicop_obj$family,
                                 rotation = bicop_obj$rotation,
                                 parameters = bicop_obj$parameters,
                                 var_types = bicop_obj$var_types),
            rvinecopulib::hbicop(u[, u_columns_2], cond_var = cond_var,
                                 family = bicop_obj$family,
                                 rotation = bicop_obj$rotation,
                                 parameters = bicop_obj$parameters,
                                 var_types = bicop_obj$var_types))
      
    } else {
      rvinecopulib::hbicop(u[, u_columns_1], cond_var = cond_var,
                           family = bicop_obj$family,
                           rotation = bicop_obj$rotation,
                           parameters = bicop_obj$parameters,
                           var_types = bicop_obj$var_types)
    }
  }
}

compute_g <- function(fit) {
  
  copula_eff <- NULL
  for (t in seq(length(fit$trees))) {
    if (!is.null(fit$trees[[t]]$edges)) {
      for (e in seq(length(fit$trees[[t]]$edges))) {
        left_node <- cond_node_key(
          fit$trees[[t]]$edges[[e]][[1]]
        )
        right_node <- cond_node_key(
          fit$trees[[t]]$edges[[e]][[2]]
        )
        
        u_0 <- weave_transformed(
          get(left_node, envir = fit$u$u_0),
          get(right_node, envir = fit$u$u_0)
        )
        u_1 <- weave_transformed(
          get(left_node, envir = fit$u$u_1),
          get(right_node, envir = fit$u$u_1)
        )
        if(is.null(copula_eff)) {
          copula_eff <- (
            log(predict(fit$trees[[t]]$pair_copulas[[e]][[2]], newdata = u_1)) - 
              log(predict(fit$trees[[t]]$pair_copulas[[e]][[1]], newdata = u_0))
          )
        } else {
          copula_eff <- copula_eff + (
            log(predict(fit$trees[[t]]$pair_copulas[[e]][[2]], newdata = u_1)) - 
              log(predict(fit$trees[[t]]$pair_copulas[[e]][[1]], newdata = u_0))
          )
        }
      }
    }
  }
  if (is.null(copula_eff)) {
    copula_eff <- 0
  }

  copula_eff

}

string_to_VineCop_num <- function(copname, rotation) {
  if (copname == "gaussian") {
    1
  } else if (copname == "gumbel") {
    if (rotation == 0) {
      4
    } else if(rotation == 90) {
      24
    } else if (rotation == 180) {
      14
    } else if (rotation == 270) {
      34
    }
  } else if (copname == "clayton") {
    if (rotation == 0) {
      3
    } else if (rotation == 90) {
      23
    } else if (rotation == 180) {
      13
    } else if (rotation == 270) {
      33
    }
  } else {
    stop(paste0("Family '", copname, "' not implemented!"))
  }
}


bicop_2_hbicop_der <- function(u, bicop_obj, cond_var = 2, deriv = "par") {
  
  #
  # Derivative of a hbicop with respect to the parameter, or the conditioning
  # variable.
  #
  
  if (any(bicop_obj$var_types == "d")) {
    stop("Derivatives for discrete pair copulas not yet implemented")
  }
  
  if ((cond_var == 2 & deriv == "u1") | (cond_var == 1 & deriv == "u2")) {
    # Return density
    predict(bicop_obj, newdata = as.matrix(cbind(u$u1, u$u2)), what = "pdf")
  } else {
    
    rot <- 1 * (bicop_obj$rotation %in% c(90, 270))
    if (cond_var == 2) {
      VineCopula::BiCopHfuncDeriv(
        u$u1, u$u2, string_to_VineCop_num(bicop_obj$family, bicop_obj$rotation), 
        (-1)**rot * bicop_obj$par, deriv = deriv
      )
    } else if (cond_var == 1 & deriv == "u1") {
      # For some reason, one cannot express the second derivative of C(u_1, u_2)
      # wrt u_1 in terms of the VineCopula::BiCopHfuncDeriv function for the given
      # rotation directly, so this annoying workaround is needed.
      if (bicop_obj$rotation == 90) {
        -VineCopula::BiCopHfuncDeriv(
          u$u2, 1 - u$u1, 
          string_to_VineCop_num(bicop_obj$family, bicop_obj$rotation) %% 10,
          bicop_obj$par, deriv = "u2"
        )
      } else if (bicop_obj$rotation == 180) {
        # d_u_1 180 degrees 
        VineCopula::BiCopHfuncDeriv(
          1 - u$u2, 1 - u$u1, 
          string_to_VineCop_num(bicop_obj$family, bicop_obj$rotation) %% 10, 
          bicop_obj$par, deriv = "u2"
        )
      } else if (bicop_obj$rotation == 270) {
        # d u_1 270 degrees
        -VineCopula::BiCopHfuncDeriv(
          1 - u$u2, u$u1, 
          string_to_VineCop_num(bicop_obj$family, bicop_obj$rotation) %% 10, 
          bicop_obj$par, deriv = "u2"
        )
      } else {
        # d_u_2 0 degrees 
        VineCopula::BiCopHfuncDeriv(
          u$u2, u$u1,  
          string_to_VineCop_num(bicop_obj$family, bicop_obj$rotation), 
          bicop_obj$par, deriv = "u2"
        )
      }
    } else if (cond_var == 1 & deriv == "par") {
      # This makes even less sense than the previous workaround, bug in 
      # VineCopula::BicopHfuncDeriv?
      if (bicop_obj$rotation == 90) {
        -VineCopula::BiCopHfuncDeriv(
          u$u2, 1 - u$u1,
          string_to_VineCop_num(bicop_obj$family, bicop_obj$rotation),
          -bicop_obj$par, deriv = deriv
        )
      } else if (bicop_obj$rotation == 180) {
        VineCopula::BiCopHfuncDeriv(
          u$u2, u$u1,
          string_to_VineCop_num(bicop_obj$family, bicop_obj$rotation),
          bicop_obj$par, deriv = deriv
        )
      } else if (bicop_obj$rotation == 270) {
        # d u_1 270 degrees
        -VineCopula::BiCopHfuncDeriv(
          1 - u$u2, 1 - u$u1,
          string_to_VineCop_num(bicop_obj$family, bicop_obj$rotation),
          -bicop_obj$par, deriv = deriv
        )
      } else {
        # d_u_2 0 degrees
        VineCopula::BiCopHfuncDeriv(
          u$u2, u$u1,
          string_to_VineCop_num(bicop_obj$family, bicop_obj$rotation),
          bicop_obj$par, deriv = deriv
        )
      }
    }
  }
}

bicop_2_logcopdens_der <- function(u, bicop_obj, deriv = "par") {
  
  #
  # Derivative of a log-copula density with respect to the parameter
  #
  
  if (any(bicop_obj$var_types == "d")) {
    stop("Derivatives for discrete pair copulas not yet implemented")
  }
  rot <- 1 * (bicop_obj$rotation %in% c(90, 270))
  if (deriv == "par") {
    VineCopula::BiCopDeriv(
      u$u1, u$u2, string_to_VineCop_num(bicop_obj$family, bicop_obj$rotation),
      (-1)**rot * bicop_obj$par, log = TRUE
    )
  } else if (deriv %in% c("u1", "u2")) {
    VineCopula::BiCopDeriv(
      u$u1, u$u2, family = string_to_VineCop_num(bicop_obj$family, bicop_obj$rotation),
      par = (-1)**rot * bicop_obj$parameters, deriv = deriv
    ) / predict(bicop_obj, newdata = weave_transformed(u$u1, u$u2))
  } else {
    stop("Erroneous 'deriv' argument in bicop_2_logcopdens_der.")
  }
}

which_edge <- function(node, edges) {
  which(
    sapply(edges, function(edge)
      equal_nodes(node, condition_node(edge[[1]], edge[[2]]))
      | equal_nodes(node, condition_node(edge[[2]], edge[[1]])))
  )
}

get_u_upair_for_edge <- function(edge, u, y_cond) {
  list(
    u1 = get(cond_node_key(edge[[1]]), envir = if(y_cond == 0) u$u_0 else u$u_1),
    u2 = get(cond_node_key(edge[[2]]), envir = if(y_cond == 0) u$u_0 else u$u_1)
  )
}
fit_pair <- function(edge, u, y, xtype, family_pair=c("gaussian", "gaussian")) {
  cop_1 <- rvinecopulib::bicop(
    weave_transformed(
      get(cond_node_key(edge[[1]]), envir = u$u_0),
      get(cond_node_key(edge[[2]]), envir = u$u_0)
    )[y == 0, ],
    var_types = c(extract_var_type(xtype[edge[[1]]$vertice]),
                  extract_var_type(xtype[edge[[2]]$vertice])),
    family_set = family_pair[1], par_method = "mle"
  )
  if(cop_1$family == "indep")
  {
    cop_1$family <- "gaussian"
    cop_1$parameters <- matrix(0,1,1)
    cop_1$npars <- 1
  }
  cop_2 <- rvinecopulib::bicop(
    weave_transformed(
      get(cond_node_key(edge[[1]]), envir = u$u_1),
      get(cond_node_key(edge[[2]]), envir = u$u_1)
    )[y == 1, ],
    var_types = c(extract_var_type(xtype[edge[[1]]$vertice]),
                  extract_var_type(xtype[edge[[2]]$vertice])),
    family_set = family_pair[2], par_method = "mle"
  )
  if(cop_2$family == "indep")
  {
    cop_2$family <- "gaussian"
    cop_2$parameters <- matrix(0,1,1)
    cop_2$npars <- 1
  }
  list(cop_1, cop_2)
}

pair_effect <- function(edge, u, copula_pair) {
  
  cop_pred_0 <- log(
    predict(copula_pair[[1]],
            weave_transformed(get(cond_node_key(edge[[1]]), envir = u$u_0),
                              get(cond_node_key(edge[[2]]), envir = u$u_0)))
  )
  cop_pred_1 <- log(
    predict(copula_pair[[2]],
            weave_transformed(get(cond_node_key(edge[[1]]), envir = u$u_1),
                              get(cond_node_key(edge[[2]]), envir = u$u_1)))
  )
  
  cop_pred_1 - cop_pred_0
  
}
