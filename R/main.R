
init_logistic_copula <- function(y, x, xtype, which_include, reg.method="glm",
                                 set_nonsig_zero=FALSE) {
  
  if (reg.method == "glm") {
    ft <- glm(y~x, family=binomial())
    if (set_nonsig_zero) {
      nonsig <- which(summary(ft)$coef[-1, 4] > .05)
      sig <- which(summary(ft)$coef[-1, 4] <= .05)
      ft_n <- glm(y~x[, -nonsig])
      ft$coef[nonsig + 1] <- 0
      ft$coeg[c(1, 1 + sig)] <- coef(ft_n)
    }
    
  } else {
    stop(paste0("reg.method = ", reg.method, " not implemented!"))
  }
  
  
  prs <- coefs_to_pars(
    y, x, xtype, coef(ft)[1], coef(ft)[-1]
  )

  trees <- TreeGraphList(length(which_include))
  u <- transform_margins_to_hash(x, xtype, prs, which_include)
  beta_0 <- parameters_to_intercept(xtype, prs)
  beta_x <- pars_to_linear_pred(x, xtype, prs) - beta_0
  copula_eff <- rep(0, length(beta_x))
  
  m_obj <- list(
    trees=trees, xtype=xtype, parameters=prs, beta_0=beta_0,
    beta_x_insample=beta_x, copula_eff_insample=copula_eff, u=u,
    which_include=which_include, beta_vec=coef(ft)
  )
  class(m_obj) <- "logistic_copula"
  attr(m_obj, "call") <- sys.call()
  
  m_obj
}

predict.logistic_copula <- function(object, new_x, ...) {
  ##' predict.logistic_copula
  ##' @name predict.logistic_copula
  ##' @aliases predict.logistic_copula
  ##' @description
  ##' Computes predicted probability of Y=1 for a logistic regression model with
  ##' a vine extension. 
  ##' @param object The model object as returned by fit_copula_interactions
  ##' @param new_x A matrix of covariate values to compute predictions for.
  ##' @param ... Not used.
  new_u <- transform_margins_to_hash(
    new_x, object$xtype, object$parameters, object$which_include
  )
  object_copy <- object
  object_copy <- set_transformed_vars(
    new_x, object_copy$which_include, object_copy
    )
  cop_eff <- compute_g(object_copy)

  plogis(cbind(1, new_x) %*% object$beta_vec + cop_eff)
}

fit_copula_interactions <- function(
    y, x, xtype, family_set=c("gaussian", "clayton", "gumbel"), 
    oos_validation=FALSE, tau=2, which_include=NULL, reg.method="glm",
    maxit_final=1000, maxit_intermediate=50, verbose=FALSE, 
    adjust_intercept=TRUE, max_t=Inf, test_x=NULL, test_y=NULL,
    set_nonsig_zero=FALSE, reltol=sqrt(.Machine$double.eps)
    ) {
  ##' fit_copula_interactions
  ##' @name fit_copula_interactions
  ##' @aliases fit_copula_interactions
  ##' @description This is the main function of the package, which
  ##' starting from an initial logistic regression model with only main effects
  ##' of each covariate, selects and fits interaction terms in the form of two
  ##' R-vine models with identical graphical structure, one for each class.
  ##' @param y A vector of n observations of the (univariate) binary outcome
  ##' variable y
  ##' @param x A (n x p) matrix of n observations of p covariates
  ##' @param xtype A vector of p characters that have to take the value
  ##' "c_a", "c_p", "d_b" or "d_b", to indicate whether each margin of the
  ##' is continuous with full support, continuous with support on the positive
  ##' real line, discrete (binary) or a counting variable.
  ##' @param family_set A vector of strings that specifies the set of
  ##' pair-copula families that the fitting algorithm chooses from. For an
  ##' overview of which values that can be specified, see the documentation for
  ##' \link[rvinecopulib]{bicop}.
  ##' @param oos_validation Whether to use an external sample for validation
  ##' instead of an in-sample likelihood based criteria. Would require that
  ##' both test_x and test_y are provided if set to TRUE. 
  ##' @param tau Parameter used when selecting the structure, where the 
  ##' the criteria is (new_likelihood - previous_likelihood - tau), 
  ##' so that an additional edge in the copulas is only accepted if it leads to
  ##' an increase in the likelihood that exceeds tau. Setting tau to NULL, has 
  ##' the same effect as -Inf. 
  ##' @param which_include The column indices of the covariates that could be
  ##' included in the copula effects.
  ##' @param reg.method The method by which the initial regression coefficients
  ##' are fitted.
  ##' @param maxit_final The maximum number of gradient optimisation iterations
  ##' to use when the full structure has been selected to refit all the
  ##' parameters. Defaults to 1000.
  ##' @param maxit_intermediate The maximum number of gradient optimisation
  ##' iterations to use when adding a newly selected component to refit the
  ##' parameters. Defaults to 10.
  ##' @param verbose Whether information about the progress should be printed 
  ##' to the console.
  ##' @param adjust_intercept Whether to intermediately refit the intercept
  ##' during the model/structure selection procedure. Defaults to true. 
  ##' @param max_t The maximum number of trees in the copula models. Defaults
  ##' to Inf, i.e., no maximum.
  ##' @param test_x Part of the optional validation set,
  ##' see @oos_validation.
  ##' @param test_y Part of the optional validation set,
  ##' see @oos_validation.
  ##' @param set_nonsig_zero If true, non-significant regression coefficients 
  ##' (in the initial glm model) will be set to zero 
  ##' @param reltol Relative convergence tolerance, see the documentation for 
  ##' \link[stats]{optim}.
  ##' @examples 
  ##' library(mclust)
  ##' data("wdbc", package="mclust")
  ##' y <- wdbc$Diagnosis == "M"
  ##' x <- as.matrix(wdbc[, 3:11])
  ##' rowss <- sample(length(y), round(length(y) * 0.5))
  ##' xtype <- rep("c_a", ncol(x))
  ##' 
  ##' md <- LogisticCopula::fit_copula_interactions(
  ##'   y[rowss], x[rowss, ], xtype, verbose = T, tau=log(length(y[rowss])),
  ##'   maxit_intermediate = 50, maxit_final = 50
  ##' )
  ##' md2 <- LogisticCopula::fit_copula_interactions(
  ##'   y[rowss], x[rowss, ], xtype, verbose = T, tau=Inf
  ##' )
  ##' 
  ##' plot(predict(md2, new_x = x[-rowss, ]),
  ##'      predict(md, new_x = x[-rowss, ]), col = y[-rowss] + 3)
  ##' @export
  
  if (is.null(which_include)) {
    which_include <- which((xtype == "c_a") | (xtype == "c_p"))
  } else {
    if(any((xtype[which_include] != "c_a") & (xtype[which_include] != "c_p"))) {
      ind <- which(
        (xtype[which_include] != "c_a") & (xtype[which_include] != "c_p")
        )
      which_include <- which_include[-ind] 
    }
  }

  if (length(family_set) > 1) {
    family_combs <- cbind(
      combn(family_set, 2),
      combn(family_set, 2)[rev(1:2), ],
      rbind(family_set,
            family_set)
    )
  } else{
    family_combs <- matrix(rep(family_set, 2), ncol=1)
  }

  if (oos_validation & (is.null(test_x) | is.null(test_y))) {
    stop("Must supply test_x and test_y for cdevriterion = 'OOS validation'!")
  }

  m_obj <- init_logistic_copula(y, x, xtype, which_include, 
                                reg.method=reg.method,
                                set_nonsig_zero=set_nonsig_zero)

  for (t in seq(min(length(which_include) - 1, max_t))) {
    # Find all valid edges
    valid_edges <- all_initially_valid_edges(
      m_obj$trees[[t]]$nodes, 
      if(!is.null(which_include)) which_include else seq(length(xtype))
    )
    
    if (length(valid_edges$edges) == 0) {
      break
    }
    
    # Fit all model for all valid edges, only gaussian
    copula_pairs <- new.env(hash=T)
    for (edge in valid_edges$edges) {
      assign(edge_key(edge), fit_pair(edge, m_obj$u, y, xtype))
    }
    
    improvement <- TRUE
    k <- 1
    while (improvement) {
      
      # Compute the copula - effects for each edge
      pair_effects <- sapply(
        1:length(valid_edges$edges),
        function(j) pair_effect(
          valid_edges$edges[[j]], m_obj$u,
          get(edge_key(valid_edges$edges[[j]]), envir=copula_pairs)
        )
      )

      # Find the best edge, whether it improves the model or not
      j_star <- choose_best_effect(
        y, m_obj$beta_0 + m_obj$beta_x + m_obj$copula_eff_insample,
        pair_effects,
        adjust_intercept=adjust_intercept
      )
      
      rm(pair_effects)
      
      best_edge <- valid_edges$edges[[j_star]]
      
      # Fit all combinations of pair - copula families to the best edge
      copula_pairs_family_comb <- lapply(
        1:ncol(family_combs),
        function(l) fit_pair(best_edge, m_obj$u, y, xtype, family_combs[, l])
      )

      # Compute the copula effects for each combination
      pair_effects_family_comb <- sapply(
        1:ncol(family_combs),
        function(l) pair_effect(
          best_edge, m_obj$u, copula_pairs_family_comb[[l]]
        )
      )

      # Find the best such combination, if there is one that improves the model
      l_star <- choose_best_effect(
        y, m_obj$beta_0 + m_obj$beta_x + m_obj$copula_eff_insample,
        pair_effects_family_comb, adjust_intercept=adjust_intercept
      )

      # Is there a combination that improves the model?¨
      m_cand <- m_obj

      # Add the edge to the tree
      m_cand$trees[[t]] <- add_edge_to_tree(
        m_cand$trees[[t]], valid_edges$edges[[j_star]]
      )
      
      m_cand$trees[[t]]$graph <- igraph::add_edges(
        m_cand$trees[[t]]$graph, valid_edges$edge_mat[, j_star]
      )
      m_cand$trees[[t]] <- add_pair_copulas_to_tree(
        m_cand$trees[[t]], copula_pairs_family_comb[[l_star]]
      )

      # Refit 
      m_cand_tmp <- try(
        fit_model(
            y, x, m_cand, maxit=min(maxit_final, maxit_intermediate),
            num_grad=FALSE, verbose=verbose, reltol=reltol
        ),
        silent=TRUE
      )

      if(length(m_cand_tmp) == 1) {
        stop("Error during fitting of model")
      } else {
        m_cand <- m_cand_tmp
        rm(m_cand_tmp)
      }

      cand_lik <- full_likelihood(
        y, x, m_cand, m_cand$beta_0, m_cand$beta_vec[-1],
        transformed_delta_vec(m_cand)
      )

      prev_lik <- full_likelihood(
        y, x, m_obj, m_obj$beta_0, m_obj$beta_vec[-1],
        transformed_delta_vec(m_obj)
      )

      if (oos_validation) {
        cand_lik <- sum(
          dbinom(test_y, 1, plogis(predict(m_cand, test_x)), TRUE)
        )
        prev_lik <- sum(
          dbinom(test_y, 1, plogis(predict(m_obj, test_x)), TRUE)
        )
        crit <- cand_lik - prev_lik
      } else {
        if (!is.null(tau)) {
          crit <- 2 * (cand_lik - prev_lik - tau)
        } else {
          crit <- Inf
        }
      }

      if (crit > 0) {
        m_obj <- m_cand
        valid_edges$edges <- valid_edges$edges[-j_star]
        valid_edges$edge_mat <- as.matrix(valid_edges$edge_mat[, -j_star])
        
        if(length(valid_edges$edges) > 0) {
          if (
            any(!sapply(seq(length(valid_edges$edges)),
                        function (j) valid_edge(m_obj$trees[[t]],
                                                valid_edges$edges[[j]],
                                                valid_edges$edge_mat[, j])))
          ) {
            rem_edg <- which(
              !sapply(seq(length(valid_edges$edges)),
                      function (j) valid_edge(
                        m_obj$trees[[t]], valid_edges$edges[[j]],
                        valid_edges$edge_mat[, j])
              )
            )
            valid_edges$edges <- valid_edges$edges[-rem_edg]
            valid_edges$edge_mat <- as.matrix(valid_edges$edge_mat[, -rem_edg])
          }
        }
      }
      
      if (length(valid_edges$edges) == 0 | crit <= 0){
        improvement <- FALSE
      }
      rm(m_cand)
      rm(pair_effects_family_comb)
      rm(copula_pairs_family_comb)
    }
    
    # ADD THE NEW TREE
    if (length(m_obj$trees[[t]]$edges) > 1){
      m_obj$trees <- append(
        m_obj$trees,
        list(TreeGraph(edges_to_next_trees_nodes(m_obj$trees[[t]]$edges)))
      )
    } else if(t < length(m_obj$trees[[1]]$nodes) - 1) {
      break
    }
    k <- k + 1
    rm(copula_pairs)
  }

  if (length(m_obj$trees[[length(m_obj$trees)]]$edges) == 0) {
    m_obj$trees <- m_obj$trees[1:(length(m_obj$trees) - 1)]
  }

  # Refit one last time
  if(maxit_final > 50) {
    if(length(m_obj$trees[[1]]$edges) > 0) {
      m_obj_tmp <- try(
        fit_model(
          y, x, m_obj, maxit=maxit_final, num_grad=FALSE, verbose=verbose,
          reltol=reltol
        ),
        silent=TRUE
      )
      if(length(m_obj_tmp) == 1) {
        print("Error during fitting of model")
        stop()
      } else {
        m_obj <- m_obj_tmp
        rm(m_obj_tmp)
      }
    }
  }

  final_lik <- full_likelihood(
    y, x, m_obj, m_obj$beta_0, m_obj$beta_vec[-1],
    transformed_delta_vec(m_obj)
  )

  rm(valid_edges)
  m_obj
}


fit_model <- function(y, x, m_obj, maxit=5, num_grad=FALSE,
                      verbose=FALSE, hessian=FALSE, reltol=sqrt(.Machine$double.eps)) {
  
  # Extract parameters :
  t_delta <- transformed_delta_vec(m_obj)

  optfit <- optim(
    fn=function(par) {
      lik <- full_likelihood(y, x, m_obj, beta_0=par[1], betas=par[2:(ncol(x) + 1)],
                             delta=par[(ncol(x) + 2):length(par)],
                             transform_delta=TRUE)
      -lik
    },
    gr=function(par) {
      if (num_grad) {
        numDeriv::grad(function(par) {
          - full_likelihood(y, x, m_obj, beta_0=par[1], betas=par[2:(ncol(x) + 1)],
                            delta=par[(ncol(x) + 2):length(par)],
                            transform_delta=TRUE)
        }, x=par)
      } else {
        -full_gradient(y, x, m_obj, beta_0=par[1], betas=par[2:(ncol(x) + 1)], 
                       par[(ncol(x) + 2):length(par)])
      }
    },
    par=c(m_obj$beta_vec, t_delta), method="BFGS",
    control=list(maxit=maxit, trace=ifelse(verbose, 100, 0), reltol=reltol),
    hessian=hessian
  )

  m_obj <- set_model(
    y, x, m_obj, optfit$par[seq(length(m_obj$beta_vec))], 
    transform_t_to_delta_vec(optfit$par[-seq(length(m_obj$beta_vec))], m_obj)
  )
  m_obj$fit_inf <- optfit
  m_obj

}


full_likelihood <- function(y, x, fit, beta_0, betas, delta, 
                            transform_delta=TRUE) {
  if (any(fit$xtype == "d_c") | any(fit$xtype == "d_b")) {
    ind_d <- which((fit$xtype == "d_c") | (fit$xtype == "d_b"))
    if (any(is.infinite(exp(betas[ind_d])))) {
      -1e100
    } else {
      if (transform_delta & !is.null(delta)) {
        t_delta <- delta
        delta <- transform_t_to_delta_vec(delta, fit)
      }
      
      fit <- set_model(y, x, fit, c(beta_0, betas), delta)
      eta <- fit$beta_0 + fit$beta_x_insample + fit$copula_eff_insample
      
      p <- plogis(eta)
      
      sum(dbinom(y, 1, p, TRUE))
    }
  }
  else {
    if (transform_delta & !is.null(delta)) {
      t_delta <- delta
      delta <- transform_t_to_delta_vec(delta, fit)
    }
    
    fit <- set_model(y, x, fit, c(beta_0, betas), delta)
    eta <- fit$beta_0 + fit$beta_x_insample + fit$copula_eff_insample
    
    p <- plogis(eta)
    
    sum(dbinom(y, 1, p, TRUE))
  }
}
