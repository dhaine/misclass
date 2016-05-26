#' Fit Bayesian 3-level logistic model to evaluate misclassification.
#'
#' Fit a Bayesian 3-level logistic model using Stan to evaluate true association,
#' quantify bias, and compare relative impact of selection vs. misclassification
#' biases.
#'
#' @useDynLib misclass, .registration = TRUE
#' @param data Data file.
#' @param iter A positive integer specifying how many iterations for each chain
#' (including warmup). The default is 1000.
#' @param warmup A positive integer specifying number of warmup (aka burnin)
#' iterations. Warmup samples should not be used for inference. The number of
#' warmup should not be larger than iter and the default is 100.
#' @param chains A positive integer specifying number of chains. Defaults to 4.
#' @param cores Number of cores to use when executing the chains in parallel (up
#' to the number of chains).
#' @param seed Positive integer. Used by \code{set.seed} to make results
#' reproducible.
#' @param nsimul Number of simulations.
#'
#' @return An object of class \code{stanfit}.
#'
#' @author Denis Haine
#'
#' @examples
#' runProgram(sim_list[[1]],
#'            iter = 200,
#'            warmup = 25,
#'            chains = 4,
#'            cores = 4,
#'            seed = 123)
#' @export
#' @import rstan
runProgram <- function(data,
                       iter = 1000,
                       warmup = 100,
                       chains = 4,
                       cores,
                       seed = 123,
                       nsimul) {
    if(warmup >= iter)
        stop("Number of warmup samples is greater than number of iterations.")
    if(cores > chains)
        stop("Number of cores should not exceed number of chains.")
    ## Model
    stanfit <- stanmodels$logistic

    ## Collect results
    out_list <- vector("list", nsimul)
    m <- matrix(NA, nrow = chains * (iter - warmup), ncol = 13)
    colnames(m) <- c("post_t", "post_s", "post_ss", "post_sm", "post_dpp", 
                     "post_dpsi", "post_dsip", "post_dss", "post_dssi", "post_dsis", 
                     "post_dps", "post_dsp", "post_tr")

    for (i in 1:nsimul) {
        ## 1. True association
        ## Data
        data2 <- data[[i]][which(data[[i]]$S1 == 0), ]
        x <- model.matrix(S2 ~ E_h, data2)
        x <- x[, 2, drop = FALSE]
        x.mean <- colMeans(x)
        data2$herd2 <- cumsum(c(TRUE, with(data2, herd[-1] != herd[-length(herd)])))
        data2$cow2 <- cumsum(c(TRUE, with(data2, cow[-1] != cow[-length(cow)])))
        stan_data <- list(N = nrow(data2),
                          y = data2$S2,
                          K = 1,
                          X = ifelse(x == 0, -x.mean, 1 - x.mean),
                          X_means = as.array(x.mean),
                          J_1 = data2$herd2,
                          N_1 = length(unique(data2$herd2)),
                          K_1 = 1,
                          Z_1 = rep(1, nrow(data2)),
                          J_2 = data2$cow2,
                          N_2 = length(unique(data2$cow2)),
                          K_2 = 1,
                          Z_2 = rep(1, nrow(data2)))
        ## Running
        fit <- sampling(stanfit,
                        data = stan_data,
                        iter = iter,
                        warmup = warmup,
                        chains = chains,
                        seed = seed,
                        cores = cores)
        m[, 1] <- as.matrix(fit, pars = "b[2]")
                
        ## 2. Single sample S1i S2i
        ## Data
        data2 <- data[[i]][which(data[[i]]$S1i == 0), ]
        x <- model.matrix(S2i ~ E_h, data2)
        x <- x[, 2, drop = FALSE]
        x.mean <- colMeans(x)
        data2$herd2 <- cumsum(c(TRUE, with(data2, herd[-1] != herd[-length(herd)])))
        data2$cow2 <- cumsum(c(TRUE, with(data2, cow[-1] != cow[-length(cow)])))
        stan_data <- list(N = nrow(data2),
                          y = data2$S2,
                          K = 1,
                          X = ifelse(x == 0, -x.mean, 1 - x.mean),
                          X_means = as.array(x.mean),
                          J_1 = data2$herd2,
                          N_1 = length(unique(data2$herd2)),
                          K_1 = 1,
                          Z_1 = rep(1, nrow(data2)),
                          J_2 = data2$cow2,
                          N_2 = length(unique(data2$cow2)),
                          K_2 = 1,
                          Z_2 = rep(1, nrow(data2)))
        ## Running
        fit <- sampling(stanfit,
                        data = stan_data,
                        iter = iter,
                        warmup = warmup,
                        chains = chains,
                        seed = seed,
                        cores = cores)
        m[, 2] <- as.matrix(fit, pars = "b[2]")
        
        ## 3. Single sample selection bias only
        ## Data
        data2 <- data[[i]][which(data[[i]]$S1i == 0), ]
        x <- model.matrix(S2 ~ E_h, data2)
        x <- x[, 2, drop = FALSE]
        x.mean <- colMeans(x)
        data2$herd2 <- cumsum(c(TRUE, with(data2, herd[-1] != herd[-length(herd)])))
        data2$cow2 <- cumsum(c(TRUE, with(data2, cow[-1] != cow[-length(cow)])))
        stan_data <- list(N = nrow(data2),
                          y = data2$S2,
                          K = 1,
                          X = ifelse(x == 0, -x.mean, 1 - x.mean),
                          X_means = as.array(x.mean),
                          J_1 = data2$herd2,
                          N_1 = length(unique(data2$herd2)),
                          K_1 = 1,
                          Z_1 = rep(1, nrow(data2)),
                          J_2 = data2$cow2,
                          N_2 = length(unique(data2$cow2)),
                          K_2 = 1,
                          Z_2 = rep(1, nrow(data2)))
        ## Running
        fit <- sampling(stanfit,
                        data = stan_data,
                        iter = iter,
                        warmup = warmup,
                        chains = chains,
                        seed = seed,
                        cores = cores)
        m[, 3] <- as.matrix(fit, pars = "b[2]")
        
        ## 4. Single sample misclassification bias only
        ## Data
        data2 <- data[[i]][which(data[[i]]$S1 == 0), ]
        x <- model.matrix(S2i ~ E_h, data2)
        x <- x[, 2, drop = FALSE]
        x.mean <- colMeans(x)
        data2$herd2 <- cumsum(c(TRUE, with(data2, herd[-1] != herd[-length(herd)])))
        data2$cow2 <- cumsum(c(TRUE, with(data2, cow[-1] != cow[-length(cow)])))
        stan_data <- list(N = nrow(data2),
                          y = data2$S2,
                          K = 1,
                          X = ifelse(x == 0, -x.mean, 1 - x.mean),
                          X_means = as.array(x.mean),
                          J_1 = data2$herd2,
                          N_1 = length(unique(data2$herd2)),
                          K_1 = 1,
                          Z_1 = rep(1, nrow(data2)),
                          J_2 = data2$cow2,
                          N_2 = length(unique(data2$cow2)),
                          K_2 = 1,
                          Z_2 = rep(1, nrow(data2)))
        ## Running
        fit <- sampling(stanfit,
                        data = stan_data,
                        iter = iter,
                        warmup = warmup,
                        chains = chains,
                        seed = seed,
                        cores = cores)
        m[, 4] <- as.matrix(fit, pars = "b[2]")
        
        ## 5. Duplicate sample parallel on S1 and S2
        ## Data
        data2 <- data[[i]][which(data[[i]]$S1_parall == 0), ]
        x <- model.matrix(S2_parall ~ E_h, data2)
        x <- x[, 2, drop = FALSE]
        x.mean <- colMeans(x)
        data2$herd2 <- cumsum(c(TRUE, with(data2, herd[-1] != herd[-length(herd)])))
        data2$cow2 <- cumsum(c(TRUE, with(data2, cow[-1] != cow[-length(cow)])))
        stan_data <- list(N = nrow(data2),
                          y = data2$S2,
                          K = 1,
                          X = ifelse(x == 0, -x.mean, 1 - x.mean),
                          X_means = as.array(x.mean),
                          J_1 = data2$herd2,
                          N_1 = length(unique(data2$herd2)),
                          K_1 = 1,
                          Z_1 = rep(1, nrow(data2)),
                          J_2 = data2$cow2,
                          N_2 = length(unique(data2$cow2)),
                          K_2 = 1,
                          Z_2 = rep(1, nrow(data2)))
        ## Running
        fit <- sampling(stanfit,
                        data = stan_data,
                        iter = iter,
                        warmup = warmup,
                        chains = chains,
                        seed = seed,
                        cores = cores)
        m[, 5] <- as.matrix(fit, pars = "b[2]")
        
        ## 6. Duplicate sample parallel on S1 and simple on S2
        ## Data
        data2 <- data[[i]][which(data[[i]]$S1_parall == 0), ]
        x <- model.matrix(S2i ~ E_h, data2)
        x <- x[, 2, drop = FALSE]
        x.mean <- colMeans(x)
        data2$herd2 <- cumsum(c(TRUE, with(data2, herd[-1] != herd[-length(herd)])))
        data2$cow2 <- cumsum(c(TRUE, with(data2, cow[-1] != cow[-length(cow)])))
        stan_data <- list(N = nrow(data2),
                          y = data2$S2,
                          K = 1,
                          X = ifelse(x == 0, -x.mean, 1 - x.mean),
                          X_means = as.array(x.mean),
                          J_1 = data2$herd2,
                          N_1 = length(unique(data2$herd2)),
                          K_1 = 1,
                          Z_1 = rep(1, nrow(data2)),
                          J_2 = data2$cow2,
                          N_2 = length(unique(data2$cow2)),
                          K_2 = 1,
                          Z_2 = rep(1, nrow(data2)))
        ## Running
        fit <- sampling(stanfit,
                        data = stan_data,
                        iter = iter,
                        warmup = warmup,
                        chains = chains,
                        seed = seed,
                        cores = cores)
        m[, 6] <- as.matrix(fit, pars = "b[2]")
        
        ## 7. Duplicate sample parallel on S2 and simple on S1
        ## Data
        data2 <- data[[i]][which(data[[i]]$S1i == 0), ]
        x <- model.matrix(S2_parall ~ E_h, data2)
        x <- x[, 2, drop = FALSE]
        x.mean <- colMeans(x)
        data2$herd2 <- cumsum(c(TRUE, with(data2, herd[-1] != herd[-length(herd)])))
        data2$cow2 <- cumsum(c(TRUE, with(data2, cow[-1] != cow[-length(cow)])))
        stan_data <- list(N = nrow(data2),
                          y = data2$S2,
                          K = 1,
                          X = ifelse(x == 0, -x.mean, 1 - x.mean),
                          X_means = as.array(x.mean),
                          J_1 = data2$herd2,
                          N_1 = length(unique(data2$herd2)),
                          K_1 = 1,
                          Z_1 = rep(1, nrow(data2)),
                          J_2 = data2$cow2,
                          N_2 = length(unique(data2$cow2)),
                          K_2 = 1,
                          Z_2 = rep(1, nrow(data2)))
        ## Running
        fit <- sampling(stanfit,
                        data = stan_data,
                        iter = iter,
                        warmup = warmup,
                        chains = chains,
                        seed = seed,
                        cores = cores)
        m[, 7] <- as.matrix(fit, pars = "b[2]")
        
        ## 8. Duplicate sample series on S1 and S2
        ## Data
        data2 <- data[[i]][which(data[[i]]$S1_series == 0), ]
        x <- model.matrix(S2_series ~ E_h, data2)
        x <- x[, 2, drop = FALSE]
        x.mean <- colMeans(x)
        data2$herd2 <- cumsum(c(TRUE, with(data2, herd[-1] != herd[-length(herd)])))
        data2$cow2 <- cumsum(c(TRUE, with(data2, cow[-1] != cow[-length(cow)])))
        stan_data <- list(N = nrow(data2),
                          y = data2$S2,
                          K = 1,
                          X = ifelse(x == 0, -x.mean, 1 - x.mean),
                          X_means = as.array(x.mean),
                          J_1 = data2$herd2,
                          N_1 = length(unique(data2$herd2)),
                          K_1 = 1,
                          Z_1 = rep(1, nrow(data2)),
                          J_2 = data2$cow2,
                          N_2 = length(unique(data2$cow2)),
                          K_2 = 1,
                          Z_2 = rep(1, nrow(data2)))
        ## Running
        fit <- sampling(stanfit,
                        data = stan_data,
                        iter = iter,
                        warmup = warmup,
                        chains = chains,
                        seed = seed,
                        cores = cores)
        m[, 8] <- as.matrix(fit, pars = "b[2]")
        
        ## 9. Duplicate sample series on S1 and simple on S2
        ## Data
        data2 <- data[[i]][which(data[[i]]$S1_series == 0), ]
        x <- model.matrix(S2i ~ E_h, data2)
        x <- x[, 2, drop = FALSE]
        x.mean <- colMeans(x)
        data2$herd2 <- cumsum(c(TRUE, with(data2, herd[-1] != herd[-length(herd)])))
        data2$cow2 <- cumsum(c(TRUE, with(data2, cow[-1] != cow[-length(cow)])))
        stan_data <- list(N = nrow(data2),
                          y = data2$S2,
                          K = 1,
                          X = ifelse(x == 0, -x.mean, 1 - x.mean),
                          X_means = as.array(x.mean),
                          J_1 = data2$herd2,
                          N_1 = length(unique(data2$herd2)),
                          K_1 = 1,
                          Z_1 = rep(1, nrow(data2)),
                          J_2 = data2$cow2,
                          N_2 = length(unique(data2$cow2)),
                          K_2 = 1,
                          Z_2 = rep(1, nrow(data2)))
        ## Running
        fit <- sampling(stanfit,
                        data = stan_data,
                        iter = iter,
                        warmup = warmup,
                        chains = chains,
                        seed = seed,
                        cores = cores)
        m[, 9] <- as.matrix(fit, pars = "b[2]")
        
        ## 10. Duplicate sample series on S2 and simple on S1
        ## Data
        data2 <- data[[i]][which(data[[i]]$S1i == 0), ]
        x <- model.matrix(S2_series ~ E_h, data2)
        x <- x[, 2, drop = FALSE]
        x.mean <- colMeans(x)
        data2$herd2 <- cumsum(c(TRUE, with(data2, herd[-1] != herd[-length(herd)])))
        data2$cow2 <- cumsum(c(TRUE, with(data2, cow[-1] != cow[-length(cow)])))
        stan_data <- list(N = nrow(data2),
                          y = data2$S2,
                          K = 1,
                          X = ifelse(x == 0, -x.mean, 1 - x.mean),
                          X_means = as.array(x.mean),
                          J_1 = data2$herd2,
                          N_1 = length(unique(data2$herd2)),
                          K_1 = 1,
                          Z_1 = rep(1, nrow(data2)),
                          J_2 = data2$cow2,
                          N_2 = length(unique(data2$cow2)),
                          K_2 = 1,
                          Z_2 = rep(1, nrow(data2)))
        ## Running
        fit <- sampling(stanfit,
                        data = stan_data,
                        iter = iter,
                        warmup = warmup,
                        chains = chains,
                        seed = seed,
                        cores = cores)
        m[, 10] <- as.matrix(fit, pars = "b[2]")
        
        ## 11. Duplicate sample parallel on S1 and series on S2
        ## Data
        data2 <- data[[i]][which(data[[i]]$S1_parall == 0), ]
        x <- model.matrix(S2_series ~ E_h, data2)
        x <- x[, 2, drop = FALSE]
        x.mean <- colMeans(x)
        data2$herd2 <- cumsum(c(TRUE, with(data2, herd[-1] != herd[-length(herd)])))
        data2$cow2 <- cumsum(c(TRUE, with(data2, cow[-1] != cow[-length(cow)])))
        stan_data <- list(N = nrow(data2),
                          y = data2$S2,
                          K = 1,
                          X = ifelse(x == 0, -x.mean, 1 - x.mean),
                          X_means = as.array(x.mean),
                          J_1 = data2$herd2,
                          N_1 = length(unique(data2$herd2)),
                          K_1 = 1,
                          Z_1 = rep(1, nrow(data2)),
                          J_2 = data2$cow2,
                          N_2 = length(unique(data2$cow2)),
                          K_2 = 1,
                          Z_2 = rep(1, nrow(data2)))
        ## Running
        fit <- sampling(stanfit,
                        data = stan_data,
                        iter = iter,
                        warmup = warmup,
                        chains = chains,
                        seed = seed,
                        cores = cores)
        m[, 11] <- as.matrix(fit, pars = "b[2]")
        
        ## 12. Duplicate sample series on S1 and parallel on S2
        ## Data
        data2 <- data[[i]][which(data[[i]]$S1_series == 0), ]
        x <- model.matrix(S2_parall ~ E_h, data2)
        x <- x[, 2, drop = FALSE]
        x.mean <- colMeans(x)
        data2$herd2 <- cumsum(c(TRUE, with(data2, herd[-1] != herd[-length(herd)])))
        data2$cow2 <- cumsum(c(TRUE, with(data2, cow[-1] != cow[-length(cow)])))
        stan_data <- list(N = nrow(data2),
                          y = data2$S2,
                          K = 1,
                          X = ifelse(x == 0, -x.mean, 1 - x.mean),
                          X_means = as.array(x.mean),
                          J_1 = data2$herd2,
                          N_1 = length(unique(data2$herd2)),
                          K_1 = 1,
                          Z_1 = rep(1, nrow(data2)),
                          J_2 = data2$cow2,
                          N_2 = length(unique(data2$cow2)),
                          K_2 = 1,
                          Z_2 = rep(1, nrow(data2)))
        ## Running
        fit <- sampling(stanfit,
                        data = stan_data,
                        iter = iter,
                        warmup = warmup,
                        chains = chains,
                        seed = seed,
                        cores = cores)
        m[, 12] <- as.matrix(fit, pars = "b[2]")
        
        ## 13. Triplicate samples
        ## Data
        data2 <- data[[i]][which(data[[i]]$S1_tri == 0), ]
        x <- model.matrix(S2_tri ~ E_h, data2)
        x <- x[, 2, drop = FALSE]
        x.mean <- colMeans(x)
        data2$herd2 <- cumsum(c(TRUE, with(data2, herd[-1] != herd[-length(herd)])))
        data2$cow2 <- cumsum(c(TRUE, with(data2, cow[-1] != cow[-length(cow)])))
        stan_data <- list(N = nrow(data2),
                          y = data2$S2,
                          K = 1,
                          X = ifelse(x == 0, -x.mean, 1 - x.mean),
                          X_means = as.array(x.mean),
                          J_1 = data2$herd2,
                          N_1 = length(unique(data2$herd2)),
                          K_1 = 1,
                          Z_1 = rep(1, nrow(data2)),
                          J_2 = data2$cow2,
                          N_2 = length(unique(data2$cow2)),
                          K_2 = 1,
                          Z_2 = rep(1, nrow(data2)))
        ## Running
        fit <- sampling(stanfit,
                        data = stan_data,
                        iter = iter,
                        warmup = warmup,
                        chains = chains,
                        seed = seed,
                        cores = cores)
        m[, 13] <- as.matrix(fit, pars = "b[2]")

        out_list[[i]] <- m
    }
    out_list
}
