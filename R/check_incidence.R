#' Fit Bayesian 3-level Poisson model to assess bias of measure of incidence
#'
#' Fit a Bayesian 3-level Poisson model using Stan to evaluate true incidence,
#' quantify bias, and compare relative impact of selection vs. misclassification
#' biases.
#'
#' @useDynLib misclass, .registration = TRUE
#' @param data Data file.
#' @param iter A positive integer specifying how many iterations for each chain
#' (including warmup). The default is 500.
#' @param warmup a positive integer specifying number of warmup (aka burnin)
#' iterations. Warmup samples should not be used for inference. The number of
#' warmup should not be larger than \code{iter} and the default is 100.
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
#' sim_list <- vector("list", 1)
#' set.seed(123)
#' sim_list <- replicate(n = 1, expr = make_data(100, 30, "saureus"), simplify = FALSE)
#' check_incidence(sim_list,
#'                 iter = 200,
#'                 warmup = 25,
#'                 chains = 1,
#'                 cores = 1,
#'                 seed = 123,
#'                 nsimul = 1)
#' @export
#' @import Rcpp
#' @importFrom rstan sampling get_sampler_params extract summary
#' @importFrom stats model.matrix
#' @importFrom utils setTxtProgressBar txtProgressBar
check_incidence <- function(data,
                            iter =500,
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
    stanfit <- stanmodels$poisson

    ## Collect results
    out_list <- vector("list", nsimul)
    m <- matrix(NA, nrow = chains * (iter - warmup), ncol = 4)
    colnames(m) <- c("post_t", "post_s", "post_ss", "post_sm")
    diag_list <- vector("list", nsimul)
    d <- matrix(NA, nrow = 4, ncol = 5)
    colnames(d) <- c("Rhat", "n_eff", "sample_size", "divergent", "treedepth")

    count_divergences <- function(fit) {
        sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
        sum(sapply(sampler_params, function(x) c(x[, 'divergent__']))[, 1])
    }
    max_treedepth <- function(fit) {
        sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
        max(sapply(sampler_params, function(x) c(x[, 'treedepth__']))[, 1])
    }

    pb <- utils::txtProgressBar(min = 0, max = nsimul, style = 3)

    for (i in 1:nsimul) {
        ## 1. True association
        ## Data
        data2 <- data[[i]][which(data[[i]]$S1 == 0), ]
        x <- stats::model.matrix(S2 ~ 1, data2)
        x <- x[, 0, drop = FALSE]
        data2$herd2 <- cumsum(c(TRUE, with(data2, herd[-1] != herd[-length(herd)])))
        data2$cow2 <- cumsum(c(TRUE, with(data2, cow[-1] != cow[-length(cow)])))
        stan_data <- list(N = nrow(data2),
                          y = data2$S2,
                          X = x,
                          J_1 = data2$herd2,
                          N_1 = length(unique(data2$herd2)),
                          K_1 = 1,
                          Z_1 = rep(1, nrow(data2)),
                          J_2 = data2$cow2,
                          N_2 = length(unique(data2$cow2)),
                          K_2 = 1,
                          Z_2 = rep(1, nrow(data2)))
        ## Running
        fit <- rstan::sampling(stanfit,
                               data = stan_data,
                               iter = iter,
                               warmup = warmup,
                               chains = chains,
                               seed = seed,
                               cores = cores,
                               control = list(adapt_delta = 0.99))
        m[, 1] <- as.matrix(fit, pars = "temp_Intercept")
        d[1, 1:5] <- c(rstan::summary(fit)$summary[1, "Rhat"],
                       rstan::summary(fit)$summary[1, "n_eff"],
                       length(rstan::extract(fit, pars = "lp__")[[1]]),
                       count_divergences(fit),
                       max_treedepth(fit))
                
        ## 2. Single samples
        ## Data
        data2 <- data[[i]][which(data[[i]]$S1i == 0), ]
        x <- stats::model.matrix(S2i ~ 1, data2)
        x <- x[, 0, drop = FALSE]
        data2$herd2 <- cumsum(c(TRUE, with(data2, herd[-1] != herd[-length(herd)])))
        data2$cow2 <- cumsum(c(TRUE, with(data2, cow[-1] != cow[-length(cow)])))
        stan_data <- list(N = nrow(data2),
                          y = data2$S2i,
                          X = x,
                          J_1 = data2$herd2,
                          N_1 = length(unique(data2$herd2)),
                          K_1 = 1,
                          Z_1 = rep(1, nrow(data2)),
                          J_2 = data2$cow2,
                          N_2 = length(unique(data2$cow2)),
                          K_2 = 1,
                          Z_2 = rep(1, nrow(data2)))
        ## Running
        fit <- rstan::sampling(stanfit,
                               data = stan_data,
                               iter = iter,
                               warmup = warmup,
                               chains = chains,
                               seed = seed,
                               cores = cores,
                               control = list(adapt_delta = 0.99))
        m[, 2] <- as.matrix(fit, pars = "temp_Intercept")
        d[2, 1:5] <- c(rstan::summary(fit)$summary[1, "Rhat"],
                       rstan::summary(fit)$summary[1, "n_eff"],
                       length(rstan::extract(fit, pars = "lp__")[[1]]),
                       count_divergences(fit),
                       max_treedepth(fit))
        
        ## 3. Single samples, selection bias only
        ## Data
        data2 <- data[[i]][which(data[[i]]$S1i == 0), ]
        x <- stats::model.matrix(S2 ~ 1, data2)
        x <- x[, 0, drop = FALSE]
        data2$herd2 <- cumsum(c(TRUE, with(data2, herd[-1] != herd[-length(herd)])))
        data2$cow2 <- cumsum(c(TRUE, with(data2, cow[-1] != cow[-length(cow)])))
        stan_data <- list(N = nrow(data2),
                          y = data2$S2,
                          X = x,
                          J_1 = data2$herd2,
                          N_1 = length(unique(data2$herd2)),
                          K_1 = 1,
                          Z_1 = rep(1, nrow(data2)),
                          J_2 = data2$cow2,
                          N_2 = length(unique(data2$cow2)),
                          K_2 = 1,
                          Z_2 = rep(1, nrow(data2)))
        ## Running
        fit <- rstan::sampling(stanfit,
                               data = stan_data,
                               iter = iter,
                               warmup = warmup,
                               chains = chains,
                               seed = seed,
                               cores = cores,
                               control = list(adapt_delta = 0.99))
        m[, 3] <- as.matrix(fit, pars = "temp_Intercept")
        d[3, 1:5] <- c(rstan::summary(fit)$summary[1, "Rhat"],
                       rstan::summary(fit)$summary[1, "n_eff"],
                       length(rstan::extract(fit, pars = "lp__")[[1]]),
                       count_divergences(fit),
                       max_treedepth(fit))
        
        ## 4. Single samples, misclassification bias only
        ## Data
        data2 <- data[[i]][which(data[[i]]$S1 == 0), ]
        x <- stats::model.matrix(S2i ~ 1, data2)
        x <- x[, 0, drop = FALSE]
        data2$herd2 <- cumsum(c(TRUE, with(data2, herd[-1] != herd[-length(herd)])))
        data2$cow2 <- cumsum(c(TRUE, with(data2, cow[-1] != cow[-length(cow)])))
        stan_data <- list(N = nrow(data2),
                          y = data2$S2i,
                          X = 0,
                          J_1 = data2$herd2,
                          N_1 = length(unique(data2$herd2)),
                          K_1 = 1,
                          Z_1 = rep(1, nrow(data2)),
                          J_2 = data2$cow2,
                          N_2 = length(unique(data2$cow2)),
                          K_2 = 1,
                          Z_2 = rep(1, nrow(data2)))
        ## Running
        fit <- rstan::sampling(stanfit,
                               data = stan_data,
                               iter = iter,
                               warmup = warmup,
                               chains = chains,
                               seed = seed,
                               cores = cores,
                               control = list(adapt_delta = 0.99))
        m[, 4] <- as.matrix(fit, pars = "temp_Intercept")
        d[4, 1:5] <- c(rstan::summary(fit)$summary[1, "Rhat"],
                       rstan::summary(fit)$summary[1, "n_eff"],
                       length(rstan::extract(fit, pars = "lp__")[[1]]),
                       count_divergences(fit),
                       max_treedepth(fit))

        out_list[[i]] <- m
        diag_list[[i]] <- d

        utils::setTxtProgressBar(pb, i)
    }
    close(pb)
    out_list <- append(out_list, diag_list)
    out_list
}
