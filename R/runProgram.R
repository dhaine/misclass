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
                       seed = 123) {
    if(warmup >= iter)
        stop("Number of warmup samples is greater than number of iterations.")
    if(cores > chains)
        stop("Number of cores should not exceed number of chains.")
    ## Model
    stanfit <- stanmodels$logistic
    ## 1. True association
    ## Data
    x <- model.matrix(S2 ~ E_h, data)
    x.mean <- colMeans(x)
    stan_data <- list(N = nrow(data),
                      y = data$S2,
                      K = 2,
                      X = x,
                      X_means = as.array(x.mean),
                      J_1 = data$cow,
                      N_1 = length(unique(data$cow)),
                      K_1 = 1,
                      Z_1 = rep(1, nrow(data)),
                      J_2 = data$herd,
                      N_2 = length(unique(data$herd)),
                      K_2 = 1,
                      Z_2 = rep(1, nrow(data)))
    ## Running
    true_fit <- sampling(stanfit,
                         data = stan_data,
                         iter = iter,
                         warmup = warmup,
                         chains = chains,
                         seed = seed,
                         cores = cores)
    true_fit
}
