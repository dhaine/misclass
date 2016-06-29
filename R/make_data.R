#' Create simulated datasets.
#'
#' Create datasets from a hypothetical cohort study. Since correlation structures
#' are often found in udder health studies, the generated datasets are deemed
#' to have occurred from the collection of two milk samples collected 1 month
#' apart from each quarters of a random sample of 30 cows per herd, from 100
#' dairy herds. The first milk sample (S1) is used to identify quarters at risk of
#' intramammary infection (IMI) at the beginning of the cohort, while the second
#' (S2) is used to identify the outcome (acquisition of a new IMI). Three
#' hypothetical exposures E_q, E_c, E_h (quarter, cow, and herd level) with known
#' strength of association (OR~3.0) are generated. As it is often the case
#' (Dufour et al., 2011, 2012), exposures are equally associated with odds of a
#' prevalent IMI on first milk sample as with odds of IMI acquisition on second
#' sample. Exposures are randomly associated with odds of eliminating an existing
#' IMI (OR=1.0). If S. aureus or CNS is chosen, default parameters are used.
#' Otherwise user has to provide his own.
#'
#' @param n_herd Number of herds.
#' @param n_cow Number of cows per herd.
#' @param bact Type of bacteria: S. aureus, CNS, or other. If other, as to provide parameters.
#' @param E_hPr Exposure distribution (0 to 1) of the binary herd-level (h) predictor. S. aureus and CNS = 0.5.
#' @param E_cPr Exposure distribution (0 to 1) of the binary cow-level (c) predictor. S. aureus and CNS = 0.5.
#' @param E_qPr Exposure distribution (0 to 1) of the binary quarter-level (q) predictor. S. aureus and CNS = 0.5.
#' @param sigma_sqhPr Herd-level variance (sigma_sq) for prevalence of intra-mammary infection (IMI). S. aureus = 0.14; CNS = 0.363.
#' @param sigma_sqcPr Cow-level variance for prevalence of IMI. S. aureus = 2.25; CNS = 0.294.
#' @param b0_Pr Intercept for IMI prevalence; aiming at a prevalence of 2.5%. S. aureus = -6.7; CNS = -2.15.
#' @param OR_hPr OR of association between herd-level variable and IMI prevalence. S. aureus and CNS = 3.
#' @param OR_cPr OR of association between cow-level variable and IMI prevalence. S. aureus and CNS = 3.
#' @param OR_qPr OR of association between observation-level variable and IMI prevalence. S. aureus and CNS = 3.
#' @param sigma_sqhI Herd-level variance for incidence of IMI. S. aureus = 0.838; CNS = 0.27.
#' @param sigma_sqcI Cow-level variance for incidence of IMI. S. aureus = 2.926; CNS = 0.256.
#' @param b0_I Intercept for IMI incidence. S. aureus = -8.3; CNS = -2.4.
#' @param OR_hI OR of association between herd-level variable and IMI incidence. S. aureus and CNS = 3.
#' @param OR_cI OR of association between cow-level variable and IMI incidence. S. aureus and CNS = 3.
#' @param OR_qI OR of association between observation-level variable and IMI incidence. S. aureus and CNS = 3.
#' @param sigma_sqhEl Herd-level variance for elimination of IMI. S. aureus = 0.15; CNS = 0.112.
#' @param sigma_sqcEl Cow-level variance for elimination of IMI. S. aureus = 2.246; CNS = 0.7.
#' @param b0_El Intercept for IMI persistency. S. aureus = -0.6; CNS = 1.6.
#' @param se_parms Vector for mode and x to determine shape parameters of Beta distribution of sensitivity (Se). S. aureus = c(0.90, 0.85); CNS = c(0.60, 0.55).
#' @param sp_parms Vector for mode and x to determine shape parameters of Beta distribution of specificity (Sp). S. aureus = 1; CNS = c(0.95, 0.90).
#' @param se_series Se improvement or loss associated with sampling strategy: duplicate series. S. aureus = -0.1; CNS = -0.25.
#' @param sp_series Sp improvement or loss associated with sampling strategy: duplicate series. S. aureus = 0; CNS = 0.05.
#' @param se_parall Se improvement or loss associated with sampling strategy: duplicate parallel. S. aureus = 0.1; CNS = 0.15.
#' @param sp_parall Sp improvement or loss associated with sampling strategy:duplicate parallel. S. aureus = 0; CNS = -0.05.
#' @param se_tri Se improvement or loss associated with sampling strategy: triplicate (2 out of 3). S. aureus = 0; CNS = 0.
#' @param sp_tri Sp improvement or loss associated with sampling strategy: triplicate (2 out of 3). S. aureus = 0; CNS = 0.10.
#'
#' @return A data frame with variables:
#' \item{herd}{Herd id.}
#' \item{cow}{Cow id.}
#' \item{quarter}{Quarter id.}
#' \item{S1}{First milk sample true status.}
#' \item{S2}{Second milk sample true status.}
#' \item{E_h}{Herd-level exposure.}
#' \item{E_c}{Cow-level exposure.}
#' \item{E_q}{Quarter-level exposure.}
#' \item{S1i}{Misclassified first milk sample.}
#' \item{S2i}{Misclassified second milk sample.}
#' \item{S1_series}{Misclassified first milk sample based on duplicate series sampling strategy.}
#' \item{S2_series}{Misclassified second milk sample based on duplicate series sampling strategy.}
#' \item{S1_parall}{Misclassified first milk sample based on duplicate parallel sampling strategy.}
#' \item{S2_parall}{Misclassified second milk sample based on duplicate parallel sampling strategy.}
#' \item{S1_tri}{Misclassified first milk sample based on triplicate sampling strategy.}
#' \item{S2_tri}{Misclassified second milk sample based on triplicate sampling strategy.}
#'
#' @references Dufour, S., Dohoo, I.R., Barkema, H.W., DesCôteaux, L., DeVries, T.J., Reyher, K.K., Roy, J.-P., Scholl, D.T., 2012 \emph{Epidemiology of coagulase-negative staphylococci intramammary infection in dairy cattle and the effect of bacteriological culture misclassification}. Journal of Dairy Science 95(6):3110-3124.
#' @examples
#' # Initiate a list to store the n data frames
#' sim_list <- vector("list", 5)
#' # Do not forget to set seed for replication
#' set.seed(123)
#' sim_list <- replicate(n = 5, expr = make_data(100, 30, "saureus"), simplify = FALSE)
#' # Or with a progress bar
#' require(pbapply)
#' sim_list <- pbreplicate(n = 5, expr = make_data(100, 30, "cns"), simplify = FALSE)
#' @export
#' @importFrom epiR epi.betabuster
#' @importFrom stats rnorm rbinom
make_data <- function(n_herd,
                      n_cow,
                      bact = c("saureus", "cns", "other"),
                      E_hPr = NULL,
                      E_cPr = NULL,
                      E_qPr = NULL,
                      sigma_sqhPr = NULL,
                      sigma_sqcPr = NULL,
                      b0_Pr = NULL,
                      OR_hPr = NULL,
                      OR_cPr = NULL,
                      OR_qPr = NULL,
                      sigma_sqhI = NULL,
                      sigma_sqcI = NULL,
                      b0_I = NULL,
                      OR_hI = NULL,
                      OR_cI = NULL,
                      OR_qI = NULL,
                      sigma_sqhEl = NULL,
                      sigma_sqcEl = NULL,
                      b0_El = NULL,
                      se_parms = NULL,
                      sp_parms = NULL,
                      se_series = NULL,
                      sp_series = NULL,
                      se_parall = NULL,
                      sp_parall = NULL,
                      se_tri = NULL,
                      sp_tri = NULL) {
    if(missing(bact))
        stop('Please choose type of bacteria.')
    if(!(bact %in% c("saureus", "cns", "other")))
        stop('Please choose between saureus, cns, or other.')
    if(bact == "other" & (is.null(E_hPr) | is.null(E_cPr) | is.null(E_qPr) |
                          is.null(sigma_sqhPr) | is.null(sigma_sqcPr) |
                          is.null(b0_Pr) | is.null(OR_hPr) | is.null(OR_cPr) |
                          is.null(OR_qPr) | is.null(sigma_sqhI) |
                          is.null(sigma_sqcI) | is.null(b0_I) | is.null(OR_hI) |
                          is.null(OR_cI) | is.null(OR_qI) | is.null(sigma_sqhEl) |
                          is.null(sigma_sqcEl) | is.null(b0_El) |
                          is.null(se_parms) | is.null(sp_parms) |
                          is.null(se_series) | is.null(sp_series) |
                          is.null(se_parall) | is.null(sp_parall) |
                          is.null(se_tri) | is.null(sp_tri)))
        stop('Missing parameter(s).')
    
    dt <- matrix(NA, nrow = n_herd*n_cow*4, ncol = 8)
    colnames(dt) <- c("herd", "cow", "quarter",
                      "S1", "S2",
                      "E_h", "E_c", "E_q")
    dt[, 1] <- rep(1:n_herd, each = n_cow*4)
    dt[, 2] <- rep(1:(n_herd*n_cow), each = 4)
    dt[, 3] <- 1:nrow(dt)

    if(bact == "saureus") {
        E_hPr = 0.5
        E_cPr = 0.5
        E_qPr = 0.5
        sigma_sqhPr = 0.14
        sigma_sqcPr = 2.25
        b0_Pr = -6.7
        OR_hPr = 3
        OR_cPr = 3
        OR_qPr = 3
        sigma_sqhI = 0.838
        sigma_sqcI = 2.926
        b0_I = -8.3
        OR_hI = 3
        OR_cI = 3
        OR_qI = 3
        sigma_sqhEl = 0.15
        sigma_sqcEl = 2.246
        b0_El = -0.6
        se_parms = c(0.90, 0.85)
        sp_parms = 1
        se_series = -0.1
        sp_series = 0
        se_parall = 0.1
        sp_parall = 0
        se_tri = 0
        sp_tri = 0
    }
    if(bact == "cns") {
        E_hPr = 0.5
        E_cPr = 0.5
        E_qPr = 0.5
        sigma_sqhPr = 0.363
        sigma_sqcPr = 0.294
        b0_Pr = -2.15
        OR_hPr = 3
        OR_cPr = 3
        OR_qPr = 3
        sigma_sqhI = 0.27
        sigma_sqcI = 0.256
        b0_I = -2.4
        OR_hI = 3
        OR_cI = 3
        OR_qI = 3
        sigma_sqhEl = 0.112
        sigma_sqcEl = 0.7
        b0_El = 1.6
        se_parms = c(0.60, 0.55)
        sp_parms = c(0.95, 0.90)
        se_series = -0.25
        sp_series = 0.05
        se_parall = 0.15
        sp_parall = -0.05
        se_tri = 0
        sp_tri = 0.10
    }
    else {
        E_hPr = E_hPr
        E_cPr = E_cPr
        E_qPr = E_qPr
        sigma_sqhPr = sigma_sqhPr
        sigma_sqcPr = sigma_sqcPr
        b0_Pr = b0_Pr
        OR_hPr = OR_hPr
        OR_cPr = OR_cPr
        OR_qPr = OR_qPr
        sigma_sqhI = sigma_sqhI
        sigma_sqcI = sigma_sqcI
        b0_I = b0_I
        OR_hI = OR_hI
        OR_cI = OR_cI
        OR_qI = OR_qI
        sigma_sqhEl = sigma_sqhEl
        sigma_sqcEl = sigma_sqcEl
        b0_El = b0_El
        se_parms = se_parms
        sp_parms = sp_parms
        se_series = se_series
        sp_series = sp_series
        se_parall = se_parall
        sp_parall = sp_parall
        se_tri = se_tri
        sp_tri = sp_tri
    }

    se <- epiR::epi.betabuster(mode = se_parms[1], conf = 0.80,
                               greaterthan = TRUE, x = se_parms[2],
                               conf.level = 0.95, max.shape1 = 100,
                               step = 0.001)
    if(length(sp_parms == 1))
        sp <- 1
    else sp <- epiR::epi.betabuster(mode = sp_parms[1], conf = 0.80,
                                    greaterthan = TRUE, x = sp_parms[2],
                                    conf.level = 0.95, max.shape1 = 100,
                                    step = 0.001)
    ## Herd level exposure
    herd_sample <- sample(unique(dt[, 1]), size = (length(unique(dt[, 1])) * E_hPr))
    dt[, 6] <- ifelse(dt[, 1] %in% herd_sample, 1, 0)
    ## Cow level exposure
    cow_sample <- sample(unique(dt[, 2]), size = (length(unique(dt[, 2])) * E_cPr))
    dt[, 7] <- ifelse(dt[, 2] %in% cow_sample, 1, 0)
    ## Positive quarters within cows
    quarter_sample <- sample(unique(dt[, 3]), size = (length(unique(dt[, 3])) * E_qPr))
    dt[, 8] <- ifelse(dt[, 3] %in% quarter_sample, 1, 0)

    ## Attribute values to S1 based on Pr distribution, Pr clustering, and Pr
    ## associations with exposures
    ## Draw values
    ##  -for each herd v_0k (level 3 random error term; between-herd variance, normally
    ##   distributed with mean of zero and constant variance sigma^2
    ##   v_k ~ N(0, sigma_sqhPr)), and
    ##  -for each cow u_0jk (level2 random error term; between-cow variance, normally
    ##   distributed with mean of zero and constant variance sigma^2
    ##   u_jk ~ N(0, sigma_sqcPr)
    cluster_herd <- matrix(ncol = 2, nrow = n_herd)
    cluster_herd[,  1] <- unique(dt[, 1])
    cluster_herd[, 2] <- rnorm(n_herd, 0, sqrt(sigma_sqhPr))
    colnames(cluster_herd) <- c("herd", "v_0k")
    dt <- merge(dt, cluster_herd, by = "herd")
    
    cluster_cow <- matrix(ncol = 2, nrow = n_herd * n_cow)
    cluster_cow[,  1] <- unique(dt[, 2])
    cluster_cow[, 2] <- rnorm(n_herd * n_cow, 0, sqrt(sigma_sqcPr))
    colnames(cluster_cow) <- c("cow", "u_0jk")
    dt <- merge(dt, cluster_cow, by = "cow")
    ## Compute probability (pi) that S1 == 1
    dt$pi <- with(dt, 1 /
                      (1 + exp(-1 * (b0_Pr + log(OR_hPr) * E_h + log(OR_cPr) *
                                     E_c + log(OR_qPr) * E_q + u_0jk + v_0k))))
    ## Draw from Bernoulli distribution with pi
    dt$S1 <- rbinom(length(dt$S1), 1, dt$pi)
    
    ## Divide the data set in quarters at risk of new IMI and quarters already infected
    atrisk <- subset(dt, S1 == 0, select = -c(v_0k, u_0jk, pi))
    notatrisk <- subset(dt, S1 == 1, select = -c(v_0k, u_0jk, pi))
    ## Attribute values to S2 for quarter at risk based on Incidence distribution,
    ## Incidence clustering, and Incidence association with exposures
    ## Draw values
    ##  -for each herd v_0k, and
    ##  -for each cow u_0jk
    cluster_herd <- matrix(ncol = 2, nrow = n_herd)
    cluster_herd[,  1] <- unique(dt$herd)
    cluster_herd[, 2] <- rnorm(n_herd, 0, sqrt(sigma_sqhI))
    colnames(cluster_herd) <- c("herd", "v_0k")
    atrisk <- merge(atrisk, cluster_herd, by = "herd")
    
    cluster_cow <- matrix(ncol = 2, nrow = n_herd * n_cow)
    cluster_cow[,  1] <- unique(dt$cow)
    cluster_cow[, 2] <- rnorm(n_herd * n_cow, 0, sqrt(sigma_sqcI))
    colnames(cluster_cow) <- c("cow", "u_0jk")
    atrisk <- merge(atrisk, cluster_cow, by = "cow")
    ## Compute probability (pi) that S2 == 1
    atrisk$pi <- with(atrisk, 1 /
                              (1 + exp(-1 * (b0_I + log(OR_hI) * E_h + log(OR_cI) *
                                             E_c + log(OR_qI) * E_q + u_0jk + v_0k))))
    ## Draw from Bernoulli distribution with pi
    atrisk$S2 <- rbinom(length(atrisk$S2), 1, atrisk$pi)
    ## Attribute values to S2 for quarter NOT at risk based on Elimination distribution
    ## and Elimination clustering
    ## Draw values
    ##  -for each herd v_0k, and
    ##  -for each cow u_0jk
    cluster_herd <- matrix(ncol = 2, nrow = n_herd)
    cluster_herd[,  1] <- unique(dt$herd)
    cluster_herd[, 2] <- rnorm(n_herd, 0, sqrt(sigma_sqhEl))
    colnames(cluster_herd) <- c("herd", "v_0k")
    notatrisk <- merge(notatrisk, cluster_herd, by = "herd")
    
    cluster_cow <- matrix(ncol = 2, nrow = n_herd * n_cow)
    cluster_cow[,  1] <- unique(dt$cow)
    cluster_cow[, 2] <- rnorm(n_herd * n_cow, 0, sqrt(sigma_sqcEl))
    colnames(cluster_cow) <- c("cow", "u_0jk")
    notatrisk <- merge(notatrisk, cluster_cow, by = "cow")
    ## Compute probability (pi) that S2 == 1
    notatrisk$pi <- with(notatrisk, 1 / (1 + exp(-1 * (b0_El + u_0jk + v_0k))))
    ## Draw from Bernoulli distribution with pi
    notatrisk$S2 <- rbinom(length(notatrisk$S2), 1, notatrisk$pi)
    
    dt <- rbind(atrisk, notatrisk)
    dt <- dt[with(dt, order(herd, cow, quarter)), ]
    ## draw Se and Sp values from beta distributions for S1i and S2i
    dt$se <- rbeta(1, se$shape1, se$shape2)
    dt$se <- with(dt, ifelse(se > 0.999999, 0.999999, se))
    if (length(sp_parms == 1))
        dt$sp <- sp_parms
    else dt$sp <- rbeta(1, sp$shape1, sp$shape2)
    if(length(sp_parms) > 1)
        dt$sp <- with(dt, ifelse(sp > 0.999999, 0.999999, sp))
    ## Compute misclassified resuls for S1i and S2i
    dt$S1i <- ifelse(dt$S1 == 1,
                     rbinom(length(dt[dt$S1==1, ]$S1), 1, dt$se),
                     rbinom(length(dt[dt$S1==1, ]$S1), 1, 1 - dt$sp))
    dt$S2i <- ifelse(dt$S2 == 1,
                     rbinom(length(dt[dt$S2==1, ]$S2), 1, dt$se),
                     rbinom(length(dt[dt$S2==1, ]$S2), 1, 1 - dt$sp))
    ## Se and Sp values for duplicate in series
    dt$se_series <- dt$se + se_series
    dt$se_series <- ifelse(dt$se_series > 0.999999, 0.999999, dt$se_series)
    dt$sp_series <- dt$sp + sp_series
    dt$sp_series <- ifelse(dt$sp_series > 0.999999 & dt$sp != 1,
                           0.999999, dt$sp_series)
    ## Compute misclassified resuls for S1i and S2i
    dt$S1_series <- ifelse(dt$S1 == 1,
                           rbinom(length(dt[dt$S1==1, ]$S1), 1, dt$se_series),
                           rbinom(length(dt[dt$S1==1, ]$S1), 1, 1 - dt$sp_series))
    dt$S2_series <- ifelse(dt$S2 == 1,
                           rbinom(length(dt[dt$S2==1, ]$S2), 1, dt$se_series),
                           rbinom(length(dt[dt$S2==1, ]$S2), 1, 1 - dt$sp_series))
    ## Se and Sp values for duplicate in parallel
    dt$se_parall <- dt$se + se_parall
    dt$se_parall <- ifelse(dt$se_parall > 0.999999, 0.999999, dt$se_parall)
    dt$sp_parall <- dt$sp + sp_parall
    dt$sp_parall <- ifelse(dt$sp_parall > 0.999999 & dt$sp != 1,
                           0.999999, dt$sp_parall)
    ## Compute misclassified resuls for S1i and S2i
    dt$S1_parall <- ifelse(dt$S1 == 1,
                           rbinom(length(dt[dt$S1==1, ]$S1), 1, dt$se_parall),
                           rbinom(length(dt[dt$S1==1, ]$S1), 1, 1 - dt$sp_parall))
    dt$S2_parall <- ifelse(dt$S2 == 1,
                           rbinom(length(dt[dt$S2==1, ]$S2), 1, dt$se_parall),
                           rbinom(length(dt[dt$S2==1, ]$S2), 1, 1 - dt$sp_parall))
    ## Se and Sp values for triplicates
    dt$se_tri <- dt$se + se_tri
    dt$se_tri <- ifelse(dt$se_tri > 0.999999, 0.999999, dt$se_tri)
    dt$sp_tri <- dt$sp + sp_tri
    dt$sp_tri <- ifelse(dt$sp_tri > 0.999999 & dt$sp != 1, 0.999999, dt$sp_tri)
    ## Compute misclassified resuls for S1i and S2i
    dt$S1_tri <- ifelse(dt$S1 == 1,
                        rbinom(length(dt[dt$S1==1, ]$S1), 1, dt$se_tri),
                        rbinom(length(dt[dt$S1==1, ]$S1), 1, 1 - dt$sp_tri))
    dt$S2_tri <- ifelse(dt$S2 == 1,
                        rbinom(length(dt[dt$S2==1, ]$S2), 1, dt$se_tri),
                        rbinom(length(dt[dt$S2==1, ]$S2), 1, 1 - dt$sp_tri))
    
    ## Clean-up
    dt <- subset(dt, select = c(herd, cow, quarter,
                                S1, S2, E_h, E_c, E_q,
                                S1i, S2i,
                                S1_series, S2_series,
                                S1_parall, S2_parall,
                                S1_tri, S2_tri))
    return(dt)
}

