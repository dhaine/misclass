\name{makeData}
\alias{makeData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makeData(n_herd, n_cow, E_hPr = 0.5, E_cPr = 0.5, E_qPr = 0.5, sigma_sqhPr = 0.14, sigma_sqcPr = 2.25, b0_Pr = -6.7, OR_hPr = 3, OR_cPr = 3, OR_qPr = 3, sigma_sqhI = 0.838, sigma_sqcI = 2.926, b0_I = -8.3, OR_hI = 3, OR_cI = 3, OR_qI = 3, sigma_sqhEl = 0.15, sigma_sqcEl = 2.246, b0_El = -0.6, se_parms = c(0.9, 0.85), sp_parms = 1, se_series = -0.1, sp_series = 0, se_parall = 0.1, sp_parall = 0, se_tri = 0, sp_tri = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n_herd}{
%%     ~~Describe \code{n_herd} here~~
}
  \item{n_cow}{
%%     ~~Describe \code{n_cow} here~~
}
  \item{E_hPr}{
%%     ~~Describe \code{E_hPr} here~~
}
  \item{E_cPr}{
%%     ~~Describe \code{E_cPr} here~~
}
  \item{E_qPr}{
%%     ~~Describe \code{E_qPr} here~~
}
  \item{sigma_sqhPr}{
%%     ~~Describe \code{sigma_sqhPr} here~~
}
  \item{sigma_sqcPr}{
%%     ~~Describe \code{sigma_sqcPr} here~~
}
  \item{b0_Pr}{
%%     ~~Describe \code{b0_Pr} here~~
}
  \item{OR_hPr}{
%%     ~~Describe \code{OR_hPr} here~~
}
  \item{OR_cPr}{
%%     ~~Describe \code{OR_cPr} here~~
}
  \item{OR_qPr}{
%%     ~~Describe \code{OR_qPr} here~~
}
  \item{sigma_sqhI}{
%%     ~~Describe \code{sigma_sqhI} here~~
}
  \item{sigma_sqcI}{
%%     ~~Describe \code{sigma_sqcI} here~~
}
  \item{b0_I}{
%%     ~~Describe \code{b0_I} here~~
}
  \item{OR_hI}{
%%     ~~Describe \code{OR_hI} here~~
}
  \item{OR_cI}{
%%     ~~Describe \code{OR_cI} here~~
}
  \item{OR_qI}{
%%     ~~Describe \code{OR_qI} here~~
}
  \item{sigma_sqhEl}{
%%     ~~Describe \code{sigma_sqhEl} here~~
}
  \item{sigma_sqcEl}{
%%     ~~Describe \code{sigma_sqcEl} here~~
}
  \item{b0_El}{
%%     ~~Describe \code{b0_El} here~~
}
  \item{se_parms}{
%%     ~~Describe \code{se_parms} here~~
}
  \item{sp_parms}{
%%     ~~Describe \code{sp_parms} here~~
}
  \item{se_series}{
%%     ~~Describe \code{se_series} here~~
}
  \item{sp_series}{
%%     ~~Describe \code{sp_series} here~~
}
  \item{se_parall}{
%%     ~~Describe \code{se_parall} here~~
}
  \item{sp_parall}{
%%     ~~Describe \code{sp_parall} here~~
}
  \item{se_tri}{
%%     ~~Describe \code{se_tri} here~~
}
  \item{sp_tri}{
%%     ~~Describe \code{sp_tri} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (n_herd, n_cow, E_hPr = 0.5, E_cPr = 0.5, E_qPr = 0.5, 
    sigma_sqhPr = 0.14, sigma_sqcPr = 2.25, b0_Pr = -6.7, OR_hPr = 3, 
    OR_cPr = 3, OR_qPr = 3, sigma_sqhI = 0.838, sigma_sqcI = 2.926, 
    b0_I = -8.3, OR_hI = 3, OR_cI = 3, OR_qI = 3, sigma_sqhEl = 0.15, 
    sigma_sqcEl = 2.246, b0_El = -0.6, se_parms = c(0.9, 0.85), 
    sp_parms = 1, se_series = -0.1, sp_series = 0, se_parall = 0.1, 
    sp_parall = 0, se_tri = 0, sp_tri = 0) 
{
    dt <- matrix(NA, nrow = n_herd * n_cow * 4, ncol = 8)
    colnames(dt) <- c("herd", "cow", "quarter", "S1", "S2", "E_h", 
        "E_c", "E_q")
    dt[, 1] <- rep(1:n_herd, each = n_cow * 4)
    dt[, 2] <- rep(1:(n_herd * n_cow), each = 4)
    dt[, 3] <- 1:nrow(dt)
    se <- epiR::epi.betabuster(mode = se_parms[1], conf = 0.8, 
        greaterthan = TRUE, x = se_parms[2], conf.level = 0.95, 
        max.shape1 = 100, step = 0.001)
    if (length(sp_parms == 1)) 
        sp <- 1
    else sp <- epiR::epi.betabuster(mode = sp_parms[1], conf = 0.8, 
        greaterthan = TRUE, x = sp_parms[2], conf.level = 0.95, 
        max.shape1 = 100, step = 0.001)
    herd_sample <- sample(unique(dt[, 1]), size = (length(unique(dt[, 
        1])) * E_hPr))
    dt[, 6] <- ifelse(dt[, 1] \%in\% herd_sample, 1, 0)
    cow_sample <- sample(unique(dt[, 2]), size = (length(unique(dt[, 
        2])) * E_cPr))
    dt[, 7] <- ifelse(dt[, 2] \%in\% cow_sample, 1, 0)
    quarter_sample <- sample(unique(dt[, 3]), size = (length(unique(dt[, 
        3])) * E_qPr))
    dt[, 8] <- ifelse(dt[, 3] \%in\% quarter_sample, 1, 0)
    cluster_herd <- matrix(ncol = 2, nrow = n_herd)
    cluster_herd[, 1] <- unique(dt[, 1])
    cluster_herd[, 2] <- rnorm(n_herd, 0, sqrt(sigma_sqhPr))
    colnames(cluster_herd) <- c("herd", "v_0k")
    dt <- merge(dt, cluster_herd, by = "herd")
    cluster_cow <- matrix(ncol = 2, nrow = n_herd * n_cow)
    cluster_cow[, 1] <- unique(dt[, 2])
    cluster_cow[, 2] <- rnorm(n_herd * n_cow, 0, sqrt(sigma_sqcPr))
    colnames(cluster_cow) <- c("cow", "u_0jk")
    dt <- merge(dt, cluster_cow, by = "cow")
    dt$pi <- with(dt, 1/(1 + exp(-1 * (b0_Pr + log(OR_hPr) * 
        E_h + log(OR_cPr) * E_c + log(OR_qPr) * E_q + u_0jk + 
        v_0k))))
    dt$S1 <- rbinom(length(dt$S1), 1, dt$pi)
    atrisk <- subset(dt, S1 == 0, select = -c(v_0k, u_0jk, pi))
    notatrisk <- subset(dt, S1 == 1, select = -c(v_0k, u_0jk, 
        pi))
    cluster_herd <- matrix(ncol = 2, nrow = n_herd)
    cluster_herd[, 1] <- unique(dt$herd)
    cluster_herd[, 2] <- rnorm(n_herd, 0, sqrt(sigma_sqhI))
    colnames(cluster_herd) <- c("herd", "v_0k")
    atrisk <- merge(atrisk, cluster_herd, by = "herd")
    cluster_cow <- matrix(ncol = 2, nrow = n_herd * n_cow)
    cluster_cow[, 1] <- unique(dt$cow)
    cluster_cow[, 2] <- rnorm(n_herd * n_cow, 0, sqrt(sigma_sqcI))
    colnames(cluster_cow) <- c("cow", "u_0jk")
    atrisk <- merge(atrisk, cluster_cow, by = "cow")
    atrisk$pi <- with(atrisk, 1/(1 + exp(-1 * (b0_I + log(OR_hI) * 
        E_h + log(OR_cI) * E_c + log(OR_qI) * E_q + u_0jk + v_0k))))
    atrisk$S2 <- rbinom(length(atrisk$S2), 1, atrisk$pi)
    cluster_herd <- matrix(ncol = 2, nrow = n_herd)
    cluster_herd[, 1] <- unique(dt$herd)
    cluster_herd[, 2] <- rnorm(n_herd, 0, sqrt(sigma_sqhEl))
    colnames(cluster_herd) <- c("herd", "v_0k")
    notatrisk <- merge(notatrisk, cluster_herd, by = "herd")
    cluster_cow <- matrix(ncol = 2, nrow = n_herd * n_cow)
    cluster_cow[, 1] <- unique(dt$cow)
    cluster_cow[, 2] <- rnorm(n_herd * n_cow, 0, sqrt(sigma_sqcEl))
    colnames(cluster_cow) <- c("cow", "u_0jk")
    notatrisk <- merge(notatrisk, cluster_cow, by = "cow")
    notatrisk$pi <- with(notatrisk, 1/(1 + exp(-1 * (b0_El + 
        u_0jk + v_0k))))
    notatrisk$S2 <- rbinom(length(notatrisk$S2), 1, notatrisk$pi)
    dt <- rbind(atrisk, notatrisk)
    dt <- dt[with(dt, order(herd, cow, quarter)), ]
    dt$se <- rbeta(1, se$shape1, se$shape2)
    dt$se <- with(dt, ifelse(se > 0.999999, 0.999999, se))
    if (length(sp_parms == 1)) 
        dt$sp <- sp_parms
    else dt$sp <- rbeta(1, sp$shape1, sp$shape2)
    if (length(sp_parms) > 1) 
        dt$sp <- with(dt, ifelse(sp > 0.999999, 0.999999, sp))
    dt$S1i <- ifelse(dt$S1 == 1, rbinom(length(dt[dt$S1 == 1, 
        ]$S1), 1, dt$se), rbinom(length(dt[dt$S1 == 1, ]$S1), 
        1, 1 - dt$sp))
    dt$S2i <- ifelse(dt$S2 == 1, rbinom(length(dt[dt$S2 == 1, 
        ]$S2), 1, dt$se), rbinom(length(dt[dt$S2 == 1, ]$S2), 
        1, 1 - dt$sp))
    dt$se_series <- dt$se + se_series
    dt$se_series <- ifelse(dt$se_series > 0.999999, 0.999999, 
        dt$se_series)
    dt$sp_series <- dt$sp + sp_series
    dt$sp_series <- ifelse(dt$sp_series > 0.999999 & dt$sp != 
        1, 0.999999, dt$sp_series)
    dt$S1_series <- ifelse(dt$S1 == 1, rbinom(length(dt[dt$S1 == 
        1, ]$S1), 1, dt$se_series), rbinom(length(dt[dt$S1 == 
        1, ]$S1), 1, 1 - dt$sp_series))
    dt$S2_series <- ifelse(dt$S2 == 1, rbinom(length(dt[dt$S2 == 
        1, ]$S2), 1, dt$se_series), rbinom(length(dt[dt$S2 == 
        1, ]$S2), 1, 1 - dt$sp_series))
    dt$se_parall <- dt$se + se_parall
    dt$se_parall <- ifelse(dt$se_parall > 0.999999, 0.999999, 
        dt$se_parall)
    dt$sp_parall <- dt$sp + sp_parall
    dt$sp_parall <- ifelse(dt$sp_parall > 0.999999 & dt$sp != 
        1, 0.999999, dt$sp_parall)
    dt$S1_parall <- ifelse(dt$S1 == 1, rbinom(length(dt[dt$S1 == 
        1, ]$S1), 1, dt$se_parall), rbinom(length(dt[dt$S1 == 
        1, ]$S1), 1, 1 - dt$sp_parall))
    dt$S2_parall <- ifelse(dt$S2 == 1, rbinom(length(dt[dt$S2 == 
        1, ]$S2), 1, dt$se_parall), rbinom(length(dt[dt$S2 == 
        1, ]$S2), 1, 1 - dt$sp_parall))
    dt$se_tri <- dt$se + se_tri
    dt$se_tri <- ifelse(dt$se_tri > 0.999999, 0.999999, dt$se_tri)
    dt$sp_tri <- dt$sp + sp_tri
    dt$sp_tri <- ifelse(dt$sp_tri > 0.999999 & dt$sp != 1, 0.999999, 
        dt$sp_tri)
    dt$S1_tri <- ifelse(dt$S1 == 1, rbinom(length(dt[dt$S1 == 
        1, ]$S1), 1, dt$se_tri), rbinom(length(dt[dt$S1 == 1, 
        ]$S1), 1, 1 - dt$sp_tri))
    dt$S2_tri <- ifelse(dt$S2 == 1, rbinom(length(dt[dt$S2 == 
        1, ]$S2), 1, dt$se_tri), rbinom(length(dt[dt$S2 == 1, 
        ]$S2), 1, 1 - dt$sp_tri))
    dt <- subset(dt, select = c(herd, cow, quarter, S1, S2, E_h, 
        E_c, E_q, S1i, S2i, S1_series, S2_series, S1_parall, 
        S2_parall, S1_tri, S2_tri))
    return(dt)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
