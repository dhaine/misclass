\name{runProgram}
\alias{runProgram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
runProgram(data, iter = 1000, warmup = 100, chains = 4, cores, seed = 123)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{iter}{
%%     ~~Describe \code{iter} here~~
}
  \item{warmup}{
%%     ~~Describe \code{warmup} here~~
}
  \item{chains}{
%%     ~~Describe \code{chains} here~~
}
  \item{cores}{
%%     ~~Describe \code{cores} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
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
function (data, iter = 1000, warmup = 100, chains = 4, cores, 
    seed = 123) 
{
    if (warmup >= iter) 
        stop("Number of warmup samples is greater than number of iterations.")
    if (cores > chains) 
        stop("Number of cores should not exceed number of chains.")
    stanfit <- stanmodels$logistic
    x <- model.matrix(S2 ~ E_h, data)
    x.mean <- colMeans(x)
    stan_data <- list(N = nrow(data), y = data$S2, K = 2, X = x, 
        X_means = as.array(x.mean), J_1 = data$cow, N_1 = length(unique(data$cow)), 
        K_1 = 1, Z_1 = rep(1, nrow(data)), J_2 = data$herd, N_2 = length(unique(data$herd)), 
        K_2 = 1, Z_2 = rep(1, nrow(data)))
    true_ass <- rstan::stan(file = stanfit, data = stan_data, 
        iter = iter, warmup = warmup, chains = chains, seed = seed, 
        cores = cores)
    return(true_ass)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
