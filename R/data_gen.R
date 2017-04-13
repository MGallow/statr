## Matt Galloway
## From STAT 8054 HW1 assignment


#' @title Normal Data Generator
#' @description True beta values are generated from p independent draws from N(0, 1/p) distribution. X_{-1} are n independent draws from (p - 1) multivariate normal N(0, Sigma) where Sigma has (j, k) entry theta^abs(j - k).
#'
#' Y is then generated using the X = (1, X_{-1}) and true beta values with an iid error term that follows distribution N(0, var). We can specify the desired number of replications (reps).
#'
#' @param n desired sample size
#' @param p desired dimension
#' @param theta parameter used to generate covariance matrix
#' @param var variance of generated y values
#' @param reps number of replications
#' @return generated design matrix (X), response values (Y)(matrix if reps > 1), true beta values
#' @export
#' @examples
#' data_gen(1000, 10, 0.5)



# we define the data generation function
data_gen = function(n, p, theta, var = 0.5, reps = 200) {

    # randomly generate betas
    betas = rnorm(p, 0, sqrt(1/p))

    # generate sigma matrix
    Sigma = matrix(rep(0, (p - 1)^2), ncol = p - 1)

    for (i in 1:(p - 1)) {
        for (j in 1:(p - 1)) {

            Sigma[i, j] = theta^abs(i - j)

        }
    }

    # use Choleskys Decomp to generate X columns
    Z = matrix(rnorm(n * (p - 1)), ncol = n)
    L = t(chol(Sigma))
    X_ = L %*% Z

    # generate full design matrix
    X1 = rep(1, n)
    X = t(rbind(X1, X_))

    # generate matrix of random noise (epsilons)(n x replications) note that we
    # generate for all replications at once
    Eps = matrix(rnorm(n * reps, 0, sqrt(var)), ncol = reps)

    # finally, generate y values
    XB = X %*% matrix(betas)
    Y = as.numeric(XB) + Eps


    returns = list(Y = Y, X = X, betas = betas)
    return(returns)

}
