## Matt Galloway


#' @title Normal Data Generator
#' @description generates random data.
#'
#' @param n desired sample size
#' @param p desired dimension
#' @param theta blah
#' @param var blah
#' @param reps blah
#' @export
#' @examples
#' data_gen(1000, 10, 0.5)



# we first define a function to generate the data note that this function will be
# used in problems 4, 5, 6
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
