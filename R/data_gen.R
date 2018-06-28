## Matt Galloway From STAT 8054 HW1 assignment problem 4


#' @title Normal Linear Data Generator
#' @description True beta values are generated from p independent draws from N(0, 1/p) distribution. X_{-1} are n independent draws from (p - 1) multivariate normal N(0, Sigma) where Sigma has (j, k) entry theta^abs(j - k).
#'
#' Y is then generated using the X = (1, X_{-1}) and true beta values with an iid error term that follows distribution N(0, var). We can specify the desired number of replications (reps).
#'
#' @param n desired sample size
#' @param p desired dimension
#' @param r number of responses
#' @param sparsity desired sparsity for beta
#' @param Sigma covariance matrix structure used to generate Y | X
#' @param s option to specify diagonal elements in Sigma
#' @param SigmaX covariance matrix structure used to generate data X
#' @param sx option to specify diagonal elements in SigmaX
#' @param var variance of generated y values
#' @param reps number of replications
#' @param ... additional arguments to pass to data generating functions
#' @return generated design matrix (X), response values (Y)(matrix if reps > 1), true beta values
#' @export

# we define the data generation function
data_gen = function(n, p, r = 1, sparsity = 0.5, Sigma = c("tridiag", "dense", "denseQR", "compound"), s = NULL, SigmaX = c("tridiag", "dense", "denseQR", "compound"), sx = NULL, ...) {

    # checks
    SigmaX = match.arg(SigmaX)
    Sigma = match.arg(Sigma)

    # randomly generate betas
    betas = matrix(rnorm(p*r, 0, sqrt(1/p)), nrow = p, ncol = r)
    betas = betas*matrix(rbinom(p*r, 1, prob = sparsity), nrow = p, ncol = r)

    # generate data X
    SigmaX = switch(SigmaX,
      "tridiag" = tridiag(p = p, n = n, ...),
      "dense" = dense(p = p, n = n, ...),
      "denseQR" = denseQR(p = p, n = n, ...),
      "compound" = compound(p = p, n = n))
    X = SigmaX$X

    SigmaX = SigmaX$S
    if (!is.null(sx)){
      diag(SigmaX) = sx
    }

    # generate sigma matrix
    Sigma = switch(Sigma,
      "tridiag" = tridiag(p = r, n = n, ...),
      "dense" = dense(p = r, n = n, ...),
      "denseQR" = denseQR(p = r, n = n, ...),
      "compound" = compound(p = r, n = n))

    Sigma = Sigma$S
    if (!is.null(s)){
      diag(Sigma) = s
    }

    # create data
    Z = matrix(rnorm(n*r), nrow = n, ncol = r)
    out = eigen(Sigma, symmetric = TRUE)
    Sigma.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
    Y = X %*% betas + Z %*% Sigma.sqrt


    returns = list(Y = Y, X = X, betas = betas, Sigma = Sigma, SigmaX = SigmaX)
    return(returns)

}



##-----------------------------------------------------



#' @title Generate tri-diagonal matrices
#' @description Generate p-dimensional matrices so that its inverse is tri-diagonal.
#' @param p desired dimension
#' @param base base multiplier
#' @param n option to generate n observations from covariance matrix S
#' @export
#' @examples
#' tridiag(p = 10, base = 0.7)


# we define the tridiag function
tridiag = function(p = 8, base = 0.7, n = NULL) {

    # generate tapered matrices
    S = matrix(0, nrow = p, ncol = p)

    for (i in 1:p) {
        for (j in 1:p) {
            S[i, j] = base^abs(i - j)
        }
    }

    # oracle
    Omega = qr.solve(S)

    # create data, if specified
    if (!is.null(n)) {

        # generate n by p matrix X with rows drawn iid N_p(0, sigma)
        Z = matrix(rnorm(n * p), nrow = n, ncol = p)
        out = eigen(S, symmetric = TRUE)
        S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
        X = Z %*% S.sqrt

        return(list(Omega = Omega, S = S, X = X))

    } else {

        return(list(Omega = Omega, S = S))

    }

}




##-----------------------------------------------------



#' @title Generate dense matrices
#' @description Generate p-dimensional matrices so that its inverse is dense.
#' @param p desired dimension
#' @param base base multiplier
#' @param n option to generate n observations from covariance matrix S
#' @export
#' @examples
#' dense(p = 10, base = 0.9)


# we define the dense function
dense = function(p = 8, base = 0.9, n = NULL) {

    # generate matrix
    S = matrix(base, nrow = p, ncol = p)
    diag(S) = 1

    # oracle
    Omega = qr.solve(S)

    # create data, if specified
    if (!is.null(n)) {

        # generate n by p matrix X with rows drawn iid N_p(0, sigma)
        Z = matrix(rnorm(n * p), nrow = n, ncol = p)
        out = eigen(S, symmetric = TRUE)
        S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
        X = Z %*% S.sqrt

        return(list(Omega = Omega, S = S, X = X))

    } else {

        return(list(Omega = Omega, S = S))

    }

}




##-----------------------------------------------------



#' @title Generate dense matrices (via spectral decomposition)
#' @description Generate p-dimensional matrices so that its inverse is dense. The matrix will be generated so its first 'num' eigen values are 1000 and the remaining are 1. The orthogonal basis is generated via QR decomposition of
#' @param p desired dimension
#' @param num number of 'large' eigen values. Note num must be less than p
#' @param n option to generate n observations from covariance matrix S
#' @export
#' @examples
#' denseQR(p = 10, num = 10)


# we define the dense function
denseQR = function(p = 8, num = 5, n = NULL) {

    # generate eigen values
    eigen = c(rep(1000, num), rep(1, p - num))

    # randomly generate orthogonal basis (via QR)
    Q = matrix(rnorm(p * p), nrow = p, ncol = p) %>% qr %>% qr.Q

    # generate matrix
    S = Q %*% diag(eigen) %*% t(Q)

    # oracle
    Omega = qr.solve(S)

    # create data, if specified
    if (!is.null(n)) {

        # generate n by p matrix X with rows drawn iid N_p(0, sigma)
        Z = matrix(rnorm(n * p), nrow = n, ncol = p)
        out = eigen(S, symmetric = TRUE)
        S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
        X = Z %*% S.sqrt

        return(list(Omega = Omega, S = S, X = X))

    } else {

        return(list(Omega = Omega, S = S))

    }

}





##-----------------------------------------------------



#' @title Generate compound symmetric matrices
#' @description Generate a p-dimensional compound symmetric matrix.
#' @param p desired dimension
#' @param n option to generate n observations from covariance matrix S
#' @export
#' @examples
#' compound(p = 10, n = 100)


# we define the dense function
compound = function(p = 8, n = NULL) {

    # generate precision matrix
    Omega = matrix(0.9, nrow = p, ncol = p)
    diag(Omega) = 1

    # generate covariance matrix
    S = qr.solve(Omega)

    # create data, if specified
    if (!is.null(n)) {

        # generate n by p matrix X with rows drawn iid N_p(0, sigma)
        Z = matrix(rnorm(n * p), nrow = n, ncol = p)
        out = eigen(S, symmetric = TRUE)
        S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
        X = Z %*% S.sqrt

        return(list(Omega = Omega, S = S, X = X))

    } else {

        return(list(Omega = Omega, S = S))

    }

}
