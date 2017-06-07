
library(statr)

X = as.matrix(dplyr::select(iris, -c(Species, Sepal.Length)))
X1 = cbind(1, X)

y = as.matrix(dplyr::select(iris, Sepal.Length))
y1 = as.matrix(dplyr::select(iris, Species))
y1 = ifelse(y1 == "setosa", 1, 0)

fit = linearr(X, y)

test_that("testing linearr(X, y)", {
  expect_identical(all(abs(fit$gradient) < 1e-5), TRUE)
})

fit = linearr(X, y, method = "MM")

test_that("testing linearr(X, y, method = MM)", {
  expect_identical(all(abs(fit$gradient) < 1e-5), TRUE)
})


fit = linearr(X, y, lam = 0.1, penalty = "ridge")

test_that("testing linearr(X, y, lam = 0.1, penalty = ridge)", {
  expect_identical(all(abs(fit$gradient) < 1e-5), TRUE)
})


fit = linearr(X, y, lam = seq(0.1, 2, 0.1), penalty = "ridge")

test_that("testing linearr(X, y, lam = seq(), penalty = ridge)", {
  expect_identical(all(abs(fit$gradient) < 1e-5), TRUE)
})


fit = linearr(X, y, lam = 0.1, alpha = 1.2, penalty = "bridge")

test_that("testing linearr(X, y, lam = 0.1, alpha = 1.2, penalty = bridge)", {
  expect_identical(all(abs(fit$gradient) < 1e-5), TRUE)
})


fit = linearr(X, y, lam = seq(0.1, 2, 0.1), alpha = seq(1.1, 1.9, 0.1), penalty = "bridge")

test_that("testing linearr(X, y, lam = seq(), alpha = seq(), penalty = bridge)", {
  expect_identical(all(abs(fit$gradient) < 1e-5), TRUE)
})


fit = linearr(X1, y)

test_that("testing linearr(X1, y)", {
  expect_identical(all(abs(fit$gradient) < 1e-5), TRUE)
})


fit = linearr(X1, y, method = "MM")

test_that("testing linearr(X1, y, method = MM)", {
  expect_identical(all(abs(fit$gradient) < 1e-5), TRUE)
})


fit = linearr(X1, y, lam = 0.1, penalty = "ridge")

test_that("testing linearr(X1, y, lam = 0.1, penalty = ridge)", {
  expect_identical(all(abs(fit$gradient) < 1e-5), TRUE)
})


fit = linearr(X1, y, lam = seq(0.1, 2, 0.1), penalty = "ridge")

test_that("testing linearr(X1, y, lam = seq(), penalty = ridge)", {
  expect_identical(all(abs(fit$gradient) < 1e-5), TRUE)
})


fit = linearr(X1, y, lam = 0.1, alpha = 1.2, penalty = "bridge")

test_that("testing linearr(X1, y, lam = 0.1, alpha = 1.2, penalty = bridge)", {
  expect_identical(all(abs(fit$gradient) < 1e-5), TRUE)
})


fit = linearr(X1, y, lam = seq(0.1, 2, 0.1), alpha = seq(1.1, 1.9, 0.1), penalty = "bridge")

test_that("testing linearr(X1, y, lam = seq(), alpha = seq(), penalty = bridge)", {
  expect_identical(all(abs(fit$gradient) < 1e-5), TRUE)
})