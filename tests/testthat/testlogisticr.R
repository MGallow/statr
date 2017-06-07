
library(statr)

X = as.matrix(dplyr::select(iris, -c(Species, Sepal.Length)))
X1 = cbind(1, X)

y = as.matrix(dplyr::select(iris, Sepal.Length))
y1 = as.matrix(dplyr::select(iris, Species))
y1 = ifelse(y1 == "setosa", 1, 0)

fit = logisticr(X, y1)

test_that("testing logisticr(X, y1)", {
  expect_identical(all(abs(fit$gradient) < 1e-5), TRUE)
})


fit = logisticr(X, y1, lam = 0.1, penalty = "ridge")

test_that("testing logisticr(X, y1, lam = 0.1, penalty = ridge)", {
  expect_identical(all(abs(fit$gradient) < 1e-5), TRUE)
})


fit = logisticr(X, y1, lam = seq(0.1, 2, 0.1), penalty = "ridge")

test_that("testing logisticr(X, y1, lam = seq(), penalty = ridge)", {
  expect_identical(all(abs(fit$gradient) < 1e-5), TRUE)
})


fit = logisticr(X, y1, lam = 0.1, alpha = 1.2, penalty = "bridge")

test_that("testing logisticr(X, y1, lam = 0.1, alpha = 1.2, penalty = bridge)", {
  expect_identical(all(abs(fit$gradient) < 1e-5), TRUE)
})


fit = logisticr(X, y1, lam = seq(0.1, 2, 0.1), alpha = seq(1.1, 1.9, 0.1), penalty = "bridge")

test_that("testing logisticr(X, y1, lam = seq(), alpha = seq(), penalty = bridge)", {
  expect_identical(all(abs(fit$gradient) < 1e-5), TRUE)
})


fit = logisticr(X1, y1)

test_that("testing logisticr(X1, y1)", {
  expect_identical(all(abs(fit$gradient) < 1e-5), TRUE)
})


fit = logisticr(X1, y1, lam = 0.1, penalty = "ridge")

test_that("testing logisticr(X1, y1, lam = 0.1, penalty = ridge)", {
  expect_identical(all(abs(fit$gradient) < 1e-5), TRUE)
})


fit = logisticr(X1, y1, lam = seq(0.1, 2, 0.1), penalty = "ridge")

test_that("testing logisticr(X1, y1, lam = seq(), penalty = ridge)", {
  expect_identical(all(abs(fit$gradient) < 1e-5), TRUE)
})


fit = logisticr(X1, y1, lam = 0.1, alpha = 1.2, penalty = "bridge")

test_that("testing logisticr(X1, y1, lam = 0.1, alpha = 1.2, penalty = bridge)", {
  expect_identical(all(abs(fit$gradient) < 1e-5), TRUE)
})


fit = logisticr(X1, y1, lam = seq(0.1, 2, 0.1), alpha = seq(1.1, 1.9, 0.1), penalty = "bridge")

test_that("testing logisticr(X1, y1, lam = seq(), alpha = seq(), penalty = bridge)", {
  expect_identical(all(abs(fit$gradient) < 1e-5), TRUE)
})