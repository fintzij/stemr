library(stemr)

context("Computing the incidence")

test_that("Incidence counts are correctly recovered", {

        set.seed(12511)

        x <- matrix(rexp(100), ncol = 10); x[,6:10] <- apply(x[,6:10], 2, cumsum)
        y <- matrix(rexp(100), ncol = 10); y[,6:10] <- apply(y[,6:10], 2, cumsum)
        z <- matrix(rexp(100), ncol = 10); z[,6:10] <- apply(z[,6:10], 2, cumsum)

        x_0 <- x; y_0 <- y; z_0 <- z

        x_0[,6:10] <- apply(x_0[,6:10], 2, function(a) c(a[1], diff(a)))
        y_0[c(2,4,6,8,10), 6:10] <- 0
        y_0[c(1,3,5,7,9), 6:10] <- apply(y_0[c(1,3,5,7,9), 6:10] , 2, function(a) c(a[1], diff(a)))
        z_0[c(1,3,5,7,9), 6:10] <- 0
        z_0[c(2,4,6,8,10), 6:10] <- apply(z_0[c(2,4,6,8,10), 6:10] , 2, function(a) c(a[1], diff(a)))

        compute_incidence(x, 5:9, rep(list(0:9), 5))                 # all rows are census rows
        compute_incidence(y, 5:9, rep(list(seq(0,9, by=2)), 5))      # census at odd rows
        compute_incidence(z, 5:9, rep(list(seq(1,9, by=2)), 5))      # census at even rows

        expect_identical(x_0, x)
        expect_identical(y_0, y)
        expect_identical(z_0, z)
})