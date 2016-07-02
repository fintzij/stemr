library(stemr)

context("C++ equivalent of findInterval")

test_that("The intervals are correctly identified for a given vector", {

        x <- c(-1, 0, 0.5, 1, 1.5)
        x_breaks <- c(0,1)

        y <- runif(10, 0, 2)
        y_breaks <- seq(0, 2, by = 0.1)

        expect_identical(find_interval(x, x_breaks, rightmost_closed = FALSE, all_inside = FALSE),
                         findInterval(x, x_breaks, rightmost.closed = FALSE, all.inside = FALSE))
        expect_identical(find_interval(x, x_breaks, rightmost_closed = TRUE, all_inside = FALSE),
                         findInterval(x, x_breaks, rightmost.closed = TRUE, all.inside = FALSE))
        expect_identical(find_interval(x, x_breaks, rightmost_closed = FALSE, all_inside = TRUE),
                         findInterval(x, x_breaks, rightmost.closed = FALSE, all.inside = TRUE))
        expect_identical(find_interval(x, x_breaks, rightmost_closed = TRUE, all_inside = TRUE),
                         findInterval(x, x_breaks, rightmost.closed = TRUE, all.inside = TRUE))

        expect_identical(find_interval(y, y_breaks, rightmost_closed = FALSE, all_inside = FALSE),
                         findInterval(y, y_breaks, rightmost.closed = FALSE, all.inside = FALSE))
        expect_identical(find_interval(y, y_breaks, rightmost_closed = TRUE, all_inside = FALSE),
                         findInterval(y, y_breaks, rightmost.closed = TRUE, all.inside = FALSE))
        expect_identical(find_interval(y, y_breaks, rightmost_closed = FALSE, all_inside = TRUE),
                         findInterval(y, y_breaks, rightmost.closed = FALSE, all.inside = TRUE))
        expect_identical(find_interval(y, y_breaks, rightmost_closed = TRUE, all_inside = TRUE),
                         findInterval(y, y_breaks, rightmost.closed = TRUE, all.inside = TRUE))
})