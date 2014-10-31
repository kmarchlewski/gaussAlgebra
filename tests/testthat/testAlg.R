library(gaussAlgebra)
context("Simple algebra tests")

test_that("Const", {
  w = const.gAlg(3.14)
  expect_that(w, is_a("gAlg"))
  tab = w@tab
  expect_equal(dim(tab),c(1,3,1))
  expect_equal(tab[1,1,1],Inf)
  expect_equal(tab[1,2,1],0)
  expect_equal(tab[1,3,1],3.14)
  expect_equal(calc(w,0),3.14)  
  expect_equal(calc(w,1),3.14)  
})

