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
  expect_equal(calc(w,2),3.14)  
})

test_that("Linear", {
  w = linear.gAlg(1)
  expect_that(w, is_a("gAlg"))
  tab = w@tab
  expect_equal(dim(tab),c(1,4,1))
  expect_equal(tab[1,1,1],Inf)
  expect_equal(tab[1,2,1],0)
  expect_equal(tab[1,3,1],0)
  expect_equal(tab[1,4,1],1)
  expect_equal(calc(w,0),0)  
  expect_equal(calc(w,1),1)  
  expect_equal(calc(w,2),2)  
})

test_that("Mult", {
  w = linear.gAlg(1)
  w = w * w
  expect_that(w, is_a("gAlg"))
  tab = w@tab
  expect_equal(dim(tab),c(1,5,1))
  expect_equal(tab[1,1,1],Inf)
  expect_equal(tab[1,2,1],0)
  expect_equal(tab[1,3,1],0)
  expect_equal(tab[1,4,1],0)
  expect_equal(tab[1,5,1],1)
  expect_equal(calc(w,0),0)  
  expect_equal(calc(w,1),1)  
  expect_equal(calc(w,2),4)  
})

test_that("Add", {
  w = linear.gAlg(1)
  v = const.gAlg(3.14)
  w = w + v
  expect_that(w, is_a("gAlg"))
  tab = w@tab
  expect_equal(dim(tab),c(1,4,2))
  i = order(tab[1,4,])
  tab[] = tab[,,i]
  expect_equal(tab[1,1,1],Inf)
  expect_equal(tab[1,2,1],0)
  expect_equal(tab[1,3,1],3.14)
  expect_equal(tab[1,4,1],0)
  expect_equal(tab[1,1,2],Inf)
  expect_equal(tab[1,2,2],0)
  expect_equal(tab[1,3,2],0)
  expect_equal(tab[1,4,2],1)
  expect_equal(calc(w,0),3.14)  
  expect_equal(calc(w,1),4.14)  
  expect_equal(calc(w,2),5.14)  
})

test_that("Complex", {
  x = linear.gAlg(1)
  l = const.gAlg(3.14)
  w = x*x*x+3*x*x-2*x+l
  v = 7*x*x+x
  expect_that(w, is_a("gAlg"))
  expect_that(v, is_a("gAlg"))
  l = 3.14
  x = seq(0,1,len=100)
  expect_equal(calc(w,x),x*x*x+3*x*x-2*x+l)  
  expect_equal(calc(v,x),7*x*x+x)
  expect_equal(calc(w+v,x),calc(w,x)+calc(v,x))  
  expect_equal(calc(w-v,x),calc(w,x)-calc(v,x))  
  expect_equal(calc(w*v,x),calc(w,x)*calc(v,x))  
})

test_that("Even dim op", {
  d = 30
  x = linear.gAlg(1,d)
  l = const.gAlg(3.14,d)
  w = x*x*x+3*x*x-2*x+l
  v = 7*x*x+x
  expect_that(w, is_a("gAlg"))
  expect_that(v, is_a("gAlg"))
  l = 3.14
  x = seq(0,1,len=100)
  X = cbind(x,matrix(runif(100*(d-1)),100,d-1))
  expect_equal(calc(w,X),x*x*x+3*x*x-2*x+l)  
  expect_equal(calc(v,X),7*x*x+x)
  expect_equal(calc(w+v,X),calc(w,X)+calc(v,X))  
  expect_equal(calc(w-v,X),calc(w,X)-calc(v,X))  
  expect_equal(calc(w*v,X),calc(w,X)*calc(v,X))  
})

