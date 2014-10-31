#' Package for efficient analitical calculations on polynomial-times-gaussian functions
#'
#' @docType package
#' @name gaussAlgebra-package
#' @rdname gaussAlgebra-package
#' @useDynLib gaussAlgebra
NULL

#' Class that represents a polynomial-times-gaussian function
#' 
#' @slot tab Array keeping all the information about the function (see below)
#' @exportClass gAlg
gAlg = setClass("gAlg", representation(tab="array"))

#' Transforms polynomial coefficients into a 1D gAlg function
#' 
#' @param ... coefficients of a polynomial
#' @export
poly = function(...) {
  A = c(Inf,0,...)
  dim(A) = c(1,length(A),1)
  gAlg(tab=A)
}

#' Generates a constant gAlg function
#' 
#' Creates a gAlg function representing: f(x1,x2,...xd) = a
#' 
#' @param a Value of the constant function
#' @param d Dimension in which the function lives
#' @export
const.gAlg = function(a=1,d=1) {
  A = cbind(Inf,0,c(a,rep(1,d-1)))
  dim(A) = c(d,3,1)
  gAlg(tab=A)
}

#' Generates a linear gAlg function
#' 
#' Creates a gAlg function representing: f(x1,x2,...xd) = xk
#' 
#' @param d Dimension in which the function lives
#' @param k Index of the dimension in which the function is linear
#' @export
linear.gAlg = function(k,d=1) {
  if ((k > d) || (k<1)) stop("linear.gAlg: a must: 1 <= k <= d")
  A = rep(0, d*4)
  dim(A) = c(d,4,1)
  A[,1,]=Inf
  A[,2,]=0
  A[,3,]=1
  A[k,3,]=0	
  A[k,4,]=1	
  gAlg(tab=A)
}


#' Generates a Gauss gAlg function
#' 
#' Creates a gAlg function representing: f(x1,x2,...xd) = Normal distribution(sd,mean)
#' 
#' @param sd Vector of standard deviations (for each dimension)
#' @param mean Vector of means (for each dimension)
#' @export
Gauss = function(sd=1,mean=0) {
  d = length(sd)
  if (length(mean) == 1) mean = rep(mean,d)
  if (length(mean) != d) stop("Gauss: size of mean and sd is different")
  B = c(sd**2,mean,dnorm(0,sd=sd))
  dim(B) = c(d,3,1)
  gAlg(tab=B)
}

calc.gAlg <- function(A,x) {
  if (is.vector(x)) x = as.matrix(x)
  if (! is.double(x)) if (is.numeric(x)) x[] = as.double(x)
  if (dim(A)[1] != dim(x)[2]) stop("dimension mismatch!");
  dims = c(dim(A)[1], dim(A)[2] - 3, dim(A)[3], dim(x)[1])
  ndims = dim(x)[1]
  ret = .C("calc",as.integer(dims),A,x,ret=double(ndims),DUP=F,NAOK=T,PACKAGE="gaussAlgebra")
  ret = ret$ret
  ret
}

#' Evaluates the gAlg function in specific points
#' 
#' @param object gAlg function
#' @param x points in which to evaluate the function
#' @exportMethod calc
setMethod("calc", signature("gAlg","numeric"), function(object,x) calc.gAlg(object@tab,x))
setMethod("calc", signature("gAlg","array"), function(object,x) calc.gAlg(object@tab,x))


mult.gAlg <- function(A,B,mult=TRUE) {
  if (dim(A)[1] != dim(B)[1]) stop("dimension mismatch!");
  dims = c(dim(A)[1], dim(A)[2] - 3, dim(A)[3], dim(B)[2]-3, dim(B)[3], ifelse(mult,1,0))
  ndims = c(dims[1],3+dims[2]+dims[4],dims[3]*dims[5])
  ret = .C(
    "conv",
    as.integer(dims),
    A,
    B,
    fg=double(ndims[1]*ndims[2]*ndims[3]),
    DUP=FALSE,
    NAOK=TRUE,
    PACKAGE="gaussAlgebra"
  )
  ret = ret$fg
  dim(ret)=ndims
  gAlg(tab=ret)
}

setMethod("%%", signature("gAlg","gAlg"), function(e1,e2) mult.gAlg(e1@tab,e2@tab,mult=FALSE) )
setMethod("*",  signature("gAlg","gAlg"), function(e1,e2) mult.gAlg(e1@tab,e2@tab,mult=TRUE) )

#setMethod("Ops", signature("gAlg","gAlg"), Ops.gvector)

plus.gAlg <- function(A,B,plus=TRUE) {
  d = dim(A)[1]
  if (d != dim(B)[1]) stop("dimension mismatch!");
  o = max(dim(A)[2],dim(B)[2])
  dims = c(d,o,dim(A)[3]+dim(B)[3]);
  AB = rep(0,prod(dims))
  dim(AB) = dims
  if (!plus) {
    B[,3:dim(B)[2],] = -B[,3:dim(B)[2],]
  }
  AB[1:d,1:dim(A)[2],1:dim(A)[3]] = A
  AB[1:d,1:dim(B)[2],1:dim(B)[3] + dim(A)[3]] = B
  gAlg(tab=AB)
}

setMethod("+", signature("gAlg","gAlg"), function(e1,e2)  plus.gAlg(e1@tab,e2@tab,plus=TRUE ) )
setMethod("-",  signature("gAlg","gAlg"), function(e1,e2) plus.gAlg(e1@tab,e2@tab,plus=FALSE) )

setMethod("lag", signature("gAlg"), function(x,dx,...) {
  if (length(dx) != dim(x@tab)[1]) stop("dimension mismatch!");
  x@tab[,2,] = x@tab[,2,] + dx
  x
})

setMethod("Arith", signature("numeric","gAlg"), function(e1,e2) {
  e1v = const.gAlg(e1,dim(e2@tab)[1])
  callGeneric(e1v,e2)
})
setMethod("Arith", signature("gAlg","numeric"), function(e1,e2) {
  e2v = const.gAlg(e2,dim(e1@tab)[1])
  callGeneric(e1,e2v)
})

setMethod("*", signature("numeric","gAlg"), function(e1,e2) {
  k = dim(e2@tab)[2]
  e2@tab[1,3:k,] = e2@tab[1,3:k,] * e1
  e2
})

setMethod("*", signature("gAlg","numeric"), function(e1,e2) {
  k = dim(e1@tab)[2]
  e1@tab[1,3:k,] = e1@tab[1,3:k,] * e2
  e1
})

as.character.gAlg = function(x,...)
{
  km = dim(x@tab)[2]-2
  d  = dim(x@tab)[1]
  n  = dim(x@tab)[3]
  a = paste0("x",1:d)
  r = outer(a,1:km-1,function(x,y) paste0(x,"^",y,"*"))
  r[,1] = ""
  if (km > 1) r[,2] = paste0(a,"*")
  
  ret = sapply(1:n,function(k) {
    ret = sapply(1:d,function(i) {
      p = x@tab[i,1:km+2,k]
      nr = paste0(r[i,],p)
      ret = paste(nr[p != 0],collapse=" + ")
      if (is.finite(x@tab[i,1,k])) {
        ret = paste0("(",ret,")*exp(-",a[i],"^2/(2*",x@tab[i,1,k],"))")
      }
      if (x@tab[i,2,k] != 0) {
        ret = paste0(ret," [ origin: ",a[i],"=",x@tab[i,2,k]," ]")
      }
      ret
    })
    paste0(ret,collapse="*\n")
  })
  paste0(ret,collapse="\n+++++++++\n")
}

#' Prints a gAlg function in a nice readable form
#' 
#' @export
print.gAlg = function(object) {
  cat(as.character(object))
}

#' @export
setMethod("show", "gAlg", print.gAlg)
