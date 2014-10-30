
setClass("NP", representation(tab="array"))
NP = function(tab) {
  new("NP", tab=tab);
}

poly = function(...) {
  A = c(Inf,0,...)
  dim(A) = c(1,length(A),1)
  NP(A)
}

const.NP = function(a=1,d=1) {
  A = cbind(Inf,0,c(a,rep(1,d-1)))
  dim(A) = c(d,3,1)
  NP(A)
}

linear.NP = function(k,d=1) {
  if ((k > d) || (k<1)) stop("linear.NP: a must: 1 <= k <= d")
  A = rep(0, d*4)
  dim(A) = c(d,4,1)
  A[,1,]=Inf
  A[,2,]=0
  A[,3,]=1
  A[k,3,]=0	
  A[k,4,]=1	
  NP(A)
}


D = function(sd=1,mean=0) {
  d = length(sd)
  if (length(mean) == 1) mean = rep(mean,d)
  B = c(sd**2,mean,dnorm(0,sd=sd))
  dim(B) = c(d,3,1)
  NP(B)
}

calc.m <- function(A,x) {
  if (is.vector(x)) x = as.matrix(x)
  if (dim(A)[1] != dim(x)[2]) stop("dimension mismatch!");
  dims = c(dim(A)[1], dim(A)[2] - 3, dim(A)[3], dim(x)[1])
  ndims = dim(x)[1]
  ret = .C("calc",as.integer(dims),A,x,ret=double(ndims),DUP=F,NAOK=T)
  ret = ret$ret
  ret
}

setMethod("calc", signature("NP","numeric"), function(object,x) calc.m(object@tab,x))
setMethod("calc", signature("NP","array"), function(object,x) calc.m(object@tab,x))


mult.NP <- function(A,B,mult=TRUE) {
  print(mult)
  if (dim(A)[1] != dim(B)[1]) stop("dimension mismatch!");
  dims = c(dim(A)[1], dim(A)[2] - 3, dim(A)[3], dim(B)[2]-3, dim(B)[3], ifelse(mult,1,0))
  ndims = c(dims[1],3+dims[2]+dims[4],dims[3]*dims[5])
  ret = .C("conv",as.integer(dims),A,B,fg=double(ndims[1]*ndims[2]*ndims[3]),DUP=F,NAOK=T)
  ret = ret$fg
  dim(ret)=ndims
  NP(ret)
}

setMethod("%%", signature("NP","NP"), function(e1,e2) mult.NP(e1@tab,e2@tab,mult=FALSE) )
setMethod("*",  signature("NP","NP"), function(e1,e2) mult.NP(e1@tab,e2@tab,mult=TRUE) )

#setMethod("Ops", signature("NP","NP"), Ops.gvector)


plus.NP <- function(A,B,plus=TRUE) {
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
  NP(AB)
}

setMethod("+", signature("NP","NP"), function(e1,e2)  plus.NP(e1@tab,e2@tab,plus=TRUE ) )
setMethod("-",  signature("NP","NP"), function(e1,e2) plus.NP(e1@tab,e2@tab,plus=FALSE) )

setMethod("lag", signature("NP"), function(x,dx,...) {
  if (length(dx) != dim(x@tab)[1]) stop("dimension mismatch!");
  x@tab[,2,] = x@tab[,2,] + dx
  x
})

setMethod("Arith", signature("numeric","NP"), function(e1,e2) {
  e1v = const.NP(e1,dim(e2@tab)[1])
  callGeneric(e1v,e2)
})
setMethod("Arith", signature("NP","numeric"), function(e1,e2) {
  e2v = const.NP(e2,dim(e1@tab)[1])
  callGeneric(e1,e2v)
})

setMethod("*", signature("numeric","NP"), function(e1,e2) {
  k = dim(e2@tab)[2]
  e2@tab[1,3:k,] = e2@tab[1,3:k,] * e1
  e2
})

setMethod("*", signature("NP","numeric"), function(e1,e2) {
  k = dim(e1@tab)[2]
  e1@tab[1,3:k,] = e1@tab[1,3:k,] * e2
  e1
})
