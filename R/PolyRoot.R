trimZeros=function(vec)
{
  nzeros=sum(cumsum(abs(rev(vec)))==0)
  if(nzeros>0) vec[-seq(length(vec), length=nzeros, by=-1L)] else vec
}

eigenPolyRoot=function(co, only.real=FALSE)
{
  co=trimZeros(co)
  N=length(co)
  mat=Matrix::Matrix(0, N-1L, N-1L, sparse=TRUE)
  if(N==3){
	mat[2L,1L]=1
  }else if(N==2L){
	ans = -co[1L]/co[2L]
	return(if(isTRUE(only.real)) ans else as.complex(ans))
  }else if(N==1L){
	return(if(isTRUE(only.real)) numeric(0) else complex(0))
  }else diag(mat[-1L,])=1
  mat[,N-1L]=-co[-N]/co[N]
  rts=eigen(mat, symm=FALSE, only=TRUE)$value ## actually eigen converts mat to a dense matrix
  rts=rts[order(abs(Im(rts)), -(Re(rts)))]
  if(isTRUE(only.real)) Re(rts[Im(rts)==0]) else rts
}


numPolyVar.default=function(e1, ...)
{# numeric or bigq
	co=sign(e1[e1!=0])
	sum(co[-1]*co[-length(co)]<0)
}
numPolyVar.list=function(e1, at, ...)
{
	e1=do.call(c, lapply(e1, numPolyEval, at=at))
	NextMethod('numPolyVar')	
}
numPolyVar=function(e1, ...) UseMethod('numPolyVar')


numPolyRealRootIso=function(e1, lower=-upper, upper=numPolyRootUBound(e1)*1.001, eps=1e-3)
{## bisection based on Sturm's theorem
	if(is.infinite(lower)) lower = -numPolyRootUBound(e1)*1.001
	if(is.infinite(upper)) upper =  numPolyRootUBound(e1)*1.001

	nroots=numPolyNRealRoot(e1, lower, upper)
	if(nroots==0) return(list())

	if(upper - lower < eps && nroots==1) return(list(c(lower, upper)))
	
	mid=.5*(lower+upper)
#	ans=list()
#	if(numPolyNRealRoot(e1, lower, mid) >=1) ans=c(ans, Recall(e1, lower, mid, eps))
#	if(numPolyNRealRoot(e1, mid, upper) >=1) ans=c(ans, Recall(e1, mid, upper, eps))
#	return(ans)
	numPolyRealRootIso.Recall=eval(match.call()[[1]], parent.frame()) ## for profiling purposes
	c(numPolyRealRootIso.Recall(e1, lower, mid, eps),numPolyRealRootIso.Recall(e1, mid, upper, eps))
}
