solve.polynomialq <-function(a, b, method='polyroot', ...)
{
	method = match.arg(method, c('polyroot', 'eigen', 'bisection'))
	
	if(!missing(b)) a <- a - b
	if(method=='eigen'){
		class(a) = 'bigq'
		a1=trimZeros(a, 'leading')
		r=numeric(length(a)-length(a1))
		a=a1

		switch(as.character(length(a)),
			   "0" =,
			   "1" = r,
			   "2" = sort(c(r,  as.numeric(- a[1L]/a[2L]))),
		   {
		   a <- rev(a)
		   a <- as.numeric( (a/a[1L])[-1L] )
		   M <- rbind( - a, cbind(diag(length(a) - 1), 0))
		   sort(c(r, eigen(M, symmetric = FALSE,
							   only.values = TRUE)$values))
		   })
	}else if(method=='polyroot')	{
		class(a) = 'bigq'
		a = as.numeric( a / sum(abs(a)) * length(a) )
		sort(polyroot(a))
	}else if(method=='bisection') {
		iso = realRootIso(a, ...)
		nroot=length(iso)
		ans=do.call('c', lapply(iso, sum))
		ans = ans /2
		attr(ans, 'intervals') = iso
		ans
	}
}

uBound = function(p, method) UseMethod('uBound')
uBound.methods=c('Cauchy', 'Lagrange', 'Kojima','Fujiwara','SumAdjRatio','Kalantari05','Deutsch',paste0('A',1:31),paste0('B',1:14))
uBound.polynomialq=function(p, method=c('Cauchy', 'Lagrange', 'Kojima','Fujiwara','SumAdjRatio','Deutsch'))
{
	e1=trimZeros(coef(p))
	method=match.arg(method, uBound.methods, several.ok=TRUE)
	nmethod=length(method)
	
	if(nmethod==0L)stop(paste("method needs to be one of", paste(dQuote(uBound.methods),collapse=',')))
	bnd=vector('list' , nmethod)
	names(bnd)=method
	
	n=length(e1)
	Co = e1/e1[n]
	
	for(m in method) bnd[[m]]=switch(m, 
		Cauchy = bq1 + max (abs(e1[-n]/e1[n])), 
		Lagrange = max(c(bq1, sum(abs(e1[-n]/e1[n])))),
		Kojima = if(any(e1==bq0)) Inf else 2 * max(abs(e1[-n]/e1[-1L])*c(.5, rep(1,max(0,n-2)))), 
		Fujiwara = 2 * max((abs(e1[-n]/e1[n])*c(.5, rep(1,max(0,n-2))))^(1/safeseq(n-1,1,by=-1))),
		SumAdjRatio = if(any(e1[-1]==0)) Inf else sum(abs(e1[-n]/e1[-1])),
		A11 	=,
		Deutsch = if(any(e1[-1]==0)) Inf else abs(Co[n-1]) + max(abs(Co[-n]/Co[-1])),
		Inf
	)
	nrslt = sapply(bnd, length)
	if(all(nrslt==0L)) {
			ans = Inf  # as.bigz(Inf) 
	}else ans=min(do.call('c', bnd))
	attr(ans, 'bounds')=bnd
	ans
}


lBound = function(p, method) UseMethod('lBound')
lBound.polynomialq=function(p, method=c('Rouche'))
{
	e1=trimZeros(coef(p))
	n=length(e1)
	method=match.arg(method, several.ok=TRUE)
	nmethod=length(method)
	bnd=vector('list', nmethod)
	names(bnd)=method
	if('Rouche'%in% method) bnd$Rouche = max(abs(e1[1L])/(abs(e1[1L]) + max(abs(e1[-1L]))), 
											 abs(e1[1L])/max(abs(e1[1L]),sum(abs(e1[-1L]))))
	nrslt = sapply(bnd, length)
	if(all(nrslt==0L)){
		ans = as.bigz(0L)
	}else ans = max(do.call('c', bnd))
	attr(ans, 'bounds')=bnd
	ans
}

sturm = function(p) UseMethod('sturm')
sturm.polynomialq = function(p)  ## this is extremely slow beyond dozens of degree
{
	e1 = trimZeros(p)
	dif = deriv(e1)
	
	ans=polyqlist(e1, dif)
	last=dif; last2=e1
	repeat{
		tmp = -(last2 %% last)
		#if(length(tmp)==0) tmp=zero
		ans=c(ans, tmp)
		if(degree(tmp)==0L) break
		#if(length(tmp)==0L) stop('not square free?')
		last2=last
		last=tmp
	}
	ans
}

squareFree =function(p, ...) UseMethod('squareFree')
squareFree.polynomialq = function(p, ...)
{
	e1=trimZeros(p)
	#if(length(e1)==1L) return(e1/e1)
	#e1=e1/e1[length(e1)]
	
	a0=GCD(e1, deriv(e1))
	b1=e1 %/% a0
	c1=deriv(e1) %/% a0
	d1=c1 - deriv(b1)
	ans=polyqlist()
	repeat{
		ai=GCD(b1, d1)
		b2=b1 %/% ai
		c2=d1 %/% ai
		ans = c(ans,  ai)
		d1=c2 - deriv(b2)
		b1=trimZeros(b2); c1=c2
		if(degree(b1)==0L && b1==1) break
	}
	class(ans) = 'polyqlist'
	ans
}

nRealRoots = function(p, lower, upper, ...) UseMethod('nRealRoots')
nRealRoots.polynomialq = function(p, lower=-upper, upper=uBound(p)*(1+as.bigq(1L,1e22)), method='Sturm', ...)
{
	stopifnot(upper>lower)
	sqFree=squareFree(p)
	sqFree=sqFree[sapply(sqFree,degree)>=1L]
	
	
	if(length(sqFree)>1L) {
		return(sum(sapply(sqFree, nRealRoots, lower,upper, method,...)))
	}
	
	e1=sqFree[[1L]]
	method=match.arg(tolower(method), 'sturm')
	if(method!='sturm') .NotYetImplemented()
	### sturm's theorem
		sturm.seq=sturm(e1)
		evals = predict(sturm.seq, c(as.bigq(lower), as.bigq(upper)), SIMPLIFY=FALSE)
		nlower=decartes(do.call('c', lapply(evals, '[[', i=1)))
		nupper=decartes(do.call('c', lapply(evals, '[[', i=2)))
		return(nlower - nupper)
}

realRootIso = function(p, ...) UseMethod('realRootIso')
realRootIso.polynomialq=function(p, lower=-upper, upper=uBound(p)*(1+as.bigq(1L,1e22)), eps=as.bigq(1, 1e3), method='bisect',...)
{
	## bisection based on Sturm's theorem
	if(is.infinite(lower)) lower = -uBound(p)*(1+as.bigq(1L,1e22))
	if(is.infinite(upper)) upper =  uBound(p)*(1+as.bigq(1L,1e22))
	eps = as.bigq(eps)
	lower = as.bigq(lower)
	upper = as.bigq(upper)

	nRealRoots=nRealRoots(p, lower, upper)
	if(nRealRoots==0) return(list())

	if(upper - lower < eps && nRealRoots==1) return(list(c(lower, upper)))
	
	mid=(lower+upper)/2
	c(realRootIso(p, lower, mid, eps,...), realRootIso(p, mid, upper, eps))
}
