solve.polynomialq <-function(a, b, method='polyroot', ...)
{
	method = match.arg(method, c('polyroot', 'eigen', 'bracket'))
	
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
	}else if(method=='bracket') {
		.NotYetImplemented()
	}
}

uBound = function(p, method) UseMethod('uBound')

uBound.polynomialq=function(p, method=c('Cauchy', 'Lagrange', 'Kojima','Fujiwara','SumAdjRatio'))
{
	e1=trimZeros(coef(p))
	n=length(e1)
	method=match.arg(method, several.ok=TRUE)
	nmethod=length(method)
	if(nmethod==0L)stop('"method needs to be one of c("Cauchy", "Lagrange", "Kojima","Fujiwara","SumAdjRatio")')
	bnd=vector('list' , nmethod)
	names(bnd)=method
	if('Cauchy'%in% method)	bnd[["Cauchy"]] = bq1 + max (abs(e1[-n]/e1[n]))
	if('Lagrange'%in% method) bnd[["Lagrange"]] = max(c(bq1, sum(abs(e1[-n]/e1[n]))))
	if('Kojima'%in% method) bnd[["Kojima"]] = if(any(e1==bq0)) Inf else 2 * max(abs(e1[-n]/e1[-1L])*c(.5, rep(1,max(0,n-2))))
	if('Fujiwara'%in% method) bnd[["Fujiwara"]] = 2 * max((abs(e1[-n]/e1[n])*c(.5, rep(1,max(0,n-2))))^(1/safeseq(n-1,1,by=-1)))
	if('SumAdjRatio'%in% method) bnd[["SumAdjRatio"]] = if(any(e1[-1]==0)) Inf else sum(abs(e1[-n]/e1[-1]))
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

