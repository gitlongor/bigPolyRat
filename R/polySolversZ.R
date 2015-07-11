solve.polynomialz <-function(a, b, method='polyroot', ...)
{
	method = match.arg(method, c('polyroot', 'eigen', 'bisection'))
	
	if(!missing(b)) a <- a - b
	if(method=='eigen'){
		class(a) = 'bigz'
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
		class(a) = 'bigz'
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

#uBound = function(p, method) UseMethod('uBound')
#uBound.methods=c('Cauchy', 'Lagrange', 'Kojima','Fujiwara','SumAdjRatio','Kalantari05','Deutsch',paste0('A',1:31),paste0('B',1:14))
uBound.polynomialz=uBound.polynomialq #=function(p, method=c('Cauchy', 'Lagrange', 'Kojima','Fujiwara','SumAdjRatio','Deutsch'))
# {
	# e1=trimZeros(coef(p))
	# method=match.arg(method, uBound.methods, several.ok=TRUE)
	# nmethod=length(method)
	
	# if(nmethod==0L)stop(paste("method needs to be one of", paste(dQuote(uBound.methods),collapse=',')))
	# bnd=vector('list' , nmethod)
	# names(bnd)=method
	
	# n=length(e1)
	# Co = e1/e1[n]
	
	# for(m in method) bnd[[m]]=switch(m, 
		# Cauchy = bq1 + max (abs(e1[-n]/e1[n])), 
		# Lagrange = max(c(bq1, sum(abs(e1[-n]/e1[n])))),
		# Kojima = if(any(e1==0)) Inf else 2 * max(abs(e1[-n]/e1[-1L])*c(.5, rep(1,max(0,n-2)))), 
		# Fujiwara = as.bigq(2 * max(as.numeric(abs(e1[-n]/e1[n])*c(.5, rep(1,max(0,n-2))))^(1/safeseq(n-1,1,by=-1)))),
		# SumAdjRatio = if(any(e1[-1]==0)) Inf else sum(abs(e1[-n]/e1[-1])),
		# A11 	=,
		# Deutsch = if(any(e1[-1]==0)) Inf else abs(Co[n-1]) + max(abs(Co[-n]/Co[-1])),
		# Inf
	# )
	# nrslt = sapply(bnd, length)
	# if(all(nrslt==0L)) {
			# ans = Inf  # as.bigz(Inf) 
	# }else ans=min(do.call('c', bnd))
	# attr(ans, 'bounds')=bnd
	# ans
# }


lBound = function(p, method) UseMethod('lBound')
lBound.polynomialz=lBound.polynomialq#=function(p, method=c('Rouche'))
# {
	# e1=trimZeros(coef(p))
	# n=length(e1)
	# method=match.arg(method, several.ok=TRUE)
	# nmethod=length(method)
	# bnd=vector('list', nmethod)
	# names(bnd)=method
	# if('Rouche'%in% method) bnd$Rouche = max(abs(e1[1L])/(abs(e1[1L]) + max(abs(e1[-1L]))), 
											 # abs(e1[1L])/max(abs(e1[1L]),sum(abs(e1[-1L]))))
	# nrslt = sapply(bnd, length)
	# if(all(nrslt==0L)){
		# ans = as.bigz(0L)
	# }else ans = max(do.call('c', bnd))
	# attr(ans, 'bounds')=bnd
	# ans
# }

#sturm = function(p) UseMethod('sturm')
sturm.polynomialz = function(p)  {
	p=as.polynomialq(p)
	.Class = c(.Class, 'polynomialq')
	NextMethod(.Generic, p)
}
















#squareFree =function(p, ...) UseMethod('squareFree')
squareFree.polynomialz = function(p, ...)  {
	p=as.polynomialq(p)
	.Class = c(.Class, 'polynomialq')
	NextMethod(.Generic, p)
}


















#nRealRoots = function(p, lower, upper, ...) UseMethod('nRealRoots')
nRealRoots.polynomialz = function(p, lower=-upper, upper=uBound(p)*(1+as.bigz(1L,1e22)), method='Sturm', ...){
	p=as.polynomialq(p)
	.Class = c(.Class, 'polynomialq')
	NextMethod(.Generic, p)
}


















#realRootIso = function(p, ...) UseMethod('realRootIso')
realRootIso.polynomialz=function(p, lower=-upper, upper=uBound(p)*(1+as.bigz(1L,1e22)), eps=as.bigq(1, 1e3), method='bisect',...){
	p=as.polynomialq(p)
	.Class = c(.Class, 'polynomialq')
	NextMethod(.Generic, p)
}