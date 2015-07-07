ratlist=function(...)UseMethod("ratlist")
ratlist.rationalz=function(...)
{
	lst=list(...)
	cls=sapply(lst, class)
	if(!all(cls=='rationalz')) {
		if(all(grepl('^rational[qz]$', cls))) {
			do.call('ratlist', lapply(lst, as.rationalq))
		}else stop('unknown object classes')
	}else structure(lst, class='ratlist')
}
ratlist.rationalq=function(...)
{
	lst=list(...)
	cls=sapply(lst, class)
	if(!all(cls=='rationalq')) {
		if(all(grepl('^rational[qz]$', cls))) {
			do.call('ratlist', lapply(lst, as.rationalq))
		}else stop('unknown object classes')
	}else structure(lst, class='ratlist')
}

ratlist.polyzlist=function(
	numerList, 
	denomList=rep(polynomialz(1), length.out=length(numerList)), 
	...)
{
	l1=length(numerList); l2=length(denomList)
	if(l1!=l2){
		L=max(l1,l2); 
		numerList=rep(numerList, length.out=L)
		denomList=rep(denomList, length.out=L)
	}
	structure(mapply(rational, numerList, denomList, SIMPLIFY=FALSE, ...), class='ratlist')
}
ratlist.polyqlist=function(
	numerList, 
	denomList=rep(polynomialq(1), length.out=length(numerList)), 
	...) NULL
	body(ratlist.polyqlist) = body(ratlist.polyzlist)


numerator.ratlist=function(x)
{
	ans = lapply(x, '[[', 'numerator')
	to.call=if(inherits(ans[[1]], 'polynomialz')) {
			'polyzlist' 
		} else if(inherits(ans[[1]], 'polynomialq')) {
			'polyqlist'
		} else stop('unknown numerator class')
	do.call(to.call, ans)
}

denominator.ratlist=function(x)
{
	ans = lapply(x, '[[', 'denominator')
	to.call=if(inherits(ans[[1]], 'polynomialz')) {
			'polyzlist' 
		} else if(inherits(ans[[1]], 'polynomialq')) {
			'polyqlist'
		} else stop('unknown denominator class')
	do.call(to.call, ans)
}

if(FALSE){
.commonDenom=function(num.polylist, denom.polylist) ## finds n[1]/d[1] + n[2]/d[2] + ... ### need to optimize for speed
{  
	Ln=length(num.polylist); Ld=length(denom.polylist)
	if(Ln!=Ld){Ld=Ln=max(Ld,Ln); num.polylist=rep(num.polylist, length.out=Ln); denom.polylist=rep(denom.polylist, length.out=Ld); }
	if(Ln==1) return(rational(num.polylist[[1L]], denom.polylist[[1L]]))
	denom.ans = LCM(denom.polylist)
	num.ans = polynomial(0)
	for(i in seq(Ln)){
		num.ans = num.ans + num.polylist[[i]] * (denom.ans / denom.polylist[[i]])
	}
	rational(num.ans, denom.ans)
}
}

.sum1ratlist=function(x) {
	if(sum(!duplicated(sapply(x, 'class')))) x = lapply(x, as.rationalq)
	ans0 = x[[1L]]
	if(length(x)==1L) return(ans0)
	Reduce("+", x[-1L], ans0)
}	
.prod1ratlist=function(x) {
	if(sum(!duplicated(sapply(x, 'class')))) x = lapply(x, as.rationalq)
	ans0 = x[[1L]]
	if(length(x)==1L) return(ans0)
	Reduce("*", x[-1L], ans0)
}	

Summary.ratlist = function(..., na.rm=FALSE)
{
    ok <- switch(.Generic, sum = , prod = TRUE, FALSE)
    if(!ok)
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.", .Generic, .Class))
	lst = c(...)
	
	switch(.Generic, 
		"sum" = {
			all.rat=lapply(list(...), .sum1ratlist)
			ans0 = all.rat[[1]]
			if(length(all.rat)==1L) return(ans0)
			Reduce("+", all.rat[-1L], ans0)	
		}, 
		"prod" = {
			all.rat=lapply(list(...), .prod1ratlist)
			ans0 = all.rat[[1]]
			if(length(all.rat)==1L) return(ans0)
			Reduce("*", all.rat[-1L], ans0)	
		}
	)
}

if(FALSE){
`+.ratlist`=function(e1, e2)  ## not elementwise addition, but sum(e1) + sum(e2)
{
	if(missing(e2)) return(.commonDenom(e1$numerator, e1$denominator))
	if(inherits(e2, 'polynomial')){
		e2 = ratlist(polylist(e2))
	}else if(inherits(e2, 'polylist')){
		e2 = ratlist(e2)
	}else if(mode(e2)=='numeric') {
		e2 = ratlist(do.call('polylist',as.list(e2)))
	}
	
	.commonDenom(c(e1$numerator, e2$numerator), 
				c(e1$denominator, e2$denominator))
}
`-.ratlist`=function(e1, e2)
{
	if(missing(e2)) return(.commonDenom(e1$numerator * (-1), e1$denominator))
	.commonDenom(c(e1$numerator, e2$numerator * -1 ), 
				 c(e1$denominator, e2$denominator))
}
}


#coef.polynomial is already defined in polynom
