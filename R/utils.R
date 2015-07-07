bz0 = as.bigz(0L)
bq0 = as.bigq(0L)
bz1 = as.bigz(1L)
bq1 = as.bigq(1L)
trimZeros = function(x, end='trailing', empty.OK=TRUE) UseMethod('trimZeros')
trimZeros.default=function(x, end='trailing', empty.OK=TRUE)
{
	end = match.arg(end, c('leading', 'trailing', 'both'))
	trailing = any(end==c('trailing','both'))
	
	if(trailing) x=rev(x)
	nzeros=sum(cumsum(abs(x))==0)
    ans = if(nzeros>0) x[-seq(nzeros)] else x
	if(!empty.OK && length(ans)==0L) ans = 0
	if(trailing) ans = rev(ans)
	if(end=='both') Recall(ans, 'leading', empty.OK) else ans
}
trimZeros.bigz=function(x, end='trailing', empty.OK=TRUE)
{
	end = match.arg(end, c('leading', 'trailing', 'both'))
	trailing = any(end==c('trailing','both'))

	if(trailing) x=rev(x)
	nzeros=sum(cumsum(abs(x))==bz0)
    ans = if(nzeros>0) x[-seq(nzeros)] else x
	if(!empty.OK && length(ans)==0L) ans = bz0
	if(trailing) ans = rev(ans)
	if(end=='both') Recall(ans, 'leading', empty.OK) else ans
}
trimZeros.bigq=function(x, end='trailing', empty.OK=TRUE)
{
	end = match.arg(end, c('leading', 'trailing', 'both'))
	trailing = any(end==c('trailing','both'))

	if(trailing) x=rev(x)
	nzeros=sum(cumsum(abs(x))==bq0)
    ans = if(nzeros>0) x[-seq(nzeros)] else x
	if(!empty.OK && length(ans)==0L) ans = bq0
	if(trailing) ans = rev(ans)
	if(end=='both') Recall(ans, 'leading', empty.OK) else ans
}

if(FALSE){
## accumulate a la Abelson and Sussman.
accumulate <- function(f, init, x, right = TRUE) {
    if(length(x) == 0)
        return(init)
    f <- match.fun(f)
    if(right)
        f(x[[1]], Recall(f, init, x[-1], right = TRUE))
    else
        Recall(f, f(init, x[[1]]), x[-1], right = FALSE)
}
}
