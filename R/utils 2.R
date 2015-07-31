safeseq = function (from = 1L, to = 1L, by = 1L, ...) 
{
  disc = by * (from - to)
  if (disc > 0) {
    vector(class(disc), 0L)
  }
  else seq(from = from, to = to, by = by, ...)
}

bz0 = as.bigz(0L)
bq0 = as.bigq(0L)
bz1 = as.bigz(1L)
bq1 = as.bigq(1L)
bzn1 = as.bigz(-1L)
bqn1 = as.bigq(-1L)
bf0 = mpfr(0L,20)
bf1 = mpfr(1L,20)
bfn1 = mpfr(-1L,20)
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

trimZeros.polynomialz = function(x, end='trailing', empty.OK=TRUE) 
  polynomialz(trimZeros(coef(x), end, empty.OK))

trimZeros.polynomialq = function(x, end='trailing', empty.OK=TRUE) 
  polynomialq(trimZeros(coef(x), end, empty.OK))


trimZeros.mpfr=function(x, end='trailing', empty.OK=TRUE)
{
  end = match.arg(end, c('leading','trailing','both'))
  trailing = any(end==c('trailing','both'))
  
  if(trailing) x=rev(x)
  nzeros=sum(cumsum(abs(x))==bf0)
  ans = if(nzeros>0) x[-seq(nzeros)] else x
  if(!empty.OK && length(ans)==0L) ans = bf0
  if(trailing) ans = rev(ans)
  if(end=='both') Recall(ans, 'leading', empty.OK) else ans
  
}
trimZeros.polynomialf = function(x, end='trailing', empty.OK=TRUE) 
  polynomialf(trimZeros(coef(x), end, empty.OK))
