Bonferroni=function(pvalue,alpha=0.05){
  N=length(pvalue)
  ps=alpha/N
  return(list(validate=which(pvalue<ps),alphaCorrige=ps))
}

Benjamini=function(pvalue,alpha=0.05){
  N=length(pvalue)
  psort=sort(pvalue)
  ben=alpha/N+(0:(N-1))*alpha/N
  ps=ben[which(psort>ben)[1]]
  print(ps)
  return(list(validate=which(pvalue<ps),alphaCorrige=ps))
}
