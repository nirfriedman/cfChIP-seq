require(preprocessCore)

QQNormalizeGenes = function(GeneDiff, CommonG = CommonGenes) {
  CommonGenes = CommonG
 # print(length(CommonGenes))
  GeneCounts.CommonGenes = GeneDiff[CommonGenes,]
  GeneCounts.CommonGenes.qq = normalize.quantiles(GeneCounts.CommonGenes)
  dimnames(GeneCounts.CommonGenes.qq) = dimnames(GeneCounts.CommonGenes)
  
  QQNorm = sapply(colnames(GeneDiff), function(s){
    X = GeneCounts.CommonGenes[,s]
    Y = GeneCounts.CommonGenes.qq[,s]
    coef(lm(Y~X-1))
  })
  
  names(QQNorm) = colnames(GeneDiff)
  return(QQNorm)
}
