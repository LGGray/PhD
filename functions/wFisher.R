#' @title wFisher
#' @description sample size-weighted Fisher's method
#' @param p A numeric vector of p-values
#' @param weight A numeric vector of weight or sample size for each experiment
#' @param is.onetail Logical. If set TRUE, p-values are combined without considering the direction of effects, and vice versa. Default: TRUE.
#' @param eff.sign A vector of signs of effect sizes. It works when is.onetail = FALSE
#' @return p : Combined p-value
#' @return overall.eff.direction : The direction of combined effects.
#' @examples wFisher(p=c(0.01,0.2,0.8), weight = c(50,60,100),is.onetail=FALSE, eff.sign=c(1,1,1))
#' @importFrom "stats" "pgamma" "qgamma"
#' @references Becker BJ (1994). “Combining significance levels.” In Cooper H, Hedges LV (eds.), A handbook of research synthesis, 215–230. Russell Sage, New York.
#' @references Fisher RA (1925). Statistical methods for research workers. Oliver and Boyd, Edinburgh.
#' @export
# https://www.nature.com/articles/s41598-021-86465-y

wFisher = function(p, weight = NULL, is.onetail = TRUE, eff.sign)
{
  if(is.null(weight)){weight = rep(1, length(p))}
  idx.na = which(is.na(p))
  if(length(idx.na)>0){
    p = p[-idx.na];
    weight = weight[-idx.na];
    if(!is.onetail)
    {
      eff.sign = eff.sign[-idx.na]
    }
  }
  NP = length(p)
  NS = length(weight)
  if(NP!=NS){stop("The length of p and weight vector must be identical.")}
  N = NS
  Ntotal = sum(weight)
  ratio = weight/Ntotal
  Ns = N*ratio
  G = c()

  if(is.onetail)
  {
    for(i in 1:length(p))
    {
      G = append(G, qgamma(p = p[i], shape = Ns[i], scale=2, lower.tail=F))
    }
    Gsum = sum(G)
    resultP = pgamma(q=Gsum, shape=N, scale=2, lower.tail=F)
  }else{
    p1 = p2 = p
    idx_pos = which(eff.sign > 0)
    idx_neg = which(eff.sign < 0)
    # positive direction
    G = c()
    p1[idx_pos] = p[idx_pos]/2
    p1[idx_neg] = 1-p[idx_neg]/2
    for(i in 1:length(p1))
    {
      G = append(G, qgamma(p = p1[i], shape = Ns[i], scale=2, lower.tail=F))
    }
    Gsum = sum(G)
    resultP1 = pgamma(q=Gsum, shape=N, scale=2, lower.tail=F)
    # negative direction
    G = c()
    p2[idx_pos] = 1-p[idx_pos]/2
    p2[idx_neg] = p[idx_neg]/2
    for(i in 1:length(p2))
    {
      G = append(G, qgamma(p = p2[i], shape = Ns[i], scale=2, lower.tail=F))
    }
    Gsum = sum(G)
    resultP2 = pgamma(q=Gsum, shape=N, scale=2, lower.tail=F)
    resultP = 2* min(resultP1, resultP2)
    if(resultP > 1.0){resultP = 1.0}
    overall.eff.direction = if(resultP1<=resultP2){"+"}else{"-"}
  }
  RES = if(is.onetail){list(p=min(1,resultP))}else{list(p=min(1,resultP), overall.eff.direction=overall.eff.direction)}
  return(RES)
}
