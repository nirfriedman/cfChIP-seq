
EstimateScale = function(X, c, b, m, v, w = NULL)  {
  if( is.null(w) )
    w = rep(1, length(X))
  
  z = (b/c)
  z2 = z**2
  
  ll = function(x) {
    m = x*m
    v = (x**2) * v
    mu = c*m+b 
    r = m**2/v+2*z*m/v+z2/v 
    -sum(w*dnbinom(X,mu=mu, size=r, log = TRUE ))
  }
  
  Y.guess = (X-b)/c
  fit = lm(Y.guess ~ m - 1)
 
#  print(Y.guess) 
#  print(fit)
  opt = optimize(ll,interval = c(1E-10, 1E10))
#  print(opt)
  c(scale = opt$minimum, ll = -opt$objective)
}

RobustEstimateScale = function(X, c, b, m, v, w = NULL) {
  if( is.null(w) )
    w = rep(1, length(X))
  xx = sapply(1:100, function(i) {
    I = sample(length(X), floor(0.75*length(X)))
    EstimateScale(X[I], c[I], b[I], m[I], v[I], w[I])[1]
  })
#  print(sort(xx))
  median(xx)
}


#
# quad fails in examples such as ExpectedPosteriorGammaPoisson(10,.05, 1, 100, 10000, quad = TRUE, K = 100)
#
ExpectedPosteriorGammaPoisson = function(Y, c, b, m, v, K = 100, quad = FALSE) {
  alpha = m**2/v
  beta = m/v
  alpha1 = alpha+Y
  beta1 = beta+c

  if( qgamma(1-1/(2*K), shape = alpha1, rate = beta1) == 0 ) {
    xs = c(0,1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1)
  } else {
    if(!quad) {
      qs = c( 
        seq(0, 0.1/K, 0.01/K)+0.005/K,
        seq(0, 1/K, 0.1/K)+0.05/K,
        seq(1/K, 0.5, 1/K)+0.5/K)
      qs = c(qs, 1-qs)
      qs = qs[qs >= 0 & qs <= 1]
      qs = unique(sort(qs))
      xs = qgamma(qs, shape = alpha1, rate = beta1)
 
    }
    if(quad) {
      gq = gauss.quad(K, "laguerre", alpha1/beta1)
      xs = gq$nodes/beta1
    }
  }
  # to avoid overshoot
  while( min(xs) > 1 )
    xs = c(min(xs)/sqrt(2), xs, max(xs)*sqrt(2))
  
  xs.mid = (xs[-1]+xs[-length(xs)])/2
  ps = pgamma(xs.mid, shape = alpha1, rate = beta1, log = TRUE)
  v1 = c(ps,0)
  v0 = c(-Inf, ps)
#  ys.old = exp(v1) - exp(v0)
  ys.log = v1+log(1 - exp(v0-v1))
  if( Y == 0 ){
    ws = exp(ys.log)
  } else {
    ws = log(1+b/(c*xs))*Y+ys.log
    ws = exp(ws-max(ws))
  }
  ws = ws/sum(ws)
  c(sum(xs*ws), sum((xs**2) * ws))
}

EstimateMeanVar = function(X, c, b, w = NULL, n0 = .1, p0 = 1, prior.w = 1) {

#  p0 = 1  # product
  q0 = n0  # sum
  w0 = 0.1
  if( is.null(w) )
    w = rep(1, length(X))
  
  z = (b/c)
  z2 = z**2
  
  ll = function(x) {
    m = x[1]
    v = x[2]
    alpha = m**2/v
    beta = m/v
    mu = c*m+b  # c*alpha/beta+b
    r = m**2/v+2*z*m/v+z2/v #alpha+2*b*beta/c+(b*beta/c)**2/alpha
    -sum(w*dnbinom(X,mu=mu, size=r, log = TRUE )) -
       prior.w*((alpha-1)*log(p0)-beta*q0-n0*lgamma(alpha)+alpha*n0*log(beta)) #+ w0*(m + 1/m)
      #0
  }
  
  test.prior = function(m,v) {
    alpha = m**2/v
    beta = m/v
    -((alpha-1)*log(p0)-beta*q0-n0*lgamma(alpha)+alpha*n0*log(beta))
  }
  
  test.prior2 = function(m) {
    opt = optimize(function(v) test.prior(m,v), c(1E-4,1E10))
    opt$objective
  }

  ll.m = function(m) {
    opt = optimize(function(v) ll(c(m,v)), c(1E-4,1E10))
    opt$objective
  }
 
  c.mean = mean(c)
  b.mean = mean(b)
  Y.guess = pmax(0,(X-b)/c)
  mu.emp = max(mean(Y.guess),0.01)
  var.emp = max(var(Y.guess), 0.01)
  opt = optim(c(mu.emp,var.emp), ll)
  
  if( opt$par[2] < 1E-3 || opt$par[1] < 1E-3 ) {
    opt = optimize(ll.m,  c(1E-5,1E10))
    m = opt$minimum
    opt = optimize(function(v) ll(c(m,v)), c(1E-4,1E10))
    v = opt$minimum
  } else {
    m = opt$par[1]
    v = opt$par[2]
  }

  if(0) {
    Y = sapply(1:length(X), function(k) ExpectedPosteriorGammaPoisson(X[k], c[k], b[k], m, v))
    meanP = mean(Y[1,])
    varP = mean(Y[2,]) - meanP**2
    #need to be fixed
    if( is.na(meanP) ) {
      meanP = m
      varP = v
    }
  } else {
    meanP = m
    varP = v
  }
    
  c(mean = m, var = v, ll = -ll(c(m,v)), meanP = meanP, varP = varP)
}

eps = 1E-8

CondExpX1 = function(Y, lambda0, mu2, var2) {
  x2 = 0:Y
  a2 = (mu2**2)/var2
  b2 = mu2/var2

  px2 = dnbinom(x2,mu = mu2, size = (mu2**2)/(var2), log = TRUE)
  px0 = dpois(x2,lambda0, log = TRUE)
  px = px2 + rev(px0)

  pp = max(px)
  px = exp(px - pp)
  ll = pp+log(sum(px))
  
  px = px / sum(px)
  EX = sum(x2 * px)
  c(EX = EX,
    ElogX = sum(log(x2+eps)*px),
    EY = (a2+EX)/(b2+1),
    ElogY = sum(digamma(a2+x2)*px) - log(b2+1),
    ll = ll)
}

CondExpX2 = function(Y, lambda0, mu1, var1, mu2, var2) {
  x2 = 0:Y
  a2 = mu2**2/var2
  b2 = mu2/var2
  
  px0 = dpois(x2,lambda0, log = TRUE)
  px1 = dnbinom(x2,mu = mu1, size = (mu1**2)/(var1), log = TRUE)
  px2 = dnbinom(x2,mu = mu2, size = (mu2**2)/(var2), log = TRUE)
  px01 = sapply(x2, function(k) { z = px1[1:(k+1)] + rev(px0[1:(k+1)]); zm = max(z); zm + log(sum(exp(z - zm)))})
  
  px = px2 + rev(px01)
  pp = max(px)
  px = exp(px - pp)
  ll = pp+log(sum(px))
  px = px / sum(px)
  
  EX = sum(x2 * px)
  c(EX = EX,
    ElogX = sum(log(x2+eps)*px),
    EY = (a2+EX)/(b2+1),
    ElogY = sum(digamma(a2+x2)*px) - log(b2+1),
    ll = ll)
}

CondExpXN = function(Y, lambda0, mu.l, var.l, mu, var) {
  x.range = 0:Y
  a = mu**2/var
  b = mu/var
  
  p.l = c(list(dpois(x.range,lambda0, log = TRUE)),
          mapply(function(m,v) dnbinom(x.range,mu = m, size = (m**2)/v, log = TRUE), mu.l, var.l, SIMPLIFY = FALSE))
  px = dnbinom(x.range,mu = mu, size = (mu**2)/(var), log = TRUE)
  
  
  Mp = rep(1:(Y+1), (Y+1):1)
  Mq = unlist(lapply((Y+1):1, function(k) 1:k))
  SubSum = function(p,q) {
    M = matrix(-Inf, nr = Y+1, nc = Y+1)
    M[lower.tri(M, diag = TRUE)] = p[Mp] + q[Mq]
    zm = rowMax(M)
    zm + log(rowSums(exp(M - zm)))
  }
  q = Reduce(SubSum, p.l)

  px = px + rev(q)
  pp = max(px)
  px = exp(px - pp)
  ll = pp+log(sum(px))
  px = px / sum(px)
  
  EX = sum(x.range * px)
  c(EX = EX,
    ElogX = sum(log(x.range+eps)*px),
    EY = (a+EX)/(b+1),
    ElogY = sum(digamma(a+x.range)*px) - log(b+1),
    ll = ll)
}

# inefficient version
CondExpXMulti = function(Y, lambda0, mu.l, var.l) {
  rows = c("EX", "EY", "ElogY")
  
  df = sapply(1:length(mu.l), 
         function(i)
           CondExpXN(Y, lambda0, mu.l[-i], var.l[-i], mu.l[i], var.l[i])[rows])
  df 
  E["EX", 1] = Y - sum(E["EX", 2:ncol(E)])
  E
}
  
if(0) {
  #test code for CondExpX1/2
  
  PlotErr = function(Os, Xs, Ms, Ms.err) {
    df = data.frame(x = Os,
                    y = (Xs/Ms), 
                    ymax = (Xs/(Ms+Ms.err)),
                    ymin = (Xs/(Ms-Ms.err)))
    p = ggplot(df, aes(x = x, y = y))
    p = p + geom_ribbon(aes(ymin = ymin, ymax=ymax), alpha = 0.25, fill = "blue", color = "blue")
    p = p + geom_point()
    p = p+ scale_y_log10(limits=c(0.95,1.05))
    p = p+ geom_hline(yintercept = 1, color = "red")
    p
  }
  
  lambda0 = 10
  mu2 = 5
  var2 = 50
  N = 1E7
  Y = rgamma(N, shape = (mu2**2)/var2, rate = mu2/var2)
  X = rpois(N, Y)
  B = rpois(N, lambda0)
  O = X+B
  ns = sapply(split(X,O), length) 
  Os = as.numeric(names(ns))
  Ms = sapply(split(X,O),mean)
  Ms.err = 1.96*sapply(split(X,O),sd)/sqrt(ns)
  Ys = sapply(split(Y,O), mean)
  Ys.err = 1.96*sapply(split(Y,O), sd)/sqrt(ns)
  LogYs = sapply(split(log(Y),O), mean)
  LogYs.err = 1.96*sapply(split(log(Y),O),sd)/sqrt(ns)
  
  Xs = sapply(Os, function(o) CondExpX1(o, lambda0, mu2, var2 ))
  qplot(Ms, Xs["EX",])
  PlotErr(Os, Xs["EX",], Ms, Ms.err) + labs(title ="X1 E[X]")
  PlotErr(Os, Xs["EY",], Ys, Ys.err) + labs(title = "X1 E[Y]")
  PlotErr(Os, Xs["ElogY",], LogYs, LogYs.err) + labs(title = "X1 E[logY]")
  
  mu1 = 10
  var1 = 5
  Y1 = rgamma(N, shape = (mu1**2)/var1, rate = mu1/var1)
  X1 = rpois(N, Y1)
  O = X+B+X1
  ns = sapply(split(X,O), length) 
  Os = as.numeric(names(ns))
  Ms = sapply(split(X,O),mean)
  Ms.err = 1.96*sapply(split(X,O),sd)/sqrt(ns)
  Ys = sapply(split(Y,O), mean)
  Ys.err = 1.96*sapply(split(Y,O), sd)/sqrt(ns)
  LogYs = sapply(split(log(Y),O), mean)
  LogYs.err = 1.96*sapply(split(log(Y),O),sd)/sqrt(ns)
#  Xs = sapply(Os, function(o) CondExpX2(o, lambda0, mu1, var1, mu2, var2 ))
  Xs = sapply(Os, function(o) CondExpXN(o, lambda0, mu1, var1, mu2, var2 ))
  PlotErr(Os, Xs["EX",], Ms, Ms.err) + labs(title ="X2 E[X]")
  PlotErr(Os, Xs["EY",], Ys, Ys.err) + labs(title = "X2 E[Y]")
  PlotErr(Os, Xs["ElogY",], LogYs, LogYs.err) + labs(title = "X2 E[logY]")
}

EstimateGammaFromStat = function(N, SY, SlogY) {
  #  print(c(N,SY, SlogY))
  z = log(SY/N) - SlogY/N
  ur = uniroot(function(a) log(a) - digamma(a) - z, c(1E-6,1E6))
  a = ur$root
  b = N*a / SY
  c(a = a, b = b)
}

if(0) {
  #code to evaluate EstimateGamma
  m = 10
  v = 1
  a = (m**2)/v
  b = m/v
  N = 1000000
  X = rgamma(N, shape = a, rate = b)
  ll = EstimateGammaFromStat(N, sum(X), sum(log(X)))
  c(m = ll["a"]/ll["b"], v = ll["a"]/(ll["b"]**2))
}

EstimateMeanVarEM = function(O, s2, BG, w = NULL, 
                             mu1 = NULL, var1 = NULL, s1 = NULL,
                             n0 = .1, p0 = 1, prior.w = 1) {

  
  Pseudo.N = n0*p0
  Pseudo.SY = n0*p0*prior.w
  Pseudo.SlogY = n0*prior.w * log(p0)
  
  if( is.null(w) )
    w = rep(1, length(X))
  
  
  EstimateY2 = function(mu2, var2) {
    MU2 = mu2 * s2
    VAR2 = var2 * (s2**2)
    if( is.null(mu1)) {
      Z = sapply(1:length(O), 
                 function(i) CondExpX1(O[i], BG[i], MU2[i], VAR2[i]))
    } else {
      MU1 = mu1*s1
      VAR1 = var1 * (s1**2)
      
      Z = sapply(1:length(O), 
                 function(i) CondExpX2(O[i], BG[i], MU1[i], VAR1[i], MU2[i], VAR2[i]))
    }
    N = length(O)
    SlogY = sum(Z["ElogY",]) #- sum(log(s2))
    SY = sum(Z["EY",]/s2)
    cc = EstimateGammaFromStat(N+Pseudo.N, SY+Pseudo.SY, SlogY+Pseudo.SlogY)
    c( cc, ll = sum(Z["ll",]))
  }

  EMIterate = function (param, Iter) {
    Stop = FALSE
    i = 1
    logl= -Inf
    m = param["m"]
    v = param["v"]
#    print(c(m,v))
    while( !Stop && i < Iter  ) {
      ll =  EstimateY2(m, v)
      if( ll["ll"] < logl+eps )
        Stop = TRUE
      m = ll["a"]/ll["b"]
      v = m/ll["b"]
      logl = ll["ll"]
      i = i+1
      #    print(c(i,ll["ll"],mu = mu2, var = var2))
    }
    c(m = as.numeric(m), v = as.numeric(v), ll = as.numeric(logl))
  }
  
  EMInit = function(c = 1) {
    Diff = O-BG
    if( !is.null(mu1) )
      Diff = Diff - mu1*s1
    Diff[Diff < 0] = 0
    Diff = Diff / s2
    I = sample(length(O), c*length(O), replace = TRUE)
    mu2 = pmax(mean(Diff[I]), 1E-5)*(1+rexp(1,2))
    var2 = pmax(var(Diff[I]), 1E-5)*(1+rexp(1,2))
    c(m = mu2, v = var2)
  }
  
  logl = -Inf
  
  RS = replicate(10,EMIterate(EMInit(), 5), simplify = FALSE)
  lls = sapply(RS, function(r) r["ll"])
  t = quantile(lls, 0.75)
  RS = lapply(RS[lls>t], function(x) EMIterate(x,10))
  lls = sapply(RS, function(r) r["ll"])
  i = which.max(lls)  
  p = EMIterate(RS[[i]],20)
  
  mu2 = p["m"]
  var2 = p["v"]
  logl = p["ll"]
  
  if(0) {
    Y = sapply(1:length(O), function(k) ExpectedPosteriorGammaPoisson(O[k], s2[k], BG[k], mu2, var2))
    meanP = mean(Y[1,])
    varP = mean(Y[2,]) - meanP**2
    
    #need to be fixed
    if( is.na(meanP) ) {
      meanP = as.numeric(mu2)
      varP = as.numeric(var2)
    }
  } else {
    meanP = as.numeric(mu2)
    varP = as.numeric(var2)
  }
  
  c(mean = as.numeric(mu2), var = as.numeric(var2), ll = logl, meanP = meanP, varP = varP)
}

# test code
if(0) {
  catn = function(...) cat(...,"\n")
  library(ggplot2)
  library(reshape2)
  library(cowplot)
  N = 1000
  b = rep(0,N)  
  b = rexp(N,1)
  c = rep(1,N)
  c = rgamma(N, 1, 1)
  w = rep(1,N)
  z = (b/c)
  z2 = z**2

  mu = 1
  var = 5
  
  X = rep(0,N)
  if(0) {
    X = rep(1,N)
  }
  if(1) {
    mu.x = mu*c
    var.x = var*c*c
    X = rpois(N, lambda = rgamma(N, shape = (mu.x**2)/var.x, rate = mu.x/var.x )) + rpois(N, lambda = b)
  }
  
  EstimateMeanVarEM(X,c,b)
  EstimateMeanVar(X,c,b)
  
  mu.test = c(eps, 0.1,0.5, 1, 1, 1, 10,40)
  var.test = c(.1, 1,1,1, 5,20, 10, 200)
  M = 100

  prior.test = data.frame(
    n0 = c( # 1, .1, 1, 1,
           2),
    p0 = c(# 1,  1, 1, 0.5,
           1),
    pw = c(#0, 1, 1, 1, 
           1),
    name =  c(#"p0_0", "p0_1", "p1_0","p1_1",
              "p2_0" ))
  
  SimpleTest = function(mu,var, lambda0 = 1, alphac = 0.5, betac= 0.5 ) {
    Est = sapply(1:M, function(i) {        
      c = rgamma(N, alphac, betac)
      b = rpois(N,lambda0)
      mu.x = mu*c
      var.x = var*c*c
      X = rpois(N, lambda = rgamma(N, shape = (mu.x**2)/var.x, rate = mu.x/var.x )) + rpois(N, lambda = b)
      w = EstimateMeanVarEM(X,c,b)
    })
    c(mu = quantile(Est[1,]), var = quantile(Est[2,]))
  }
  
  if(0) {
  catn("mu =",1,"var = ",5, "lam0 = ",0, "a/b c = ",100,100)
  SimpleTest(1,5, 0, 100,100)
  catn("mu =",1,"var = ",5, "lam0 = ",1, "a/b c = ",100,100)
  SimpleTest(1,5,1, 100,100)
  catn("mu =",1,"var = ",5, "lam0 = ",10, "a/b c = ",100,100)
  SimpleTest(1,5,10, 100,100)
  catn("mu =",1,"var = ",5, "lam0 = ",1, "a/b c = ",10,10)
  SimpleTest(1,5,10, 10,10)
  catn("mu =",1,"var = ",5, "lam0 = ",1, "a/b c = ",2,2)
  SimpleTest(1,5,10, 2,2)
  
  catn("mu =",1,"var = ",5, "lam0 = ",1, "a/b c = ",1,1)
  SimpleTest(1,5,10, 1,1)
  }
  
  for( i in 1:length(mu.test) ) {
    mu = mu.test[i]
    var = var.test[i]
    est = sapply(1:M, function(i) {
#      c = rgamma(N, .5, .5)
      c = 2**runif(N, min = -1, max= 1)
#      b = rgamma(N, 5,1)
      b = rep(1, N)
      mu.x = mu*c
      var.x = var*c*c
      X = rpois(N, lambda = rgamma(N, shape = (mu.x**2)/var.x, rate = mu.x/var.x )) + rpois(N, lambda = b)
      unlist(
        c( lapply(1:nrow(prior.test),
                  function(j) {
                    w = EstimateMeanVarEM(X,c,b, 
                                          n0 = prior.test[j,"n0"], 
                                          p0 = prior.test[j,"p0"], 
                                          prior.w = prior.test[j,"pw"])
                    names(w) = paste0("EM-",prior.test[j, "name"], ".", names(w))
                    w
                  }),
           lapply(1:nrow(prior.test),
                  function(j) {
                    w = EstimateMeanVar(X,c,b, 
                                        n0 = prior.test[j,"n0"], 
                                        p0 = prior.test[j,"p0"], 
                                        prior.w = prior.test[j,"pw"])
                    names(w) = paste0(prior.test[j, "name"], ".", names(w))
                    w
                 })))
    })
    
    df = melt(est, as.is = TRUE)
    df$test = sapply(strsplit(df$Var1, "\\."), function(x) x[1])
    df$type = sapply(strsplit(df$Var1, "\\."), function(x) x[2])
    p = plot_grid(
      ggplot(df[grepl("mean",df$type),], aes(x = test, y = value ))+
        geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = mu, color = "blue")+scale_y_log10() + facet_grid( ~ type),
      ggplot(df[grepl("var", df$type),], aes(x = test, y = value ))+
        geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = var, color = "blue")+scale_y_log10()+facet_grid(.~type),
      labels = c(paste("mean(", mu, ")"), paste("var(",var,")")),
      ncol = 1
    )
    print(p)
    ggsave(p, filename = paste0("~/Temp/test-", mu, "-", var, ".png"))
  }
}


LRTMeanVar = function(X,c,b,Pos, Neg, Weights = NULL) {
  S = c(Pos,Neg)
  if(is.null(Weights)) {
    Weights = rep(1, length(S))
    names(Weights) = S
  }
  l12 = EstimateMeanVar(X[S],c[S],b[S], Weights[S])
  l1 = EstimateMeanVar(X[Pos], c[Pos], b[Pos], Weights[Pos])
  l2 = EstimateMeanVar(X[Neg], c[Neg], b[Neg], Weights[Neg])
  LR = -2*(l12["ll"] - l1["ll"] - l2["ll"])
  pVal = pchisq(LR, 2,  lower.tail = FALSE)
  list(pVal = as.numeric(pVal), LR = LR, l.both = l12, l.Pos = l1, l.Neg = l2)
}


if(0) {
  
  library(ggplot2)
  library(cowplot)
  N = 10000
  Range = 0:80

  
  m = 20
  c = 1
  p.list = do.call("c",
                   lapply(c(.5,1, 4), function(c) 
                     do.call("c",
                             lapply(c(10,30,60), function(v) 
                               lapply(c(0,10,15), function(b) {
                                 m1 = m - b
                                 alpha = m1**2/v
                                 beta = m1/v
                                 mu = c*alpha/beta+b
                                 r = alpha+2*b*beta/c+(b*beta/c)**2/alpha
                                 
                                 Y = rpois(N, rgamma(N,shape = alpha, rate = beta)*c+b)
                                 Z = dnbinom(Range,mu=mu, size=r,)
                                 
                                 p = ggplot(data = data.frame(Y = Y), aes(Y)) + theme_classic()
                                 p = p + geom_histogram(binwidth = 1)
                                 p = p + geom_step(data = data.frame(x = Range-0.5, y = Z*N), aes(x,y), color = "red")
                                 p = p + labs(subtitle = paste("m =", m, "var =", v, "c =", c, "bg =",b))
                                 p
                               })
                             ))))  
  
  p = plot_grid(plotlist = p.list, nc = 3)
  p
  
  
  alpha = 100
  beta = 1
  n = 100000
  c = rgamma(n, shape = 5, rate = .5)/10
  b = runif(n, 0,100)
#  b = rep(0,n)
  
  Y = rgamma(n, shape = alpha, rate = beta)
  X = rpois(n, c*Y+b)
  X.max = max(X)
  df = data.frame(x = X, y = Y)
  p = ggplot(data=df, aes(x))+theme_classic()
  p = p + geom_histogram(binwidth = 1)
  p
  
  df = data.frame(x = (X-b)/c)
  p = ggplot(data=df, aes(x))+theme_classic()
  p = p + geom_histogram(binwidth = 1)
  cat("mean(X) =", mean(X), "\n")
  mu = (X-b)/c
  cat("mu =", mean(mu), "\n")
  cat("var(X) =", var(X), "\n")
  r = alpha+2*b*beta/c+(b*(beta*c))**2*alpha
  cat("sigma^2 =", mean(mu+(mu**2)*1/r), "\n")
  
  p = p + geom_step(data=data.frame(x = -0.5+0:X.max, y = dnbinom(0:X.max, size  = mean(r), mu=mean(mu))*n), aes(x,y), color = "blue")
  p
  
  
  Ns = exp(seq(log(100),log(n), 0.1))
  Est = lapply(Ns, function(m) EstimateMeanVar(X[1:m], c[1:m], b[1:m]))
  Est.means = sapply(Est, function(x) x["mean"])
  Est.vars =     sapply(Est, function(x) x["var"])      
  p = ggplot(data.frame(Ns = Ns, mean = Est.means, var = Est.vars), aes(x = Ns)) 
  p = p+geom_line(aes(y = mean), color = "blue")
  p = p+geom_point(aes(y = mean), color = "blue")
  p = p+geom_line(aes(y = var), color = "red")
  p = p+geom_point(aes(y = var), color = "red")
  p = p + ylim(75,125)
  p
}
