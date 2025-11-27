#UNIFORM DISTRIBUTION
#f(x) = P(X=k) = 1/k
UPx<- function(k){
  return (1/k)
}

Umean<- function(x){
  return(sum(x)*(1/length(x)))
}

Uvar<- function(x){
  return((sum(x-Umean(x))**2)*1/length(x))
}

#BERNOULLI DISTRIBUTION

#BINOMIAL DISTRIBUTION
Bfx<- function(n,p,x){
  return(choose(n,x)*(p**x)*((1-p)**(n-x)))
}

Bmean<-function(n,p){
  return (n*p)
}

Bvar<- function(n,p){
  return(n*p*(1-p))
}

#GEOMETRIC DISTRIBUTION

#POISSON DISTRIBUTION

PA<- function(A,x){
  return(exp(-A)*((A**x)/factorial(x)))
}

Amean<- function(A){
  return (A)
}

Avar<- function(A){
  return (A)
}

#Approximations of the Binomial Distribution
#It is generally accepted that:
#  If n > 30 and np < 5 then B(n, p) ≈ P(np)
#If n > 30 and nq < 5 then B(n, q) ≈ P(nq)
#If n > 30, np ≥ 5 and nq ≥ 5 then B(n, p) ≈ N(np, √npq)