#DENSITY FUNCTTION
#The probability that something happens has a function f(x) whose integral between -infinite and infinite is 1
#P(a <= X <= b) = integral from a to b of f(x)

#DISTRIBUTION FUNCTION
#The distribution function of a continuous random variable X is a function F (x) that associates each value xi with the probability that the variable takes a value less than or equal to that value.
#F(xi) = P(X <= xi) = integral of f(x) from -infinite to xi

#You can calculate the probabilities as areas

#mean or expected value = integral of xf(x)

#variance = integral of x²f(x) - mean²


#CONTINUOUS UNIFORM DISTRIBUTION
U<-function(a, b){
  return(1/(b-a))
}

Umean<-function(a,b){
  return ((a+b)/2)
}

Uvar<-function(a,b){
  return(((b-a)**2)/12)
}

#Density function is constant with linear growth

#NORMAL DISTRIBUTION
#u = mean (determines where it is centered.)
#o = standart deviation (determines its width)
N<-function(x,u,o){
  return((1/(o*sqrt(2*pi)))*exp(((x-u)**2)/2*(o**2)))
}

#Standard Normal Distribution -> u=0 and o=1
#For example:
#P(Z > 0.52) = 1 − P(Z ≤ 0.52) = 1 − F (0.52) = 1 − 0.6985 = 0.3015

#The theorem for standardization is:
#X ∼ N(μ, σ) ⇒ Z = X − μ /σ ∼ N(0, 1).

#CHI SQUARE DISTRIBUTION
