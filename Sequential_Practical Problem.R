

# Question 2 --------------------------------------------------------------

# Given parameters
# Under Null Hypothesis
theta_0 = 7
# Under Alternative Hypothesis
theta_1 = 10

# Given Strength
alpha = 0.2
beta = 0.2

# Approximate values of boundary A and B
A = (1-beta)/alpha
B = beta/(1-alpha)
# Logarithm of the boundary values
a = log(A)
b = log(B)

# Number of observations
n=1:12

# Function for plotting the lower limit
low = function(n)
{
  return((b-n*(theta_0-theta_1))/(log(theta_1/theta_0)))
}

# Function for plotting the upper limit
upper = function(n)
{
  return((a-n*(theta_0-theta_1))/(log(theta_1/theta_0)))
}

# Sequence of observations - To determine the continuation/acceptance/rejection
# of the 12 th observation
X = c(7,15,11,3,4,5,7,9,12,14,0,11)
S_n = cumsum(X)
limit_mat = matrix(c(low(n),S_n,upper(n)),ncol=3,nrow=12,byrow=F)

# Function for plotting continuation region of Wald's SPRT
n = 1:100
plot_cont_reg = function(l,u){
  plot(n,low(n),type = "l",xlab = "n",ylab = bquote(S[n]),xlim=c(6,20),ylim=c(50,100),
       main=bquote("Continuation Region for WALD's SPRT-"~S[n]^c),
       lwd=2,lty=2)
  lines(n,upper(n),type="l",lty=2,lwd=2)
  polygon(c(n, rev(n)), c(upper(n), rev(low(n))),
          col = "#EC7063")
  #points(12,98)
  #text(13,98,"(12,98)")
  text(12,65,"Acceptance Region")
  text(8,90,"Rejection Region")
}
# Continuation region graph
plot_cont_reg(low,upper)

# Sequence for different values of t0
t_0 = seq(-2,2,0.01)#c(-.5,-.3,-.2,0,.2,.5,.7)

# Function for generating OC Function
OC = function(t_0)
{
  oc = array(0)
  for(i in 1:length(t_0))
  {
    oc[i] = ((A^t_0[i])-1)/(((A^t_0[i])-(B^t_0[i])))
    if(t_0[i] == 0)
    {
      oc[i] = log(A)/(log(A/B))
    }
  }
  return(oc)
}

# OC Function values at different values of t0
oc = OC(t_0)

#Function for generating theta
theta = function(t_0)
{
  th = array(0)
  for(i in 1:length(t_0))
  {
    th[i] = t_0[i]*(theta_1-theta_0)/(((theta_1/theta_0)^t_0[i])-1)
    if(t_0[i] == 0)
    {
      th[i] = (theta_1-theta_0)/log(theta_1/theta_0)
    }
  }
  return(th)
}

# Theta values for different values of t0
th = theta(t_0)

# Expectation of Z1
E_Z1 = (theta_0-theta_1)+th*log(theta_1/theta_0)

# Variance of Z1
v_Z1 = th*(log(theta_1/theta_0))^2

#Function generating ASN
ASN = function(oc,e,v){
  asn = c()
  for(i in 1:length(t_0)){
    asn[i] = (b*oc[i]+a*(1-oc[i]))/E_Z1[i]
    if(t_0[i] == 0){
      asn[i] = ((b^2)*oc[i]+(a^2)*(1-oc[i]))/v_Z1[i]
    }
  }
  return(asn)
}

# ASN Function values for different values of OC
asn = ASN(oc,E_Z1,v_Z1)


## Wald's Approximate OC and ASN Curves

par(mfrow=c(1,2))

plot(th,oc,main=bquote("Wald's Approximate OC Curve"~L(theta)),
     ylab=bquote(L(theta)),xlab=bquote(theta))
plot(th,asn,main=bquote("Wald's Approximate ASN Curve"~E[theta](N)),
     ylab=bquote(E[theta](N)),xlab=bquote(theta))
abline(h=1.653)
abline(v=c(7,10),col="blue",lty=2)


# Comparison of Wald's SPRT with Fixed Sample Size Test

mp_size=function(alpha,beta){
  n = array(0)
  k = array(0)
  for(i in 1:40)
  {
    for(j in 1:100)
    {
      null_hyp = ppois(j,i*7,lower.tail = FALSE)
      alt_hyp = ppois(j,i*10,lower.tail = FALSE)
      if(null_hyp<=alpha && alt_hyp>(1-beta))
      {
        n[i] = i
        k[j] = j
      }
    }
  }
  n=n[!is.na(n)]
  k=k[!is.na(k)] # We get the value of quantile k-1
  h0 = ppois(k[1],n[1]*7,lower.tail = FALSE)
  h1 = ppois(k[1],n[1]*10,lower.tail = FALSE)
  z = c(n[1],k[1]+1,h0,h1) # Adding 1 to k-1, to get the quantile k
  names(z) = c("n","k","Size","Power")
  return(z)
}
mp_size(0.2,0.2)


# Question 3 --------------------------------------------------------------

# Given parameters
# Under Null Hypothesis
theta_0 = 0
# Under Alternative Hypothesis
theta_1 = 1

# Given Strength
alpha = 0.01
beta = 0.01

# Approximate values of boundary A and B
A = (1-beta)/alpha
B = beta/(1-alpha)
# Logarithm of the boundary values
a = log(A)
b = log(B)

# Number of observations
n=1:12

# Function for plotting the lower limit
low = function(n)
{
  return(0.50 - (4.595/n))
}

# Function for plotting the upper limit
upper = function(n)
{
  return(0.50 + (4.595/n))
}

# Sequence of observations - To determine the continuation/acceptance/rejection
# of the 12 th observation
#X = c(7,15,11,3,4,5,7,9,12,14,0,11)
#S_n = cumsum(X)
#limit_mat = matrix(c(low(n),S_n,upper(n)),ncol=3,nrow=12,byrow=F)

# Function for plotting continuation region of Wald's SPRT
n = 1:100
plot_cont_reg = function(l,u){
  plot(n,low(n),type = "l",xlab = "n",ylab = bquote(S[n]),
       main=bquote("Continuation Region for WALD's SPRT-"~S[n]^c),
       lwd=2,lty=2,ylim = c(-5,5))
  lines(n,upper(n),type="l",lty=2,lwd=2)
  polygon(c(n, rev(n)), c(upper(n), rev(low(n))),
          col = "#EC7063")
  #points(12,98)
  text(60,-2,"Acceptance Region")
  text(60,2,"Rejection Region")
  #text(8,90,"Rejection Region")
}
# Continuation region graph
plot_cont_reg(low,upper)

# Sequence for different values of t0
t_0 = seq(-2,2,0.01)#c(-.5,-.3,-.2,0,.2,.5,.7)

# Function for generating OC Function
OC = function(t_0)
{
  oc = array(0)
  for(i in 1:length(t_0))
  {
    oc[i] = ((A^t_0[i])-1)/(((A^t_0[i])-(B^t_0[i])))
    if(t_0[i] == 0)
    {
      oc[i] = log(A)/(log(A/B))
    }
  }
  return(oc)
}

# OC Function values at different values of t0
oc = OC(t_0)

#Function for generating theta
theta = function(t_0)
{
  th = array(0)
  for(i in 1:length(t_0))
  {
    #th[i] = t_0[i]*(theta_1-theta_0)/(((theta_1/theta_0)^t_0[i])-1)
    #if(t_0[i] == 0)
    #{
      th[i] = ((theta_1+theta_0) - ((theta_1-theta_0)*t_0[i]))/2
    #}
  }
  return(th)
}

# Theta values for different values of t0
th = theta(t_0)

# Expectation of Z1
E_Z1 = (theta_1-theta_0)*(th - 0.5*(theta_0+theta_1))

# Variance of Z1
v_Z1 = (theta_1-theta_0)^2

#Function generating ASN
ASN = function(oc,e,v)
{
  asn = c()
  for(i in 1:length(t_0))
  {
    if(e[i] != 0)
      asn[i] = (b*oc[i]+a*(1-oc[i]))/e[i]
    else
      asn[i] = ((b^2)*oc[i]+(a^2)*(1-oc[i]))/v[i]
  }
  return(asn)
}

# ASN Function values for different values of OC
asn = ASN(oc,E_Z1,v_Z1)


## Wald's Approximate OC and ASN Curves

par(mfrow=c(1,2))

plot(th,oc,main=bquote("Wald's Approximate OC Curve"~L(theta)),
     ylab=bquote(L(theta)),xlab=bquote(theta))
plot(th,asn,main=bquote("Wald's Approximate ASN Curve"~E[theta](N)),
     ylab=bquote(E[theta](N)),xlab=bquote(theta))
#abline(h=1.653)
#abline(v=c(7,10),col="blue",lty=2)


# Comparison of Wald's SPRT with Fixed Sample Size Test

#mp_size=function(alpha,beta){
#  n = array(0)
#  k = array(0)
#  for(i in 1:40)
#  {
#    for(j in 1:100)
#    {
#      null_hyp = ppois(j,i*7,lower.tail = FALSE)
#      alt_hyp = ppois(j,i*10,lower.tail = FALSE)
#      if(null_hyp<=alpha && alt_hyp>(1-beta))
#      {
#        n[i] = i
#        k[j] = j
#      }
#    }
#  }
#  n=n[!is.na(n)]
#  k=k[!is.na(k)] # We get the value of quantile k-1
#  h0 = ppois(k[1],n[1]*7,lower.tail = FALSE)
#  h1 = ppois(k[1],n[1]*10,lower.tail = FALSE)
#  z = c(n[1],k[1]+1,h0,h1) # Adding 1 to k-1, to get the quantile k
#  names(z) = c("n","k","Size","Power")
#  return(z)
#}
#mp_size(0.2,0.2)


# Question 4 --------------------------------------------------------------
rm(list=ls())
# Given parameters
# Under Null Hypothesis
theta_0 = 4
# Under Alternative Hypothesis
theta_1 = 5

# Given Strength
alpha = 0.01
beta = 0.01

# Approximate values of boundary A and B
A = (1-beta)/alpha
B = beta/(1-alpha)
# Logarithm of the boundary values
a = log(A)
b = log(B)

# Number of observations
n=1:12

# Function for plotting the lower limit
low = function(n)
{
  return((-(4.595/n) - log(0.8))/(1/20))
}

# Function for plotting the upper limit
upper = function(n)
{
  return(((4.595/n) - log(0.8))/(1/20))
}

# Sequence of observations - To determine the continuation/acceptance/rejection
# of the 12 th observation
#X = c(7,15,11,3,4,5,7,9,12,14,0,11)
#S_n = cumsum(X)
#limit_mat = matrix(c(low(n),S_n,upper(n)),ncol=3,nrow=12,byrow=F)

# Function for plotting continuation region of Wald's SPRT
n = 1:100
plot_cont_reg = function(l,u){
  plot(n,low(n),type = "l",xlab = "n",ylab = bquote(S[n]),
       main=bquote("Continuation Region for WALD's SPRT-"~S[n]^c),
       lwd=2,lty=2,ylim=c(-50,50))
  lines(n,upper(n),type="l",lty=2,lwd=2)
  polygon(c(n, rev(n)), c(upper(n), rev(low(n))),
          col = "#EC7063")
  #points(12,98)
  text(20,-40,"Acceptance Region")
  text(20,40,"Rejection Region")
  #text(8,90,"Rejection Region")
}
# Continuation region graph
plot_cont_reg(low,upper)

# Sequence for different values of t0
t_0 = seq(0.01,2,0.01)

# Function for generating OC Function
OC = function(t_0)
{
  oc = array(0)
  for(i in 1:length(t_0))
  {
    oc[i] = ((A^t_0[i])-1)/(((A^t_0[i])-(B^t_0[i])))
    if(t_0[i] == 0)
    {
      oc[i] = log(A)/(log(A/B))
    }
  }
  return(oc)
}

# OC Function values at different values of t0
oc = OC(t_0)

#Function for generating theta
theta = function(t_0)
{
  th = array(0)
  for(i in 1:length(t_0))
  {
    #th[i] = t_0[i]*(theta_1-theta_0)/(((theta_1/theta_0)^t_0[i])-1)
    #if(t_0[i] == 0)
    #{
    th[i] = 1/((t_0[i]*(1/20))/(1-(0.8^t_0[i])))
    #}
  }
  return(th)
}

# Theta values for different values of t0
th = theta(t_0)

# Expectation of Z1
E_Z1 = log(0.8) + th*(1/20)

# Variance of Z1
v_Z1 = (th*(1/20))^2

#Function generating ASN
ASN = function(oc,e,v)
{
  asn = c()
  for(i in 1:length(t_0))
  {
    if(e[i] != 0)
      asn[i] = (b*oc[i]+a*(1-oc[i]))/e[i]
    else
      asn[i] = ((b^2)*oc[i]+(a^2)*(1-oc[i]))/v[i]
  }
  return(asn)
}

# ASN Function values for different values of OC
asn = ASN(oc,E_Z1,v_Z1)


## Wald's Approximate OC and ASN Curves

par(mfrow=c(1,2))

plot(th,oc,main=bquote("Wald's Approximate OC Curve"~L(theta)),
     ylab=bquote(L(theta)),xlab=bquote(theta))
plot(th,asn,main=bquote("Wald's Approximate ASN Curve"~E[theta](N)),
     ylab=bquote(E[theta](N)),xlab=bquote(theta))
#abline(h=1.653)
#abline(v=c(7,10),col="blue",lty=2)


# Comparison of Wald's SPRT with Fixed Sample Size Test

#mp_size=function(alpha,beta){
#  n = array(0)
#  k = array(0)
#  for(i in 1:40)
#  {
#    for(j in 1:100)
#    {
#      null_hyp = ppois(j,i*7,lower.tail = FALSE)
#      alt_hyp = ppois(j,i*10,lower.tail = FALSE)
#      if(null_hyp<=alpha && alt_hyp>(1-beta))
#      {
#        n[i] = i
#        k[j] = j
#      }
#    }
#  }
#  n=n[!is.na(n)]
#  k=k[!is.na(k)] # We get the value of quantile k-1
#  h0 = ppois(k[1],n[1]*7,lower.tail = FALSE)
#  h1 = ppois(k[1],n[1]*10,lower.tail = FALSE)
#  z = c(n[1],k[1]+1,h0,h1) # Adding 1 to k-1, to get the quantile k
#  names(z) = c("n","k","Size","Power")
#  return(z)
#}
#mp_size(0.2,0.2)


# Question 1 --------------------------------------------------------------

# Given parameters
# Under Null Hypothesis
theta_0 = 0.5
# Under Alternative Hypothesis
theta_1 = 0.8

# Given Strength
alpha = 0.2
beta = 0.2

# Approximate values of boundary A and B
A = (1-beta)/alpha
B = beta/(1-alpha)
# Logarithm of the boundary values
a = log(A)
b = log(B)

# Number of observations
n=1:9

# Function for plotting the lower limit
low = function(n)
{
  (-1.386 - n*log(0.4))/log(4)
}

# Function for plotting the upper limit
upper = function(n)
{
  (1.386 - n*log(.4))/log(4)
}

# Sequence of observations - To determine the continuation/acceptance/rejection
# of the 12 th observation
X = c(1,0,0,1,1,1,1,1,1)
S_n = cumsum(X)
limit_mat = matrix(c(low(n),S_n,upper(n)),ncol=3,nrow=9,byrow=F)

# Function for plotting continuation region of Wald's SPRT
n = 1:100
plot_cont_reg = function(l,u){
  plot(n,low(n),type = "l",xlab = "n",ylab = bquote(S[n]),
       main=bquote("Continuation Region for WALD's SPRT-"~S[n]^c),
       lwd=2,lty=2)
  lines(n,upper(n),type="l",lty=2,lwd=2)
  polygon(c(n, rev(n)), c(upper(n), rev(low(n))),
          col = "#EC7063")
  #points(12,98)
  text(60,10,"Acceptance region")
  #text(12,65,"Acceptance Region")
  text(20,40,"Rejection Region")
}
# Continuation region graph
plot_cont_reg(low,upper)

# Sequence for different values of t0
t_0 = seq(-2,2,0.01)#c(-.5,-.3,-.2,0.1,.2,.5,.7)
t_0 = t_0[-201]

# Function for generating OC Function
OC = function(t_0)
{
  oc = array(0)
  for(i in 1:length(t_0))
  {
    oc[i] = ((A^t_0[i])-1)/(((A^t_0[i])-(B^t_0[i])))
    if(t_0[i] == 0)
    {
      oc[i] = log(A)/(log(A/B))
    }
  }
  return(oc)
}

# OC Function values at different values of t0
oc = OC(t_0)

#Function for generating theta
theta = function(t_0)
{
  th = array(0)
  for(i in 1:length(t_0))
  {
    #th[i] = t_0[i]*(theta_1-theta_0)/(((theta_1/theta_0)^t_0[i])-1)
    #if(t_0[i] == 0)
    #{
      th[i] = (1-(0.5/0.2)^t_0[i])/(1-4^(t_0[i]))
    #}
  }
  return(th)
}

# Theta values for different values of t0
th = theta(t_0)

# Expectation of Z1
E_Z1 = th*log(4) + log(0.4)

# Variance of Z1
v_Z1 = th*(1-th)*(log(4))^2

#Function generating ASN
ASN = function(oc,e,v){
  asn = c()
  for(i in 1:length(t_0)){
    asn[i] = (b*oc[i]+a*(1-oc[i]))/E_Z1[i]
    if(t_0[i] == 0){
      asn[i] = ((b^2)*oc[i]+(a^2)*(1-oc[i]))/v_Z1[i]
    }
  }
  return(asn)
}

# ASN Function values for different values of OC
asn = ASN(oc,E_Z1,v_Z1)


## Wald's Approximate OC and ASN Curves

par(mfrow=c(1,2))

plot(th,oc,main=bquote("Wald's Approximate OC Curve"~L(theta)),
     ylab=bquote(L(theta)),xlab=bquote(theta))
plot(th,asn,main=bquote("Wald's Approximate ASN Curve"~E[theta](N)),
     ylab=bquote(E[theta](N)),xlab=bquote(theta))
abline(h=1.653)
abline(v=c(7,10),col="blue",lty=2)
