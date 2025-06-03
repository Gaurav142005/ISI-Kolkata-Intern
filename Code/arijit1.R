m=1000
K=c(100,200,500,1000,2000,5000)


#r=matrix(c(1/2,1/10,1/3,1/9,1/4,1/7,1/5,1/6),nrow=2)		#  Weibull distribution Mean=1/rate; shape=(2,2,5,5,8,8,11,11)
#r=matrix(c(2,2,1/2,1/2,1/5,1/5,1/8,1/8),nrow=2)		#  Weibull distribution Mean=Gamma(1+1/shape)/rate; shape=(1/2,1/2,1/5,1/5,1/8,1/8,1/11,1/11)
r=matrix(c(.25,.25, .5,.5, 1,1, 2,2),nrow=2)			#  Exponential distribution Mean=1/rate; Rate=(.25,.25, .5,.5, 1,1, 2,2)
#r=matrix(c(12,12,9,9,6,6,3,3 ),nrow=2)			#  Gamma Distribution with scale=1,  centering at 12, 9, 6,3 (12,12,9,9,6,6,3,3 )
#r=matrix(c(1,50,2,49,3,48,4,47),nrow=2)			# Mixture of two Gamma Distributions with scale=1,  centering at 5, (1,9,2,8,3,7,4,6) 
#r=matrix(c(4,10,1,4, 2,4, 3,4),nrow=2)			# Aymmetric Mixture of Beta (4,10;10,4); (1,4; 4,1); (2,4; 4,2); (3,4; 4,3)
#r=matrix(c(4,4,3,3,2,2,1,1),nrow=2)			# Symmetric Mixture of Beta (4,4,3,3,2,2,1,1)
#r=matrix(c(3,-3,2,-2,1,-1,0,0),nrow=2)			# Normal two Mixture (3,-3,2,-2,1,-1,0,0)
#r=matrix(c(4,-4,3,-3,2,-2,1,-1),nrow=2)   		# Two Normal distributions with means (3,-3,2,-2,1,-1,0,0) along with one standard normal
#r=matrix(c(5,-5, 4,-4,3,-3,2,-2),nrow=2)   		# t distribution with  df 6, 10 and 20 added in the loop, location (mean) parameter are (5,-5, 4,-4,3,-3,2,-2)
count=matrix(0,nrow=length(K),ncol=ncol(r))
for (di in 1:4)   
{
for (j in 1:length(K))
{n=K[j]
for (i in 1:m)
{
#x=rbeta(n, 2, 1, 0)
#x=rbeta(n, 15, 4, 0)

########################< mixture of exponential distribution>#########
a=r[1,di]
b= r[2,di]
x1=rexp(n, a)
x2=rexp(n, b)
u=runif(n,0,1)
x= (1-floor(2*u))*x1 + floor(2*u)*x2

#########################< mixture of Beta distribution>#########
#a=r[1,di]
#b= r[2,di]
#x1=rbeta(n, a, b,0)
#x2=rbeta(n, b,a,0)
#u=runif(n,0,1)
#x= (1-floor(2*u))*x1 + floor(2*u)*x2

#########################< mixture of Gamma distribution>#########
#s1=r[1,di]
#s2= r[2,di]
#x1=rgamma(n, s1, 1)
#x2=rgamma(n, s2, 1)
#u=runif(n,0,1)
#x= (1-floor(2*u))*x1 + floor(2*u)*x2

#########################< mixture of Weibull distribution>#########
#s1=r[1,di]
#s2= r[2,di]
#x1=rweibull(n, s1, 1)
#x2=rweibull(n, s2, 1)
#u=runif(n,0,1)
#x= (1-floor(2*u))*x1 + floor(2*u)*x2

#######################<mixture of t-distrn >################
#x1= r[1,di]+rt(n,20,0)
#x2= r[2,di]+rt(n,20,0)
#u=runif(n,0,1)
#x= (1-floor(2*u))*x1 + floor(2*u)*x2
#######################################
#
####################< mixture of Normal distribution >################
#
#x0=rnorm(n)
#x1= r[1,di]+rnorm(n)
#x2= r[2,di]+rnorm(n)
#u=runif(n,0,1)              #### We may not need to change mixing coefficient everytime #######
#x= (1-floor(2*u))*x1 + floor(2*u)*x2									     # mixture of 2 normal; either X0 or X1 or X2 is chosen		
#x= (1-floor(3*u))*(2-floor(3*u))*0.5*x0 + (2-floor(3*u))*floor(3*u)*x1 + (floor(3*u)-1)*floor(3*u)*0.5*x2    # mixture of 3 normal; either X0 or X1 or X2 is chosen

#########################################################
mu1=mean(x)
mu2=(sum((x - mu1)^2))/n
mu3=(sum((x - mu1)^3))/n
mu4=(sum((x - mu1)^4))/n
mu5=(sum((x - mu1)^5))/n
mu6=(sum((x - mu1)^6))/n
mu7=(sum((x - mu1)^7))/n
mu8=(sum((x - mu1)^8))/n
mu9=(sum((x - mu1)^9))/n
mu10=(sum((x - mu1)^10))/n                                                                                                                                                                                                       
mu11=(sum((x - mu1)^11))/n
mu12=(sum((x - mu1)^12))/n
####################
KURT=mu4/mu2^2
Vhat=(mu6 - mu3^2)/mu2^3
###################
alp0=c(1, -(6*mu5 - 3*mu2*mu3), -mu3, -3*(mu2^2)*Vhat)
GAM=matrix(c(mu12 - mu6^2, mu7, mu9 - mu3*mu6, mu8 - mu2*mu6, mu7, mu2, mu4, mu3, mu9 - mu3*mu6, mu4, mu6 - mu3^2,
mu5 - mu2*mu3, mu8 - mu2*mu6, mu3, mu5 - mu2*mu3, mu4 - mu2^2), ncol=4, byrow=T)
####################
SIG2=(t(alp0)%*%GAM%*%alp0)/mu2^6
uppbd=9 + qnorm(0.95)*sqrt(SIG2/n)
lowbd=9 - qnorm(0.95)*sqrt(SIG2/n)
if (Vhat <= lowbd) 
{count[j,di]=count[j,di]+1
}
}
}
}
