
library(pracma)

S_0<- 50 # Initial Price
K<-50.00 #Strike Price
r_0<- 0.07
sigma<-0.13 #Volatility sigma= 13%
MT<- 1 # Maturity Time T=1 year
d<- 12 #Discretize T year into d subintervals -----Start with 10^4 for d = 12,and n = 1000 for d = 52
delta<-MT/d
epsilon<-.05 # error tolerance
n<-10^4 #Number of sample paths. Start with 10^4 for d = 12,and n = 1000 for d = 52
tic()

#..~..~..~..~..~..~..~..~..~..~..~..~..~..~..~..~..~..~..~..~..~..~..~..~..~..
#generate n*d random normal numbers----------------
r<- matrix(numeric(n*d),nrow=n) #just the frame 
z<- matrix(rnorm(n*d),nrow=n) # this fills the frame

#--------------generate sample paths of interest rate 
r[1,]<- r_0+0.18*(0.086-r_0)*delta+0.02*sqrt(delta)*z[1,]# generates the samples of 1st column
for (i in 2:d){   
  r[i,]<-r[i-1,]+0.18*(0.086-r[i-1,])*delta+0.02*sqrt(delta)*z[i,] #second column 
}
#--------------------------------------------------Asset Price
#generate random normal numbers----------------

S<- matrix(numeric(n*d),nrow=n)
x<- matrix(rnorm(n*d),nrow=n)

#-generate sample path of Asset Price

S[,1]<-S_0*exp((r_0-sigma^2/2)*delta+sigma*sqrt(delta)*x[,1])
for (i in 2:d){
  S[,i]<-S[,i-1]*exp((r[,i-1]-sigma^2/2)*delta+sigma*sqrt(delta)*x[,i])
}

#------------------------------------------------
#generate the instances of discounted payoff( this is wrong, see notes)
EuroCallPayoff<-pmax(S[,d]-K,0)*exp(-apply(r,1,sum)*delta) 
EuroPutPayoff<-pmax(K-S[,d],0)*exp(-apply(r,1,sum)*delta)

#4a) estimate the option price using the sample mean of discounted payoff
EuroCallPrice<-mean(EuroCallPayoff)
EuroPutPrice<-mean(EuroPutPayoff)

#------------------------------------------------
#calculate the estimation error based on 99% confidence level
est_Callerror<-2.58*sd(EuroCallPayoff)/sqrt(n)
#est_Puterror<-2.58*sd(EuroPutPayoff)/sqrt(n)

#99% confidence interval for the option prices
CI_Call <- c(EuroCallPrice-est_Callerror,EuroCallPrice+est_Callerror)
CI_Put <- c(EuroPutPrice-est_Puterror,EuroPutPrice+est_Puterror)

EuroCallPrice
CI_Call



r_BS<-0.07

#exact option prices using Black-Scholes formula
ExactEuroCall<-S_0*pnorm((log(S_0/K)+(r_BS+sigma^2/2)*MT)/(sigma*sqrt(MT)))-K*exp(-r_BS*MT)*pnorm((log(S_0/K)+(r_BS-sigma^2/2)*MT)/(sigma*sqrt(MT)))
ExactEuroPut<-K*exp(-r_BS*MT)*pnorm((log(K/S_0)-(r_BS-sigma^2/2)*MT)/(sigma*sqrt(MT)))-S_0*pnorm((log(K/S_0)-(r_BS+sigma^2/2)*MT)/(sigma*sqrt(MT)))


ExactEuroCall
ExactEuroPut

#true relative errors of the estimated option prices
error_EuroPut<-abs(EuroPutPrice-ExactEuroPut)/ExactEuroPut
error_EuroCall<-abs(EuroCallPrice-ExactEuroCall)/ExactEuroCall

#print out the results
cat("With S(0)=",S_0, ", r=",r_BS,", volatility=",sigma, ", T=", MT, ", K=",K, ", monitoring frequency=",d, ", and sample paths n=",n, ". The estimated European call option price is" , EuroCallPrice, ", with relative error as", error_EuroCall, ".")
#------------------------------------------------
------------------------------------------------
toc() 
#Required Sample Size
v<-rnorm(n)
hat_sigma<-sd(v) #calculate the sample standard deviation
C <- 1.1 #amplifying constant
N<-ceiling((2.58*C*hat_sigma/epsilon)^2)


#do the estimation using another independent sample with N instances
y<-rnorm(N)
est<-mean(y)
error<-2.58*sd(y)/sqrt(N)

N
est
error
