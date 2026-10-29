# Kanamycin schedules for the noise experiment: 30 daily transfers, 6 patches per regime
# White noise: truncated normal series (mean 2, sd 1, lower bound 0)
# Red noise: the same values reordered to a lag-1 autocorrelation of 0.8
# Writes the series and the daily kanamycin volumes (µl of the 1 mg/ml stock) to schedules/
# run from the repository folder

require(truncnorm)
library(tidyverse)

find_tr_mu_sd = function(mu_desired, sd_desired, a) {
  # mu_desired: desired mean
  # sd_desired: desired sd
  # a: lower bound
  
  optfun = function(par) {
    mu_parametric = par[1]
    sd_parametric = par[2]
    
    alpha = (a - mu_parametric)/sd_parametric
    Z = 1 - pnorm(alpha)
    
    mu_realised = mu_parametric + sd_parametric*dnorm(alpha)/Z
    
    sd_realised = sqrt(sd_parametric^2*(1+alpha*dnorm(alpha)/Z-(dnorm(alpha)/Z)^2))
    
    sum(c((mu_realised-mu_desired)^2, (sd_realised-sd_desired)^2))
  }
  
  optout = optim(par = c(mu_desired, sd_desired), optfun)
  optout$par
}


# generate random variable with desired distribution
tr_pars = find_tr_mu_sd(mu_desired = 2, sd_desired = 1, a = 0) 

#x = rtruncnorm(n = 30, a = 0,
#               mean = tr_pars[1], sd = tr_pars[2])
#mean(x)
#sd(x)
#hist(x)

set.seed(2023)
white_noise <- as.data.frame(replicate(n = 6, rtruncnorm(n = 30, a = 0,
                            mean = tr_pars[1], sd = tr_pars[2]) , simplify = T))

par(mfrow = c(2, 3)) 
for (i in 1:ncol(white_noise)){
  plot(white_noise[,i], type = 'l')
}

apply(white_noise, 2, mean)
apply(white_noise, 2, sd)

write.table(white_noise, "schedules/white_noise_mean2_sd1.txt", sep="\t", dec= ",") 

par(mfrow = c(2, 3)) 
for (i in 1:ncol(white_noise)){
  plot(hist(white_noise[,i]))
}

#Once you've done that, you can restructure x to add whatever 
#temporal autocorrelation you want. Below are two functions - 
#one that lets you tune autocorrelation, and one that should 
#automatically let you set any level of autocorrelation that you want.


# add temporal autocorrelation to the resulting variable
redify_noise = function(x, weight = 2) {
  # x: a random variable
  # weight: how strong should autocorrelation be?#
  # higher weight gives stronger autocorr
  
  xred = numeric(length(x))
  xused = rep(0, length(x))
  
  startpoint = sample(length(x), 1)
  xred[1] = x[startpoint]
  xused[startpoint] = 1
  
  for(i in 2:length(x)) {
    wt = exp(-abs(weight) * abs(xred[i-1]-x[xused == 0]))
    wt = wt/sum(wt)
    
    if(length(wt) > 1) {
      ps = sample(which(xused == 0), 1, prob = wt)
    } else {
      ps = which(xused == 0)
    }
    
    xred[i] = x[ps]
    xused[ps] = 1
  }
  
  xred
}

#redx = redify_noise(x, weight = 2)
#mean(redx)
#sd(redx)
#hist(redx)
# unchanged mean and sd

#plot(redx, type = 'l')
#cor(redx[-1], redx[-length(redx)])
#cor(x[-1], x[-length(x)])
# but with added autocorrelation

#redify the 6 series 
set.seed(2023)
redified <- as.data.frame(apply(white_noise,2,redify_noise))
write.table(redified, "schedules/redified_mean2_sd1_weight2.txt", sep="\t", dec= ",") 

par(mfrow = c(2, 3)) 
for (i in 1:ncol(redified)){
  plot(redified[,i], type = 'l')
}

par(mfrow = c(2, 3)) 
for (i in 1:ncol(redified)){
  hist(redified[,i])
}

apply(redified, 2, mean)
apply(redified, 2, sd)

# generate series with fixed autocorrelation strength
fixedred = function(x, cor_desired = 0.8, weight0 = 1, maxiter = 1000, delta = 0.05, stepsize = 0.1) {
  # x: random variable
  # cor_desired: what is the desired temporal autocorrelation? must be between 0 and 1
  # weight0: starting point for weight - defaults to 1
  # maxiter: maximum number of iterations
  # delta: tolerance (i.e. how close to cor_desired to we need to get?)
  # stepsize: size of steps in weights used by optimizer
  
  weightuse = weight0
  n = 0
  while(n < maxiter) {
    xred = redify_noise(x, weight = weightuse)
    xcor = cor(xred[-1], xred[-length(xred)])
    
    if(abs((xcor-cor_desired)/cor_desired) < delta) {
      n = maxiter
    } else {
      weightuse = weightuse + sign(cor_desired - xcor)*stepsize
    }
  }
  
  return(list(xred=xred, xcor = xcor, weight = weightuse))
}


#red_out = fixedred(x, cor_desired = 0.9)
#plot(red_out$xred, type = 'l')
#red_out$xcor
#hist(red_out$xred)


#redify the 6 series 
set.seed(2023)
fixed    <- lapply(white_noise, fixedred, cor_desired = 0.8)
fixed_corr_red <- cbind.data.frame(V1=fixed$V1$xred, V2=fixed$V2$xred, V3=fixed$V3$xred, 
                   V4=fixed$V4$xred, V5=fixed$V5$xred, V6=fixed$V6$xred)

write.table(fixed_corr_red, "schedules/fixed_corr_red_mean2_sd1_corr0.8.txt", sep="\t", dec= ",") 

par(mfrow = c(2, 3)) 
for (i in 1:ncol(fixed_corr_red)){
  plot(fixed_corr_red[,i], type = 'l')
}

apply(fixed_corr_red, 2, mean)
apply(fixed_corr_red, 2, sd)

par(mfrow = c(2, 3)) 
for (i in 1:ncol(fixed_corr_red)){
  hist(fixed_corr_red[,i])
}  
  
#calculate the amount will be taken from the stock solution 1mg/ml
  
white.kan.mean2.sd1    <- round((white_noise * 6), 1)
red.kan.mean2.sd1.w2   <- round((redified * 6), 1)
red.kan.mean2.sd1.c0.8 <- round((fixed_corr_red * 6), 1)
  
write.table(white.kan.mean2.sd1, "schedules/white.kan.mean2.sd1.txt", sep="\t", dec= ",") 
write.table(red.kan.mean2.sd1.w2, "schedules/red.kan.mean2.sd1.w2.txt", sep="\t", dec= ",") 
write.table(red.kan.mean2.sd1.c0.8, "schedules/red.kan.mean2.sd1.c0.8.txt", sep="\t", dec= ",")  

#I'll use white and fixedred series

#Constant treatment
constant.d <- summary(rbind(gather(white_noise), gather(fixed_corr_red)))
  
#control 1-6
round(0.2558 * 6, 1)
round(1.4679 * 6, 1)
round(2.1644 * 6, 1)
round(2.1540 * 6, 1)
round(2.7264 * 6, 1)
round(4.9004 * 6, 1)
  
