## power simulation code 
## Laurel Brehm, originally 06/06/19 but updated  3/14/24

## to run all of these functions, you'll need the following libraries
library(lme4)  ## for mixed models
library(faux)  ## for faux factorial data
library(boot)  ## for bootstrapping and other resampling/transformations
library(dae)   ## for latin square generator
### install what you need like this-- only do once, not repeatedly...
##### install.packages('faux')
##### install.packages('dae') 

## there's currently an issue in my office computer version of R with the package Matrix-- lme4 depends on a specific version of it
## if you get the following error, function 'cholmod_factor_ldetA' not provided by package 'Matrix''
## run the following code to get the correct dependencies of lme4, then rerun the library line
## install.packages("lme4", type = "source")
## library(lme4)

## primer on a few key functions we're using today
## these are our data generation functions:
rnorm(10,m=1,sd=1)
rep(10,2)
seq(1:20)
## embedding anything in this makes it repeat n times
replicate(3,rep(10,2))
## these are our stats functions:
?t.test()
?lmer()

## here are some power simulation functions, starting simple and becoming more complicated.
powerFunc2Gr<- function(m1,sd1,m2,sd2,n){
  ## sample data according to inputs
  group1 <- rnorm(n/2,m1,sd1)
  group2 <- rnorm(n/2,m2,sd2)
  ## test differences between samples: is it significant?
  ts <- t.test(group1,group2)
  ps <- ts[3]  ## the p-value of the test
  ## the function will output whatever is on the last line-- here, output the pvalue
  ps  
}

## try running the function! what happens if you run it several times with the same values?
powerFunc2Gr(1,1,2,1,10)

## because we're sampling, we get different results!
## so, embed in function 'replicate' to run repeatedly.
out <- replicate(1000,powerFunc2Gr(m1=1,sd1=1,m2=2,sd2=1,n=100))
power01 <- mean(out<0.05)
## show the output
power01

## can also get more elaborate: loop over multiple sample sizes!
## look at increasing samples by 10 over 10 steps
## make a small data frame to put this in:
power02 <- rbind(seq(1:10)*10,rep(0,10))
rownames(power02)<- c("totalN","power")

## now, we loop over values s for subjects. we could loop over any of the input variables.
for(s in 1:10){
  out <- replicate(1000,powerFunc2Gr(1,1,2,1,s*10))  ##the default is that the inputs will be assigned to variables in the same order as defined when you wrote the function
  power02[2,s] = mean(out<0.05)
}
## show the output
power02


## same exercise, but looping over effect size, two dimensionally-- and we can plot it
## we will create a plot of effect size in terms of cohen's d = (m1-m2)/sd
## (if m1 and sd are always 1, then we can get different cohen's d values by increasing m2)
## https://en.wikipedia.org/wiki/Effect_size#Cohen's_d
power03=matrix(nrow=50,ncol=9)
rownames(power03)=paste0('subs=',1:50*10)
colnames(power03)=paste0('d=0.',1:9)

## now, we loop over values s for subjects, and d for cohen's d value (effect size)
## this one's a bit slow- need replicates= 5000 for a really nice plot. run smaller to see effects yourself
## 
for(d in 1:9){
  for(s in 1:50){
    out <- replicate(500,powerFunc2Gr(1,1,1+(d*.1),1,s*10))  ##the default is that the inputs will be assigned to variables in the same order as defined when you wrote the function
    power03[s,d] = mean(out<0.05)
}}
## show the output
power03

contour(x=(1:50*10),y=(1:9*.1),z=power03,
        xlab='sample size',ylab='cohen d',col=c(rep('gray',8),rep('black',3)))
### realistic social science effect sizes are .5 and below: you need more observations than you think you do


### 2x2 fully-crossed, Latin-square design with two repeated measures
## here's what that means:
## each person sees each item in one of four conditions of Variable1(a,b) and Variable2 (a,b)-- these are our 4 'groups'
## this gives us a total of s*i observations, distributed evenly across g (here, 4) groups
## let's say we have 12 people and 24 items-- so that's (12*24)/4 observations per group
g=4
s=12
i=24
gn= (i*s)/g
group1 <- rnorm(gn,m=1,sd=1)
group2 <- rnorm(gn,2,1)
group3 <- rnorm(gn,3,1)
group4 <- rnorm(gn,4,1)
rs <- c(group1,group2,group3,group4)

## stack together the group identifiers to create a 2x2 design
## in so doing, code contrasts for each variable: .5 and -.5, corresponding to effects coding
## groups 1 and 2 belong to v1, level A= .5, groups 3 and 4 belong to v1, level B=-.5
## groups 1 and 3 belong to v2, level A= .5, groups 2 and 4 belong to v2, level B=-.5
v1 <- c(rep(.5,gn*2), rep(-.5,gn*2))
v2 <- rep(c(rep(.5,gn), rep(-.5,gn)),2)

interaction.plot(x.factor=v1,trace.factor=v2,response=rs)


### embed the same simulated data set into a function
powerFuncLMB<- function(s,i){
  ## here s=n subjects and i= n items
  ## these are my two repeated measures-- subjects see multiple items, and items are seen by multiple subjects
  # !! make sure that s and i are divisible by number of groups in the design !! 
  
  ## observations per group, per person = i/groups
  sg <- i/4
  ## observations per group, per item = s/groups
  ig <- s/4
  ## observations per group, combining all people = (s*i)/groups
  g <- s*i/4
    
  ## set up dependent measure per group
  ##in these data, I put in 2 main effects-- go up 1 unit on response value r between each group
  ## concretely: the difference between 1 and 2 is the same as 3 and 4, and the difference between 1 and 3 is the same as 2 and 4
  ## but! there's a ton of variability at the level of the observations
  group1 <- rnorm(g,m=1,sd=1)
  group2 <- rnorm(g,2,1)
  group3 <- rnorm(g,3,1)
  group4 <- rnorm(g,4,1)
  rs <- c(group1,group2,group3,group4)

  ## stack together the group identifiers to create a 2x2 design
  ## in so doing, code contrasts for each variable: .5 and -.5, corresponding to effects coding
  ## groups 1 and 2 belong to v1, level A= .5, groups 3 and 4 belong to v1, level B=-.5
  ## groups 1 and 3 belong to v2, level A= .5, groups 2 and 4 belong to v2, level B=-.5
  v1 <- c(rep(.5,g*2), rep(-.5,g*2))
  v2 <- rep(c(rep(.5,g), rep(-.5,g)),2)
  
  ## assign participants to observations
  ## repeat sequence of subjects sg times, and then repeat that for all 4 groups
  ## and append 's' to the numeric for ease of tracking
  ss <- paste0('s',rep(rep(seq(1:s),sg),4))
  
  ## now paste together 
  ds <- as.data.frame(cbind(rs,v1,v2,ss))
  
  ##update types
  ds$rs <- as.numeric(as.character(ds$rs))
  ds$v1 <- as.numeric(as.character(ds$v1))
  ds$v2 <- as.numeric(as.character(ds$v2))
  ds$ss <- as.factor(ds$ss)
  
  ## sort by subject, by condition
  ds <- ds[order(ss),]
  
  ##assign items to observations
  ## in this design, each person sees each item once
  ## so, we want to make sure we rotate through as in a latin square
  ## what we are setting up here is like having 4 lists in your experiment.
  ##participant 1 has items assigned to conditions that start at number 1
  ## participant 2 starts condition 1 at ig+1, with 1:ig assigned to condition 4
  ## participant 3 starts condition 1 at 2ig+1, participant 4 starts condition 1 at 3ig+1
  ## repeat this ig times to fill out across participants and append 'i' to beginning
  ii <- paste0('i',rep(c(1:i, (ig+1):i, 1:ig, (2*ig+1):i, 1:(2*ig), (3*ig+1):i, 1:(3*ig)),ig))
  
  ##stick this column on, update type
  ds$ii <- ii
  ds$ii <- as.factor(ds$ii)
  
  ## adjust our response values based upon the subject and item they came from
  ## these are 'random effects'
  ## since they're an 'adjustment', the mean of the distribution is always 0
  ### unless these are quite large, you'll get a 'boundary (singular) fit' warning, so we will suppress that 
  sr <- rnorm(s,m=0,sd=.1)
  ir <- rnorm(i,m=0,sd=.1)
  
  ## find the place in each vector that corresponds to the subject and item, and add those adjustments
  ds$rs <- ds$rs + sr[as.numeric(ds$ss)] + ir[as.numeric(ds$ii)] 

  ## now we can test a model!
  lmer1 <- suppressMessages(lmer(rs~v1*v2 + (1+v1+v2|ii) + (1+v1+v2|ss),ds,REML=F,
                                 control = lmerControl(calc.derivs = FALSE)))
  ### suppress the warning messages since I know we'll get a lot of singular fit warnings
  ### suppress calc.derivs to speed up performance
  
  ## to test p-values from a mixed effect model, compare to a model without the effect of interest
  ## examine that interaction: it should have power = .05 (our alpha/significance level)
  lmer2 <- suppressMessages(lmer(rs~v1+v2 + (1+v1+v2|ii) + (1+v1+v2|ss),ds,REML=F,
                                 control = lmerControl(calc.derivs = FALSE)))
  ps <- anova(lmer1,lmer2)[2,8]
  ps  
}

## I'd usually run this at 1000 replicates, but for speed here, run only 100:
out <- replicate(100,powerFuncLMB(s=12,i=24))
summary(out)
power04<- mean(out<0.05)   ###power is at just about 0.05-- what we expect!
power04

## run again with fewer samples (fewer s, fewer i)
out <- replicate(100,powerFuncLMB(s=4,i=4))
power05 <- mean(out<0.05) 
power05

## with few subjects, few items, the power level can be pretty divergent from what we expect
## this is why underpowered designs can reflect false positives
## it's the same reason why you need to run these simulations many times to see a 'true' result'
## make sure to run enough replicates in R and in real life so that you're estimating the central tendency of the effect, not the extremes



## here's a version of a power analysis for a real study
## this one uses a binomial (0/1) DV
## also create the data set using 'faux' package
## I'm going to use this to calculate power for Bethany's thesis study on memory and prosody--
## it is a replication of Fraundorf et al 2010, in the visual world paradigm.  
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2935187/
##
##items are like:
####Both the British and the French biologists had been searching Malaysia and Indonesia for
#### the endangered monkeys. Finally, the BRITISH spotted one of the monkeys in Malaysia and 
#####planted a radio tag on it.
#### dv = memory for items based upon prosody, 2x2 design
## accented = word is pitch-accented or not (levels: LHStar, Hstar, contrasts .5, -.5)
## position =  (levels: First, Second, contrasts .5, -.5)
## we can reconstruct means for each beta using the model table in the paper
## we can consider these a high-end estimate of the effect size 
## (an effect the same size as before, which could be inflated)
intercept = 1.92
accented = 1.21 * c(-.5,0.5)
position = 0.32 * c(-.5,0.5)
interaction = 0.7* c(.25,-.25,-.25,.25)
betas = c(intercept + accented[1] + position[1] + interaction[1],
          intercept + accented[1] + position[2] + interaction[2],
          intercept + accented[2] + position[1] + interaction[3],
          intercept + accented[2] + position[2] + interaction[4])
probs1 = inv.logit(betas)
### then we can also create a smaller, but still reasonable effect size
## this one has a difference that's 2/3 of the original size (10%, not 15%)
## it would still be an effect size we'd care about, so let's make sure we have power to observe it
## the technical term here is 'SESOI'-- smallest effect size of interest
## do this by adding .05 to each of the first two (smaller) numbers
probs2 = probs1 + c(.05,.05,0,0)

powerFuncBMB<- function(s,i,probs){
    ## first, create the data set with faux
    within = list(accented=c("HStar","LHStar"),position=c("First","Second"))
    groups = 4 ## the combo of all within vars
    rps = i/4  ## replicates-- number of items per condition
    between = list(ss = paste0('s',1:s))
    ## create means from probabilities
    means =rep(probs,s)
    ## generate design matrix
    ds <- sim_design(n=(i/groups), within, between, mu=means, sd=.05,long=T,plot=F)
    ## we want to get i/conditions repeats per person-- the conditions are the 4 combos of accented and position
    ## mu is our response rates per condition on scale of probability, put into our 4 item sets
    ## sd is a 'reasonable' rate of 5% in probability space-- I'm putting this in as a guess since I don't actually know the sds around the means
    ## this is an ok thing to do as long as you're conservative!
    ## order the design by participants
    ds <- ds[order(ds$ss),]
    
    ## now create an item vector to add on top of this-- 
    ## we are going to use a function to give us a latin square, and repeat it as many times as we need, then take of that what we need
    iin = paste0("ii",rep(designLatinSqrSys(i),ceiling(s/i)))
    ds$ii = as.factor(iin[1:(s*i)])
    
    ## generate random effects for participants
    ir  <- rep(rnorm((i/rps),m=0,sd=.05),each=s)
    sr  <- rep(rnorm(s,m=0,sd=.05),each=(i/rps))
    
    ## find the place in each vector that corresponds to the subject and item, and add those adjustments
    ds$y <- ds$y + sr[as.numeric(ds$ss)] + ir[as.numeric(ds$ii)] 
    
    ## reassign any too big and too small values, if there are any
    if(dim(ds[ds$y > 1,])[1]!=0){
       ds[ds$y > 1,]$y <- 1} else{}
    if(dim(ds[ds$y < 0 ,])[1]!=0){
       ds[ds$y < 0,]$y <- 0} else{}
    
    ## take the probability values and turn them into a binomial
    ds$y_binom <- rbinom(s*i,size=1,prob=ds$y)
    
    ##reassign the first value if there's complete separation of values as ones 
    ##(would need to do similar in the other direction if proportions are quite close to zero)
    ifelse(mean(ds[ds$accented == "HStar" & ds$position=="First",]$y_binom) == 1,
        ds[ds$accented == "HStar" & ds$position=="First",]$y_binom <- c(0,rep(1,(i*s/g-1) )),'ok')
    ifelse(mean(ds[ds$accented == "HStar" & ds$position=="Second",]$y_binom) == 1,
        ds[ds$accented == "HStar" & ds$position=="Second",]$y_binom <- c(0,rep(1,(i*s/g-1) )),'ok')
    ifelse(mean(ds[ds$accented == "LHStar" & ds$position=="First",]$y_binom) == 1,
        ds[ds$accented == "LHStar" & ds$position=="First",]$y_binom <- c(0,rep(1,(i*s/g-1) )),'ok')
    ifelse(mean(ds[ds$accented == "LHStar" & ds$position=="Second",]$y_binom) == 1, 
        ds[ds$accented == "LHStar" & ds$position=="Second",]$y_binom <- c(0,rep(1,(i*s/g-1) )),'ok')
    
    ## set contrasts
    contrasts(ds$accented) <- c(-.5,0.5)   
    contrasts(ds$position) <- c(-.5,0.5)

    ## run a model
    glmer1 <- glmer(y_binom ~ accented*position + (1|ii) + (1|ss),ds,family='binomial',
                    control = glmerControl(calc.derivs = FALSE))

    ## extract the pval of interest  -- it will be (n rows) * (n completed rows)+(position of item in last row) - here (4*3)+2=14
    unlist(summary(glmer1)[10])[14]
}

## make 2 small data frames to put the results in:
bpower1 <- rbind(seq(1:6)*8,rep(0,6))
rownames(bpower1)<- c("totalN","power")
bpower2 <- bpower1

## here, we loop over s for subjects. we could loop over any of the input variables.
## do for both effect sizes-- the original, and the 67%
for(ss in 1:6){
  out <- replicate(100,powerFuncBMB(s=ss*8,i=20,probs=probs1))
  bpower1[2,ss] = mean(out<0.05)}
bpower1

for(ss in 1:6){
  out <- replicate(100,powerFuncBMB(s=ss*8,i=20,probs=probs2))
  bpower2[2,ss] = mean(out<0.05)}
bpower2
