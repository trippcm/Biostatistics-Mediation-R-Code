### The following R code implements the examples used in the paper
### "A Classical Regression Framework for Mediation Analysis: 
### Fitting One Model to Estimate Mediation Effects"
### by Christina T. Saunders and Jeffrey D. Blume 
### updated August 5, 2017

### Although the BRAIN-ICU data (as described in the manuscript) cannot be made publically available,
### the following code is used to implement the examples in the manuscript. 

###########################################
####### CODE TO IMPLEMENT EXAMPLES ########
###########################################


### Load libraries
library(MASS) # for mvrnorm() function
library(xtable) # for creating tables of results

set.seed(12345) # set seed for reproducibility
nboot <- 10000 # number of bootstrap replications
mcs <- 10000 # number of monte carlo simulations

### function to estimate the indirect effect (and its variance) 
### for a unit change in the exposure for continuous x, m, y 
### from a *simple mediation model*
### using the Essential Mediation Components (EMCs) formulae, 
### Sobel's formula, bootstrapping, and Monte Carlo methods.
### This function is used in Example 1 Part 1 and Example 2.


### NOTE: When (x-x*) = 1, the EMC from the simple mediation model
### is equal to the indirect effect and the portion eliminated.

mediation_fxn <- function(x,m,y){
  #### Fit simple mediation models
  full <-lm(y ~ x + m)
  sub <- lm(y ~ x)
  pathA <- lm(m ~ x)
  
  # functionals from full model for using formula
  v.xm <- vcov(full)["x","m"]
  v.mx <- vcov(full)["m","x"]
  v.mm <- vcov(full)["m","m"]
  
  # fitted values and residuals for residual bootstrap
  f.fit <- full$fitted.values
  f.resid <- full$residuals
  
  # extract estimated coefficients
  ahat <- pathA$coef["x"]
  bhat <- full$coef["m"]
  chat <- sub$coef["x"]
  cprimehat <- full$coef["x"]
  
  # extract sample variances
  svar.c <- vcov(sub)["x","x"]
  svar.cprime <- vcov(full)["x","x"]
  svar.a <- vcov(pathA)["x","x"]
  svar.b <- vcov(full)["m","m"]
  
  # create var-cov matrix for Monte Carlo simulation
  # (this specification assumes cov=0)
  Sigma.diff <- matrix(c(svar.c,0,0,svar.cprime),byrow=T,ncol=2)
  Sigma.prod <- matrix(c(svar.a,0,0,svar.b),byrow=T,ncol=2)
  
  ### Formula for the EMC and its variance ####
  emc <- v.xm %*% solve(v.mm) * (-1*bhat)
  var.emc <- v.xm * solve(v.mm) * v.mx
  var.bm <- v.mm
  
  ### Calculate IDE using other methods ###
  ide.diff <- chat- cprimehat
  ide.prod <- ahat*bhat
  
  ### Calculate Var(IDE) using delta method approximations
  var.sobel <- ahat^2*svar.b + bhat^2*svar.a
  var.exact <- ahat^2*svar.b + bhat^2*svar.a + svar.a*svar.b
  var.unbiased <- ahat^2*svar.b + bhat^2*svar.a - svar.a*svar.b
  
  ############################
  ######  Bootstrap   ########
  ############################
  
  ## bootstrap cases
  boot.ide <- boot.c<- boot.cprime<- boot.a<- boot.b<- rep(NA,nboot)
  for(k in 1:nboot){
    # sample with replacement from the rows of (y,x,m)
    dat <- cbind(y,x,m)
    tmp <- sample(1:nrow(dat),nrow(dat),replace=TRUE)
    dat.tmp <- as.data.frame(dat[tmp, ])
    # models
    mf <- lm(y ~ x + m, data=dat.tmp)
    ms <- lm(y ~ x, data=dat.tmp)
    patha <- lm(m ~ x, data=dat.tmp)
    # store results
    boot.ide[k] <- ms$coef["x"] - mf$coef["x"]
    boot.c[k] <- ms$coef["x"]
    boot.cprime[k] <- mf$coef["x"]
    boot.a[k] <- patha$coef["x"]
    boot.b[k] <- mf$coef["m"]
  }
  
  ## bootstrapped covariances of model coefficients from bootstrapping cases
  bootcov.diff <- cov(boot.c, boot.cprime)
  bootcov.prod <- cov(boot.a, boot.b)
  bootrho.diff <- bootcov.diff/(sd(boot.c)*sd(boot.cprime))
  bootrho.prod <- bootcov.prod/(sd(boot.a)*sd(boot.b))
  
  #### bootstrap residuals 
  bootr.ide<- bootr.c<- bootr.cprime<- bootr.a<- bootr.b <- rep(NA,nboot)
  for(k in 1:nboot){
    f.resid.dat <- f.resid
    # sample with replacement from the residuals from the full model to get E*.full
    f.tmp <- sample(1:length(f.resid.dat),length(f.resid.dat),replace=TRUE)
    f.resid.dat.tmp <- (f.resid.dat[f.tmp])
    # bootstrapped Y values: Y*.full = Yhat.full + E*.full
    ystar.f <- f.fit + f.resid.dat.tmp
    
    # regress bootstrapped Y* on fixed design matrix to obtain bootstrap regression coefficients
    mf <- lm(ystar.f ~ x + m)
    ms <- lm(ystar.f ~ x)
    patha <- lm(m~x)
    bootr.ide[k] <- ms$coef["x"] - mf$coef["x"]
    bootr.c[k] <- ms$coef["x"]
    bootr.cprime[k] <- mf$coef["x"]
    bootr.a[k] <- patha$coef["x"]
    bootr.b[k] <- mf$coef["m"]
  }
  
  ###################
  ### Monte Carlo ###
  ###################
  Sigma.diff.cov <- matrix(c(svar.c,bootcov.diff,bootcov.diff,svar.cprime),byrow=T,ncol=2)
  Sigma.prod.cov <- matrix(c(svar.a,bootcov.prod,bootcov.prod,svar.b),byrow=T,ncol=2)
  mc.diff <- mc.prod<- mc.diff.cov<- mc.prod.cov<- rep(NA,mcs)
  for(l in 1:mcs){
    # assuming off-diagonals of Sigma are 0
    draw.diff <- mvrnorm(n=1,mu=c(chat,cprimehat),Sigma=Sigma.diff)
    mc.diff[l] <- draw.diff[1]-draw.diff[2]
    
    draw.prod <- mvrnorm(n=1,mu=c(ahat,bhat),Sigma=Sigma.prod)
    mc.prod[l] <- draw.prod[1]*draw.prod[2]
    
    ## to use bootstrapped covariances:
    # draw.diff.cov <- mvrnorm(n=1, mu=c(chat,cprimehat), Sigma=Sigma.diff.cov)
    # mc.diff.cov[l] <- draw.diff.cov[1] - draw.diff.cov[2]
    
    # draw.prod.cov<- mvrnorm(n=1,mu=c(ahat,bhat),Sigma=Sigma.prod.cov)
    # mc.prod.cov[l]<- draw.prod.cov[1]*draw.prod.cov[2]
  }
  
  ###################
  ##### 95% CIs #####
  ###################
  
  ## function to output confidence intervals
  pretty95ci <- function(lower,upper){paste("(",lower,", ",upper,")",sep="")}
  
  ##### Sobel Intervals
  sobel.ci <- c(ide.prod - qnorm(.975,0,1)*sqrt(var.sobel),
               ide.prod + qnorm(.975,0,1)*sqrt(var.sobel))
  sobel.ci.out <- pretty95ci(round(sobel.ci[1],3),round(sobel.ci[2],3))
  
  exact.ci <- c(ide.prod - qnorm(.975,0,1)*sqrt(var.exact),
               ide.prod + qnorm(.975,0,1)*sqrt(var.exact))
  exact.ci.out <- pretty95ci(round(exact.ci[1],3), round(exact.ci[2],3))
  
  unbiased.ci <- c(ide.prod - qnorm(.975,0,1)*sqrt(var.unbiased),
                  ide.prod + qnorm(.975,0,1)*sqrt(var.unbiased))
  unbiased.ci.out <- pretty95ci(round(unbiased.ci[1],3),round(unbiased.ci[2],3))
  
  #### Resampling Intervals
  boot.ci <- quantile(boot.ide,c(0.025,0.975))
  bootr.ci <- quantile(bootr.ide, c(0.025, 0.975))
  
  #### MC with covariance = 0
  mcdiff.ci <- quantile(mc.diff,c(0.025,0.975))
  mcprod.ci <- quantile(mc.prod,c(0.025,0.975))
  
  # using EMC formula and t quantiles
  p <- 3 # number of parameters (intercept, coefficient for x, coefficient for m)
  N <- length(y)
  emc.ci <- c(ide.diff + v.xm * solve(v.mm) * qt(.975,df = N-p)*sqrt(var.bm),
             ide.diff - v.xm * solve(v.mm) * qt(.975,df = N-p)*sqrt(var.bm))
  
  
  #########################
  ##### Store results #####
  #########################
  
  # IDE results
  ide.results <- rbind(emc,
                       ide.diff,
                       ide.prod,
                       mean(boot.ide),
                       mean(bootr.ide),
                       mean(mc.prod),
                       mean(mc.diff))
  rownames(ide.results) <- c("emc",
                             "ide.diff",
                             "ide.prod",
                             "bootcase.ide",
                             "bootresid.ide",
                             "mc.prod",
                             "mc.diff")
  # Varhat(IDEhat) results
  var.results <- rbind(var.emc,
                       var.sobel,
                       var.exact,
                       var.unbiased,
                       var(boot.ide),
                       var(bootr.ide),
                       var(mc.prod),
                       var(mc.diff))
  rownames(var.results) <- c("var.emc",
                             "var.sobel",
                             "var.exact",
                             "var.unbiased",
                             "var.boot.cases",
                             "var.bootresid",
                             "var.mcprod",
                             "var.mcdiff")
  # SE(IDE) results
  sdtab <- round(rbind(sqrt(var.emc),
                      sqrt(var.sobel),
                      sqrt(var.exact),
                      sqrt(var.unbiased),
                      sd(boot.ide),
                      sd(bootr.ide),
                      sd(mc.prod),
                      sd(mc.diff)),2)
  rownames(sdtab) <- c("EMC",
                      "Sobel",
                      "Exact",
                      "Unbiased",
                      "Boot cases",
                      "Boot resid",
                      "MC ab",
                      "MC c-c'")
  
  # CI results 
  citab <- round(rbind("EMC CI:"=emc.ci,
                      "Sobel CI:"=sobel.ci,
                      "Exact CI:"=exact.ci,
                      "Unbiased CI:"=unbiased.ci,
                      "Boot case CI:"= boot.ci,
                      "Boot resid CI:" = bootr.ci,
                      "MC ab:"= mcprod.ci,
                      "MC c-c'"= mcdiff.ci),2)
  
  list("IDE results" = ide.results,
       "VAR results" = var.results,
       "Boot cases coefs" = list("a" = boot.a,"b" = boot.b,"c" = boot.c,"cprime" = boot.cprime),
       "Boot residual coefs" = list("a" = bootr.a, "b" = bootr.b, "c" = bootr.c, "cprime" = bootr.cprime),
       "SDs" = sdtab,
       "CI tab" = citab)
}

#######################################
## Example 1: Simple mediation model ##
#######################################

## Log baseline creatinine to SOFA to S100B (N=121)

# load data and assign variables
load("example1.Rdata")
y <- s100b # outcome y
m <- sofa # mediator m
x <- bl.cr # predictor x

## call mediation_fxn() function
example1<- mediation_fxn(x,m,y)

# output results
results<- example1
xtable(results[["IDE results"]],digits=4, caption = "IDEhat")
xtable(results[["SDs"]], caption="Standard errors")
xtable(results[["VAR results"]],digits=4, caption = "VARhat(IDEhat)")
xtable(results[["CI tab"]], caption="Confidence intervals")

#####################################################
## Example 1 part 2: mediation with X and X^2 term ##
#####################################################

xsquared<- x^2
ptm <- proc.time() # start clock
full <- lm(y ~ x + xsquared + m)
v.xm <- vcov(full)[c("x","xsquared"),"m"]
v.m <- vcov(full)["m","m"]
beta.m <- full$coefficients["m"]
emc <- -1 * v.xm %*% solve(v.m) %*% beta.m ## EMCs
emc[1] + emc[2] ## IDEhat
var.emc <- v.xm %*% solve(v.m) %*% t(v.xm) ## Var(EMCs)
sqrt(diag(var.emc)) ## SDs of linear and quadratic components of IDEhat
sqrt(var.emc[1,1] + var.emc[2,2] + 2*var.emc[1,2]) # Var(IDEhat)
proc.time() - ptm # Stop the clock

# how to implement with mediate fxn
library(mediation)
library(mgcv)
ptm <- proc.time() # start clock
model.m <- gam(m ~ x + I(x^2))
model.y <- gam(y ~ x + I(x^2) + m)
med <- mediate(model.m = model.m, model.y = model.y, 
               sims=5000, 
               treat=c("x","I(x^2)"),
               mediator="m", boot=TRUE)
proc.time() - ptm # Stop the clock

summary(med)
sd(med$d1.sims) ## SE(IDE)

## compare total effect to fitted model -- they match!
# med$tau.coef # 232.8977
# lm(y ~ x + I(x^2)) # 81.83 + 151.07

#######################################
## Example 2: Simple mediation model ##
#######################################

# assign variables
x <- bl.cr
m <- daily.sofa
y <- rbans.global.score

example2 <- mediation_fxn(x,m,y)

# output results
results<- example2
xtable(results[["IDE results"]],digits=4, caption = "IDEhat")
xtable(results[["SDs"]], caption="Standard errors")
xtable(results[["VAR results"]],digits=4, caption = "VARhat(IDEhat)")
xtable(results[["CI tab"]], caption="Confidence intervals")

###########################################################################
## Example 3: Model with confounders and exposure-confounder interaction ##
###########################################################################

# assign variables
x <- bl.cr # exposure
m <- daily.sofa # mediator
y <- rbans.global.score # outcome
conf <- charlson.score # confounder

full <- lm(y ~ x + m + conf + x:conf)
m.mod <- lm(m ~ x + conf + x:conf)
sub <- lm(y ~ x + conf + x:conf)

## formula using EMCs
emc <- -vcov(full)[c("x","x:conf"),"m"] %*% solve(vcov(full)["m","m"]) %*% full$coef["m"]
x.treat <- 1 # x
x.ctrl <- 0 # x*
(emc[1] + emc[2]*mean(conf))*(x.treat-x.ctrl) ## portion eliminated, which equals the NIE

## mediation formula
(m.mod$coef["x"] + m.mod$coef["x:conf"]*mean(conf))*(full$coef["m"])*(x.treat - x.ctrl)

## compare time using EMC formula and mediate function

## using EMC formula...
x.treat <- 1
x.ctrl <- 0
ptm <- proc.time() # start clock
full <- lm(y ~ x + m + conf + x:conf)
sub <- lm(y ~ x + m + conf + x:conf)
emc <- -vcov(full)[c("x","x:conf"),"m"]%*% solve(vcov(full)["m","m"]) %*% full$coef["m"]
(ide <- (emc[1] + emc[2]*mean(conf))*(x.treat-x.ctrl))
(cde <- full$coef["x"] + full$coef["x"]*mean(conf))
(te <- sub$coef["x"] + sub$coef["x"]*mean(conf))
# variance
v.emc <- vcov(full)[c("x","x:conf"),"m"]%*% solve(vcov(full)["m","m"]) %*% vcov(full)["m",c("x","x:conf")]
v.ide <- (v.emc[1,1]) + (mean(conf)^2 * v.emc[2,2]) + (2*mean(conf)*v.emc[1,2])
proc.time() - ptm # Stop the clock

### look at the 75th percentile of X vs the mean
x.treat <- quantile(x, p=.75)
x.ctrl <- quantile(x, p=.50)
ptm <- proc.time() # start clock
(emc[1] + emc[2]*mean(conf))*(x.treat-x.ctrl)
proc.time() - ptm

### repeat using the mediate function...
full <- lm(y ~ x + m + conf + x:conf)
m.mod <- lm(m ~ x + conf + x:conf)

ptm <- proc.time() # start clock
med <- mediate(model.m = m.mod, 
               model.y = full, 
               treat = "x", 
               mediator="m",
               treat.value = 1,
               control.value = 0,
               sims = 5000)
proc.time() - ptm

sd(med$d1.sims[1, ])## SE(IDE)
summary(med)

## compare 75th to 50th percentiles of exposure
x.treat <- quantile(x,p=.75)
x.ctrl <- mean(x)
system.time(
  med2 <- mediate(model.m = m.mod, 
                  model.y = full, 
                  treat = "x", 
                  mediator="m",
                  treat.value = x.treat,
                  control.value = x.ctrl,
                  sims = 5000)
)

summary(med2)

########################################
## Example 4: Multiple mediator model ##
########################################

## Example with 2 continuous mediators

# assign variables
y <- rbans.global.score
x <- bl.cr
m <- daily.sofa
m2 <- tot.benz

## single-model approach: fit full model
full <- lm(y ~ x + m + m2)

## estimate total indirect effect through set of mediators (m,m2)
x.treat <- 1
x.ctrl <- 0
v.xm <- vcov(full)["x",c("m","m2")]
v.mm <- vcov(full)[c("m","m2"), c("m","m2")]
v.mx <- vcov(full)[c("m","m2"), "x"]
beta.m <- c(full$coefficients["m"], full$coefficients["m2"])
emc <- -1 * v.xm %*% solve(v.mm) %*% beta.m
total.IDE <- emc %*% (x.treat - x.ctrl)
## check that this equals difference of coefficients
# sub<- lm(y ~ x)
# sub$coef["x"] - full$coefficients["x"] ## total IDE

## variance of total IDE
v.xm <- vcov(full)["x",c("m","m2")]
v.mm <- vcov(full)[c("m","m2"), c("m","m2")]
v.mx <- vcov(full)[c("m","m2"), "x"]
var.emc <- v.xm %*% solve(v.mm) %*% v.mx
sqrt(var.emc) # SE(EMC)

## specific IDE through m
(-1 * vcov(full)["x","m"] * solve(vcov(full)["m","m"]) * full$coefficients["m"])*(x.treat - x.ctrl)
## check that this equals difference of coefficients
# sub3<- lm(y ~ x + m2) # model excluding m
# sub3$coefficients["x"] - full$coefficients["x"] # IDE of M

## variance of IDE through m
vcov(full)["x","m"] * solve(vcov(full)["m","m"]) * vcov(full)["m","x"]

## specific IDE through m2
-1 * vcov(full)["x","m2"] * solve(vcov(full)["m2","m2"]) * full$coefficients["m2"]
## check that this equals difference of coefficients
# sub2<- lm(y ~ x + m) # model excluding m2
# sub2$coef["x"] - full$coefficients["x"] # IDE of M2

## variance of IDE through m2
vcov(full)["x","m2"] * solve(vcov(full)["m2","m2"]) * vcov(full)["m2","x"]

###################################################
## MacKinnon's parallel (aka Hayes' single-step) ##
## and VanderWeele's regression-based approach   ##
###################################################

## fit full model and a separate model for each mediator
full<- lm(y ~ x + m + m2)
pathA1<- lm(m ~ x)
pathA2<- lm(m2 ~ x)
# IDE through m estimated using ...
a1b1<- pathA1$coefficients["x"]*full$coefficients["m"]
a2b2<- pathA2$coefficients["x"]*full$coefficients["m2"]
a1b1 + a2b2 ## these sum to total IDE

# extract model coefficients
a1<- pathA1$coefficients["x"]
a2<- pathA2$coefficients["x"]
b1<- full$coefficients["m"]
b2<- full$coefficients["m2"]
# extract sample variances of coefficients
s2_a1<- vcov(pathA1)["x","x"]
s2_a2<- vcov(pathA2)["x","x"]
s2_b1<- vcov(full)["m","m"]
s2_b2<- vcov(full)["m2","m2"]
sb1b2<- vcov(full)["m","m2"] # cov between b1 and b2

## variance of IDE through M1: a1^2 * s2_b1 + b1^2 * s2_a1 (from MacKinnon's book pg 107)
(a1^2 * s2_b1) + (b1^2 * s2_a1)

## variance of IDE through M2: a2^2 * s2_b2 + b2^2 * s2_a2
(a2^2 * s2_b2) + (b2^2 * s2_a2)

## variance of total IDE through set M1 and M2:
(s2_a1*b1^2) + (s2_b1*a1^2) + (s2_a2 * b2^2) + (s2_b2 * a2^2) + (2*a1*a2*sb1b2)


###############################################
## Hayes' serial model / sequential approach ##
###############################################

# assuming X --> M --> M2 --> Y
full<- lm(y ~ x + m + m2)
# m2 depends on M and X
sub1<- lm(m2 ~ x + m)
# M depends on X
sub2<- lm(m ~ x)

## IDE of X through M to Y =  alpha1*beta1
alpha1<- sub2$coefficients["x"]
s2_alpha1<- vcov(sub2)["x","x"]

beta1<- full$coefficients["m"]
s2_beta1<- vcov(full)["m","m"]

alpha1*beta1

## IDE of X through M2 to Y = alpha2*beta2
alpha2<- sub1$coef["x"]
s2_alpha2<- vcov(sub1)["x","x"]

beta2<- full$coefficients["m2"]
s2_beta2<- vcov(full)["m2","m2"]

alpha2*beta2

## IDE of X through M1 to M2 to Y = alpha_1*d21*beta2
d21<- sub1$coefficients["m"]
s2_d21 <- vcov(sub1)["m","m"]
alpha1 * d21 * beta2

(alpha1*beta1) +  (alpha2*beta2) + (alpha1 * d21 * beta2) ## these sum to total IDE

## variance of alpha1*beta1 is a1^2 * s2_b1 + b1^2 * s2_a1
v_a1b1 <- (alpha1^2 * s2_beta1) + (beta1^2 * s2_alpha1)
v_a1b1
(alpha1*beta1) + c(-1,1)*qnorm(.975)*sqrt(v_a1b1)

## variance of alpha2*beta2 is a2^2 * s2_b2 + b2^2 * s2_a2
v_a2b2 <- (alpha2^2 * s2_beta2) + (beta2^2 * s2_alpha2)
sqrt(v_a2b2)
alpha2*beta2 + c(-1,1)*qnorm(.975)*sqrt(v_a2b2)

## variance of alpha1*d21*beta2:
v_a1_d21_b2<- (alpha1^2 * d21^2 * s2_beta2) + (alpha1^2 * beta2^2 * s2_d21) + (d21^2 * beta2^2 + s2_alpha1)
sqrt(v_a1_d21_b2)
alpha1*d21*beta2 + c(-1,1)*qnorm(.975)*sqrt(v_a1_d21_b2)

## switch temporal order of mediators: X --> M2 --> M --> Y
full<- lm(y ~ x + m2 + m)
# M depends on M2 and X
sub1<- lm(m ~ x + m2)
# M2 depends on X
sub2<- lm(m2 ~ x)

## IDE of X through M2 to Y: alpha1*beta1
sub2$coef["x"]*full$coefficients["m2"]
## IDE of X through M to Y: alpha2*beta2
sub1$coef["x"]*full$coefficients["m"]
## IDE of X through M2 to M1 to Y: alpha_1 *d21 *Beta2
sub2$coefficients["x"] * sub1$coef["m2"] * full$coefficients["m"]
# these specific IDEs again sum to total IDE
