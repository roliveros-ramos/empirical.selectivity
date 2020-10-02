library(fields)
library(mgcv)
library(r4ss)

# remotes::install_github("roliveros-ramos/fks")
# remotes::install_github("roliveros-ramos/r4ss.selectivity")

library(fks)
library(r4ss.selectivity)

# Computing Empirical Selectivity -----------------------------------------

# path = "../../SS_models/BigEye_2020/" # path to SS outputs
# rep = SS_output2(path)
# rep = YellowtailRockfish_2017 # data demo 1
# rep = BigSkate_2019           # data demo 2
rep = BigEye_2020             # data demo 3
es = empirical_selectivity(rep, fleet=2, by="length")
plot(es)

# subsetting length/age use "range", yr is an indexation argument as usual
es1 = es[yr=50:160] # extract years 50 to 160
es2 = es[length=c(60,190)] # extract sizes within [60, 190], age=c(2,20) for age selex
plot(es1)
plot(es2)

# Model Fitting -----------------------------------------------------------

# One block
dat = weighted.mean(es2, w="catch") # weight empirical selectivity by catch
# dat = weighted.mean(es, w="equal") # weight empirical selectivity by catch
plot(dat)

# error functions using for fitting the empirical selectivity to a model
d_mul  = function(x, y) -dnorm(y, mean=x, sd=x*(1-x), log = TRUE)
lnorm  = function(x, y) -dlnorm(y, mean=log(x), sd=sd(log(x/y)),log = TRUE)
lnorm2 = function(x, y) -dlnorm(y, mean=log(x), log = TRUE)

patt = 1 # 2-parameter logistic
ss_01a = fit_selectivity(dat, pattern=patt)
ss_01b = fit_selectivity(dat, pattern=patt, FUN=d_mul)
ss_01c = fit_selectivity(dat, pattern=patt, FUN=lnorm)
ss_01d = fit_selectivity(dat, pattern=patt, FUN=lnorm2)

plot(dat)
lines(ss_01a, col=1, lwd=2)
lines(ss_01b, col=2, lwd=2)
lines(ss_01c, col=3, lwd=2)
lines(ss_01d, col=4, lwd=2)

ss_01 = ss_01a

patt = 24 # 6-parameter double normal
ss_24a = fit_selectivity(dat, pattern=patt)
ss_24b = fit_selectivity(dat, pattern=patt, FUN=d_mul)
ss_24c = fit_selectivity(dat, pattern=patt, FUN=lnorm)
ss_24d = fit_selectivity(dat, pattern=patt, FUN=lnorm2)

plot(dat)
lines(ss_24a, col=1, lwd=2)
lines(ss_24b, col=2, lwd=2)
lines(ss_24c, col=3, lwd=2)
lines(ss_24d, col=4, lwd=2)

ss_24 = ss_24a

patt = 27 # splines with knots as parameters
# the number k includes the 2 external knots, k=5 would be
# equivalent to k=3 in the freeknotspline package
ss_27 = fit_selectivity(dat, pattern=patt, k=5) # default, use genetic algorithm
ss_27a = fit_selectivity(dat, pattern=patt, k=5, control=list(optimizer="golden"))

plot(dat)
lines(ss_01, col="black", lwd=2)
lines(ss_24, col="blue", lwd=2)
lines(ss_27, col="red", lwd=2)

summary(ss_01)
summary(ss_24)
summary(ss_27)

## Using blocks

# User defined blocks
ss_block0 = fit_selectivity(es2, pattern=24, blocks=c(17, 61, 180),
                            w="catch")
plot(ss_block0)
summary(ss_block0)

# equidistant split of time in n blocks
ss_block1 = fit_selectivity(es2, pattern=24, blocks = 3, w="catch")
plot(ss_block1)
summary(ss_block1)

# Fitting multiple models (i.e. selectivity patterns)
ss_mult = fit_selectivity(es2, pattern=c(1, 24, 27), blocks=3, k=5, w="catch")
plot(ss_mult)
summary(ss_mult)

# optimization of blocks positions (method="cluster")
# blocks indicates the maximum number of blocks to fit
ss_block2 = fit_selectivity(es2, pattern=27, blocks = 6, w="catch",
                            method="cluster", FUN=lnorm,
                            control=list(min_block_size=10))
plot(ss_block2)
summary(ss_block2)

# write the config for the ctl file
SS_writeselec(ss_01)
SS_writeselec(ss_24)
SS_writeselec(ss_27)
SS_writeselec(ss_mult, t=1) # config for block 1
SS_writeselec(ss_block0, t=3) # config for block 3

