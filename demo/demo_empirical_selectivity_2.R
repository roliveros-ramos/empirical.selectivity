library(fields)
library(mgcv)
library(r4ss)
packageVersion("r4ss")
#[1] '1.43.0'

# remotes::install_github("roliveros-ramos/fks")
# remotes::install_github("roliveros-ramos/r4ss.selectivity")
library(fks)
library(empirical.selectivity)
packageVersion("empirical.selectivity")
#[1] '0.1.5'


# Computing Empirical Selectivity -----------------------------------------

# path = "../../SS_models/BigEye_2020/" # path to SS outputs
# rep = SS_output2(path)
#rep = YellowtailRockfish_2017 # data demo 1
#rep = BigSkate_2019           # data demo 2
rep = BigEye_2020             # data demo 3
es = empirical_selectivity(rep, fleet=20, by="length")
plot(es)

es = empirical_selectivity(rep, fleet=1, by="length")
# es = empirical_selectivity(rep, fleet=1, by="age")#Yellowtail

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
#default is  FUN = function(x, y) (log(x) - log(y))^2
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
lines(ss_27a, col="green", lwd=2) #same results as expected

summary(ss_01)
summary(ss_24)
summary(ss_27)

# Fitting multiple models (i.e. selectivity patterns)
ss_mult = fit_selectivity(es2, pattern=c(1, 24, 27),  k=5, w="catch")
plot(ss_mult)
summary(ss_mult)

ss_mult = fit_selectivity(es2, pattern=c(1, 24, 27),  k=10, w="catch")
plot(ss_mult)
summary(ss_mult)

#compare different number of knots
ss_mult = fit_selectivity(es2, pattern=27,  k=c(3,10), w="catch")
#ss_mult
plot(ss_mult) # from 3 to 10 knots are tried, only the best fitting curve is shown
summary(ss_mult)
SS_writeselec(ss_mult,phase=2)


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
#pattern 27: splines with k=3, k=4 and k5 + pattern 1 and 24 are compared
ss_mult = fit_selectivity(es2, pattern=c(1, 24, 27), blocks=3, k=5, w="catch")
plot(ss_mult)
summary(ss_mult)



# optimization of blocks positions (method="cluster")
# kmeans with number of centers = number of blocks
# as it can be seen in
# r4ss.selectivity:::.clusterBlocks
# blocks indicates the maximum number of blocks to fit
ss_block2 = fit_selectivity(es2, pattern=27, blocks = 6, w="catch",
                            method="cluster", FUN=lnorm,
                            control=list(min_block_size=10))
#method = NULL (is equally spaced)
# Not implemented yet, future work:
#if (method == "optim")
#  breaks = .optimBlocks(object, breaks, w)

plot(ss_block2)
summary(ss_block2)
#warning sometimes it could fail, try lowering the initial size
es3 = es[length=c(55,190)] # extract sizes within [60, 190], age=c(2,20) for age selex
ss_block3 = fit_selectivity(es3, pattern=27, blocks = 6, w="catch",
                            method="cluster", FUN=lnorm,
                            control=list(min_block_size=10))
#
plot(ss_block3)
summary(ss_block3)

# manually work on the first block
es4 = es[yr=17:48, length=c(60,190)] # extract years 50 to 160
# extract sizes within [60, 190], age=c(2,20) for age selex
ss_first_block = fit_selectivity(es4, pattern=27,  k=c(3,10), w="catch")
plot(ss_first_block)
summary(ss_first_block)

# write the config for the ctl file
SS_writeselec(ss_01)
SS_writeselec(ss_24)
SS_writeselec(ss_27)
SS_writeselec(ss_mult, t=1) # config for block 1
SS_writeselec(ss_block0, t=3) # config for block 3
SS_writeselec(ss_first_block, t=1) # config for block 3


# Using an external matrix ------------------------------------------------
# how to use it with an external N matrix if you have your own model

rep = oneFleet  # external matrix
ww = catch # catch weighted

years = seq(from=1980, length=nrow(rep), by=0.25) # quarter
marks = seq(from=5, length=ncol(rep))
es = empirical_selectivity(rep, fleet="LL", by="length",
                           years=years, bins=marks)

plot(es)

# subsetting lengths
es1 = es[yr=2000:2015] # extract years 50 to 160
es2 = es[length=c(30,120)] # extract sizes within [60, 190], age=c(2,20) for age selex
plot(es1)
plot(es2)

dat = weighted.mean(es2, w=ww)
plot(dat)

lnorm  = function(x, y) -dlnorm(y, mean=log(x), sd=sd(log(x/y)),log = TRUE)
ss_24 = fit_selectivity(dat, pattern=24, FUN = lnorm)

plot(ss_24)
lines(ss_01a, col=1, lwd=2)
lines(ss_01b, col=2, lwd=2)
lines(ss_01c, col=3, lwd=2)
lines(ss_01d, col=4, lwd=2)



