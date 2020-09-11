library(r4ss)
library(fields)
library(mgcv)
library(fks)
library(r4ss.selectivity)

# remotes::install_github("roliveros-ramos/fks")
# remotes::install_github("roliveros-ramos/r4ss.selectivity")

# Computing Empirical Selectivity -----------------------------------------

# path = "../../SS_models/BigEye_2020/"
# rep = SS_output2(path)
# rep = YellowtailRockfish_2017
# rep = BigSkate_2019
rep = BigEye_2020
es = empirical_selectivity(rep, fleet=2, by="length")
plot(es)
# subseting length/age use "range", yr the argument as usual
es1 = es[yr=50:160] # extract years 50 to 160
es2 = es[length=c(60,190)] # extract sizes within [60, 190], age=c(2,20) for age selex
plot(es1)
plot(es2)

# Model Fitting -----------------------------------------------------------

# One block
dat = weighted.mean(es2, w="catch")
plot(dat)
# the number k includes the 2 external knots, k=5 would be
# equivalent to k=3 in the freeknotspline package

d_mul = function(x, y) -dnorm(y, mean=x, sd=x*(1-x), log = TRUE)
lnorm = function(x, y) -dlnorm(y, mean=log(x), sd=sd(log(x/y)),log = TRUE)
lnorm2 = function(x, y) -dlnorm(y, mean=log(x), log = TRUE)

patt = 24
ss_01a = fit_selectivity(dat, pattern=patt)
ss_01b = fit_selectivity(dat, pattern=patt, FUN=d_mul)
ss_01c = fit_selectivity(dat, pattern=patt, FUN=lnorm)
ss_01d = fit_selectivity(dat, pattern=patt, FUN=lnorm2)

plot(dat)
lines(ss_01a, col=1, lwd=2)
lines(ss_01b, col=2, lwd=2)
lines(ss_01c, col=3, lwd=2)
lines(ss_01d, col=4, lwd=2)

ss_24 = fit_selectivity(dat, pattern=24)
ss_27 = fit_selectivity(dat, pattern=27, k=5)
ss_27a = fit_selectivity(dat, pattern=27, k=5, control=list(optimizer="golden"))


lines(ss_24, col="black", lwd=2)
lines(ss_27, col="red", lwd=2)

summary(ss_01)
summary(ss_24)
summary(ss_27)

# Using blocks

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
ss_block2 = fit_selectivity(es2, pattern=24, blocks = 6, w="catch",
                            method="cluster", FUN=lnorm,
                            control=list(min_block_size=20))
plot(ss_block2)
summary(ss_block2)

clus = attr(attr(ss_block2, "breaks"), "cluster")
par(mfrow=c(2,1), mar=c(3,3,1,1))
plot(clus$year, clus$cluster, col=clus$cluster, pch=19, cex=0.7, las=1)
plot(clus$year, clus$block, col=clus$cluster, pch=19, cex=0.7, las=1)

# write the config for the ctl file
SS_writeselec(ss_01)
SS_writeselec(ss_24)
SS_writeselec(ss_27)
SS_writeselec(ss_mult, t=1) # config for block 1
SS_writeselec(ss_block0, t=3) # config for block 3

x = 1:100
y = cumsum(runif(100))

mod = fks(x=x, y=y, k=5)

# fixing the seed (reproducibility)
# blocks is "maximum number of clusters"
# fitting issues related to max empirical selectivity value
