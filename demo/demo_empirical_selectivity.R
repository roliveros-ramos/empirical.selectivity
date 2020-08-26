library(r4ss)
library(fields)
library(mgcv)
library(fks)
library(r4ss.selectivity)


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
# the number k includes the 2 external nodes, k=5 would be
# equivalent to k=3 in the freeknotspline package
ss_27 = fit_selectivity(dat, pattern=27, k=5)
ss_24 = fit_selectivity(dat, pattern=24)
ss_01 = fit_selectivity(dat, pattern=1)

lines(ss_27, col="red", lwd=2)
lines(ss_24, col="black", lwd=2)
lines(ss_01, col="green", lwd=2)

# Using blocks

# User defined blocks
ss_block0 = fit_selectivity(es2, pattern=24, blocks=c(17, 100, 130, 180),
                            w="catch")
plot(ss_block0)

# equidistant split of time in n blocks
ss_block1 = fit_selectivity(es2, pattern=24, blocks = 3, w="catch")
plot(ss_block1)

# Fitting multiple models (i.e. selectivity patterns)
ss_mult = fit_selectivity(es2, pattern=c(1, 24, 27), blocks=3, k=5, w="catch")

# optimization of blocks positions (method="optim")
# if block positions are provided, are used as start search point
ss_block2 = fit_selectivity(es2, pattern=24, blocks = 3, w="catch", method="optim")
plot(ss_block2)

# write the config for the ctl file
SS_writeselec(ss_01)
SS_writeselec(ss_24)
SS_writeselec(ss_27)
SS_writeselec(ss_mult, t=1) # config for block 1
SS_writeselec(ss_block0, t=3) # config for block 3

