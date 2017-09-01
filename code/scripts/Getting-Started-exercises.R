############################################
####  INTRODUCTION: GETTING STARTED     ####
####           EXERCISES                ####
############################################

#1 Normal Percentiles
# The qnorm function returns the percentiles
# (quantiles) of a normal distribution. Use the 
# qnorm function to find the 95th percentile of
# the standard normal distribution. Then, use 
# the qnorm function to find the quartiles of 
# the standard normal distribution (the quartiles 
# are the 25th, 50th, and 75th percentiles).
# Note you can use c(.25, .5, .75) as the first
# argument to qnorm().

#2 Chi-square Density Curve
# Use the curve function to display the graph
# of the CHI-SQ(1) density. The chi-square density 
# function is dchisq().

#3 Gamma Densities 
# Use the curve function to display the graph 
# of the gamma density with shape parameter 1 
# and rate parameter 1. The gamma density 
# function is dgamma. Consult the help file
# ?dgamma to see how to specify the
# parameters.

# Then use the curve function with add=TRUE
# to display the graphs of the gamma density
# with shape parameter k and rate 1 
# for k = 2,3, all in same graphics window.

#4 Binomial Probabilities
# Let X by the number of "ones" obtained in 
# 12 rolls of a fair die. Then X has a 
# Binomial(n=12,p=1/3) distribution. Compute
# a table of binomial probabilities for x =
# 0,1,...,12 by two methods: (1) See binomial
# probability density formula and use 0:12 for
# sequence of x values and choose() function to
# compute the binomial coefficients in the formula;
# (2) use dbinom function. Compare results

#5 Binomial CDF
# Let X by the number of "ones" obtained in 
# 12 rolls of a fair die. Then X has a 
# Binomial(n=12,p=1/3) distribution. Compute
# a table of cumulative binomial probabilities
# (the CDF) for x = 0,1,...,12 by two methods: 
# (1) Using cumsum() and first result from
# previous exercise; (2) use pbinom function.
# Compare results. What is P(X >7)?

#6 President's Heights
# Refer to the data about the heights of the
# United States Presidents as compared with their 
# main opponent in the presidential election.
# Create a scatterplot of the loser's height vs 
# the winner's height using the plot() function.

# Presidential Height in Centimeters
winner = c(185, 182, 182, 188, 188, 188, 185, 185, 
           177, 182, 182, 193, 183, 179, 179, 175)
opponent=c(175, 193, 185, 187, 188, 173, 180, 177,
           183, 185, 180, 180, 182, 178, 178, 173)
