############################################
####  INTRODUCTION: GETTING STARTED     ####
############################################

## Preliminaries

x = c(109, 65, 22, 3, 1)
z <- c(109, 65, 22, 3, 1)

y = rpois(200, lambda=.61)
w <- rpois(200, lambda=.61)

## Basic Operations

temps = c(51.9, 51.8, 51.9, 53)

temps

temps - 32

(5/9) * (temps - 32)

CT = c(48, 48.2, 48, 48.7)
temps - CT

# Presidential Height in Centimeters
winner = c(185, 182, 182, 188, 188, 188, 185, 185, 177,
           182, 182, 193, 183, 179, 179, 175)
opponent = c(175, 193, 185, 187, 188, 173, 180, 177, 183,
             185, 180, 180, 182, 178, 178, 173)

length(winner)

year = seq(from=2008, to=1948, by=-4)

years = seq(2008, 1948, -4)

winner[4] = 189
winner[5] = 189

winner[4:5] = 190

winner

mean(winner)

mean(opponent)

difference = winner - opponent

data.frame(year, winner, opponent, difference)

taller.won = winner > opponent
taller.won

table(taller.won)

barplot(rev(difference), xlab="Election years 1948 to 2008",
        ylab="Height difference in cm")

plot(winner, opponent)

# Fatalities due to Prussian calvary horsekicks
k = c(0, 1, 2, 3, 4)
x = c(109, 65, 22, 3, 1)

barplot(x, names.arg=k)

p = x / sum(x)
p

r = sum(p * k)
r

v = sum(x * (k - r)^2) / 199
v

f = r^k * exp(- r) / factorial(k)
f

f = dpois(k, r)
f

floor(200*f) #expected counts

x

cbind(k, p, f)

## R Scripts

# Horsekicks.R
# Prussian horsekick data

k = c(0, 1, 2, 3, 4)
x = c(109, 65, 22, 3, 1)
p = x / sum(x) #relative frequencies
print(p)
r = sum(k * p) #mean
v = sum(x * (k - r)^2) / 199 #variance
print(r)
print(v)
f = dpois(k, r)
print(cbind(k, p, f))

getwd()
source("horsekicks.R")

y = rpois(200, lambda=.61)
kicks = table(y) #table of sample frequencies
kicks

Theoretical = dpois(0:3, lambda=.61)
Sample = kicks / 200
cbind(Theoretical, Sample)

mean(y)

var(y)

example(mean)

## Functions

var.n = function(x){
  v = var(x)
  n = NROW(x)
  v * (n - 1) / n
}

temps = c(51.9, 51.8, 51.9, 53)

var(temps)

var.n(temps)

# Many R functions require functions as arguments. For
# example, the integrate function (numerical integration)
# is one where the integrand must be supplied as an
# argument

# This function returns the integrand of the Beta function
# evaluated at a given point x. Arguments a and b
# specify the exponents:
f = function(x, a=1, b=1)
  x^(a-1) * (1-x)^(b-1)

# f function can be used to evaluate the integrand
# along a sequence of x values: here in a
# sequence from 0 to 1 with steps of 0.2
x = seq(0, 1, .2) 

# we call f with x as a vector, it returns a vector
# of same length as x
f(x, a=2, b=2)

# now we obtain numerical integration result for
# a = b = 2 by calling the R integrate function
# with function f as the first argument;
integrate(f, lower=0, upper=1, a=2, b=2)

# R has function beta to compute this integral
beta(2, 2)

# Here we graph a function using curve()
# We are going to graph the function 
# f(x) = (x^(a-1))((1-x)^(b-1)) for
# a = b = 2 which is (x*(1-x)) or the 
# integrand in our previous example
curve(x*(1-x), from=0, to=1, ylab="f(x)")

## Vectors and Matrices

## Class mobility example

# Vectors are one-dimensional and heterogeneous
probs = c(.45, .05, .01, .48, .70, .50, .07, .25, .49)

# A matrix is two-dimensional and heterogeneous
P = matrix(probs, nrow=3, ncol=3)

# values entered column by column
P

# can assign row, column names
# are identical here
rownames(P) <- colnames(P) <- c("lower", "middle", "upper")
P

# rows should and do sum to 1
rowSums(P)

# can use apply to compute at margin
apply(P, MARGIN=1, FUN=sum)

# transition probabilities at two generations
# are product of the matrix
P2 = P %*% P
P2

# to extract element values
P2[1, 3]

# just first row
P2[1, ]

P4 = P2 %*% P2
P8 = P4 %*% P4
P8

# here we compute a matrix of constants
# note are doing it "byrow"
Q = matrix(c(  0.45, 0.48, 0.07,
               0.05, 0.70, 0.25,
               0.01, 0.50, 0.49), nrow=3, ncol=3)
Q

## Dataframes

## First six records
head(USArrests)

# Sample size and dimension
NROW(USArrests)
dim(USArrests) #dimension

# Names of variables
names(USArrests)

# Structure of Data
str(USArrests)
is.data.frame(USArrests)

arrests = as.matrix(USArrests)
str(arrests)

# Missing values
any(is.na(USArrests))

# Working with dataframes
summary(USArrests)

# Extracting data from a dataframe
USArrests["California", "Murder"]
USArrests["California", ]

USArrests$Assault

# Histograms
# We previously saw that the distribution
# of assaults may be positively skewed because 
# the sample mean is larger than the sample 
# median. A histogram of the data helps to 
# visualize the shape of the distribution.
hist(USArrests$Assault)

# The skewness is easier to observe here:
library(MASS)
truehist(USArrests$Assault)

# The two histogram functions hist and 
# truehist have different default methods
# for determining the bin width, and the 
# truehist function by default produces
# a probability histogram. The optional 
# arguments of hist below match the defaults
# of truehist.
hist(USArrests$Assault, prob=TRUE, 
     breaks="scott")

# Attaching a dataframe

# Is we attach a dataframe, the variables can
# be referenced directly:
attach(USArrests)

# So can compute the percent of crimes
# that are murders using vectorized operations:
murder.pct = 100 * Murder / (Murder + Assault + Rape)
head(murder.pct)

# should detach() df when done with it

# alternative to attach is to use with().
# Useful for displaying plots or summary
# statistics, but variables created using
# with() are local to with(), ie the
# expr block in this case:
with(USArrests, expr={
  murder.pct = 100 * Murder / (Murder + Assault + Rape)
})
murder.pct

# Scatterplots and correlations

# With numeric data, is interesting to display
# scatterplots of different pairs of data to 
# look for possible relations between variables.
# plot() will display a single scatterplot
# because USArrests is still attached:
plot(UrbanPop, Murder)


# Does not show a strong relation between
# murders and percent population.

# pairs() function will display an array
# of scatterplots:
pairs(USArrests)

# Shows much more info about the data.

# cor() shows degree of linear association
# between two variables. Is Pearson corr coef:
cor(UrbanPop, Murder)

# cor() against the dataframe gives cor matrix:
# Confirms what we saw in scatterplot
cor(USArrests)

##### Importing Data

# Entering data manually

# A table of gas mileages on four new models
# of Japanese luxury cars is given. Do the
# models have the same gas mileage, on average?
y1 <- c(22, 26);y1    
y2 = c(28, 24, 29)
y3 = c(29, 32, 28);y3
y4 = c(23, 24);y4 
y = c(y1, y2, y3, y4);y
Model = c(rep("A", 2), rep("B", 3), 
          rep("C", 3), rep("D", 2))
Model
mileages = data.frame(y, Model)

str(mileages)
mileages

# Importing data from a text file
# often one needs to import text from 
# an external (ASCII) file. The data is
# often delimited by a space, tab or comma.

?read.table # note sep, header arguments

getwd()
file.exists("lunatics.txt")
lunatics = read.table("lunatics.txt", 
                      header=TRUE)

str(lunatics)

# lunatics small, so can simply use
# print() function to view result of
# read.table() call used to import the data.
lunatics

file.exists("college.txt")
dat = read.table("college.txt",
                 # \t is for tab delimited:
                 header=TRUE, sep="\t")

# Data available on the Internet. We use
# read.table() with complete URL and skip=60
# to read the data beginning on line 61
pidigits = read.table(
  "http://www.itl.nist.gov/div898/strd/univ/data/PiDigits.dat",
  skip=60)
head(pidigits)

# even though only one variable, it is a 
# dataframe

# are the digits of pi uniformly distributed?
# We summarize in a table:
table(pidigits)

# proportions easier to understand. This is
# another example of vectorized operations:
prop = table(pidigits) / 5000  #proportions
prop

# Recall that variance of a sample proportion
# is p(1-p)/n. So if true proportion is 0.1
# for every digit, then the standard error is:
sqrt(.1 * .9 / 5000)

# But if true proportions are unknown, we use
# sample estimates of proportions, gives us
# slightly different estimates for se.

# We replace constants 0.1 and 0.9 by vectors
# of length 10, so result is also a vector of
# length 10 rather than a scalar:
se.hat = sqrt(prop * (1-prop) / 5000)

# We use rbind() to collect results together
# into a matrix. We round result and include
# estimate of standard error in display:
round(rbind(prop, se.hat, prop-2*se.hat, 
            prop+2*se.hat), 4)

# Do any of the sample proportions fall
# outside of interval 0.1 + or - 2(se)? No.

# This helps to visualize with a horizontal
# reference line at 0.1
barplot(prop, xlab="digit", ylab="proportion")
abline(h = .1)

# Packages

# R functions are grouped into packages.
# A number of recommended packages are 
# included in the R distribution.
# examples are boot (bootstrap), MASS, and 
# lattice (graphics). In addition, thousands 
# of contributed packages are available to
# install. Type the command
library()

# to show installed packages.
# Suppose that we are interested
# in the examples related to the "law" data 
# in the bootstrap package. Typing 'law'
# at the prompt produces the following error 
# because no object named "law" is found
law

# get a similar warning with
data(law)

# must install bootstrap and then
install.packages("bootstrap")

# can query data in bootstrap:
data(package="bootstrap")

# to use objects in a package, 
# must first load it:
library(bootstrap)

# another useful feature of library
# function is to get help for a function.
# Either of these will display
# a summary of the package:
library(help=bootstrap)
help(package=bootstrap)

# now law is loaded and accessible
# in the R workspace:
law

# suppose we want to compute sample means
# or the correlation between LSAT and GPA:
mean(law)

# R wants this:
colMeans(law)

# or this:
sapply(law,mean)

# here is the correlation:
cor(law$LSAT, law$GPA)

?cor
example(cor)