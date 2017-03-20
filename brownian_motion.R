#adapted from a phytools.org tutorial
library(phytools)


#### part 1: brownian motion of a trait over time, following a single lineage

t <- 0:100  # length of time for simulation. Can be thought of as generations.
sig2 <- 0.01 # this is the instantaneous rate of change. This determines the relative size of jumps that tend to occur. when sig2 (sigma squared) is small, character state changes tend to be smaller. when its big (maximum of 1), individual state changes tend to be larger.
## first, simulate a set of random deviates. in other words, this is the series of character state changes through time (t)
x <- rnorm(n = length(t) - 1, sd = sqrt(sig2))
x #look at the values of x to get an idea of this. the character state starts at 0 then changes by the amount (x). Imagine this a size of some trait, like leaf length. when x is negative, it means it got smaller, and when positive, it got bigger
## now compute their cumulative sum
x <- c(0, cumsum(x)) #we want cumulative sum because the state a character is in at any time (t) is the sum of all its past transitions. ie, it started at "0", got bigger by this much, smaller by that much, bigger again, etc. you add those all up to see the state at teh end! 
plot(t, x, type = "l", ylim = c(-2, 2))

#Repeat these simulations several times. Then, repeat some more, changing t and sig2, one at a time. 

#1) How does the size of t affect how close the final character state (phenotype) is to the initial character state (which is always 0)? 

#2) How does the size of sig2 affect how close the final character state (phenotype) is to the initial character state? 


#### part 1: brownian motion of a trait over time, following multiple single lineages

t <- 0:100
sig2 <- 0.01
nsim <- 100
X <- matrix(rnorm(n = nsim * (length(t) - 1), sd = sqrt(sig2)), nsim, length(t) - 
              1)
X <- cbind(rep(0, nsim), t(apply(X, 1, cumsum)))
plot(t, X[1, ], xlab = "time", ylab = "phenotype", ylim = c(-2, 2), type = "l")
apply(X[2:nsim, ], 1, function(x, t) lines(t, x), t = t)

#You should have seen that this block of code does thes same thing we did before but plotting several independent lineages all at once. Try a couple of times, again changing sig2 and t. You might need to change the y axis to fit your character state range! (change the numbers in 'ylim = c(-2, 2)' to whatever you want. 

#3) How does this simulation exercise affect your answers to questions 1 and 2, above? Which simulation is more effective at demonstrating the effects of sig2 and t? Why?




#carrie: below is from phytools.org. the first part is just changing sig2 but i already included that! look through it to see what you wanna use. i think it would be good to get this formatted in an juppyter or markdown notebook!   













X <- matrix(rnorm(n = nsim * (length(t) - 1), sd = sqrt(sig2/10)), nsim, length(t) - 
              1)
X <- cbind(rep(0, nsim), t(apply(X, 1, cumsum)))
plot(t, X[1, ], xlab = "time", ylab = "phenotype", ylim = c(-2, 2), type = "l")
apply(X[2:nsim, ], 1, function(x, t) lines(t, x), t = t)

X <- matrix(0, nsim, length(t))
for (i in 1:nsim) X[i, ] <- c(0, cumsum(rnorm(n = length(t) - 1, sd = sqrt(sig2))))
plot(t, X[1, ], xlab = "time", ylab = "phenotype", ylim = c(-2, 2), type = "l")
for (i in 1:nsim) lines(t, X[i, ])

nsim <- 10000
X <- matrix(rnorm(n = nsim * (length(t) - 1), sd = sqrt(sig2)), nsim, length(t) - 
              1)
X <- cbind(rep(0, nsim), t(apply(X, 1, cumsum)))
v <- apply(X, 2, var)
plot(t, v, type = "l", xlab = "time", ylab = "variance among simulations")
