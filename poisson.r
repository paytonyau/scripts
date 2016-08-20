## This script calculates the probability of >1 event occurring given different
## average event frequencies of occurrence (using a Poisson PDF).
## The application is to understand the impact of MOI (ratio of viruses:bacteria)
## on the probability of co-infection.

## First Graphic: a Poisson Bar Plot
# Variables
x_max <- 10
y_max <- 0.3
number_viruses <- 0:x_max
color_lines <- ifelse(number_viruses<2, 'black', 'maroon3')
viral_mois <- c(0.1,0.5,1:4) # six values please

# Make figure
pdf("MOI_coinf_all.pdf", 15, 10)
par(mfrow=c(2,3))

for (i in viral_mois) {
  plot(number_viruses, dpois(number_viruses,i), type='h', lwd=10, col=color_lines, ylim=c(0,y_max), xlab='', ylab='')
  text(x=x_max/2, y=y_max-0.01, labels=paste("MOI:",i), font=2, cex=2)
  # note that lower tail as False implies P[X > x] (see ?ppois)
  text(x=x_max*0.93, y=y_max/3, labels=paste(100*signif(ppois(1, i, lower.tail=F), digits=3), "%", sep=""), col="maroon3", cex=2)
  title(xlab="Phage per Bacterium", ylab="Probability")
  }

dev.off()

# The above makes an important assumption: the mean number of infections per
# bacterium is assumed to be equal to the MOI (and is therefore the mean (λ) of
# the Poisson distribution). In other words:
# - all bacteria are equally likely to be infected
# - every virus present will go on to infect a bacterium
#
# These assumptions are unlikely to be strictly true and remember that if a
# bacterium is infected by 2 or more phage, they are not necessarily different
# genotypes. This depends on the proportions of different phage as well.
#
# More precisely: if there are two types of phage (A and B) in equal proportion
# (50:50) and there are n phage infecting a bacterium then the probabilities are:
# - P(all type A) = (1/2) ^ n
# - P(all type B) = (1/2) ^ n
# - P(all type A or B) = 2 * (1/2) ^ n
# - P(not all one type) = 1 - ( 2 * (1/2) ^ n )
#
# = 0, 0.5, 0.75, 0.875, 0.9375... for n = 1, 2, 3, 4, 5
#
# The next graph multiplies these probabilities by the Poisson ones and sums
# (up to arbitrarily large values of n). Yes, I know some maths could probably
# be done but this works too.

## Second Graphic: a Bar Plot
# Variables
viral_mois <- 1:10 # re-declared
n <- 1:10000
diff_probs <- c()

# Maths for different infections
for (i in viral_mois) {
  diff_probs <- c( diff_probs, sum(dpois(n, i) * (1 - ( 2 * (1/2) ^ n ))) )
  }

# Make figure
pdf("MOI_coinf_2types.pdf", 7, 5)
barX <- barplot(diff_probs, ylim=c(0,1), names.arg=viral_mois)
title(main="Diverse Infections Assuming 50:50 Mix", xlab="Multiplicity of Infection", ylab="Proportion of Diverse Infections")
text(x=barX, y=diff_probs+0.03, paste(signif(100*diff_probs, 3), "%", sep=""), cex=0.8, xpd=TRUE)
dev.off()
