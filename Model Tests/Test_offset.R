
## 0. Assumption:
##  A count model whose outcomes are constrained through predictors.
##  In the real case the possible outcomes are constrained through the rates, the priors, and the prior initial state.

n <- 500
x <- runif(-1, 1)
y_log <- rnorm(n, x * 2, 1) ## positive response with a few values close to zero
hist(exp(y_log))

## Simulate and fit a model for one pair of offset_true and offset_zero
simulateOffset <- function(o_t = offset_true[1], o_0 = offset_zero[1], y_l = y_log, pred = x) {
  
  count <- rpois(y, exp(y_log * o_t))
  
  is0 <- count == 0
  
  f <- glm(count ~ pred + 0, offset = ifelse(is0, log(o_0), log(o_t)), family = poisson(link = "log"))
  return(f)
}


## Variable pairs of offsets
offset_true <- rep(seq(0.1, 1, by = 0.2), 200) # length 20 * 5 = 100
l <- length(offset_true)
l_split <- l/2
offset_zero <- c(offset_true[1:l_split], runif(l_split, 0, 1))
diff_offset <- offset_zero - offset_true ## becomes negative, when offset_zero is smaller that the true offset
isequeal_offset <- offset_zero == offset_true
unequalisred <- (-isequeal_offset) + 2

## Do the fits
fits <- mapply(FUN = simulateOffset, offset_true, offset_zero, SIMPLIFY = F)


## Question 1:
## - How does setting the offset for zero bias the estimate?
##   1.1. Are estimates for populations with more zeroes biased in a direction?
##   1.2 How does the relative difference of the false to the true offset bias the estimate?
##   1.3 What increases uncertainty of the estimate?


## 1.1. The estimate is mildly affected by a false offset for zero: estimate is spread out more. There are also more smaller, estimates, but this is due to the wider spread.
## (run a few times)
beta <- sapply(fits, function(f) f$coefficients["pred"])
plot(density((beta[isequeal_offset]))) ## 100 estimates with the correct offset
lines(density(beta[!isequeal_offset]), col = 2) ## 100 estimates with a wrong offset for zeroes

table(beta[isequeal_offset] > beta[!isequeal_offset])
median(beta[isequeal_offset])
median(beta[!isequeal_offset])

## 1.2 The estimate is systematically biased: offset smaller, greater estimate
plot(beta ~ offset_zero, col = unequalisred) # black are the esimates where zero has a true offset
plot(beta ~ diff_offset, col = unequalisred) # black are the esimates where zero has a true offset


## 1.3 
## 1. The estimate is more affected (in terms of uncertainty) by the size of the offset
plot(beta ~ offset_zero, col = unequalisred) # black are the esimates where zero has a true offset

## The uncertainty of the estimate is inflated by the number of zeroes, but this is seemingly unaffected by whether zeroes have the true offset
Count <- sapply(fits, function(f) f$y)
n_zeroes <- apply(Count, 2, function(x) sum(x == 0))
plot(beta ~ n_zeroes, col = unequalisred)


## Question 2:
## - How does setting the wrong offset for zero bias the predictions?
##  2.1 The predictions with zero observations are obviously biased.
##  2.2 Are the predictions that are non-zero observations biased?
##  2.3 How does the relative difference of the false to the true offset bias the predictions at zeroes?
##  2.4 How does the relative difference of the false to the true offset bias the predictions at non-zeroes?


## Predictions
Predictions <- mapply(function(f, o) predict(f), fits, offset_true)
Residuals <- apply(Predictions, 2, function(p) y_log - p)

matplot(y_log, Residuals[,1:5])
matplot(y_log, Residuals[,l_split:(l_split+5)])
matplot(y_log, Predictions[,1:5])


