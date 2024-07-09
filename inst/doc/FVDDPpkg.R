## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
set.seed(2023)

## ----setup, echo = T, results = 'hide', message = F---------------------------
library(FVDDPpkg)

## -----------------------------------------------------------------------------
FVDDP = initialize(theta = 1.28, sampling.f = function(x) rpois(x, 5),
                   density.f = function(x) dpois(x, 5), TRUE)

## -----------------------------------------------------------------------------
FVDDP

## ----results ='hide'----------------------------------------------------------
FVDDP = update(fvddp = FVDDP, y.new = c(7, 4, 9, 7))

## -----------------------------------------------------------------------------
FVDDP

## ----results='hide'-----------------------------------------------------------
FVDDP = propagate(fvddp = FVDDP, delta.t = 0.6)

## -----------------------------------------------------------------------------
FVDDP

## ----hide=T-------------------------------------------------------------------
FVDDP = update(fvddp = FVDDP, y.new = c(4, 7, 7, 10, 10, 5))

## -----------------------------------------------------------------------------
FVDDP

## ----results ='hide'----------------------------------------------------------
FVDDP_NONATOMIC = initialize(theta = 0.7, sampling.f = function(x) rbeta(x, 4, 7),
                          density.f = function(x) dbeta(x, 4, 7), atomic = FALSE)
FVDDP_PAST_NONATOMIC = update(fvddp = FVDDP_NONATOMIC, y.new = c(0.210, 0.635, .541))
FVDDP_FUTURE_NONATOMIC = update(fvddp = FVDDP_NONATOMIC, y.new = c(0.210))
FVDDP_FUTURE_NONATOMIC = propagate(fvddp = FVDDP_FUTURE_NONATOMIC, delta.t = 0.4)
FVDDP_FUTURE_NONATOMIC = update(fvddp = FVDDP_FUTURE_NONATOMIC, y.new = c(.635))

## ----results='hide'-----------------------------------------------------------
FVDDP_SMOOTH_NONATOMIC = smooth(fvddp.past = FVDDP_PAST_NONATOMIC, fvddp.future = FVDDP_FUTURE_NONATOMIC,
                                t.past = 0.75, t.future = 0.3, y.new = c(0.210, 0.635, 0.479))

## -----------------------------------------------------------------------------
FVDDP_SMOOTH_NONATOMIC

## ----results ='hide'----------------------------------------------------------
FVDDP_ATOMIC = initialize(theta = 0.7, sampling.f = function(x) rbeta(x, 10, 0.6),
                          density.f = function(x) dbinom(x, 10, 0.6), atomic = TRUE)
FVDDP_PAST_ATOMIC = update(fvddp = FVDDP_ATOMIC, y.new = c(2, 6, 5))
FVDDP_FUTURE_ATOMIC = update(fvddp = FVDDP_ATOMIC, y.new = c(2))
FVDDP_FUTURE_ATOMIC = propagate(fvddp = FVDDP_FUTURE_ATOMIC, delta.t = 0.4)
FVDDP_FUTURE_ATOMIC = update(fvddp = FVDDP_FUTURE_ATOMIC, y.new = c(6))

## ----results='hide'-----------------------------------------------------------
FVDDP_SMOOTH_ATOMIC = smooth(fvddp.past = FVDDP_PAST_ATOMIC, fvddp.future = FVDDP_FUTURE_ATOMIC,
                             t.past = 0.75, t.future = 0.3, y.new = c(2, 6, 4))

## -----------------------------------------------------------------------------
FVDDP_SMOOTH_ATOMIC

## ----results = 'hide'---------------------------------------------------------
FVDDP =initialize(theta = 3, sampling.f= function(x) rnorm(x, -1, 3),
                  density.f = function(x) dnorm(x, -1, 3), atomic = FALSE)
FVDDP = update(fvddp = FVDDP, y.new = c(-1.145, 0.553, 0.553, 0.553))

## ----results='hide'-----------------------------------------------------------
FVDDP_APPR_PROP = approx.propagate(fvddp = FVDDP, delta.t = 0.45, N = 20000)

## -----------------------------------------------------------------------------
FVDDP_APPR_PROP

## ----results='hide'-----------------------------------------------------------
FVDDP_EXACT_PROP = propagate(fvddp = FVDDP, delta.t = 0.45)

## -----------------------------------------------------------------------------
error.estimate(FVDDP_EXACT_PROP, FVDDP_APPR_PROP)

## ----results='hide'-----------------------------------------------------------
FVDDP_SMOOTH_APPR = approx.smooth(fvddp.past = FVDDP_PAST_ATOMIC, fvddp.future = FVDDP_FUTURE_ATOMIC,
                             t.past = 0.75, t.future = 0.3, y.new = c(2, 6, 4), N = 50000)

## -----------------------------------------------------------------------------
FVDDP_SMOOTH_APPR

## -----------------------------------------------------------------------------
error.estimate(FVDDP_SMOOTH_ATOMIC, FVDDP_SMOOTH_APPR)

## -----------------------------------------------------------------------------
PRUNED = prune(fvddp = FVDDP_SMOOTH_ATOMIC, eps = 1e-02)

## -----------------------------------------------------------------------------
PRUNED

## ----results = 'hide'---------------------------------------------------------
y = posterior.sample(fvddp = FVDDP_EXACT_PROP, N = 100)

## -----------------------------------------------------------------------------
table(round(y, 3))

## ----results = 'hide'---------------------------------------------------------
y = predictive.struct(fvddp = FVDDP_EXACT_PROP, N = 100)

## -----------------------------------------------------------------------------
table(round(y, 3))

