## Installation
Install this packge via `remotes::install_github("pchiroque/betabart")`.

## Usage 
The response `y` must be a `data.frame` with columns `y1`,`yb` and `y0`. 
This can be acquired by applying the function `prepare.response` on a vector
with $y\in[0,1]$.

The covariates must be a ´data.frame` and have the same numbers of the rows as the response.

### Example 

```R
set.seed(4)
y <- sample(c(seq(0.2,0.8,length=6),rep(0,5),rep(1,4)))

y.test <- prepare.response(y)

sampleX <- data.frame(age=seq(20,30,length=15),height=seq(150,170,length=15))

#### Fitting ####
time <- Sys.time()
fit <- beta_bart(x = sampleX, y = y.test,
                 mcmc = list(burnin = 1, sample = 5))
print(duration <- Sys.time() - time)

```




