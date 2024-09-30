#' @export

betabart.predict <- function(fit){

y1index <- grep("y1",names(fit$fit))
ybindex <- grep("yb",names(fit$fit))
y0index <- grep("y0",names(fit$fit))

# probability of each case to be one
y1.predict <- map_df(fit$fit[y1index[-1]],pnorm)%>%
  map_df(function(x) sapply(x,function(p)rbinom(1,1,p)))

# The expected value of the yb
yb.predict <- map_df(fit$fit[ybindex],function(a)sapply(a,function(a)rbeta(1,exp(a),1)))

# probability of each case to be zero
y0.predict <- map_df(fit$fit[y0index],pnorm)%>%
  map_df(function(x) sapply(x,function(p)rbinom(1,1,p)))

y.predict <- do.call(cbind,lapply(1:length(y),function(i)y1.predict[[i]]+(y1.predict[[i]]==0)*(
  ifelse(y0.predict[[i]]==0,yb.predict[[i]],0))))
y.predict
}

