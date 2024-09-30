#' @export
prepare.response <- function(y){
  y1 <- 1*(y==1)
  yb <- ifelse(!y%in%c(0,1),y,NA)
  y0 <- ifelse(y==1,NA,1*(y==0))
  data.frame(y1,yb,y0)
}
