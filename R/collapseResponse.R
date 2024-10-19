#' @export
collapse.response <- function(y){
  do.call(cbind,lapply(1:nrow(y),function(i)ifelse(y$y1[[i]]==0,(y$y1[[i]]==0)*(
    ifelse(y$y0[[i]]==0,y$yb[[i]],0)),1)))
}
