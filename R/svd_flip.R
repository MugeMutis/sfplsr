svd_flip <- function(x,y){
  abs_idx <- which.max(abs(x))
  sgn <- sign((x[abs_idx]))
  x <- x*sgn
  y <- y*sgn
  
  return(list(x = x, y = y))
}
