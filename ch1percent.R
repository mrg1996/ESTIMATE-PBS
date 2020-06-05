ch1percent <- function(dt){
  require(dplyr)
  a <- round(nrow(dt) * 0.001,0)
  dt1 <- dt %>% arrange(desc(dt$PBS))
  dt2 <- dt1[1:a,]
  return(dt2)
}