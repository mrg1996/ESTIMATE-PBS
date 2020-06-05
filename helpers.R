applyWindow <- function(dt){
  f <- c()
  for (m in 1:22){
    dtemp <- dt[dt$CHR == m,]
    a <- 0
    while (a < max(dtemp$POS)){
      b <- a+ 100000
      dtf1<- dtemp[dtemp$POS >= a & dtemp$POS <= b,]
      if (nrow(dtf1) != 0){
        f <- c(f,dtf1$POS[dtf1$PBS == max(dtf1$PBS)])
      }
      a <- a + 100000
    }
  }
  dtf <- dt[dt$POS %in% f,]
  return(dtf)
} 


