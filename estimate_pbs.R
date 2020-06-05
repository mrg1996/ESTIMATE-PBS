## Estimate the PBS values 
require(ggplot2)
source("helpers.R")
source("ch1percent.R")
source("data2014.R")
require(ggpubr)
library("png")
library("grid")
library("gridExtra")
# List files

fl <- c("/home/vu/test/nenet/CEU-NENET.fst", 
        "/home/vu/test/nenet/VIET-NENET.fst",
        "/home/vu/test/nenet/CEU-VIET.fst")
# Load data
fstr <-lapply(fl, function(x) {
  data.table::fread(x)
  })
# Get the common SNPs in file ".fst" 
common <- intersect(intersect(fstr[[1]]$SNP,fstr[[2]]$SNP),fstr[[3]]$SNP)
fst <-lapply(fstr , function(x){
  x <- x[SNP %in% common,] 
})
# Estimate PBS
t <- lapply(fst, function(x) {
  -log10(1-x$FST)
})
out <- data.frame(CHROM = fst[[1]]$CHR,
           SNP = fst[[1]]$SNP,  
           POS = fst[[1]]$POS,
           PBS = 1/2 * (t[[1]] + t[[2]] -t[[3]]))

# Clean
out <- out[complete.cases(out), ]
out <- out[out$PBS != "Inf" & out$PBS != "-Inf", ]
summary(out)


# Apply windows

outw <- applyWindow(out)
outw1 <- outw[order(outw$CHR),]
outw2 <- outw1[!duplicated(outw1$PBS),]

#Choose 0.1% highest PBS
outper<- ch1percent(outw2)
summary(outper)
c <- 1 :nrow(outper)
out2 <- data.frame(Chr = outper$CHROM, 
                   start = (outper$POS %/% 100000) *100000, 
                   end = (outper$POS %/%100000)*100000 +100000,
                   rank = c)

out3 <- out2[order(out2$Chr),]
# Visualize
p <- vector("list")
p[[1]] <- khanty
p[[2]] <- Mansi
p[[3]] <- NENET
qqman::manhattan(khanty,
                 chr = "CHROM", 
                 bp = "POS", p = "PBS",
                 snp = "SNP", logp = F, 
                 col = c("orange3", "blue4"),
                 ylab = "PBS", 
                 main = "Khanty")
outwchr11 <- outper1[out3$CHR == 11, ]
summary(chr11)
ggplot(chr11) + geom_line(aes(x = POS, y = PBS)) 

#print multiple qqman to png
png("/home/vu/test/comp.png")
for(i in 1:3){
  qqman::manhattan(p[[i]],
                   chr = "CHROM", 
                   bp = "POS", p = "PBS",
                   snp = "SNP", logp = F, 
                   col = c("orange3", "blue4"),
                   ylab = "PBS")
}
dev.off()

Mansi <- out 
write.table(out3,file = "/home/vu/test/nenet/nenet.bed", sep = "\t",
            row.names = FALSE,
            col.names = F)






