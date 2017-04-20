h1 <- c(10,20,10,20,10,20,10,20)
h2 <- c(20,10,20,10,20,10,20,10)

l1 <- c(1,3,1,3,1,3,1,3)
l2 <- c(3,1,3,1,3,1,3,1)

mat <- rbind(h1,h2,l1,l2)

par(mar=c(4,4,1,1))
plot(1:8,rep(0,8), ylim=c(0,35), pch="", xlab="Time", ylab="Gene Expression")

for (i in 1:nrow(mat)) {
  lines(1:8,mat[i,], lwd=3, col=i)
}

legend(1,35,rownames(mat), 1:4, cex=0.7)