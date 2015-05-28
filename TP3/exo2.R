donn <- read.table("donnees-tp3/Synth1-40.txt", header=F)
X <- donn[,1:2]
z <- donn[,3]
A = matrix(   c(1, 0, 0, 1), nrow=2,ncol=2,byrow = TRUE)
M1 = c(-2,1)
M2 = c(1,1)
B = M2-M1
A%*%B
milieu = (M1+M2)/2
png("frontiere_Synth1-40.png")
plot(donn[,1:2], col=(donn[,3]))
abline(v=milieu[1])
dev.off()

donn <- read.table("donnees-tp3/Synth1-100.txt", header=F)
X <- donn[,1:2]
z <- donn[,3]
A = matrix(   c(1, 0, 0, 1), nrow=2,ncol=2,byrow = TRUE)
M1 = c(-2,1)
M2 = c(1,1)
B = M2-M1
A%*%B
milieu = (M1+M2)/2
png("frontiere_Synth1-100.png")
plot(donn[,1:2], col=(donn[,3]))
abline(v=milieu[1])
dev.off()

donn <- read.table("donnees-tp3/Synth1-500.txt", header=F)
X <- donn[,1:2]
z <- donn[,3]
A = matrix(   c(1, 0, 0, 1), nrow=2,ncol=2,byrow = TRUE)
M1 = c(-2,1)
M2 = c(1,1)
B = M2-M1
A%*%B
milieu = (M1+M2)/2
png("frontiere_Synth1-500.png")
plot(donn[,1:2], col=(donn[,3]))
abline(v=milieu[1])
dev.off()

donn <- read.table("donnees-tp3/Synth1-1000.txt", header=F)
X <- donn[,1:2]
z <- donn[,3]
A = matrix(   c(1, 0, 0, 1), nrow=2,ncol=2,byrow = TRUE)
M1 = c(-2,1)
M2 = c(1,1)
B = M2-M1
A%*%B
milieu = (M1+M2)/2
png("frontiere_Synth1-1000.png")
plot(donn[,1:2], col=(donn[,3]))
abline(v=milieu[1])
dev.off()
