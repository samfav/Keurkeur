distanceEuclidienne <- function(x1,x2,y1,y2){
	return sqrt((x1-x2)^2 + (y1-y2)^2)
}

nbFamille <- length(unique(zapp))
napp <- dim(Xapp)[1]
nbDimension <- dim(Xapp)


ceuc.app <- function(Xapp,zapp){
x<- NULL 
mu <- NULL
n <- NULL
for( i in 1:nbDimension){
	n[i] <- length(zapp[zapp==i])
	for( j in 1:napp){
		ifelse(zapp[j]==1, z <- 1, z <- 0)
		x[i][i] <- x[i][i] + z*Xapp[j,i]
	}
	x[i][] <- x[i][] * (1/n[i])
	}
	return x[][]
}


mu<-ceuc.app(Xapp,zapp)

napp <- dim(Xtst)[1]
ceuc.val <- function(mu,Xtst){
for( i in 1:napp){
	for( j in 1:nbclass){
		dist[i][j] <- distanceEuclidienne(Xapp[i][j], mu[i][j])
	}
	}
}

kppv.app <- function(Xapp,zapp,Xval,zval,nppv){
}

kppv.val <- function(Xapp,zapp,L,Xtst){
}

%Test de nos fonctions

donn <- read.table("data/Synth1-40.txt",header=F)
X <- donn[,1:2]
z<- donn[,3]

Xapp <- X[c(1:15,21:35),]
zapp <- z[c(1:15,21:35)]
Xtst <- X[c(16:20,36:40),]
ztst <- z[c(16:20,36:40)]

ceuc.val(ceuc.app(Xapp,zapp),Xtst)
