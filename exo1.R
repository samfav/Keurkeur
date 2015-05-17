distanceEuclidienne <- function(x1,xt,y1,yt){
	return(sqrt((x1-xt)^2 + (y1-yt)^2))
}

nbFamille <- length(unique(zapp))

ceuc.app <- function(Xapp,zapp){
napp <- dim(Xapp)[1]
nbDimension <- dim(Xapp)[2]
x<- NULL
n <- NULL
x <- rbind(rep(0, nbFamille),rep(0,nbDimension))
for(z in 1:nbFamille){
	for(i in 1:napp){
		if(zapp[i]==z){
			n[z] <- length(zapp[zapp==z])
			for(j in 1:nbDimension){
				x[z,j] <- x[z,j] + Xapp[i,j]
			}
		}
		}
		for(j in 1:nbDimension){
				x[z,j] <- x[z,j] * (1/n[z])
		}	
	}
	return(x)
}
			

mu<-rbind(ceuc.app(Xapp,zapp))

ceuc.val <- function(mu,Xtst){
napp <- dim(Xtst)[1]
min <- rbind(rep(99999, napp))
nbDimension <- dim(Xtst)[2]
etiquette <- rbind(rep(0,napp))
distance <- rbind(rep(0, nbFamille),rep(0,napp))
for(i in 1:napp){
	for(z in 1:nbFamille){
		distance[z,i] <- distanceEuclidienne(Xtst[i,1],mu[z,1],Xtst[i,2],mu[z,2])
		if(distance[z,i]< min[i]){
			min[i] <- distance[z,i]
			etiquette[i]<-z
		}
	}
}
	return(etiquette)
}


kppv.app <- function(Xapp,zapp,Xval,zval,nppv){
}

kppv.val <- function(Xapp,zapp,L,Xtst){
}

#Test de nos fonctions

donn <- read.table("data/Synth1-40.txt",header=F)
X <- donn[,1:2]
z<- donn[,3]

Xapp <- X[c(1:15,21:35),]
zapp <- z[c(1:15,21:35)]
Xtst <- X[c(16:20,36:40),]
ztst <- z[c(16:20,36:40)]

ceuc.val(ceuc.app(Xapp,zapp),Xtst)
