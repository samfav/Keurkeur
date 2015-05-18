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

front.ceuc <- function(ceuc.val, mu, Xapp, zapp)
{
	minX <- min(Xapp[,1])-1
	maxX <- max(Xapp[,1])+1
	minY <- min(Xapp[,2])-1
	maxY <- max(Xapp[,2])+1
	# grille d'affichage 
	grilleX <- seq(from=minX,to=maxX,by=0.01)
	naffX <- length(grilleX)
	grilleY <- seq(from=minY,to=maxY,by=0.01)
	naffY <- length(grilleY)

	grille <- cbind(rep.int(grilleX,times=rep(naffY,naffX)),rep(grilleY,naffX))

	# calcul des valeurs de la fonction 
	valf <- ceuc.val(mu, grille)
	plot(Xapp, col=c("red","green","blue","magenta","orange")[zapp])
	contour(grilleX, grilleY, matrix(valf,nrow=naffX,byrow=T), add=T, drawlabels=FALSE, levels=1.5)
}


nppv <- rbind(c(1,2,3,5))

getVoisins <- function(Xapp,zapp,K,Xtst){
	napp <- dim(Xapp)[1]
	ntst <- dim(Xtst)[1]
	distance <- matrix(0, ntst, napp)
	voisins <- matrix(0,ntst,K)
	for(i in 1:ntst){
		min <- matrix(99999,K,2)
		for(j in 1:napp){
		distance[i,j] <- distanceEuclidienne(Xtst[i,1],Xapp[j,1],Xtst[i,2],Xapp[j,2])
		if(j <= K ){
			min[j,1] <- distance[i,j]
			min[j,2] <- j
			min <- min[order(min[,1]),]
		}
		else if(distance[i,j]< min[K,1]){
			min[K,1] <- distance[i,j]
			min[K,2] <- j
			min <- min[order(min[,1]),]
		}
		}
		voisins[i,] <- min[,2]
	}
	return(voisins)
}

kppv.app <- function(Xapp,zapp,Xval,zval,nppv){
}

kppv.val <- function(Xapp,zapp,K,Xtst){
	voisins <- getVoisins(Xapp,zapp,K,Xtst)
	return(voisins)
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
front.ceuc(ceuc.val,ceuc.app(Xapp,zapp),Xapp,zapp)
