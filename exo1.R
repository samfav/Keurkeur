#1.1
distanceEuclidienne <- function(x1,xt,y1,yt){
	return(sqrt((x1-xt)^2 + (y1-yt)^2))
}


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
			if(K!=1){
				min <- min[order(min[,1]),]
			}
		}
		else if(distance[i,j]< min[K,1]){
			min[K,1] <- distance[i,j]
			min[K,2] <- j
			if(K!=1){
				min <- min[order(min[,1]),]
			}
		}
		}
		voisins[i,] <- min[,2]
	}
	return(voisins)
}

kppv.app <- function(Xapp,zapp,Xval,zval,nppv){
erreur <- 1
for(i in nppv){
    famille_K <- kppv.val(Xapp, zapp, i, Xval)
    erreur_K <- sum((famille_K == zval)==TRUE)/length(zval)
    print(erreur_K)
    if(erreur > erreur_K){
      erreur <- erreur_K
      K <- i
    }
  }
  return(K)
}

kppv.val <- function(Xapp,zapp,K,Xtst){
	ntst <- dim(Xtst)[1]
	etiquette <- rbind(rep(0,ntst))
	voisins <- getVoisins(Xapp,zapp,K,Xtst)
	for(i in 1:ntst){
		count1 <- 0
		count2 <- 0
		for(k in 1:K){
			if(zapp[voisins[i,k]]==1){
				count1 <- count1 + 1

			}
			else{
				count2 <- count2 + 1
			}
		}
		if(count1 > count2){
			etiquette[i]<-1
		}
		else{
			etiquette[i]<-2
		}
	}
	return(etiquette)
}

#Test de nos fonctions

donn <- read.table("data/Synth1-40.txt",header=F)
X <- donn[,1:2]
z<- donn[,3]

Xapp <- X[c(1:15,21:35),]
zapp <- z[c(1:15,21:35)]
Xtst <- X[c(16:20,36:40),]
ztst <- z[c(16:20,36:40)]

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

front.kppv <- function(kppv.val, K, Xapp, zapp)
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
	valf <- kppv.val(Xapp, zapp, K, grille)
	plot(Xapp, col=c("red","green","blue","magenta","orange")[zapp])
	contour(grilleX, grilleY, matrix(valf,nrow=naffX,byrow=T), add=T, drawlabels=FALSE, levels=1.5)
}

nbFamille <- length(unique(zapp))
nppv <- c(3,5,7,9,11,13,15,17)

K<- kppv.app(Xapp,zapp,Xtst,ztst,nppv)


#ceuc.val(ceuc.app(Xapp,zapp),Xtst)
#front.ceuc(ceuc.val,ceuc.app(Xapp,zapp),Xapp,zapp)
#front.kppv(kppv.val, K, Xapp, zapp)

#1.2
#Lecture des données

donn <- read.table("data/Synth1-40.txt", header=F)
X <- donn[,1:2]
z <- donn[,3]
Xapp40 <- X[c(1:40),]
zapp40 <- z[c(1:40)]


donn <- read.table("data/Synth1-100.txt", header=F)
X <- donn[,1:2]
z <- donn[,3]
Xapp100 <- X[c(1:100),]
zapp100 <- z[c(1:100)]


donn <- read.table("data/Synth1-500.txt", header=F)
X <- donn[,1:2]
z <- donn[,3]
Xapp500 <- X[c(1:500),]
zapp500 <- z[c(1:500)]

donn <- read.table("data/Synth1-1000.txt", header=F)
X <- donn[,1:2]
z <- donn[,3]
Xapp1000 <- X[c(1:1000),]
zapp1000 <- z[c(1:1000)]

#on prend quadratique

calcul_pik <- function(X,Z){
	nbFamille <- length(unique(Z))
	pi <- NULL
	n <- NULL
	n <- dim(X)[1]
	for(z in 1:nbFamille){
		n[z] <- length(Z[Z==z])
		pi[z]=n[z]/length(Z)
	}
	return(pi)
}

calcul_vk <- function(X, Z){
	x<- NULL
	n <- NULL
	napp <- dim(X)[1]
	nbDimension <- dim(X)[2]
	mu <- ceuc.app(X,Z)
	Vk <- rbind(rep(0, nbFamille),rep(0,nbDimension))
	for(z in 1:nbFamille){
		for(i in 1:napp){
			if(Z[i]==z){
				n[z] <- length(Z[Z==z])
				for(j in 1:nbDimension){
					temp <- (distanceEuclidienne(X[i,1], mu[j,1], X[i,2], mu[j,2]))^2
					Vk[z,j] <- Vk[z,j] + temp
					}
				}
			}
			for(j in 1:nbDimension){
				Vk[z,j] <- Vk[z,j] * (1/n[z])
				}	
		}
	return(Vk)
}

mu40 <-ceuc.app(Xapp40,zapp40)
mu100 <-ceuc.app(Xapp100,zapp100)
mu500 <- ceuc.app(Xapp500,zapp500)
mu1000 <- ceuc.app(Xapp1000,zapp1000)

pi40  <- calcul_pik(Xapp40,zapp40)
pi100 <- calcul_pik(Xapp100,zapp100)
pi500 <- calcul_pik(Xapp500,zapp500)
pi1000 <- calcul_pik(Xapp1000,zapp1000)

vk40  <- calcul_vk(Xapp40,zapp40)
vk100 <- calcul_vk(Xapp100,zapp100)
vk500 <- calcul_vk(Xapp500,zapp500)
vk1000 <- calcul_vk(Xapp1000,zapp1000)


#1.2.2
separ1 <- function(X, z) {
	g <- max(z)

	Xapp <- NULL
	zapp <- NULL
	Xtst <- NULL
	ztst <- NULL

	for (k in 1:g) {
	    indk <- which(z==k)
    	ntot <- length(indk)
	    napp <- round(ntot*2/3)
    	ntst <- ntot-napp

	    itot <- sample(indk)
    	iapp <- itot[1:napp]
	    itst <- itot[(napp+1):ntot]

    	Xapp <- rbind(Xapp, X[iapp,])
	    zapp <- c(zapp, z[iapp])
    	Xtst <- rbind(Xtst, X[itst,])
	    ztst <- c(ztst, z[itst])
	}

	res <- NULL
	res$Xapp <- Xapp
	res$zapp <- zapp
	res$Xtst <- Xtst
	res$ztst <- ztst

	res
}

#Calcul de E
calcul_E <- function(zappDonnees,zappCalcule) {
	m = length(zappDonnees)
	count1 <- 0
	for(i in 1:m){
			if(zappDonnees[i]!= zappCalcule[i]){
				count1 <- count1 + 1
				}
			}
			return(count1/m)
}

calcul_E(ztst,ceuc.val(ceuc.app(Xapp,zapp),Xtst))

separation_aleatoire <- function (Xdata, zdata) {
	t_erreur <- rbind(rep(0,20))

	t_test <- 0
	for (i in 1:20) {
		donn.sep <- separ1(Xdata,zdata)
		Xapp <- donn.sep$Xapp
		zapp <- donn.sep$zapp
		Xtst <- donn.sep$Xtst
		ztst <- donn.sep$ztst
		t_erreur[i] <- calcul_E(ztst,ceuc.val(ceuc.app(Xapp,zapp),Xtst))
	}
	return (t_apprentissage)
}

erreur40 <- separation_aleatoire(Xapp40, zapp40)
erreur100 <- separation_aleatoire(Xapp100, zapp100)
erreur500 <- separation_aleatoire(Xapp500, zapp500)
erreur1000 <- separation_aleatoire(Xapp1000, zapp1000)

mean(erreur40)
mean(erreur100)
mean(erreur500)
mean(erreur1000)