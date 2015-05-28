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

calcul_vk_quadratique <- function(X, Z){
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



donn <- read.table("data/Synth1-1000.txt", header=F)
X <- donn[,1:2]
z <- donn[,3]
Xapp1 <- X[c(1:1000),]
zapp1 <- z[c(1:1000)]

donn <- read.table("data/Synth2-1000.txt", header=F)
X <- donn[,1:2]
z <- donn[,3]
Xapp2 <- X[c(1:1000),]
zapp2 <- z[c(1:1000)]

donn <- read.table("data/Synth3-1000.txt", header=F)
X <- donn[,1:2]
z <- donn[,3]
Xapp3 <- X[c(1:1000),]
zapp3 <- z[c(1:1000)]


nbFamille <- 2

%Question 1.1.1

%Analyse discriminante quadratique 

adq.app <- function(Xapp,zapp){
	mu <- ceuc.app(Xapp,zapp)
	pi  <- calcul_pik(Xapp,zapp)
	sigma <- var(Xapp)
	return(rbind(c(mu,pi,sigma)))
}

parametresADQ1 <- adq.app(Xapp1,zapp3)
parametresADQ2 <- adq.app(Xapp2,zapp3)
parametresADQ3 <- adq.app(Xapp3,zapp3)


%Analyse discrimante linéaire 

adl.app <- function(Xapp,zapp){
	mu <- ceuc.app(Xapp,zapp)
	pi  <- calcul_pik(Xapp,zapp)
	tempsigma <- var(Xapp)
	element1 <- sum(tempsigma[1,1],tempsigma[2,2])/2
	element2 <- sum(tempsigma[1,2],tempsigma[2,1])/2
	sigma <- matrix(c(element1,element2),2,2)
	return(rbind(c(mu,pi,sigma)))
}

parametresADL1 <- adl.app(Xapp1,zapp1)
parametresADL2 <- adl.app(Xapp2,zapp2)
parametresADL3 <- adl.app(Xapp3,zapp3)

%Classifieur Bayesien naif

nba.app <- function(Xapp,zapp){
	mu <- ceuc.app(Xapp,zapp)
	pi  <- calcul_pik(Xapp,zapp)
	tempsigma <- var(Xapp)
	element1 <- sum(tempsigma[1,1],tempsigma[2,2])/2
	element2 <- sum(tempsigma[1,2],tempsigma[2,1])/2
	sigmatemp <- matrix(c(element1,element2),2,2)
	sigma <- diag(x=2) * sigmatemp
	return(rbind(c(mu,pi,sigma)))
}

parametresNBA1 <- nba.app(Xapp1,zapp1)
parametresNBA2 <- nba.app(Xapp2,zapp2)
parametresNBA3 <- nba.app(Xapp3,zapp3)

ad.val <- function(parametres,Xtst){
	p1 = alpha * mvdnorm(parametres[1:2],parametres[6:8]) * parametres[4]
	p2 = alpha * mvdnorm(parametres[3:4],parametres[9:10]) * (1-parametres[4])
	print(p1)
	print(p2)
}



%On calcule les propbabilité en utilisant de Bayes.
%Densite de x sachant 1 -> deduire probabilité et p(1|x) = alpha * f(x|1) 
[Densite normal] * pi
%affecter famille ayant la plus grande probabilité. 
%p(2|x) alpha * (1-pi)
% Le calcul de la densite normal doit etre prit a partir de librairie 






