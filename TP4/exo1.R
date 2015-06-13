library(MASS)

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

donn <- read.table("data/Synth1-1000.txt", header=F)
Xdata <- donn[,1:2]
zdata <- donn[,3]
donn.sep <- separ1(Xdata,zdata)
Xapp1 <- donn.sep$Xapp
zapp1 <- donn.sep$zapp
Xtst1 <- donn.sep$Xtst
ztst1 <- donn.sep$ztst


donn <- read.table("data/Synth2-1000.txt", header=F)
Xdata <- donn[,1:2]
zdata <- donn[,3]
donn.sep <- separ1(Xdata,zdata)
Xapp2 <- donn.sep$Xapp
zapp2 <- donn.sep$zapp
Xtst2 <- donn.sep$Xtst
ztst2 <- donn.sep$ztst

donn <- read.table("data/Synth3-1000.txt", header=F)
Xdata <- donn[,1:2]
zdata <- donn[,3]
donn.sep <- separ1(Xdata,zdata)
Xapp3 <- donn.sep$Xapp
zapp3 <- donn.sep$zapp
Xtst3 <- donn.sep$Xtst
ztst3 <- donn.sep$ztst

nbFamille <- 2

#Question 1.1.1

#Analyse discriminante quadratique 

adq.app <- function(Xapp,zapp){
	parametres <- list() 
	pi  <- calcul_pik(Xapp,zapp)
	for(i in 1:2){
		Xappi = Xapp[zapp==i,]
		mu = colMeans(Xappi)
		sigma <- var(Xappi) * nrow(Xappi)/(nrow(Xappi)-1)
		parametres[[i]] = list(pi=pi[i],mu=mu,sigma=sigma)
	}
	return(parametres)
}

parametresADQ1 <- adq.app(Xapp1,zapp3)
parametresADQ2 <- adq.app(Xapp2,zapp3)
parametresADQ3 <- adq.app(Xapp3,zapp3)


#Analyse discrimante linéaire 

adl.app <- function(Xapp,zapp){
	parametres <- list() 
	pi  <- calcul_pik(Xapp,zapp)
	sigma <- 0
	for(i in 1:2){
        Xappi = Xapp[zapp==i,]
        mu = colMeans(Xappi)
        sigma =(sigma + var(Xappi) * nrow(Xappi)) / (nrow(Xapp) - 2)
        parametres[[i]] = list(pi=pi[i],mu=mu,sigma=sigma)
    }
    return(parametres)
}

parametresADL1 <- adl.app(Xapp1,zapp1)
parametresADL2 <- adl.app(Xapp2,zapp2)
parametresADL3 <- adl.app(Xapp3,zapp3)

#Classifieur Bayesien naif

nba.app <- function(Xapp,zapp){
	parametres <- list() 
	pi  <- calcul_pik(Xapp,zapp)
	for(i in 1:2){
		Xappi = Xapp[zapp==i,]
		mu <- colMeans(Xappi)
	    sigma = diag(diag(var(Xappi) * nrow(Xappi)/(nrow(Xappi)-1)))
        parametres[[i]] = list(pi=pi[i],mu=mu,sigma=sigma)
    }
	return(parametres)
}

parametresNBA1 <- nba.app(Xapp1,zapp1)
parametresNBA2 <- nba.app(Xapp2,zapp2)
parametresNBA3 <- nba.app(Xapp3,zapp3)

mvdnorm <- function(X, mu, Sigma)
{
	X <- as.matrix(X)
	n <- dim(X)[1]
	p <- dim(X)[2]
	B <- chol(Sigma)
	U <- (X-matrix(rep(mu,n),nrow=n,byrow=T))%*%ginv(B)
	dens <- exp(-rowSums(U*U)/2) * (2*pi)^(-p/2) / det(B)
}


ad.val <- function(parametres,Xtst){
	ztst <- matrix(data = 99, nrow = dim(Xtst)[1], ncol = 3)
	for (i in 1:dim(Xtst)[1]){
		densite1 = mvdnorm(Xtst[i,],parametres[[1]]$mu,parametres[[1]]$sigma)
		densite2 = mvdnorm(Xtst[i,],parametres[[2]]$mu,parametres[[2]]$sigma)
		densite = densite1  * parametres[[1]]$pi + densite2 * parametres[[2]]$pi
		p1 = densite1 * parametres[[1]]$pi / densite
		p2 =  densite2 * parametres[[2]]$pi / densite
		if(p1 > p2){
			ztst[i,1]<-p1
			ztst[i,2]<-p2
			ztst[i,3]<-1
		}
		else{
			ztst[i,1]<-p1
			ztst[i,2]<-p2
			ztst[i,3]<-2
		}
	}
	return(ztst)
} 

probNBA1 <- ad.val(parametresNBA1,Xtst1)
probNBA2 <- ad.val(parametresNBA2,Xtst2)
probNBA3 <- ad.val(parametresNBA3,Xtst3)

probADQ1 <- ad.val(parametresADQ1,Xtst1)
probADQ2 <- ad.val(parametresADQ2,Xtst2)
probADQ3 <- ad.val(parametresADQ3,Xtst3)

probADL1 <- ad.val(parametresADL1,Xtst1)
probADL2 <- ad.val(parametresADL2,Xtst2)
probADL3 <- ad.val(parametresADL3,Xtst3)


prob.ad <- function(param, X, z, niveaux)
{
    discretisation=50
    deltaX <- (max(X[,1]) -min(X[,1]))/discretisation
    deltaY <- (max(X[,2]) -min(X[,2]))/discretisation
    minX <- min(X[,1])-deltaX
    maxX <- max(X[,1])+deltaX
    minY <- min(X[,2])-deltaY
    maxY <- max(X[,2])+deltaY
    
    # grille d'affichage 
    grilleX <- seq(from=minX,to=maxX,by=deltaX)
    naffX <- length(grilleX)
    grilleY <- seq(from=minY,to=maxY,by=deltaY)
    naffY <- length(grilleY)
    grille <- cbind(rep.int(grilleX,times=rep(naffY,naffX)),rep(grilleY,naffX))

    # calcul des valeurs de la fonction 
    valf <- ad.val(param, grille)$prob[,1]
    plot(X, col=c("red","green","blue","magenta","orange")[z])
    contour(grilleX, grilleY, matrix(valf,nrow=naffX,byrow=T), add=T, drawlabels=FALSE, levels=niveaux)
}

#prob.ad(parametresADL1,Xtst1,ztstADL1,2)

#Test sur donnes simulees

separ2 <- function(X, z)
{
	g <- max(z)

	Xapp <- NULL
	zapp <- NULL
	Xval <- NULL
	zval <- NULL
	Xtst <- NULL
	ztst <- NULL

	for (k in 1:g)
	{
	    indk <- which(z==k)
    	ntot <- length(indk)
	    napp <- round(ntot/2)
		nval <- round(ntot/4)
    	ntst <- ntot-napp-nval

	    itot <- sample(indk)
    	iapp <- itot[1:napp]
    	ival <- itot[(napp+1):(napp+nval)]
	    itst <- itot[(napp+nval+1):ntot]

    	Xapp <- rbind(Xapp, X[iapp,])
	    zapp <- c(zapp, z[iapp])
    	Xval <- rbind(Xval, X[ival,])
	    zval <- c(zval, z[ival])
    	Xtst <- rbind(Xtst, X[itst,])
	    ztst <- c(ztst, z[itst])
	}

	res <- NULL
	res$Xapp <- Xapp
	res$zapp <- zapp
	res$Xval <- Xval
	res$zval <- zval
	res$Xtst <- Xtst
	res$ztst <- ztst

	res
}

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

separation_aleatoire2 <- function (Xdata, zdata) {
	t_erreurADL <- rbind(rep(0,20))
	t_erreurADQ <- rbind(rep(0,20))
	t_erreurNBA <- rbind(rep(0,20))
	t_test <- 0
	for (i in 1:20) {
		donn.sep <- separ2(Xdata, zdata)
		Xapp <- donn.sep$Xapp
		zapp <- donn.sep$zapp
		Xval <- donn.sep$Xval
		zval <- donn.sep$zval
		Xtst <- donn.sep$Xtst
		ztst <- donn.sep$ztst
		parametresADL <- adl.app(Xapp,zapp)
		parametresADQ <- adq.app(Xapp,zapp)
		parametresNBA <- nba.app(Xapp,zapp)
		t_erreurADL[i] <- calcul_E(ztst,ad.val(parametresADL,Xtst))
		t_erreurADQ[i] <- calcul_E(ztst,ad.val(parametresADQ,Xtst))
		t_erreurNBA[i] <- calcul_E(ztst,ad.val(parametresNBA,Xtst))
	}
	return (t_erreurADL)
}

erreur1 <- separation_aleatoire2(Xapp1, zapp1)
erreur2 <- separation_aleatoire2(Xapp2, zapp2)
erreur3 <- separation_aleatoire2(Xapp3, zapp3)

