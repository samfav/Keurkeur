setwd("/Users/MarieChatelin/GitHub/Keurkeur/TP4/data/")
library(MASS)

#On choisit presque toujours d'ajouter un intercept qui est le beta(0) pour rendre notre fonction indépendantes des variations des varibles

init <- function () {
	donn <- read.table("Synth1-1000.txt", header=F)
	X <- donn[,1:2]
	z <- donn[,3]
	Xapp1 <- X[c(1:1000),]
	zapp1 <- z[c(1:1000)]
}

calcul_cout <- function(Xapp, zapp, nbDimension) {
	t_cout = matrix(sample(0), nbDimension, 1)
	napp = length(zapp)
	for (i in 1:napp) {
		if (zapp[i] == 1) {
			t_cout[i] = 1
		}
		else t_cout[i] = 0
	}
	return (t_cout)
}

calcul_log <- function(beta, Xapp, t_cout, nbDimension) {
	n = dim(Xapp)[1]
	logL = 0
	for (i in 1:n) {
		for (j in 1:nbDimension) {
			X = c(X, Xapp[i,j])
			#Utiliser cbind ?
		}
		logL = logL + t_cout[i] * t(beta) %*% X - log(exp(t(beta) %*% X + 1))
	}
	return (logL)
}

calcul_gradient <- function(beta, Xapp, t_cout, nbDimension) {
	n = dim(Xapp)[1]
	result = matrix(sample(0),nbDimension, 1)
	X = matrix(sample(0),n, nbDimension)
	for (i in 1:n) {
		for (j in 1:nbDimension) {
			X = c(X, Xapp[i,j])
			#Utiliser cbind ?
		}
		result = result + t_cout[i] * X - X * exp(t(beta) %*% X) / (1 + exp(t(beta) %*% X))
	}
	return (result)
}

calcul_H <- function(beta, Xapp, zapp, nbDimension){
    n = dim(Xapp)[1]
    W = matrix(sample(0), n,n)
    for(i in 1:n){
		for (j in 1:nbDimension) {
			X = c(X, Xapp[i,j])
			#Utiliser cbind ?
		}
        W[i,i] = exp(t(beta) %*% X) / (1 + exp(t(beta) %*% X))^2
    }
    H = -as.matrix(t(Xapp)) %*% as.matrix(W) %*% as.matrix(Xapp)
    return (H)
}

log.app <- function(Xapp, zapp, intr, epsi) {
	n = dim(Xapp)[1]
	nbDimension = dim(Xapp)[2]
    if (intr == 1) {
        nbDimension = nbDimension + 1
    }
    beta = matrix(sample(0), nbDimension, 1)
    t_cout = calcul_cout(Xapp, zapp, nbDimension)
    logL = calcul_log(beta, Xapp, t_cout, nbDimensioin)
	H = matrix(sample(1), nbDimension, nbDimension)	
	ecart = 1
    for (i in 1:(nbDimension - 1)) {
		print (beta)  
		print(epsi)
		print(ecart)
		if (ecart > epsi) {
			print(3)
           	H = calcul_H(Xapp, zapp, beta, nbDimension)
			print (H)
           	grad = calcul_gradient(beta, Xapp, t_cout, nbDimension)
           	beta = beta - ginv(H)%*%grad
        	ecart = abs(beta[i +1] - beta[i]) 
			print (beta)
		 }	
    }
    return(beta)
}


#2.1.2
La fonction log.val, permettant d’évaluer un ensemble de test, prendra comme arguments d’entrée le tableau de données Xtst à classer et la matrice beta des paramètres déterminés par la fonction log.app. Cette fonction fera un test sur les dimensions de la matrice beta, pour déterminer si une ordonnée à l’origine doit être ou non ajoutée à la matrice d’exemples Xtst à évaluer. Elle devra retourner une structure contenant la matrice prob des probabilités a posteriori estimées et le vecteur des classements associés. Ces fonctions pourront s’appuyer sur une fonction calculant les probabilités a posteriori à partir d’une matrice de paramètres beta et d’un tableau de données X.

log.val <- function(Xtst, beta) {
	nbDimension = dim(Xtst)[2]
	n = dim(Xtst)[1]
	nbDimension_beta = dim(beta)[1]
	if (nbDimension != nbDimension_beta) {
		nbDimension = nbDimesion_beta
	} 
	#DIMENSION ?!
	Xtst = cbind(Xtst, beta[1])
	result = matrix(sample(0), nbDimension, nbDimension + 1)
	proba = matrix(sample(0), nbDimension, nbDimension)
	classement = matrix(sample(0), n, 1)
	
	for (i in 1:n) {
		for (j in 1:nbDimension) {
			X = c(X, Xtst[i,j])
			#Utiliser cbind ?
		}
		proba = proba + exp(t(beta) %*% X) / (exp(t(beta) %*% X) + 1)
	}
	result = cbind(proba, classement)
	return (result)
}

#2.1.3
Il est possible de généraliser le modèle de régression logistique de manière très simple. La stratégie consiste à transformer les données dans un espace plus complexe, dans lequel les classes peuvent être séparées par un hyperplan. La régression logistique est alors effectuée dans cet espace. Par exemple, dans le cas où les individus sont décrits par les variables X1 , X2 et X3 , la régression logistique quadratique consiste à effectuer la régression logistique classique dans l’espace correspondant aux variables X1 , X2 , X3 , X1X2 , X1X3 , X2X3 , (X1 ) 2 , (X2 ) 2 et (X3 ) 2 , que l’on notera X 2 , plutôt que dans l’espace X = {X1 , X2 , X3}. Le modèle ainsi défini est donc plus flexible ; mais le nombre de paramètres à estimer étant plus important (p(p+ 3)/2 au lieu de p), il peut également s’avérer moins robuste que le modèle classique déterminé dans l’espace des caractéristiques initiales.

#Ici on est en quadratique

#Il faut écrire une fonction qui transforme X 2 en Xapp afin de pouvoir réappliquer les fonctions au dessus.


#2.2 Test sur données simulées On souhaite comparer les performances des deux modèles de régression logistique sur les trois jeux de données simulées étudiés précédemment. Pour ce faire, on utilisera le même protocole expérimental que pour l’analyse discriminante : pour chaque jeu de données, calculer le taux d’erreur (de test) moyen sur les N = 20 séparations effectuées. Qu’observe-t-on ? Comment interpréter ces résultats ? On pourra utiliser les fonctions prob.log et prob.log2 1 , disponible sur le site de l’UV, pour visualiser les courbes de niveau ou les frontières de décision obtenues respectivement par régression logistique et par régression logistique quadratique.
