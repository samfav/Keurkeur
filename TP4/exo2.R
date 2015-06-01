calcul_cout <- function(Xapp, zapp) {
	t_cout = matrix(sample(0, dim_p, 1, replace = TRUE), dim_p, 1)
	napp = length(zapp)
	for (i in 1:napp) {
		if (zapp[i] == 1) {
			t_cout[i] = 1
		}
		else t_cout[i] = 0
	}
	return (t_cout)
}

calcul_gradient <- function(beta, Xapp, t_cout) {
	n = length(t_cout)
	nbDimension <- dim(Xapp)[2]
	result = matrix(c(sample(0),n), nbDimension)
	for (j in 1:nbDimension) {
		for (i in 1:n) {
			result = result + t_cout[i] * Xapp[i][j] - exp(t(beta) * Xapp[i][j]) / (1 + exp(t(beta) * Xapp[i][j]))
		}
	}
	return (result)
}

#2.1.1

calcul_log <- function(beta, Xapp, t_cout) {
	n = length(t_cout)
	nbDimension <- dim(Xapp)[2]
	logL = matrix(c(sample(0),n), nbDimension)
	for (j in 1:nbDimension) {
		for (i in 1:n) {
			logL = logL + t_cout[i] * t(beta) * Xapp[i][j] - log(exp(t(beta) * Xapp[i][j]) + 1)
		}
	}
	return (logL)
}

log.app <- function(Xapp, zapp, intr, epsi) {
	nbDimension = dim(Xapp)[2]
	n = dim(Xapp)[1]
	if (intr == 1) {
		nbDimension = nbDimension + 1
	}
	beta = matrix(sample(0, nbDimension, 1, replace = TRUE), nbDimension, 1)
	logL = calcul_log(beta, Xapp, t_cout)
	for (i in 1:n) {
		ecart = abs(beta[i +1] - beta[i])	
		if (ecart > epsi) {
			%on change beta
			logL = calcul_log(beta, Xapp, t_cout)
		}
	}
	gradient = calcul_gradient(beta, Xapp, t_cout)	
	return(beta)
}

#2.1.2
La fonction log.val, permettant d’évaluer un ensemble de test, prendra comme arguments d’entrée le tableau de données Xtst à classer et la matrice beta des paramètres déterminés par la fonction log.app. Cette fonction fera un test sur les dimensions de la matrice beta, pour déterminer si une ordonnée à l’origine doit être ou non ajoutée à la matrice d’exemples Xtst à évaluer. Elle devra retourner une structure contenant la matrice prob des probabilités a posteriori estimées et le vecteur des classements associés. Ces fonctions pourront s’appuyer sur une fonction calculant les probabilités a posteriori à partir d’une matrice de paramètres beta et d’un tableau de données X.

log.val <- function(Xtst, beta) {
	result = matrix(sample(0, nbDimension, nbDimension + 1, replace = TRUE), nbDimension, nbDimension + 1)
	proba = matrix(sample(0, nbDimension, nbDimension, replace = TRUE), nbDimension, nbDimension)
	nbDimension = dim(Xtst)[2]
	n = dim(Xtst)[1]
	ordonnee = dim(beta)[1]
	
	if (ordonnee != n) {
		Xtst = ajout_ordonnee(Xtst, ordonnee)
	}
	
	for (i in 1:n) {
		for (j in 1:nbDimension) {
			proba[i][j] = exp(t(beta) * Xtst[i][j]) / (exp(t(beta) * Xtst[i][j]) + 1)
		}
	}
	
	return (result)
}

#2.1.3
Il est possible de généraliser le modèle de régression logistique de manière très simple. La stratégie consiste à transformer les données dans un espace plus complexe, dans lequel les classes peuvent être séparées par un hyperplan. La régression logistique est alors effectuée dans cet espace. Par exemple, dans le cas où les individus sont décrits par les variables X1 , X2 et X3 , la régression logistique quadratique consiste à effectuer la régression logistique classique dans l’espace correspondant aux variables X1 , X2 , X3 , X1X2 , X1X3 , X2X3 , (X1 ) 2 , (X2 ) 2 et (X3 ) 2 , que l’on notera X 2 , plutôt que dans l’espace X = {X1 , X2 , X3}. Le modèle ainsi défini est donc plus flexible ; mais le nombre de paramètres à estimer étant plus important (p(p+ 3)/2 au lieu de p), il peut également s’avérer moins robuste que le modèle classique déterminé dans l’espace des caractéristiques initiales.
