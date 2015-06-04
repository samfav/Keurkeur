On souhaite appliquer les trois modèles d’analyse discriminante et les deux modèles de régression logistique à la prédiction du diabète chez les individus d’une population d’amérindiens. On pourra charger les données au moyen du code suivant :

Donn <- read.csv("Pima.csv", header=T)
X <- Donn[,1:7]
z <- Donn[,8]
Xapp1 <- X[c(1:1000),]
zapp1 <- z[c(1:1000)]

Au cours d’une expérience 2 , on pourra utiliser le code suivant pour créer les données sur lesquelles apprendre le modèle de régression logistique quadratique :

Xapp2 <- Xapp
Xtst2 <- Xtst
for (p in 1:(dim(Xapp)[2]-1)) {
	for (q in (p+1):dim(Xapp)[2]) {
		Xapp2 <- cbind(Xapp2, Xapp[,p]*Xapp[,q])
		Xtst2 <- cbind(Xtst2, Xtst[,p]*Xtst[,q])
	}
}
for (p in 1:(dim(Xapp)[2]-1)) {
	Xapp2 <- cbind(Xapp2, Xapp[,p]^2)
	Xtst2 <- cbind(Xtst2, Xtst[,p]^2)
}

On utilisera ensuite le même protocole expérimental que pour les tests sur données simulées, en répétant l’expérience N = 100 fois. Calculer les taux moyens d’erreur de test pour chacun des cinq modèles étudiés. Que constate-t-on ? Comment expliquez-vous ces résultats ?

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
