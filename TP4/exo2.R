setwd("/Users/MarieChatelin/GitHub/Keurkeur/TP4/data")
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

calcul_log <- function(beta, Xapp, zapp) {
    p = exp(as.matrix(Xapp)%*%beta)/(1+exp(as.matrix(Xapp)%*%beta))
    logL<- t(as.matrix(Xapp))%*%(zapp-p)
    return (logL)
}

calcul_gradient <- function(beta, Xapp, t_cout, nbDimension) {
    n = dim(Xapp)[1]
    result = matrix(sample(0),nbDimension, 1)
    X = matrix(sample(0),n, nbDimension)
    for (i in 1:n) {
        if (nbDimension == 3) {
            X = c(beta[1], Xapp[i,1], Xapp[i,2])
        }
        else {
                X = c(Xapp[i,1], Xapp[i,2])
             }
        result = result + t_cout[i] * X - X * exp(t(beta) %*% X) / (1 + exp(t(beta) %*% X))
    }
    return (result)
}

calcul_H <- function(beta, Xapp, zapp, nbDimension){
    n = dim(Xapp)[1]
    W = matrix(sample(0), n,n)
    for(i in 1:n){
        if (nbDimension == 3) {
            X = c(beta[1], Xapp[i,1], Xapp[i,2])
        }
        else {
                X = c(Xapp[i,1], Xapp[i,2])
             }
        W[i,i] = exp(t(beta) %*% X) / (1 + exp(t(beta) %*% X))^2
    }
    H = -as.matrix(t(Xapp)) %*% as.matrix(W) %*% as.matrix(Xapp)
    return (H)
}
calcul_hessienne<-function(beta, Xapp, nbDimension){
    n = dim(Xapp)[1]
    W = matrix(sample(0), n,n)
    prod = as.matrix(Xapp)%*%beta
    W = exp(prod)/((1+exp(prod))^2)
    H = -1 * as.matrix(t(Xapp)) %*% as.matrix(diag(as.numeric(W))) %*% as.matrix(Xapp)
    return (H)
}

calcul_p <-function(beta, Xapp){
    n = dim(Xapp)[1]
    p = matrix(sample(0),n, 1)
    prod = as.matrix(Xapp) %*% beta
    p = exp(prod)/(1+exp(prod))
    return (p)
}

log.app2 <- function(Xapp, zapp, intr, epsi){
    n = dim(Xapp)[1]
    nbDimension = dim(Xapp)[2]
    if (intr == 1) {
        nbDimension = nbDimension + 1
        unit =  matrix(sample(1), n, 1)
        Xapp = cbind(as.matrix(unit), as.matrix(Xapp))
    }
    beta = matrix(sample(0), nbDimension, 1)
    t_cout = calcul_cout(Xapp, zapp, nbDimension)
    ecart = 1
    i<-0
    while (ecart > epsi){
        H = -1 * calcul_hessienne(beta, Xapp, nbDimension)
        p = calcul_p(beta, Xapp)
        beta_n = beta + solve(H)%*%t(Xapp)%*%(t_cout - p) # solve ou ginv
        #ecart = dist(rbind(t(beta), t(beta_n)), "euclidian")
        ecart = sqrt(sum((beta_n-beta)^2))
        beta = beta_n
        i = i+1
    }  
    logL = calcul_log(beta, Xapp, zapp)
    print(logL)
    print (i)
    return(beta)
}

log.app <- function(Xapp, zapp, intr, epsi) {
    n = dim(Xapp)[1]
    nbDimension = dim(Xapp)[2]
    if (intr == 1) {
        nbDimension = nbDimension + 1
        unit =  matrix(sample(1), n, 1)
        Xapp = cbind(as.matrix(unit), as.matrix(Xapp))
    }
    beta = matrix(sample(0), nbDimension, 1)
    t_cout = calcul_cout(Xapp, zapp, nbDimension)
    H = matrix(sample(1), nbDimension, nbDimension) 
    ecart = 1
    i<-0
    while (ecart > epsi){
        H = calcul_H(beta, Xapp, zapp, nbDimension)
        grad = calcul_gradient(beta, Xapp, t_cout, nbDimension)
        beta_n = beta - solve(H)%*%grad # solve ou ginv
        #ecart = dist(rbind(t(beta), t(beta_n)), "euclidian")
        ecart = sqrt(sum((beta_n-beta)^2))
        beta = beta_n
        i = i+1
    }  
    logL = calcul_log(beta, Xapp, zapp)
    print(logL)
    print (i)
    return(beta)
}


#2.1.2
#La fonction log.val, permettant d’évaluer un ensemble de test, prendra comme arguments d’entrée le tableau de données Xtst à classer et la matrice beta des paramètres déterminés par la fonction log.app. Cette fonction fera un test sur les dimensions de la matrice beta, pour déterminer si une ordonnée à l’origine doit être ou non ajoutée à la matrice d’exemples Xtst à évaluer. Elle devra retourner une structure contenant la matrice prob des probabilités a posteriori estimées et le vecteur des classements associés. Ces fonctions pourront s’appuyer sur une fonction calculant les probabilités a posteriori à partir d’une matrice de paramètres beta et d’un tableau de données X.

log.val <- function(Xtst, beta) {
    nbDimension = dim(Xtst)[2]
    n = dim(Xtst)[1]
    result = NULL
    nbDimension_beta = dim(beta)[1]
    if (nbDimension != nbDimension_beta) {
        nbDimension = nbDimension_beta
        unit =  matrix(sample(1), n, 1)
        Xtst = cbind(as.matrix(unit), as.matrix(Xtst))
    } 
    proba = calcul_p(beta, Xtst)
    classement = matrix(sample(0), n,1)
    Xtst = as.matrix(Xtst)
    beta = as.matrix(beta)
    prob1 <- exp(Xtst%*%beta)/(1 + exp(Xtst%*%beta))
    prob2 <- 1/(1 + exp(Xtst%*%beta))

    proba <- cbind(prob1,prob2)
    classement<-max.col(proba)

    result$prob = proba
    result$classement = classement
    return (result)
}

#####
#tests
val = log.val(Xapp1, beta)
new_data = cbind(Xapp1, as.matrix(val$classement))
new_data = cbind(donnees$Xtst, donnees$ztst, as.matrix(val$classement))
plot(new_data[,1:2], pch=new_data[,3])

#2.1.3
#Il est possible de généraliser le modèle de régression logistique de manière très simple.
# La stratégie consiste à transformer les données dans un espace plus complexe, dans lequel les classes peuvent 
# être séparées par un hyperplan. La régression logistique est alors effectuée dans cet espace. 
# Par exemple, dans le cas où les individus sont décrits par les variables X1 , X2 et X3 , la régression 
# logistique quadratique consiste à effectuer la régression logistique classique dans l’espace correspondant 
# aux variables X1 , X2 , X3 , X1X2 , X1X3 , X2X3 , (X1 ) 2 , (X2 ) 2 et (X3 ) 2 , que l’on notera X 2 , plutôt 
# que dans l’espace X = {X1 , X2 , X3}. Le modèle ainsi défini est donc plus flexible ; mais le nombre de
#  paramètres à estimer étant plus important (p(p+ 3)/2 au lieu de p), il peut également s’avérer moins 
#  robuste que le modèle classique déterminé dans l’espace des caractéristiques initiales.

#Ici on est en quadratique

XtoXapp <-function(X){ # valable uniquement pour X =2
    X = as.matrix(X)
    nbDimension = dim(X)[2]
    n = dim(X)[1]
    X1 = as.matrix(X[,1])
    X2 = as.matrix(X[,2])
    X1X2 = as.matrix(X1 * X2)
    X1_2 = as.matrix(X1^2)
    X2_2 = as.matrix(X2^2)
    Xapp = cbind(X1,X2,X1X2,X1_2,X2_2)
    return (Xapp)
}


#2.2 Test sur données simulées On souhaite comparer les performances des deux modèles de régression logistique
# sur les trois jeux de données simulées étudiés précédemment. Pour ce faire, on utilisera le même protocole 
# expérimental que pour l’analyse discriminante : pour chaque jeu de données, calculer le taux d’erreur (de test)
#  moyen sur les N = 20 séparations effectuées. Qu’observe-t-on ? Comment interpréter ces résultats ? On pourra 
#  utiliser les fonctions prob.log et prob.log2 1 , disponible sur le site de l’UV, pour visualiser les courbes 
#  de niveau ou les frontières de décision obtenues respectivement par régression logistique et par régression 
#  logistique quadratique.

separ1 <- function(X, z) {
    g <- max(z)

    Xapp <- NULL
    zapp <- NULL
    Xtst <- NULL
    ztst <- NULL

    for (k in 1:g)
    {
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

calcul_taux_erreur <- function(ztst, classement){
    taux = 0
    taux = sum(ztst-classement)/dim(as.matrix(ztst))[1]*100
    return (taux)
}

estim_erreur <- function(donn, quad) {
    X <- donn[,1:2]
    z <- donn[,3]
    estim = 0
    if (quad == 1) {
        X= XtoXapp(X)
    }
    for (i in 1:20) {
        donnees = separ1(X,z)
        beta = log.app2(donnees$Xapp, donnees$zapp, 0, 0.00001)
        val = log.val(donnees$Xtst, beta)
        estim = estim + calcul_taux_erreur(donnees$ztst, val$classement)
    }
    
    return (estim/20)
}

prob.log <- function(param, X, z, niveaux)
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

    grille <- as.matrix(grille)

    # calcul des valeurs de la fonction 
    valf <- log.val(param, grille)$prob[,1]
    plot(X, col=c("red","green","blue","magenta","orange")[z])
    contour(grilleX, grilleY, matrix(valf,nrow=naffX,byrow=T), add=T, drawlabels=FALSE, levels=niveaux)
}

prob.log2 <- function(param, X, z, niveaux)
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
    grille <- cbind(grille, grille[,1]*grille[,2], grille[,1]^2, grille[,2]^2)

    grille <- as.matrix(grille)

    # calcul des valeurs de la fonction 
    valf <- log.val(param, grille)$prob[,1]
    plot(X, col=c("red","green","blue","magenta","orange")[z])
    contour(grilleX, grilleY, matrix(valf,nrow=naffX,byrow=T), add=T, drawlabels=FALSE, levels=niveaux)
}
















#
