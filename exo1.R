distanceEuclidienne <- function(x1,x2,y1,y2){
	return sqrt((x1-x2)^2 + (y1-y2)^2)
}

nbFamille <- length(unique(zapp))
napp <- dim(Xapp)[1]
nbDimension <- dim(Xapp)[2]


#On parcourt sur le nombre de famille puis
#on parcourt une fois le vecteur Xapp sur Y, 
#puis dans cette boucle on itere sur le nombre de dimension. 


ceuc.app <- function(Xapp,zapp){
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
