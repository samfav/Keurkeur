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

