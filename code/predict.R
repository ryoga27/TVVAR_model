predict = function(y, result, t, p = 1){
    d = length(result$nu)
    nu = matrix(result$nu, nrow = d, ncol = 1)
    Phi = result$Phi
    Zt = make_Zt(y, t, p = p)
    out = nu + as.matrix(Phi[, , t], nrow = d, ncol = d)%*%Zt
    return(c(out))
}
