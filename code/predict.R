predict = function(y, result, t){
    nu = result$nu
    Phi = result$Phi
    Zt = make_Zt(y, t, p)
    out = nu + Phi[, , t]%*%Zt
    return(c(out))
}
