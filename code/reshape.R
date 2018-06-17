vec = function(mat){
    n = nrow(mat)
    m = ncol(mat)
    out = matrix(0, nrow = n*m, ncol = 1)
    out[, 1] = c(mat)
    return(out)
}

### sample ###
# n = 5
# m = 3
#
# mat_test = matrix(1:(n*m), nrow = n, ncol = m)
# vec_test = vec(mat_test)
# is.matrix(vec_test)
# nrow(vec_test) == n*m
# ncol(vec_test) == 1
### end ###

make_yt = function(y, t){
    T = nrow(y)
    d = ncol(y)
    out = matrix(0, nrow = d, ncol = 1)
    out[1:d] = y[t, 1:d]
    return(out)
}

### sample code ###
# T_test = 5
# d_test = 2
#
# y_test = matrix(0, nrow = T_test, ncol = d_test)
# y_test[1:length(y_test)] = 1:length(y_test)
# print(y_test)
#
# t_test = 3
#
# yt_test = make_yt(y = y_test, t = t_test)
# print(yt_test)
# is.matrix(yt_test)
# nrow(yt_test) == d_test
### end ###

make_Y = function(y, p){
    T = nrow(y)
    d = ncol(y)
    out = rep(0, d*(T-p))
    s = 1
    for(t in (p+1):T){
        out[1:d + (s-1)*d] = make_yt(y, t)
        s = s + 1
    }
    out = matrix(out, nrow = d*(T-p), ncol = 1)
    return(out)
}

### sample code ###
# T_test = 5
# d_test = 2
#
# y_test = matrix(0, nrow = T_test, ncol = d_test)
# y_test[1:length(y_test)] = 1:length(y_test)
# print(y_test)
#
# p_test = 1
#
# Y_test = make_Y(y = y_test, p = p_test)
# print(Y_test)
#
# p_test = 2
# Y_test = make_Y(y = y_test, p = p_test)
# print(Y_test)
# is.matrix(y_test)
# nrow(Y_test) == d_test*(T_test - p_test)
### end ###

make_Zt = function(y, t, p){
    T = nrow(y)
    d = ncol(y)
    out = matrix(0, nrow = d*p, ncol = 1)
    for(l in 1:p){
        out[1:d + d*(l-1), ] = make_yt(y, t-l)
    }
    return(out)
}

### sample code ###
# T_test = 5
# d_test = 2
# p_test = 1
#
# y_test = matrix(0, nrow = T_test, ncol = d_test)
# y_test[1:length(y_test)] = 1:length(y_test)
# print(y_test)
#
# t_test = 3
#
# Zt_test = make_Zt(y = y_test, t = t_test, p = p_test)
# print(Zt_test)
#
# p_test = 2
# Zt_test = make_Zt(y = y_test, t = t_test, p = p_test)
# print(Zt_test)
# is.matrix(Zt_test)
# nrow(Zt_test) == d_test*p_test
# ncol(Zt_test) == 1
### end ###

make_Z = function(y, p){
    T = nrow(y)
    d = ncol(y)
    out = matrix(0, nrow = d*(T-p), ncol = d + d^2*p*(T-p))
    Id = diag(d)
    for(t in 1:(T-p)){
        out[1:d + d*(t-p), 1:d] = Id
    }
    s = 1
    for(t in (p+1):T){
        out[1:d + d*(s-1), 1:(d^2*p) + d + (d^2*p)*(s-1)] = kronecker(t(make_Zt(y, t = t, p = p)), Id)
        s = s + 1
    }
    return(out)
}

### sample code ###
# T_test = 5
# d_test = 2
# p_test = 1
#
# y_test = matrix(0, nrow = T_test, ncol = d_test)
# y_test[1:length(y_test)] = 1:length(y_test)
# print(y_test)
#
# t_test = 3
#
# Zt_test = make_Zt(y = y_test, t = t_test, p = p_test)
# print(Zt_test)
#
# Z_test = make_Z(y = y_test, p = p_test)
# print(Z_test)
# is.matrix(Z_test)
# nrow(Z_test) == d_test*(T_test - p_test)
# ncol(Z_test) == d_test + d_test^2*p_test*(T_test - p_test)
#
# T_test = 10
# d_test = 2
# p_test = 1
#
# y_test = matrix(0, nrow = T_test, ncol = d_test)
# y_test[1:length(y_test)] = 1:length(y_test)
# print(y_test)
#
# t_test = 3
#
# Zt_test = make_Zt(y = y_test, t = t_test, p = p_test)
# print(Zt_test)
#
# Z_test = make_Z(y = y_test, p = p_test)
# print(Z_test)
# is.matrix(Z_test)
# nrow(Z_test) == d_test*(T_test - p_test)
# ncol(Z_test) == d_test + d_test^2*p_test*(T_test - p_test)
#
# p_test = 2
#
# Zt_test = make_Zt(y = y_test, t = t_test, p = p_test)
# print(Zt_test)
# Z_test = make_Z(y = y_test, p = p_test)
# print(Z_test)
# is.matrix(Z_test)
# nrow(Z_test) == d_test*(T_test - p_test)
# ncol(Z_test) == d_test + d_test^2*p_test*(T_test - p_test)
### end ###

make_Gamma = function(y, p, Phip_init = NA){
    if(any(is.na(Phip_init))){
        Phip_init = matrix(0, nrow = ncol(y), ncol = ncol(y)*p)
    }
    T = nrow(y)
    d = ncol(y)
    betap = vec(Phip_init)
    out = matrix(0, nrow = d^2*p*(T-p), ncol = 1)
    out[1:(d^2*p), 1] = - Phip_init
    return(out)
}

### sample code ###
# T_test = 5
# d_test = 2
# p_test = 1
# Phip_init_test = matrix(1:(d_test^2*p_test), nrow = d_test, ncol = d_test*p_test)
#
# y_test = matrix(0, nrow = T_test, ncol = d_test)
# y_test[1:length(y_test)] = 1:length(y_test)
# print(y_test)
#
# test_Gamma = make_Gamma(y = y_test, p = p_test)
# print(test_Gamma)
#
# test_Gamma = make_Gamma(y = y_test, p = p_test, Phip_init = Phip_init_test)
# print(test_Gamma)
# is.matrix(test_Gamma)
# nrow(test_Gamma) == d_test^2*p_test*(T_test - p_test)
#
# p_test = 2
# Phip_init_test = matrix(1:d_test^2*p_test, nrow = d_test, ncol = d_test*p_test)
# test_Gamma = make_Gamma(y = y_test, p = p_test, Phip_init = Phip_init_test)
# print(test_Gamma)
# is.matrix(test_Gamma)
# nrow(test_Gamma) == d_test^2*p_test*(T_test - p_test)
### end ###

make_W = function(y, p){
    T = nrow(y)
    d = ncol(y)
    out = matrix(0, nrow = d^2*p*(T-p), ncol = d + d^2*p*(T-p))
    Id = diag(d^2*p)
    for(t in 1:(T-p)){
        out[(1 + d^2*p*(t-1)):(d^2*p + d^2*p*(t-1)), (d + 1 + (d^2*p*(t-1))):(d + d^2*p +  (d^2*p*(t-1)))] = -Id
    }
    for(t in 1:(T-p-1)){
        out[(d^2*p + 1 + d^2*p*(t-1)):(d^2*p + d^2*p + d^2*p*(t-1)), (d + 1 + (d^2*p*(t-1))):(d + d^2*p +  (d^2*p*(t-1)))] = Id
    }
    return(out)
}

### sample code ###
# T_test = 5
# d_test = 2
# p_test = 1
#
# y_test = matrix(0, nrow = T_test, ncol = d_test)
# y_test[1:length(y_test)] = 1:length(y_test)
# print(y_test)
#
# test_W = make_W(y = y_test, p = p_test)
# print(test_W)
# is.matrix(test_W)
# nrow(test_W) == d_test^2*p_test*(T_test - p_test)
# ncol(test_W) == d_test + d_test^2*p_test*(T_test - p_test)
#
# T_test = 5
# p_test = 2
#
# test_W = make_W(y = y_test, p = p_test)
# print(test_W)
# is.matrix(test_W)
# nrow(test_W) == d_test^2*p_test*(T_test - p_test)
# ncol(test_W) == d_test + d_test^2*p_test*(T_test - p_test)
### end ###

make_psi = function(y, p, Phip_init = NA){
    if(any(!is.na(Phip_init))){
        Phip_init = Phip_init
    }
    Y = make_Y(y = y, p = p)
    Gamma = make_Gamma(y = y, p = p, Phip_init = Phip_init)
    out = rbind(Y, Gamma)
    return(out)
}

### sample code ###
# T_test = 5
# d_test = 2
# p_test = 1
#
# y_test = matrix(0, nrow = T_test, ncol = d_test)
# y_test[1:length(y_test)] = 1:length(y_test)
# print(y_test)
#
# psi_test = make_psi(y = y_test, p = p_test)
# print(psi_test)
# is.matrix(psi_test)
# nrow(psi_test) == d_test*(T_test - p_test) + d_test^2*p_test*(T_test - p_test)
# ncol(psi_test) == 1
#
# p_test = 2
# psi_test = make_psi(y = y_test, p = p_test)
# print(psi_test)
# is.matrix(psi_test)
# nrow(psi_test) == d_test*(T_test - p_test) + d_test^2*p_test*(T_test - p_test)
# ncol(psi_test) == 1
### end ###

make_zeta = function(y, p){
    Z = make_Z(y, p)
    W = make_W(y, p)
    out = rbind(Z, W)
    return(out)
}

### sample code ###
# T_test = 5
# d_test = 2
# p_test = 1
#
# y_test = matrix(0, nrow = T_test, ncol = d_test)
# y_test[1:length(y_test)] = 1:length(y_test)
# print(y_test)
#
# zeta_test = make_zeta(y = y_test, p = p_test)
# print(zeta_test)
# is.matrix(zeta_test)
# nrow(zeta_test) == (d_test + d_test^2*p_test)*(T_test - p_test)
# ncol(zeta_test) == (d_test + d_test^2*p_test*(T_test - p_test))
### end ###

make_nu = function(y, p, beta){
    T = nrow(y)
    d = ncol(y)
    out = matrix(0, nrow = d, ncol = 1)
    out[1:d, 1] = beta[1:d, ]
    return(out)
}

### sample code ###
# T_test = 5
# d_test = 2
# p_test = 1
#
# y_test = matrix(0, nrow = T_test, ncol = d_test)
# y_test[1:length(y_test)] = 1:length(y_test)
# print(y_test)
#
# beta_test = matrix(1:(d_test + (d_test^2*p_test)*(T_test - p_test)), nrow = d_test + (d_test^2*p_test)*(T_test - p_test), ncol = 1)
# print(beta_test)
#
# nu_test = make_nu(y = y_test, p = p_test, beta = beta_test)
# is.matrix(nu_test)
# nrow(nu_test) == d_test
### end ###

make_Phi = function(y, p, beta){
    T = nrow(y)
    d = ncol(y)
    out = array(0, dim = c(d, d*p, T-p))
    for(t in 1:(T-p)){
        for(j in 1:d){
            out[, j, t] = beta[1:d + d + d*(j-1) + d*d*(t-1)]
        }
    }
    return(out)
}

### sample code ###
# T_test = 5
# d_test = 2
# p_test = 1
#
# y_test = matrix(0, nrow = T_test, ncol = d_test)
# y_test[1:length(y_test)] = 1:length(y_test)
# print(y_test)
#
# beta_test = matrix(1:(d_test + (d_test^2*p_test)*(T_test - p_test)), nrow = d_test + (d_test^2*p_test)*(T_test - p_test), ncol = 1)
# print(beta_test)
#
# Phi_test = make_Phi(y = y_test, p = p_test, beta = beta_test)
# print(Phi_test)
# dim(Phi_test) == c(d_test, d_test, T_test - p_test)
#
# p_test = 2
#
# beta_test = matrix(1:(d_test + (d_test^2*p_test)*(T_test - p_test)), nrow = d_test + (d_test^2*p_test)*(T_test - p_test), ncol = 1)
# print(beta_test)
#
# Phi_test = make_Phi(y = y_test, p = p_test, beta = beta_test)
# print(Phi_test)
# dim(Phi_test) == c(d_test, d_test*p_test, T_test - p_test)
### end ###

make_beta = function(y, p, nu, Phi){
    T = nrow(y)
    d = ncol(y)
    out = matrix(0, nrow = d + d^2*p*(T-p), ncol = 1)
    out[1:d, 1] = nu
    for(t in 1:(T-p)){
        for(j in 1:d){
            out[1:d + d*(j-1) + (T-p)*(t-1) + d, 1] = Phi[, j, t]
        }
    }
    return(out)
}

### sample code ###
# T_test = 5
# d_test = 2
# p_test = 1
#
# y_test = matrix(0, nrow = T_test, ncol = d_test)
# y_test[1:length(y_test)] = 1:length(y_test)
# print(y_test)
#
# nu_test = matrix(1:d_test, nrow = d_test, ncol = 1)
# Phi_test = array(1:(d_test^2*p_test*(T_test - p_test)) + d_test, c(d_test, d_test*p_test, T_test - p_test))
#
# beta_test = make_beta(y = y_test, p = p_test, nu = nu_test, Phi = Phi_test)
# print(beta_test)
# is.matrix(beta_test)
# nrow(beta_test) == d_test^2*p_test*(T_test - p_test) +  d_test
# ncol(beta_test) == 1
### end ###
