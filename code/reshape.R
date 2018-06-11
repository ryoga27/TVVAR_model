if(require("glmnet") == FALSE){
    install.packages("glmnet")
}
library(glmnet)

make_yt = function(y, t){
    T = nrow(y)
    d = ncol(y)
    out = rep(0, length = d)
    out = as.matrix(y[t, 1:d], ncol = 1)
    return(out)
}

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

make_Zt = function(y, t, p){
    T = nrow(y)
    d = ncol(y)
    out = rep(0, d*p)
    s = 1
    for(l in 1:p){
        out[1:d + d*(s-1)] = make_yt(y, t-l)
        s = s + 1
    }

    out = matrix(out, nrow = d*p, ncol = 1)
    return(out)
}

make_Z = function(y, p){
    T = nrow(y)
    d = ncol(y)
    out = matrix(0, nrow = d*(T-p), ncol = d + d^2*p*(T-p))
    Id = diag(d)
    for(s in 1:(T-p)){
        out[1:d + d*(s-1), 1:d] = Id
    }
    s = 1
    for(t in (p+1):T){
        out[1:d + d*(s-1), 1:(d^2*p) + d + (d^2*p)*(s-1)] = kronecker(t(make_Zt(y, t = t, p = p)), Id)
        s = s + 1
    }
    return(out)
}

make_Gamma = function(y, p){
    T = nrow(y)
    d = ncol(y)
    out = matrix(0, nrow = d^2*p*(T-p), ncol = 1)
    return(out)
}

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

make_nu = function(beta, y, p){
    T = nrow(y)
    d = ncol(y)
    out = beta[1:d, ]
    return(out)
}

make_Phi = function(beta, y, p){
    T = nrow(y)
    d = ncol(y)
    out = array(0, dim = c(d, d*p, T-p))
    for(t in 1:(T-p)){
        for(j in 1:d){
            out[, j, t] = beta[(d + 1 + d*(j-1) + d*d*(t-1)):(d + d + d*(j-1) + d*d*(t-1))]
        }
    }
    return(out)
}

make_beta = function(nu, Phi, y, p){
    T = nrow(y)
    d = ncol(y)
    out = matrix(0, nrow = d + d^2*p*(T-p), ncol = 1)
    out[1:d] = nu
    for(t in 1:(T-p)){
        for(j in 1:d){
            out[1:d + d] = Phi[, j, t]
        }
    }
    return(out)
}
