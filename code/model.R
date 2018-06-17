# if(require("glmnet") == FALSE){
#     install.packages("glmnet")
# }
# library(glmnet)

TVVAR_model = function(y, p){
    if(!is.matrix(y)){
        stop("y must be matrix")
    }
    zeta = make_zeta(y = y, p = p)
    psi = make_psi(y = y, p = p)
    beta = solve(t(zeta)%*%zeta)%*%t(zeta)%*%psi
    nu = make_nu(y = y, p = p, beta = beta)
    Phi = make_Phi(y = y, p = p, beta = beta)
    out = list(
        nu = nu,
        Phi = Phi
    )
    return(out)
}

### sample code ###
# source("./reshape.R")
# source("./plot.R")
#
# data_path = "../data/sample_data.csv"
# data_set = read.csv(data_path)
# n = nrow(data_set)
# y = log(data_set[-n, ]) - log(data_set[-1, ])
# y = as.matrix(y)
#
# result = TVVAR_model(y = y, p = 1)
# plot_Phi(result$Phi)
### end ###

# TVVAR_model_lasso = function(y, p, lambda){
#     if(!is.matrix(y)){
#         stop("y must be matrix")
#     }
#     T = nrow(y)
#     d = ncol(y)
#     zeta = make_zeta(y = y, p = p)
#     psi = make_psi(y = y, p = p)
#     penalty_factor = rep(1, ncol(zeta))
#     penalty_factor[1:d] = 0
#     beta = coef(glmnet(x = zeta, y = psi, intercept = FALSE, alpha = 1, lambda = lambda, penalty.factor = penalty_factor))
#     beta = matrix(beta[-1], nrow = d + d^2*(T-p), ncol = 1)
#
#     nu = make_nu(y = y, p = p, beta = beta)
#     Phi = make_Phi(y = y, p = p, beta)
#     out = list(
#         nu = nu,
#         Phi = Phi
#     )
#     return(out)
# }

### sample code ###
source("./reshape.R")
source("./plot.R")

# data_path = "../data/sample_data.csv"
# data_set = read.csv(data_path)
# n = nrow(data_set)
# y = log(data_set[-1, ]) - log(data_set[-n, ])
# y = as.matrix(y)
#
# lambda = 0.00001
# result_lasso = TVVAR_model_lasso(y = y, p = 1, lambda = lambda)
# plot_Phi(result_lasso$Phi)
### end ###

# cv_TVVAR_model_lasso = function(y, p, lambdas){
#     T = nrow(y)
#     d = ncol(y)
#     L = length(lambdas)
#     mse  = rep(0, L)
#     out = list()
#     for(l in 1:L){
#         cat("lambda", lambdas[l], round(l/L, 5)*100, "%", "\n")
#         val = 0
#         for(t in 1:(T - trunc(T/2))){
#             if(t%%10 == 0){
#                 cat("\t", "time", t, round(t/(T - trunc(T/2)), 5)*100, "%", "\n")
#             }
#             y_train = y[1:(trunc(T/2) + t), ]
#             result = TVVAR_model_lasso(y_train, p = p, lambda = lambdas[l])
#             y_pred = predict(y_train, result, t = (trunc(T/2) + t - p))
#             epsilon = y[trunc(T/2) + 1, ] - y_pred
#             val = val + epsilon^2
#         }
#         mse[l] = (1/(T - trunc(T/2)))*t(val)%*%val
#     }
#     lambda_opt = lambdas[mse == min(mse)]
#     out$lambdas = lambdas
#     out$mse = mse
#     out$lambda_opt = lambda_opt
#     result = TVVAR_model_lasso(y_train, p = p, lambda = lambda_opt)
#     out$result = result
#     return(out)
# }

### sample code ###
# source("./reshape.R")
# source("./predict.R")
# source("./plot.R")
#
# data_path = "../data/sample_data.csv"
# data_set = read.csv(data_path)
# n = nrow(data_set)
# y = log(data_set[-n, ]) - log(data_set[-1, ])
# y = as.matrix(y)
#
# lambdas = 10^c(-10:10)
# result_cv_lasso = cv_TVVAR_model_lasso(y = y, p = 1, lambdas = lambdas)
# result_cv_lasso$result$nu
# plot_Phi(result_cv_lasso$result$Phi)
### end ###

sparse_TVVAR_model = function(y, p, lambda, penalty = "lasso"){
    if(!is.matrix(y)){
        stop("y must be matrix")
    }
    T = nrow(y)
    d = ncol(y)
    zeta = make_zeta(y = y, p = p)
    psi = make_psi(y = y, p = p)
    result = coordinate_descent_algorithm(Y = psi, X = zeta, T = T, d = d, p = p, lambda = lambda, eps = 10^-6, penalty = penalty)
    beta = result$coefficients
    beta = matrix(beta, nrow = d + d^2*(T-p), ncol = 1)

    nu = make_nu(y = y, p = p, beta = beta)
    Phi = make_Phi(y = y, p = p, beta)
    out = list(
        nu = nu,
        Phi = Phi
    )
    return(out)
}


### sample code ###
# source("./reshape.R")
# source("./penalties.R")
# source("./_coordinate_descent_algorithm.R")
# source("./plot.R")
#
# data_path = "../data/sample_data.csv"
# data_set = read.csv(data_path)
# n = nrow(data_set)
# y = log(data_set[-n, ]) - log(data_set[-1, ])
# y = as.matrix(y)
#
# lambda = 0.00001
# result_lasso = sparse_TVVAR_model(y = y, p = 1, lambda = lambda)
# plot_Phi(result_lasso$Phi)
### end ###

cv_sparse_TVVAR_model = function(y, p, lambdas, penalty = "lasso", display = TRUE, save = FALSE){
    T = nrow(y)
    d = ncol(y)
    L = length(lambdas)
    mse  = rep(0, L)
    out = list()
    for(l in 1:L){
        if(display == TRUE){
        }
        val = 0
        if(display == TRUE){
            threshold = 0
        }
        for(t in 1:(T - trunc(T/2))){
            if(display == TRUE){
                cat("lambda:", lambdas[l], "(", round(l/L, 2)*100, "%", ")", "\n")
                cat("rolling...", t, "(", round(t/(T - trunc(T/2)), 3)*100, "%", ")", "\n")
            }
            # if(display == TRUE){
            #     if(threshold <= t/(T - trunc(T/2))){
            #         cat("\t", round(t/(T - trunc(T/2)), 2)*100, "%")
            #         threshold = threshold + 0.1
            #     }
            # }
            y_train = y[1:(trunc(T/2) + t), ]
            result = sparse_TVVAR_model(y_train, p = p, lambda = lambdas[l], penalty = penalty)
            if(save == TRUE){
                if(t == (T - trunc(T/2))){
                    path_result = paste0("../results/", penalty, "_", lambdas[l], ".Rdata")
                    save(list = ls(all = TRUE), file = path_result)
                }
            }
            y_pred = predict(y_train, result, t = (trunc(T/2) + t - p))
            epsilon = y[trunc(T/2) + 1, ] - y_pred
            val = val + epsilon^2
        }
        mse[l] = (1/(T - trunc(T/2)))*t(val)%*%val
    }
    lambda_opt = lambdas[mse == min(mse)]
    if(length(mse) != 1){
        lambda_opt = min(lambda_opt)
    }
    out$lambdas = lambdas
    out$mse = mse
    out$lambda_opt = lambda_opt
    result = sparse_TVVAR_model(y_train, p = p, lambda = lambda_opt, penalty = penalty)
    out$result = result
    return(out)
}

# T_test = 5
# d_test = 2
# p_test = 1
# y = matrix(0, nrow = T_test, ncol = d_test)


# penalty = "lasso"
#
# lambdas = c(10^c(-10:10))
# result_lasso = cv_sparse_TVVAR_model(y = y, p = 1, lambda = lambdas, penalty = penalty, save = TRUE)
