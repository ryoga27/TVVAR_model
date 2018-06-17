coordinate_descent_algorithm = function(
    Y, X, T, d, p,
    lambda = 0.000001,
    eps = 10^-7,
    iter_max = 100,
    penalty = "lasso",
    display = TRUE
){
    m = T-p

    nu_list = list()
    beta_list = list()

    Y_mean = mean(Y)
    y = Y - Y_mean

    y_mat = matrix(0, nrow = m, ncol = d)
    for(j in 1:d){
        for(t in 1:m){
            y_mat[t, j] = y[j + d*(t-1)]
        }
    }

    y_mean = matrix(0, nrow = d, ncol = 1)
    for(j in 1:d){
        y_mean[j, ] = mean(y_mat[, j])
    }

    x = X[, -c(1:d)]
    x = scale(x)
    X_mean = attr(x, "scaled:center")
    X_sd   = attr(x, "scaled:scale")

    n = nrow(x)
    k = ncol(x)

    nu = y_mean
    nu_list[[1]] = nu

    mu = matrix(0, nrow = n, ncol = 1)
    mu[1:(d*(T-p)), ] = rep(nu, m)

    # beta = rep(0, length = k)
    beta = solve(t(x)%*%x)%*%t(x)%*%(y - mu)

    obj = rep(0, iter_max)

    obj[1] = (1/n)*sum((y - mu - x%*%beta)^2) + penalty_function(beta, lambda, penalty = penalty)

    for(s in 1:iter_max){
        if(display == TRUE){
            cat("iteratoin times:", s, "\n")
            cat("\t", "calculate beta...", "\n")
            threshold = 0
        }
        for(j in 1:k){
            if(display == TRUE){
                if(threshold <= j/k){
                    cat("\t", round(j/k, 2)*100, "%")
                    threshold = threshold + 0.1
                }
            }
            r = y - mu - x[, -j]%*%beta[-j]
            z = (1/n)*x[, j]%*%r
            beta[j] = penalty_solution(z, lambda = lambda, penalty = penalty)
        }
        if(display == TRUE){
            cat("\n")
            cat("\t", "beta:", beta[1:5], "...", "\n")
        }

        beta_list[[s+1]] = beta

        if(display == TRUE){
            cat("\t", "calculate nu...", "\n")
        }
        nu = rep(0, length = d)
        val = 0
        for(j in 1:d){
            val = 0
            for(iota in 1:m){
                val = val + y[j + d*(iota-1), ] - x[j + d*(iota-1), ]%*%beta
            }
            nu[j] = (1/m)*val
        }
        if(display == TRUE){
            cat("\t", "nu:", nu, "\n")
        }

        nu_list[[s+1]] = nu
        beta_list[[s+1]] = beta

        obj[s + 1] = (1/n)*sum((y - mu - x%*%beta)^2) + penalty_function(beta, lambda, penalty = penalty)

        if(display == TRUE){
            cat("\t", "object function:", obj[s+1], "\n")
        }
        mu = matrix(0, nrow = n, ncol = 1)
        mu[1:(d*(T-p)), ] = rep(nu, m)

        convergence = (abs(obj[s + 1] - obj[s]) < eps)
        if(convergence == TRUE){
            if(display == TRUE){
                cat("convergenced", "\n")
            }
            process_list = list(
                nu = nu_list,
                beta = beta_list,
                obj = obj[1:(s+1)],
                times_iter = s
            )
            break
        }
    }

    nu_hat = nu + Y_mean
    beta_hat = beta/X_sd

    coefficients = matrix(c(nu_hat, beta_hat), nrow = d + k, ncol = 1)
    args_list = list(
        coefficients = coefficients,
        lambda = lambda
    )
}


### sample code ###
# source("./reshape.R")
# source("./predict.R")
# source("./plot.R")
# source("./penalties.R")
#
# data_path = "../data/sample_data.csv"
# data_set = read.csv(data_path)
# n = nrow(data_set)
# y = log(data_set[-n, ]) - log(data_set[-1, ])
# y = as.matrix(y)
#
# p = 1
# T = nrow(y)
# d = ncol(y)
# zeta = make_zeta(y = y, p = p)
# psi = make_psi(y = y, p = p)
# lambda = 0.00001
# result = coordinate_descent_algorithm(Y = psi, X = zeta, T = T, d = d, p = p, lambda = lambda, eps = 10^-6)
# save(list = ls(), file="../results/result_test.Rdata")
# load("../results/result_test.Rdata")
# pdf("../figures/figure_test.pdf", width = 15, height = 10)
#     plot_Phi(make_Phi(y = y, p = 1, beta = result$coef))
# dev.off()
#
# penalty = "scad"
# result = coordinate_descent_algorithm(Y = psi, X = zeta, T = T, d = d, p = p, lambda = lambda, penalty = penalty, eps = 10^-6)
# path_result = paste0("../results/", "sample_data", "_", penalty,"_", lambda, ".Rdata")
# save(list = ls(), file=path_result)
# path_figrue = paste0("../figures/", "sample_data", "_", penalty,"_", lambda, ".pdf")
# pdf(path_figrue, width = 15, height = 10)
#     plot_Phi(make_Phi(y = y, p = 1, beta = result$coef))
# dev.off()
# scp ryoga@sl:Files/TVVAR_model/figures/*.pdf .
