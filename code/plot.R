plot_Phi = function(Phi, time = 1:dim(Phi)[3], time_label = time, ylim = NA){
    T = dim(Phi)[3]
    d = dim(Phi)[1]

    out = matrix(0, nrow = T, ncol = d*d)

    for(i in 1:d){
        for(j in 1:d){
            for(t in 1:T){
                out[t, i + d*(j-1)] = Phi[i, j, t]
            }
        }
    }

    par(mfrow = c(d, d))
    if(is.na(sum(ylim))){
        ylim = c(-max(abs(out)), max(abs(out)))
    }
    for(s in 1:(d*d)){
        plot(time, out[, s], type = "n", ylim = ylim, ylab = "", xaxt = "n", xlab = "time")
        if(any(out[, s] == 0)){
            time0 = time[out[, s] == 0]
            left = time0 - 1/2 - 1/10
            right = time0 + 1/2 + 1/10
            for(t0 in 1:length(time0)){
                polygon(c(left[t0], right[t0], right[t0], left[t0]), c(2*ylim[1], 2*ylim[1], 2*ylim[2], 2*ylim[2]), col = "gray", border = NA)
            }
        }
        axis(side = 1, at = time, labels = time_label, cex.axis = 1)
        grid()
        abline(h = 0)
        points(time, out[, s], type = "l", ylim = ylim, lwd = 1)
    }
    par(mfrow = c(1, 1))
}


plot_multi_Phi = function(Phi, time = 1:dim(Phi[[1]])[3], ylim = NA){
    T = dim(Phi[[1]])[3]
    d = dim(Phi[[1]])[1]

    out_list = list()
    n_Phi = length(Phi)

    Phi_max = 0

    for(i_Phi in 1:n_Phi){
        out = matrix(0, nrow = T, ncol = d*d)
        for(i in 1:d){
            for(j in 1:d){
                for(t in 1:T){
                    out[t, i + d*(j-1)] = Phi[[i_Phi]][i, j, t]
                }
            }
        }
        out_list[[i_Phi]] = out
        out_max = max(abs(out))
        if(Phi_max < out_max){
            Phi_max = out_max
        }
    }


    par(mfrow = c(d, d))
    if(is.na(sum(ylim))){
        ylim = c(-Phi_max, Phi_max)
    }
    for(s in 1:(d*d)){
        plot(time, out_list[[1]][, s], type = "n", ylim = ylim, ylab = "", xaxt = "n", xlab = "time")
        # if(any(out[, s] == 0)){
        #     time0 = time[out[, s] == 0]
        #     left = time0 - 1/2 - 1/10
        #     right = time0 + 1/2 + 1/10
        #     for(t0 in 1:length(time0)){
        #         polygon(c(left[t0], right[t0], right[t0], left[t0]), c(2*ylim[1], 2*ylim[1], 2*ylim[2], 2*ylim[2]), col = "gray", border = NA)
        #     }
        # }
        axis(side = 1, at = time, label = time, cex.axis = 1)
        grid()
        abline(h = 0)
        for(i_Phi in 1:n_Phi){
            points(time, out_list[[i_Phi]][, s], type = "l", ylim = ylim, lwd = 1, col = i_Phi)
        }
    }
    par(mfrow = c(1, 1))
}
