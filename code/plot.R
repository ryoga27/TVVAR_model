plot_Phi = function(Phi, time = 1:dim(Phi)[3]){
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
    ylim = c(-max(abs(out)), max(abs(out)))
    for(s in 1:(d*d)){
        s = 1
        plot(time, out[, s], type = "n", ylim = ylim, ylab = "", xaxt = "n", xlab = "time")
        time0 = time[out[, s] == 0]
        left = time0 - 1/2 - 1/10
        right = time0 + 1/2 + 1/10
        for(t0 in 1:length(time0)){
            polygon(c(left[t0], right[t0], right[t0], left[t0]), c(2*ylim[1], 2*ylim[1], 2*ylim[2], 2*ylim[2]), col = "gray", border = NA)
        }
        axis(side = 1, at = time, label = time, cex.axis = 1)
        grid()
        abline(h = 0)
        points(time, out[, s], type = "l", ylim = ylim, lwd = 1)
    }
    par(mfrow = c(1, 1))
}
