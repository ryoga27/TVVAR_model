plot_Phi = function(result){
    Phi = result$Phi
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
        plot(1:T, out[, s], type = "n", ylim = ylim, ylab = "", xlab = "time")
        grid()
        points(1:T, out[, s], type = "l", ylim = ylim)
        abline(h = 0)
    }
    par(mfrow = c(1, 1))
}
