## function to plot distributions of dissimilarities for both observed and simulated data
## x is the result of the 'effSize.calc' function
## use groups to select a subset of results from the object; can be numeric, indicating rows, or character, indicating level in 'trt'
## to do: added ability to plot theta, I, S_total, S_sample, J to plot.dist (use the params argument)
plot.dist <- function(x, groups=NULL, params=NULL){
        par(mar=c(5, 5, 1, 2))
        if(is.null(groups)) groups <- 1:nrow(x)
        for(i in groups) {
                if(is.numeric(i)) main <- paste('group =', as.character(x[i, 1])) else main <- paste('group =', i)
                x.obs <- attr(x, 'observed distances')[[i]]
                x.dens <- attr(x, 'simulations density')[[i]]
                if(length(x.obs) > 1) {
                        dens.max <- max(c(max(density(x.obs)$y), max(x.dens[['y']])))
                        if(attr(x, 'distance index') == 'bray') {
                                plot(seq(0, 1.2, length.out=10), seq(0, 1.2 * dens.max, length.out=10), type='n', 
                                     xlab='Bray-Curtis dissimilarity', ylab='Density', main=main, 
                                     cex.lab=1.25, cex.axis=1.25, bty='n', xaxs='i', xaxp=c(0, 1, 4))
                        }
                        if(attr(x, 'distance index') == 'sorensen') {
                                plot(seq(0, 1.2, length.out=10), seq(0, 1.2 * dens.max, length.out=10), type='n', 
                                     xlab='Sorensen dissimilarity', ylab='Density', main=main, 
                                     cex.lab=1.25, cex.axis=1.25, bty='n', xaxs='i', xaxp=c(0, 1, 4))
                        }
                        if(attr(x, 'distance index') == 'raup') {
                                plot(seq(0, 1.2, length.out=10), seq(0, 1.2 * dens.max, length.out=10), type='n', 
                                     xlab='Raup-Crick dissimilarity', ylab='Density', main=main, 
                                     cex.lab=1.25, cex.axis=1.25, bty='n', xaxs='i', xaxp=c(0, 1, 4))
                        }
                        for(j in 1:nrow(x.dens[['x']])) lines(x.dens[['x']][j, ], x.dens[['y']][j, ], lwd=0.1, col='grey')   # work on this to plot simulation densities
                        lines(density(x.obs), col='black', lwd=2)
                        legend('topleft', legend=c('observed', 'expected (neutral)'), lty='solid', lwd=c(2, 1), col=c('black', 'grey'), bty='n')
                } else {
                        dens.max <- max(density(x.dens$x)$y)
                        if(attr(x, 'distance index') == 'bray') {
                                plot(seq(0, 1.2, length.out=10), seq(0, 1.2 * dens.max, length.out=10), type='n', 
                                     xlab='Bray-Curtis dissimilarity', ylab='Density', main=main, 
                                     cex.lab=1.25, cex.axis=1.25, bty='n', xaxs='i', xaxp=c(0, 1, 4))
                        }
                        if(attr(x, 'distance index') == 'sorensen') {
                                plot(seq(0, 1.2, length.out=10), seq(0, 1.2 * dens.max, length.out=10), type='n', 
                                     xlab='Sorensen dissimilarity', ylab='Density', main=main, 
                                     sub='Only two communities, density plot of one simulated comparison',
                                     cex.lab=1.25, cex.axis=1.25, bty='n', xaxs='i', xaxp=c(0, 1, 4))
                        }
                        lines(density(x.dens$x), col='black', lty='dashed', lwd=2)
                        abline(v=x.obs, col='black', lty='solid', lwd=2)
                        legend('topleft', legend=c('observed', 'expected (neutral)'), lty=c('solid', 'dashed'), lwd=2, col='black', bty='n')     
                }
        }
}
