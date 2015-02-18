#!/usr/bin/env Rscript

library('getopt')
opt = getopt(matrix(c(
	'envfile', 'e', 1, 'character'
),ncol=4,byrow=T))

if( is.null(opt$envfile) ) opt$envfile = ''

# adapted from Serdar's script
# This file creates motif pngs to import into Drupal
# create motif pngs
plot.motifs <- function(env){
	dir.create('cmonkey_motifs')
	library("seqLogo")
	for(bcl in 1:env$k.clust){
		bicluster.plot = env$plotClust(bcl, dont.plot=T)
		if(!is.null(bicluster.plot))
			if (bicluster.plot$nrow >1){
				if(!is.na(bicluster.plot$e.val[[1]])){
					pssms.1 = bicluster.plot$`upstream meme`$motif.out$pssms[[1]]
					e.val.1 = bicluster.plot$`upstream meme`$e.val[[1]]
					pwm.1 = makePWM(t(pssms.1))
				}
				if(!is.na(bicluster.plot$e.val[[1]])){
					pssms.2 = bicluster.plot$`upstream meme`$motif.out$pssms[[2]]
					e.val.2 = bicluster.plot$`upstream meme`$e.val[[2]]
					pwm.2 = makePWM(t(pssms.2))
				}

				mySeqLogo = seqLogo::seqLogo
				bad = sapply( body(mySeqLogo), "==", "grid.newpage()") | sapply( body(mySeqLogo), "==", "par(ask = FALSE)")
				body(mySeqLogo)[bad] = NULL
				norm = function(x) scale(x, center=FALSE, scale=colSums(x))

				matrix.1 = c(pwm.1, pwm.2)
				mot.ind = c("green", "red")

				for(i in 1:2){
					cat('bicluster',bcl,'motif',i,'\n')
					pwm = matrix.1[[i]]
					png.name = paste("cmonkey", "_motifs/", sprintf("motif_%04d", bcl), "_", i, ".png", sep="")
					png(filename = png.name, bg = "transparent")
					pushViewport(viewport(x=0.50, y=0.50))
					mySeqLogo(pwm, xfontsize=12, yaxis=T)
					e.value = c(e.val.1, e.val.2)
					if(!is.na(e.value[[i]])){
						grid.text(sprintf("PSSM #%d; e=%.3g", i, e.value[[i]]), x=0.5, y=1, hjust=0.5, vjust=1, gp=gpar(fontsize=20,col= mot.ind[[i]]))
					}
					popViewport()
					dev.off()
				}
		}
	}
}

load(opt$envfile)
if(!exists('env')){
	cat('problem, no env object found!\n\n\n')
}
plot.motifs(env)
