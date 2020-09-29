## This script contains all the functions that are needed to process data, calculate off-rates and create figures for the DIVORSEQ method

## This functions extract the summits and makes a new GRanges object with the area around the summits. By default +/- 50 bp (101 in total) 
make.summits.object <- function(narrowPeak, plus=50, minus=50) {
	## extract the summits
	summits <- narrowPeak
	start(summits) <- start(narrowPeak) + narrowPeak$peak - minus
	end(summits) <- start(narrowPeak) + narrowPeak$peak + plus
	summits$peak <- minus   # this way you can find the summit position easily again
	return(summits)
}

## 

#########################################
##                                                                    		#
##         Reading in bw files  / calculate coverage            	#
##                                                                    		#
#########################################

## Convert the bigWig file to a list containing the coverage per chr, so we can make the coverage per genes out of that
Cover.bw <- function(bw) {			
  result <- list()													#Makes a list called results where the results will be stored
  for(chr in levels(seqnames(bw))) {								#For each chromosome it:
    subs = bw[seqnames(bw)== chr,]									#takes a subset of the data of only the chromosome
    cover = coverage(subs, weight= values(subs)$score, method="hash" )		#calculates a per-basepair coverage (Rle)
	result[[chr]] <- as.numeric(cover[[chr]])						#convert to numeric vector (non-Rle)
	}
  result
}


## calculate the coverage under a specific peak
## when using this function with sequences (i.e. a matrix full of A,C,G T), set extraLevelSeqMat=TRUE
align.data.peaks2 <- function (features, data,matrix.size=NULL, chrs = NULL,extraLevelSeqMat=FALSE) {
	if( is.null(chrs)) chrs <- as.vector(unique(seqnames(features))) # use this rather than levels() to keep the order chrI - chrXVI intact rather than alphabetical
	
	if(is.null(names(features))) stop("This function uses the names of the features, please make sure the the features are named!")
	
	feat.size <- max(width(features)) 
	if(!is.null(matrix.size)) {
		if(matrix.size %% 2 ==0) stop("please make sure the matrix.size is odd, otherwise the function will crash...")
		feat.size <- matrix.size # if the matrix size is specified use that, otherwise use the biggest size
	}
	
	## make a matrix that can be filled with the data points
	m <- matrix(NA, nrow=length(features), ncol=feat.size)
	colnames(m) <- -((feat.size -1)/2):((feat.size -1)/2)# for now assume symmetry
    rownames(m) <- names(features)
	
	## where the magic happens
	for(chro in chrs) { 
		peaks <- features[ seqnames(features)==chro ]
		for (g in names(peaks)) { 
			from <- start(peaks[g])
			to <- end(peaks[g]) 
			
			## fill the matrix. the left point is the (negative) position from the peak, the right most coordinate is the width minus the peak
			columns <- which(colnames(m) == as.character(-peaks[g]$peak)) : which(colnames(m) == (width(peaks[g]) -peaks[g]$peak -1)) 
			if(extraLevelSeqMat){
			m[ g,columns ]  <- data[[chro]][[1]][ from:to ] # the sequence matrix has a list per chromosomes and chrom is also a list. So we need to select that
			} else {
			m[ g,columns ]  <- data[[chro]][ from:to ] 
			}
		}
	}	
	m
}

## calculate the total coverage per peak per time points and put in a simple data frame
calc.coverage.per.peak <- function(coverage.list, sample.names, func, rm.na=F) {
	## check if all the coverage matrices have the same rownames!
	rownames.list <- list()
	for(i in names(coverage.list)) rownames.list[[i]] <- rownames(coverage.list[[i]]) # put all the rownames of each coverage matrix in a list of vectors
	rownames.check <-	sapply( rownames.list, FUN = identical, rownames.list[[1]] ) # compare all to the first item of the list i.e. the first sample
	if(! all(rownames.check)) stop("Not all rownames are identical. Please check")
	
	## create the dataframe that will be filled with the summary statistic of choice per peak
	cov.peak <- matrix(NA,nrow=nrow(coverage.list[[1]]), ncol=length(sample.names)) #first make an empty matrix to fill in
	colnames(cov.peak) <- sample.names; rownames(cov.peak) <-  rownames.list[[1]] 
	
	## fill in the matrix with the summary statistic of choice
	
	for( i in sample.names) cov.peak[,i] <- apply(coverage.list[[i]], 1,func)
	if(any( is.na(cov.peak))){
		na.peaks <- apply(is.na(cov.peak),1, FUN= any)
		na.peaks.names <- rownames(cov.peak)[na.peaks]
		if(!rm.na) warning("There were NA's found in the following peaks: \n",  na.peaks.names,"\nYou can re-run with rm.na=T, but this means the values for these peaks may not be comparable! ")
	}
	if(rm.na) for( i in sample.names) cov.peak[,i] <- apply(coverage.list[[i]], 1,func, na.rm=rm.na)
	
	## to fit the exponential decay curves it is best to have the samples as rows and each this you want to fit (peak) as columns
	cov.peak.df <-data.frame(cov.peak)
	if(exists("na.peaks")) cov.peak.df <- cbind(na.peaks,cov.peak.df) # add a row with a logical that shows which peak contained NA's
	
	## for downstream analysis is is better to transpose the DF (to have the peaks in the columns)
	## but in that case the logical vector (na.peaks) is converted to numeric (0's and 1's) so for now
	## the peaks are in the rows. They still need to be transposed for curve fitting!
	
	return(cov.peak.df) 
}



#########################################
##                                                                    		#
##        					Finding motifs       				     	#
##                                                                    		#
#########################################

## This function allows for selecting the strongest motif within a certain distance from each other 
## if two motif are closer together than the distance cutoff, take the highest motif instead of the closest
closest.motif.hi <- function(fw, rv,size=201, dist.cutoff=1){
	stopifnot(identical(names(fw), names(rv)))
	
	MAT <- matrix(NA,ncol=2,nrow=length(fw)); colnames(MAT) <- c("fw","rev")
	rownames(MAT) <- names(fw)
	
	size.conv <- (-(size -1)/2) : ((size -1)/2)
	size.mid <- (size-1)/2+1		# what is the midpoints of the motif?
	
	## function to find the closest motif (that can be applied to both the forward and the reverse motif
	find.high <- function(mat, list.motif, column,dist.cutoff){
		for(i in names(list.motif)){
			if(length(list.motif[[i]]) > 0){
				midpoints <- (start(list.motif[[i]])+end(list.motif[[i]]) )/2
				close.to.mid <- which( abs((midpoints - size.mid)) - min(abs(midpoints - size.mid)) < dist.cutoff) # which motif is closest to the midpoint?
				if(length( close.to.mid) >1) { ## sometimes 2 motifs are equally close
					high.motif <- which(mcols(list.motif[[i]])$score[close.to.mid] == max(mcols(list.motif[[i]])$score[close.to.mid])) ## which of the two is the highest
					close.to.mid <- close.to.mid[high.motif]
					## if the motifs have the same score, just pick the one that is actually closest:
					if(length(close.to.mid) >1) close.to.mid <- close.to.mid[which(midpoints[close.to.mid] == min(midpoints[close.to.mid]))]
				}
				mat[i,column] <- size.conv[midpoints[close.to.mid]]
			}
		}
		return(mat)
	}
	
	MAT <- find.high(mat=MAT, list.motif=fw,column="fw",dist.cutoff=dist.cutoff)
	MAT <- find.high(mat=MAT, list.motif=rv,column="rev",dist.cutoff=dist.cutoff)
	
	## which ones are closer then the distance cutoff? Find them and only keep the one with the highest score!
	equidist <- which(abs(abs(MAT[,"fw"]) - abs(MAT[,"rev"])) < dist.cutoff)

	if(length(equidist) > 0){
		MAT.equi <- data.frame(MAT[names(equidist),],closest=NA) # make a separate df which is filled in 
		for(i in names(equidist)){
			## get the scores and distance for the forward motifs
			dist.fw <- (( start(fw[[i]]) + end(fw[[i]]) )/2 ) - size.mid
			scores.fw <- mcols(fw[[i]])$score
			
			## get the scores and distances for the reverse motifs
			dist.rv <- (( start(rv[[i]]) + end(rv[[i]]) )/2 ) - size.mid
			scores.rv <- mcols(rv[[i]])$score
			
			## which fw and rev motif are equally closest? and which of these two has the highest motif?
			close.fw.rv <- c(which(dist.fw == MAT[i,"fw"]), which(dist.rv == MAT[i,"rev"])) ; names(close.fw.rv) <- c("fw","rev")
			scores.fw.rv <- c(scores.fw[close.fw.rv["fw"]],scores.rv[close.fw.rv["rev"]]); names(scores.fw.rv) <- c("fw","rev")
			highest.fw.rv <- which(scores.fw.rv == max(scores.fw.rv))
			
			MAT.equi[i,"closest"] <- names(highest.fw.rv)
	
		}
		## fill NA's in the MAT, just so we can fill it back in later
		MAT[names(equidist),] <- NA
	}
	## not sure why this has to be so complex, but I leave it for now for compatibility reasons
	return.min <- function(x){
		if(length(unique(x))==1){
			xx <- unique(x)		
		} else {
			m <- min(abs(x),na.rm=T)
			xx <- which(abs(x) == m)
		}
	return(xx)
	}
	
	##find the closest motif
	min.fw.rev <- apply(MAT,1,return.min)
	MAT <- data.frame(MAT,closest=c("fw","rev")[min.fw.rev],stringsAsFactors=FALSE)
	
	## fill them back in!
	if(length(equidist) >0) MAT[names(equidist),] <- MAT.equi
	
	return(MAT)
}


## a function that takes a list of peak names and the forward/reverse found motifs and aggregates them!
find.double.motifs <- function(dub.mot,fw,rv){
  DNA <- c(A="T",C="G",G="C",T="A") # needed to convert the reverse base
  double.motif.base <- list()
  
  ## loop over the peaks with more than 1 motif
  for(i in dub.mot) {
    temp.list <- list()
    FW <- fw[[i]]
    if(length(FW) > 0) {
	## you can extract the sequence that was used with subject(XSTRINGVIEWSOBJECT)
       temp.list[["fw"]]  <- unlist(strsplit(as.character(subject(FW)[start(FW)-1]),split="")) # we want the base just upstream of the motif
	   names(temp.list[["fw"]]) <- round(mcols(FW)$score,3)
    } else {
      temp.list[["fw"]] <- NA
    }
    RV <-  rv[[i]]
    if(length(RV) > 0) { # for the reverse use the "DNA" vector to take the complementary base
       temp.list[["rev"]]  <- DNA[unlist(strsplit(as.character(subject(RV)[end(RV)+1]),split=""))] # we want the base just downstream of the motif
	   names(temp.list[["rev"]]) <- round(mcols(RV)$score,3)
    } else {
      temp.list[["rev"]] <- NA
    }
    double.motif.base[[i]] <- temp.list
  }
  D <- lapply(double.motif.base,unlist) # we don't care if the peak was found forward or reverse
  return(D)
}

#############################################
##                                                                    				#
## 	Reorienting coverage matrices based on motif orientation   	#
##                                                                    				#
#############################################

## function to re-orient the Abf1 binding/coverage matrix, needs an object created by the closest.motif.hi function!
## if there is no motif, the whole peak has the value NA
## The function can deal with characters (DNA) or numeric (coverage) data. Set DNA.input to TRUE or "auto" for DNA sequences.
reorient.matrix <- function(occupancy.list ,closest.motifs,trim.cols=FALSE,DNA.input=FALSE){
	stopifnot(colnames(closest.motifs) %in% c("fw","rev","closest"))
	stopifnot(is.list(occupancy.list))
	if(is.null(names(occupancy.list))) stop("Please provide a list with names!")
	
	## what is the extra size to add?
	extra.size <- max(abs(closest.motifs[,1:2]),na.rm=TRUE)
	
	## select the peaks that are not NA, that is the ones where at least 1 motif was found
	not.NA <- closest.motifs[!is.na(closest.motifs[,3]),]
	
	reoriented.list <- list()
	
	## check if the input is DNA or not
	if( is.logical(DNA.input) && !is.na(DNA.input)){
		DNA <- DNA.input 
	} else if(!is.logical(DNA.input) && grepl("auto", tolower(DNA.input))) { # for now assume all are DNA 
		m <- occupancy.list[[1]]
		m <- m[!is.na(m)] # select all non-NA values
		if(all(grepl("A|C|G|T|N|a|c|g|t|n",m))){ # detect bases, for now also N is detected
			DNA <- TRUE
			warning("DNA sequence detected! If this was indeed the intention, please ignore this warning")
		} else {
			DNA <- FALSE
			warning("NO DNA sequence detected!! If this was indeed the intention, please ignore this warning")
		}
	} else stop("Please provide a valid DNA.input parameter. This can be either 'auto', TRUE or FALSE '")
	
	##  vector for the reverse complement
	RC <- c(a="t",c="g",g="c",t="a",A="T",C="G",G="C",T="A")
	
	## fill in the matrix for each component of the list
	for(i in names(occupancy.list)){
		## make a new matrix that can be filled in 
		new.mat <- matrix(NA,ncol=ncol(occupancy.list[[i]]) + extra.size*2 , nrow=nrow(occupancy.list[[i]])) 
		rownames(new.mat) <- rownames(occupancy.list[[i]])
		cols <- as.numeric(colnames(occupancy.list[[i]]))
		colnames(new.mat) <- seq(min(cols)-extra.size,max(cols)+extra.size,by=1)
		
		## if the input is DNA the fw is the same, but the reverse needs to be reverse complemented!
		if(DNA){
			for(j in rownames(not.NA)){
				fw.or.rev <- not.NA[j,"closest"]
				if(fw.or.rev == "fw"){
					columns <- seq( min(cols) - not.NA[j,fw.or.rev]  , max(cols) - not.NA[j,fw.or.rev]  ,  by=1 )
					new.mat[j,as.character(columns)] <- occupancy.list[[i]][j,]
				} else if(fw.or.rev == "rev"){ # only for the revs you need to take the reverse complement
					columns <- seq(max(cols)+not.NA[j,fw.or.rev], min(cols) + not.NA[j,fw.or.rev]  )
					new.mat[j, as.character(columns)] <- RC[occupancy.list[[i]][j,]]			
				} else {
					stop(" Could not find if the motif if forward or reverse and the function did not figure it out itself. Please investigate!")
				}
			} 
		} else {
			## fill it in per row, shift all to have the peak motif in the center
			## reorient the ones that have the peak in reverse orientation! 
			for(j in rownames(not.NA)){
				fw.or.rev <- not.NA[j,"closest"]
				if(fw.or.rev == "fw"){
					columns <- seq( min(cols) - not.NA[j,fw.or.rev]  , max(cols) - not.NA[j,fw.or.rev]  ,  by=1 )
					new.mat[j,as.character(columns)] <- occupancy.list[[i]][j,]
				} else if(fw.or.rev == "rev"){
					columns <- seq(max(cols)+not.NA[j,fw.or.rev], min(cols) + not.NA[j,fw.or.rev]  )
					new.mat[j, as.character(columns)] <- occupancy.list[[i]][j,]
				} else {
					stop(" Could not find if the motif if forward or reverse and the function did not figure it out itself. Please investigate!")
				}
			}
		}
		
		
		## put the new matrix in the list
		if(trim.cols) new.mat <- new.mat[,colnames(occupancy.list[[i]])]
		reoriented.list[[i]] <- new.mat
	}
	
	return(reoriented.list)

}

#########################################
##                                                                    		#
##          		Fitting exponential decay functions           	#
##                                                                    		#
#########################################



## now fit the models on all at the same time, you need to use lapply!
## you need to use tryCatch, because if one of the items fails, all fail!
nop <- function(...)invisible(NULL)

## Based on this blog: http://douglas-watson.github.io/post/2018-09_exponential_curve_fitting/
## it expects a DF with 2 columns, one "ChIP" and the other "time". You can use apply as well on this (see scripts that use the function)
try.exp.fit <- function(ChIP, Data){ # if it can't fit for whatever reason it will return NULL
  fit <- NULL 
  tryCatch( {fit <- nls(ChIP ~ SSasymp(time, yf, y0, log_alpha), data = Data)},
			error=nop,
            warning=nop,
            finally=nop)
	return(fit)
}

## predict the curve based on the fit. Needs a dataframe with the new time points to predict
try.pred <- function(fit, newData){ # if it can't fit for whatever reason it will return NULL
  pred <- NULL 
  tryCatch( {pred <- predict(fit, newdata = newData)},
			error=nop,
            warning=nop,
            finally=nop)
	return(pred)
}

## calculate the pseudo-R squared
R.squared <- function(model, y) {
   RSS <- sum(residuals(model)^2)
   TSS <- sum((y - mean(y))^2)
   (1 - (RSS/TSS))
}


## extract the coefficients  of the fit. (Useful to compare all the peaks)
get.coefficients <- function(data.fit,data.points, peak.names){
	
	coefficients.all <- matrix(NA, ncol=12, nrow=length(peak.names))
	rownames(coefficients.all) <- peak.names; colnames(coefficients.all) <- c("yf","y0","log_alpha","alpha","R2","mean_life","t50","y50%","t90","y90%","reads_t00","reads_t90")
	
	try.coef <- function(...){ # if it can't fit for whatever reason it will return NULL
		coefs <- vector(length=3)
		coefs[1:3] <- NA 
		if(!is.null(...)){
			tryCatch( {coefs <- coef(...)},
				error=nop,
				warning=nop,
				finally=nop)
		}
		return(coefs)
	}
		
	coefficients.all[peak.names,c("yf","y0","log_alpha")] <- t(sapply(data.fit[peak.names], try.coef))
	coefficients.all[,"alpha"] <- exp(coefficients.all[,"log_alpha"])
	
	colnames(coefficients.all[,c("log_alpha","alpha","mean_life")]) <- c("log_off_rate","off_rate","res_time")
	
	## calculate the pseudo-R squared
	coefficients.all[,"R2"] <- sapply(peak.names, function(x) { if(!is.na(coefficients.all[x,1]))  R.squared(data.fit[[x]],data.points[,x]) else NA} )	
	
	## the mean lifetime can be detemined by the formula tau=1/k(off) or mean_life=1/alpha
	coefficients.all[,"mean_life"]  <- 1/coefficients.all[,"alpha" ]
	
	## how to decuct time from the function:
	## time =  log((y - yf) / (y0 - yf)) / -alpha 
	## y in this case is "90% depletion" = (y0 - yf) * 0.9, and that needs to be subtracted from y0
	## depl.90 = y0 - (y0 - yf) * 0.9
	## t1/2 (or t50) can also be calculated as: log(2)/alpha
	
	## make a function that will calculate the time at which 90 percent is reached
	#y.90perc.depl <- y0 - (y0 - yf) * 0.9    # this is function slightly more readable
	#t.90perc.depl =  log((y.90perc.depl - yf) / (y0 - yf)) / -alpha
	calc.perc.depl <- function(coefs.mat, depl.perc){
		y.perc <- coefs.mat[,"y0"] - (coefs.mat[,"y0"]  - coefs.mat[,"yf"] ) * depl.perc # y0 - (y0 - yf) * 0.9 or 0.5
		t.perc <-  log((y.perc - coefs.mat[,"yf"]) / (coefs.mat[,"y0"] - coefs.mat[,"yf"])) / -coefs.mat[,"alpha"]  # log((y.90perc.depl - yf) / (y0 - yf)) / -alpha 
		ty.perc <- cbind(t.perc,y.perc)
		return(ty.perc)
	}
	
	coefficients.all[peak.names,c("t50","y50%")] <- calc.perc.depl(coefs.mat = coefficients.all[peak.names,],depl.perc=0.5)
	coefficients.all[peak.names,c("t90","y90%")] <- calc.perc.depl(coefs.mat = coefficients.all[peak.names,],depl.perc=0.9)
	
	## rename a few columns with more sensible columns names
	colnames(coefficients.all)[which(colnames(coefficients.all) %in% c("log_alpha","alpha","mean_life"))] <- c("log_off_rate","off_rate","mean_res_time")
		
	## fill the reads_t00 and reads_t90 cloumns with the average binding of the triplicates
	t00.idx <- which(data.points[,"time"] == "0") ; t90.idx <- which(data.points[,"time"] == "90")
	coefficients.all[peak.names, "reads_t00"] <- sapply(data.points[t00.idx, peak.names],mean)
	coefficients.all[peak.names, "reads_t90"]  <- sapply(data.points[t90.idx, peak.names],mean)

	return(coefficients.all)

}


## this function will plot data + the fit for a single peak
plot.fit.figure <- function(data.points, data.fit, data.pred, peak.name, Col=NULL, div.y=TRUE, ylim=NULL,ylab=NULL, pred.x=NULL,pch=19, perc.depl=90,
									leg.args=TRUE,mgp=c(2.5,1,0),leg.cex=1,abl=TRUE,...){ # 
	
	if(is.null(Col)) Col <- c(rep("black",11),rep("#E69F00",11),rep("#56B4E9",11)) 
	
	if(!is.null(ylim)){
		if(length(ylim) ==1 && ylim == 0) {
			ylim=c(0,  max(data.points[,peak.name]/if(div.y) 1000 else 1) )*1.05	# max data and starting at 0, give it 5% extra on the top
		} else if(length(ylim) == 1 && ylim == 1) {
			ylim=range(data.points[,peak.name]/if(div.y) 1000 else 1 ) 			#take the whole range, same as ylim=NULL actually
		} else {
			ylim=ylim									# take the ylim that was provided
		}
	}
	
	if(is.null(ylab)){
		if(div.y) {
			ylab <- expression(paste("ChIP signal (reads x10"^3,")"))	
		} else {
			ylab <- "ChIP signal (reads)"
		}
	}
	
	dots <- list(...)
	if(is.null(pred.x)) pred.x <- data.frame(time=0:90)
	
	if(is.null(Col)) Col <-  c(rep("black",11),rep("#E69F00",11),rep("#56B4E9",11))
		
	if(div.y) y <- data.points[,peak.name]/1000 else y <- data.points[,peak.name]
	
	plot(x = data.points$time, y=y, ylim=ylim, ylab=ylab,type="n",  las= if(div.y) 1 else NULL,...)
		
	## plot the fit, but only if the peak could be fitted at all, predict/fit != NULL
	if(!is.null(data.pred[[peak.name]]) && !is.null(data.fit[[peak.name]])){	
		if(div.y) y.pred <- data.pred[[peak.name]]/1000 else y.pred <- data.pred[[peak.name]]
		lines(x= pred.x$time, y=y.pred)
		
		## take out the coefficients (yf and alpha)
		coefs <- coef(data.fit[[peak.name]])
		alpha <- exp(coefs["log_alpha"])
		y0 <- coefs["y0"]
		yf <- coefs["yf"]
		
		## how to decuct time from the function:
		## time =  log((y - yf) / (y0 - yf)) / -alpha 
		## y in this case is "90% depletion" = (y0 - yf) * 0.9, and that needs to be subtracted from y0
		## depl.90 = y0 - (y0 - yf) * 0.9
		
		if(grepl("res|mean", tolower(perc.depl))) {
			t.perc.depl <- 1/alpha
			y.perc.depl <- yf + (y0 - yf)*exp(-alpha * t.perc.depl )
		} else if(is.numeric(perc.depl)){
			if(perc.depl > 1) perc.depl <- perc.depl/100
			y.perc.depl <- y0 - (y0 - yf) * perc.depl
			t.perc.depl =  log((y.perc.depl - yf) / (y0 - yf)) / -alpha 
		} else {
			perc.depl <- NULL
		}
				
		# for now draw ablines at the site of the 90% depletion to show that is is correct
		if(!is.null(perc.depl) && abl) abline(h=y.perc.depl/if(div.y) 1000 else 1,v=t.perc.depl,lty=2)	
		
		r.val = R.squared(data.fit[[peak.name]],data.points[,peak.name]) # function from Philip, must be loaded in!
		leg.text = vector('expression',if(!is.null(perc.depl)) 3 else 2)
		leg.text[1] = substitute(expression(italic(R^2) == MYVALUE), 
						list(MYVALUE = format(r.val, dig=2)))[2]
		leg.text[2] = substitute(expression(paste(italic(k[off]) == MYVALUE, " min"^-1)), 
						list(MYVALUE = sprintf("%0.3f",alpha)))[2]
		if(is.numeric(perc.depl)){ # if the percentage depletion was provided, add a line saying at what timepoint the given percentage was reached
			leg.text[3] = substitute(expression(MYVALUE), 
						list(MYVALUE = sprintf("%0.0f%% depletion at %0.0f min",perc.depl*100,t.perc.depl)))[2]
		} else if(grepl("res|mean", tolower(perc.depl))) { # or if the residence time was asked:
			leg.text[3] = substitute(expression(MYVALUE), 
						list(MYVALUE = sprintf("Residence time = %0.0f min", t.perc.depl)))[2]
		}
		legend("top", legend= leg.text, bty="n",cex=leg.cex)
	} else {
		legend("top", legend=c("Fit could not be made"), bty="n",cex=leg.cex)
	}
	points(x = data.points$time, y=y, pch=pch,col= Col, ...)
		 
	## create the legend for now hard-coded (not so nice)
	if(is.logical(leg.args) && leg.args){
		legend("topright",legend=c("A9","A16","B5"),bty="n", text.col=c("black","#E69F00","#56B4E9"),text.font=2,cex=leg.cex)
	} else if(is.list(leg.args)) {
		do.call(legend, leg.args)
	} else if(!is.null(leg.args)){
		warning("Did not call legend, need a list of arguments to pass to the legend function!")
	}
}





