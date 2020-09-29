## This script is used to process all data files for the paper:
## "Genome-wide off-rates reveal how DNA binding dynamics shape transcription factor function"
##
## The data files that are generated with this code can subsequently be used to create the figures with the script "DIVORSEQ_figures.R"
##
## All the code was run using R version 3.2.2
##

library(rtracklayer)	# v1.30.4, Loads in packages (among others): GenomicRanges_1.22.4, IRanges_2.4.8, S4Vectors_0.8.11
library (seqinr)  	# needed for read.fasta
library(Biostrings) # needed for matchPWM and DNAStringset etc


wd <- "/path/to/working/directory/" # put here your working directory, please end with a "/"
path_to_scripts <- "/path/to/scripts/directory/" # put here the path to your directory containing the scripts, please end with a "/"
path_to_files <- "/path/to/data/files/" # put here the path to your directory containing all data files, please end with a "/" 
output_path <- "/path/to/output/files/" # put here the path that you will use to output all your file, please end with a "/"
setwd(wd)

## Source the functions
source(paste0(path_to_scripts,"DIVORSEQ_functions.R" )) # source the script containing all required scripts

### Create a few objects that will be used later
##  The timepoints that were used
times <- c(0,5,10,15,20,30,40,50,60,75,90)
times.no10 <- c(0,5,15,20,30,40,50,60,75,90)

## the new data points used for the prediction
prd <- data.frame(time=0:90)

## load in all the peaks that were called. These were called by pooling the 3 replicate t00 timepoints, using a FE cutoff of 2.

## this is what the columns mean in the macs2 output of the narrowPeak file
## columns "pValue" and "qValue" are actually "-log10pValue" and "-log10qValue" but GRanges does not allow these names.
extraCols.narrowPeak <- c(fold_enrichment = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")

all.t00.peaks.FE2 <- import(paste0(path_to_files,"All_t00_combined_FE2_peaks.narrowPeak") , format="BED", extraCols = extraCols.narrowPeak)
names(all.t00.peaks.FE2) <- gsub("All_|_combined","",all.t00.peaks.FE2$name)

## find the summits of the undepleted sample in 101 bps window
all.t00.summits.FE2 <- make.summits.object(all.t00.peaks.FE2) 

## set a strict filter for the peaks that will be used eventually. Use the peaks with a fold enrichment cutoff of > 4
strict.peaks <- names(all.t00.peaks.FE2[all.t00.peaks.FE2$fold_enrichment > 4]) # 421 peaks




######
###
###		Calculating coverage
###
######

## first get the names of all the files of the Abf1 depletion:
Abf1.bw.names <- list.files(path_to_files)
Abf1.bw.idx <- grep("Abf1_.*bw",Abf1.bw.names)
Abf1.bw.names <- Abf1.bw.names[Abf1.bw.idx ] # should be 32 bw files, 11 time points in 3 replicates except for the 10 min time point

## create a matrix that will be used to check the total coverage. (Samples were scaled to 1 million reads, and all reads were 101 bp long, total coverage should be 101 * 10^6)
total.coverage.cen101  <- matrix(NA, ncol=16, nrow=length(Abf1.bw.names)) # 16 columns because yeast has 16 chromosomes
colnames(total.coverage.cen101) <- levels(seqnames(all.t00.summits.FE2)); rownames(total.coverage.cen101) <- gsub("Abf1_|.bw","",Abf1.bw.names)

## read in the bw files and calculate the coverage
## run this all, we also would like to know what the total coverage is per replicate (for later)
## this will take some time to run, but it makes sure only 1 bw file is read in at the same time
## and it will prevent that the R sessions uses too much memory and crashes.
## read in the reads centered and smoothed by 101 bp

coverage.all.cen101 <- list()

for( i in Abf1.bw.names ) {
  bw <- import( con=paste0(path_to_files,i), format="bigWig")
  coverage.all.cen101[[i]] <- Cover.bw(bw) 
  total.coverage.cen101[gsub("Abf1_|.bw","",i),] <-  sapply( coverage.all.cen101[[i]],sum)
}

## check if it worked:
apply(total.coverage.cen101,1, sum) # all samples should have a total coverage of 101*10^6 

## make an object that contains the coverage per peak for each sample in a matrix per basepair.
## the size is not specified, so it will use the size of biggest feature. Since all features are 101 bp, that size will be used.
## this takes quite a while to run
coverage.matrix.cen101 <- list() # only to this the first time you run this
for(i in names(coverage.all.cen101)) coverage.matrix.cen101[[i]] <- align.data.peaks2(features=all.t00.summits.FE2,data=coverage.all.cen101[[i]])

## in addition calculate the total coverage at all peaks, this will be used to normalize the each sample 

## calculate the number of reads in each peak in each sample:
cov.peak.samp <- matrix(NA,ncol=length(coverage.all.cen101),nrow=length(all.t00.peaks.FE2))
colnames(cov.peak.samp) <- names(coverage.all.cen101)
rownames(cov.peak.samp) <- names(all.t00.peaks.FE2)

## and fill in the matrix (takes a bit of times since it loops over 948 peaks)
for(peak in names(all.t00.peaks.FE2)){
	CHROM <- as.character(seqnames(all.t00.peaks.FE2[peak]))
	START <-  start(all.t00.peaks.FE2[peak])
	END  <- end(all.t00.peaks.FE2[peak])
		
	for(samp in names(coverage.all.cen101)){
		cov.peak.samp[peak,samp] <- sum(coverage.all.cen101[[samp]][[CHROM]][START:END])
	}
}


## since it takes a while to create these objects, save them.
save(coverage.matrix.cen101, file=paste0(output_path,"coverage.matrix.cen101.RData"))
save(total.coverage.cen101, file=paste0(output_path,"total.coverage.cen101.RData"))
save(cov.peak.samp,file=paste0(output_path,"cov.peak.samp.RData"))

######
###
###		Normalization
###
######

######!!	Only run this the first time!!!		################################

## load the coverage matrix, the total coverage per sample and the coverage per peak per sample
load(paste0(output_path,"coverage.matrix.cen101.RData"))
load(paste0(output_path,"total.coverage.cen101.RData"))
load(paste0(output_path,"cov.peak.samp.RData"))

## calculate the total coverage and the coverage in the called peaks per sample:
cov.each.samp <- apply(cov.peak.samp,2,sum)
cov.total <- apply(total.coverage.cen101,1, sum)

## calculate the percentage of the coverage that is not covered by the called peaks.
percent.outside.peaks <- 1-cov.each.samp/cov.total 

##  calculate the coverage per peak
cov.per.peak <- calc.coverage.per.peak(coverage.list = coverage.matrix.cen101, sample.names = Abf1.bw.names , func=sum, rm.na=T)

## resructure the data.frame
cov.per.peak2 <- t(cov.per.peak) 
cov.per.peak2 <- data.frame(time=  as.numeric(gsub("Abf1_|_min.*","",rownames(cov.per.peak2))), cov.per.peak2) # add a column with the time points
cov.per.peak2 <- cov.per.peak2[order(cov.per.peak2[,"time"]) ,] # reorder
 
## normalize all peaks by the percentage of background reads
## this equalizes the number of background reads for all samples
all.peaks.norm.background <- cov.per.peak2
all.peaks.norm.background[,-1] <- all.peaks.norm.background[,-1]  /percent.outside.peaks[rownames(cov.per.peak2)]

## save this
save(all.peaks.norm.background, file=paste0(output_path,"all.peaks.norm.background.RData"))

##############################################################

######
###
###		Peak filtering
###
######

######!!	Only run this the first time!!!		################################

## only peaks with a G or C at -8 bp from the motif can be trusted, because only at these motifs Abf1 can be crosslinked (see manuscript)
## these are the same as the peaks with a G or C at +8 bp in Rossi et al., 2018 (doi: 10.1101/gr.229518.117 )

## first  find all motifs close to the the peak summit, for this the genome sacCer3 (v64, 2011) is needed
## this can be downloaded here: http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-1-1_20110203.tgz
sequences <- read.fasta(paste0(path_to_files,"S288C_reference_sequence_R64-1-1_20110203.fsa"))
names(sequences) <- c(paste0("chr", as.roman(1:16)), "chrmt" )  # rename the sequences, they should be in this order

## make a a summits object with 1001 bp windows
all.t00.summits.FE2.1001bp <- make.summits.object(all.t00.peaks.FE2, plus=500, minus=500) 

## only to this the first time you run this
## get a matrix with sequences around each peak
seq.mat <- align.data.peaks2(features=all.t00.summits.FE2.1001bp ,data=sequences)

## and save this:
save(seq.mat,file=paste0(output_path,"seq.mat.RData"))
##############################################################

load(paste0(output_path,"seq.mat.RData"))

## search for motifs, use the motif from MacIsaac et al., Bioinformatics, 2006: http://yetfasco.ccbr.utoronto.ca/showPFM.php?mot=YKL112W_684.0
Abf1.MacIsaac.pfm <- as.matrix(read.table(paste0(path_to_files,"ABF1_motif_YETFASCO_McIsaac2006.txt" ),stringsAsFactors=FALSE,row.names=1))
Abf1.MacIsaac.pfm <- Abf1.MacIsaac.pfm[sort(rownames(Abf1.MacIsaac.pfm)),]

## however, this formatting is not a real PFM as the program expects it, therefore it needs to be converted:
Abf1.MacIsaac.pfm2<- round(Abf1.MacIsaac.pfm*1000) 
## this gives a rounding problem: 
#apply(Abf1.MacIsaac.pfm2,2,sum) # column V9 sums up to 1001 while the others do to 1000 (as expected).
Abf1.MacIsaac.pfm2["T","V9"] <- Abf1.MacIsaac.pfm2["T","V9"] -1 # the "T" has a high number of occurences (less influence of subtracting) and .55 in decimals (so "closest to 0")
class(Abf1.MacIsaac.pfm2) <- "integer" # it must be an integer matrix

## we need the priors (normally all bases are assumed to be equal, but they are not!)
priors <- table(seq.mat[,as.character(-200:200)]) # take the priors from the area surrounding the peak summit
priors  <- priors / sum(priors)
priors <- c(A= sum(priors[c("a","t")])/2 , C=sum(priors[c("c","g")])/2, G = sum(priors[c("c","g")])/2, T=sum(priors[c("a","t")])/2)

## and convert to pwm
Abf1.MacIsaac.pwm <- PWM(Abf1.MacIsaac.pfm2, prior.params = priors)
Abf1.MacIsaac.pwm.rev <- reverseComplement(Abf1.MacIsaac.pwm)

## convert the sequences that are NA around the peak to "N"
## this happens because some peaks are close to the edge of the chromsome, this has little effect on finding the motif, but it crashes the function
seq.mat.noNA <- seq.mat
seq.mat.noNA[is.na(seq.mat.noNA)] <- "N" 

## A DNA stringset object is needed to search for motifs. 201bp windows centered on the peak summit are used:
seq.mat.DNA.201bp <- DNAStringSet(apply(seq.mat.noNA[,as.character(-100:100)],1,paste,collapse=""))

## The + and - strand need to be matched separately!!
motif.loc.201bp <- lapply(seq.mat.DNA.201bp,FUN=matchPWM, pwm=Abf1.MacIsaac.pwm , with.score=TRUE,min.score="85%") #80% is the default
motif.loc.201bp.rev <- lapply(seq.mat.DNA.201bp,FUN=matchPWM, pwm=Abf1.MacIsaac.pwm.rev, with.score=TRUE,min.score="85%") #80% is the default

## This allows for counting the number of motifs that are found!
nr.motif <- list()
nr.motif[["pos_nr"]] <- sapply(motif.loc.201bp,length)
nr.motif[["rev_nr"]] <- sapply(motif.loc.201bp.rev,length)
nr.motif[["both_nr"]]   <- nr.motif[["rev_nr"]] + nr.motif[["pos_nr"]]

## Search for the motif that has the strongest motif score i.e. is closest to the consensus
## note: The function that is used can also use a distance cutoff, and search for the closest motif. Leave at 100 for 201 windows to ignore the distance to the peak
## IMPORTANT: the name might be a bit misleading, but this output a data.frame, the third column is called "closest" but in
## this case just gives the "strongest" motif (the motif that is closest to the consensus).
Abf1.highest <- closest.motif.hi(fw = motif.loc.201bp, rv = motif.loc.201bp.rev,size=201,dist.cutoff=100) # this is equivalent to taking the strongest motif in the 201 bp window

## now reorient the sequence matrix on the highest motif, such that the sequences are centered on the midpoint of the highest motif
## and the all the sequences are reoriented in the same orientation
seq.mat.reorient.hi <- reorient.matrix(occupancy.list  = list(sequences=seq.mat), closest.motifs = Abf1.highest,trim.cols=TRUE,DNA.input="auto" )[[1]] # rewritten to deal with DNA properly
#seq.mat.reorient.hi[1:30,as.character(-6:6)] # this should show that all are oriented on the motif: RTCRYnnnnnACG # NA means not motif found

## now check for each peak what the base is at position -8 (G/C or A/T)
pos.8.motif.strict <- list() 
for(base in c("A","C","G","T")) pos.8.motif.strict[[base]] <-  names(which(seq.mat.reorient.hi[strict.peaks,"-8"] == tolower(base)))

GC.peaks.strict <- intersect(strict.peaks,unlist(pos.8.motif.strict[c("C","G")]))

## however, for sites with more than 1 motif, a very strict cutoff was used such that both motifs need to have a G/C at position -8:

## find the peaks with more than 1 motif 
double.motif.seqs <- find.double.motifs(names(nr.motif[["both_nr"]][nr.motif[["both_nr"]] > 1]) ,motif.loc.201bp , motif.loc.201bp.rev )
double.motif.seqs <- double.motif.seqs[intersect(names(double.motif.seqs ), strict.peaks)] # 71/421 peaks

## find all the peaks with any motif with an A/T at position -8:
AT.idx <- sapply(lapply(double.motif.seqs, grepl,pattern="A|T" ),any) ; AT.idx.name <- names(AT.idx[AT.idx])

strict.peaks.GC <- setdiff(GC.peaks.strict ,AT.idx.name) # 195 peaks

######
###
###		Fitting Exponential Decay
###
######

load(paste0(output_path,"all.peaks.norm.background.RData"))

## fit the exponential decay function on the peaks with FE > 4 and only motifs with a G/C at -8 bp from the motif midpoint:
GC.peaks.norm.background.fit  <- lapply(all.peaks.norm.background[,strict.peaks.GC], FUN = try.exp.fit, Data=all.peaks.norm.background ) # you don't want to fit it on the time column
#all.peaks.norm.background.fit  <- lapply(all.peaks.norm.background[,-1], FUN = try.exp.fit, Data=all.peaks.norm.background ) # in case you want to fit all peaks, for Abf1 not recommended!
GC.peaks.norm.background.pred <- lapply(GC.peaks.norm.background.fit, FUN = try.pred, newData=prd)

## calculate some coefficients that are interesting
coefs.strict.GC <- get.coefficients(data.fit=GC.peaks.norm.background.fit, data.points =all.peaks.norm.background, peak.names=strict.peaks.GC)

######
###
###		Add motif scores and number of motifs to coefficients object
###
######

## add the number of motifs:
coefs.strict.GC <- data.frame(coefs.strict.GC,nr_motif=nr.motif[["both_nr"]][rownames(coefs.strict.GC)])

## calculate the maximum score of each motif for each peak (so peaks  with 2 motifs get the score of the maximum motif)
all.fw.peaks.max.score <- unlist(lapply(lapply(sapply(motif.loc.201bp , mcols),`[[`,"score"),FUN=function(x) {if(length(x) >1) max(x) else x}))
all.rv.peaks.max.score <- unlist(lapply(lapply(sapply(motif.loc.201bp.rev , mcols),`[[`,"score"),FUN=function(x) {if(length(x) >1) max(x) else x}))

## put these values of the fw and reverse motifs in a matrix and calculate the maximum of the two
all.peaks.max <- matrix(0,ncol=2,nrow=length(motif.loc.201bp)) ; rownames(all.peaks.max) <- names(motif.loc.201bp);colnames(all.peaks.max) <- c("fw","rev")
all.peaks.max[names(all.fw.peaks.max.score),"fw"] <- all.fw.peaks.max.score
all.peaks.max[names(all.rv.peaks.max.score),"rev"] <- all.rv.peaks.max.score

all.peaks.max <- cbind(all.peaks.max,pmax=pmax(all.peaks.max[,"fw"],all.peaks.max[,"rev"]))
all.peaks.max[ all.peaks.max[,"pmax"] == 0,"pmax"] <- NA

## and add them:
coefs.strict.GC$motif_score <-  all.peaks.max[rownames(coefs.strict.GC),"pmax"]

######
###
###		Determine the residence time quartiles
###
######

## before the residence time quartiles are calculated, first the peaks that overlap with telomeres are filtered out
## this is done because they behave differently compare to the other peaks and therefore it might not be fair to compare
## differences between the quartiles since the telomeric peaks are enriched in the quartile with the longes residence time

## get the annotation of the telomeres, this one can be used: http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff.gz
## but it is advisable to remove the FASTA sequences of the genome at the end.
all.SGD.features <- import(paste0(path_to_files,"saccharomyces_cerevisiae_nofasta.gff"))
telomeres <- all.SGD.features[all.SGD.features$type=="telomere"]

## get the strict peaks (FE > 4) and see which overlap with telomeres
all.t00.peaks.FE2.strict <- all.t00.peaks.FE2[strict.peaks]
peaks.telo.strict <-  all.t00.peaks.FE2.strict [overlapsAny(all.t00.peaks.FE2.strict , telomeres)] # 28 peaks strict peaks overlap with telomeres

## don't use the telomeric peaks to determine the residence time quartiles
GC.peaks.no.telo <- setdiff(rownames(coefs.strict.GC),names(peaks.telo.strict)) # 191 peaks remain

## determine the residence time quantiles
quantile.GC <- quantile(coefs.strict.GC[GC.peaks.no.telo,"off_rate"])

## and add them to the coeffficents object
coefs.strict.GC[GC.peaks.no.telo,"res_time_quartile"] <-  as.character(cut(coefs.strict.GC[GC.peaks.no.telo ,"off_rate"],breaks=quantile.GC,labels= c("longest","long","short","shortest"),include.lowest=T))
coefs.strict.GC[intersect(rownames(coefs.strict.GC),names(peaks.telo.strict)),"res_time_quartile"] <- "telo" # and label the telomeric peaks


######
###
###		creation of dataset EV1:
###
######

dataset_EV1 <- coefs.strict.GC[,c("reads_t00","reads_t90","yf","y0","off_rate","mean_res_time","nr_motif","motif_score","res_time_quartile")]

dataset_EV1 <- cbind(dataset_EV1,chr=seqnames(all.t00.summits.FE2[rownames(dataset_EV1)]),
                     start= start(all.t00.summits.FE2[rownames(dataset_EV1)]),
                     end=end(all.t00.summits.FE2[rownames(dataset_EV1)]))
dataset_EV1 <- dataset_EV1[,c(10:12,1:9)]

dataset_EV1[,-c(1,2,3,8,10,12)] <- round(dataset_EV1[,-c(1,2,3,8,10,12)],3)
dataset_EV1[,"off_rate"] <- round(dataset_EV1[,"off_rate"],4)
colnames(dataset_EV1)[which(colnames(dataset_EV1) == "reads_t00")] <- "binding_t00"
colnames(dataset_EV1)[which(colnames(dataset_EV1) == "reads_t90")] <- "binding_t90"

write.table(dataset_EV1,file=paste0(output_path,"dataset_EV1.txt"),quote=FALSE,sep="\t")