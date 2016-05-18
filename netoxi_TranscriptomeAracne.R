########################################################################
# Transcriptome data preparation:                                      #
#	- read in count data (exported from subread)		       #
#	- QC                                                           #
# 		- keep only samples with at least 1 Mio. reads         #  
#               - filter out loci with less than 2 counts per million  #
#                 in at least 3 samples                                #
#	- normalize within experiment                                  #
#	- normalize between experiments (quantile norm.)               #
#	- average replicates                                           #
#	- sample correlation (correlogram)                             #
#	- prepare input for bash aracne                                #
#		- filter out transcript with low variation             #
#		- genefilter IQR 0.5                                   #
#	- map transcription factor IDs to the latest chlamy anno       #
#	- prepare bash aracne output for Cytoscape                     #
#	                                                               #
# Diana Coman Schmid                                                   #
# Eawag 2014-2015                                                      #
# diana.comanschmid@eawag.ch                                           #
########################################################################

####################
### RNA-seq data ###
####################



# source("https://bioconductor.org/biocLite.R")
# biocLite("edgeR")

library("limma")
library("edgeR")


### experiments with replicates ###

# read in the sample description file (columns: Exp[SRP002], ExpFull[SRP002284], Sample[SRR039670_1.fastq.gz_subread.sam_counts.txt]) for experiments that have replicates
rseq.s <- read.table("/media/dianacs/data1/Express_data/RNAseq/SRP/samplesWithReps.txt",header=TRUE,sep="\t")



###zebrafish embryo rnaseq counts
rseq.s <- read.csv("E://ZebrafishEmbryo//subreadCounts//finalSubreadCounts//Archive//featureCounts_table.csv",header=TRUE,row.names=1)
dim(rseq.s)
rseq.s[1:5,1:5]

barplot(apply(rseq.s,2,mean))
barplot(colSums(rseq.s))

# keep only samples with at least 1 Mio. reads
counts.df.n <- rseq.s[,colSums(rseq.s) > 1e+06]
dim(counts.df.n)

# filter out loci with less than 2 counts per million in at least 3 samples
isexpr <- rowSums(cpm(counts.df.n)>2) >= 3
counts.ex <- counts.df.n[isexpr,]
dim(counts.ex)
counts.ex[1:5,1:5]


# TMM normalization (within experiment)
x <- calcNormFactors(counts.ex)

libSize <- colSums(counts.ex)
normalizedCounts <- log2(10e+06 *(counts.ex / (libSize * x))+1)

# plot TMM normalized expression values distribution
plotDensities(normalizedCounts)



# normalize between experiments (quant. norm)
merge.norm.rseq.qnorm <- normalizeQuantiles(normalizedCounts)
merge.norm.rseq.qnorm[1:5,1:5]
class(merge.norm.rseq.qnorm)
range(as.matrix(merge.norm.rseq.qnorm))

# plot expression levels distribution before and after quantile normalization
plot(density(as.matrix(merge.norm.rseq.qnorm)),col="blue",xlab="mRNA levels (log2)")
lines(density(as.matrix(normalizedCounts)),col="magenta")

write.csv(merge.norm.rseq.qnorm,"E://ZebrafishEmbryo//subreadCounts//finalSubreadCounts//Archive//normZebrafishEmb.csv")

#######################################################################################
#######################################################################################


# assign replicates to their corresponding experiment
# make list with unique experiment names (ExpFull) as elements
rseq.l <- as.list(unique(rseq.s$ExpFull))
names(rseq.l) <- unique(rseq.s$ExpFull)


# to each element of the list (experiment) assign the corresponding replicate samples
for (l in 1:length(names(rseq.l))){
  rseq.l[[l]] <- as.character(rseq.s[rseq.s$ExpFull == names(rseq.l)[l],"Sample"])
}


# count and plot the sample number per experiment
ex.no <-cbind(as.numeric(summary(rseq.l)[,"Length"]),as.numeric(summary(rseq.l)[,"Length"]))
colnames(ex.no) <- c("SampleNo","same")
rownames(ex.no) <- rownames(summary(rseq.l))

dev.off()
pdf(file="/media/dianacs/data1/Express_data/RNAseq/SRP/RNAseq_part1/subread_aligned/counts_part1and2/samplesPerExp_barplot.pdf",width=11, height=8.5,pointsize=12, paper='special')
par(mar=c(8,8,4,4))
xLabLocs <- barplot(t(ex.no[,1]), space=0.4, xaxt='n', ann=FALSE,ylab='Number of RNA-seq samples per Experiment',col='gold',border=NA,width = 1)
axis(1, cex.axis=0.75, las=2, at=xLabLocs, labels=rownames(ex.no))
text(x= xLabLocs, y= ex.no[,1]+1, labels=as.character(ex.no[,1]), xpd=TRUE,cex=0.75,col="magenta")
dev.off()

# specify the working directory (RNAseq count files exported from subread)
wd.counts <- "/media/dianacs/data1/Express_data/RNAseq/SRP/RNAseq_part1/subread_aligned/counts_part1and2/counts/"

# load the transcript IDs (are the same and in the same orders in each file; default export from subread)
ids <- as.matrix(read.table(file.path(wd.counts,"SRR1622094_1.fastq.gz_subread.sam_counts.txt"),header=TRUE, sep='\t')[,1])

# normalize within experiment and plot QC for each experiment
# specify the working directory for each experiment [contains the .sra files (not used here) and the annotation file (manually curated and assembeled meta-data)]

wds <- c("/media/dianacs/data1/Express_data/RNAseq/SRP/EERAD124/EERAD124/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP002/SRP002284/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP003/SRP003630/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP005/SRP005483/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP009/SRP009466/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP010/SRP010062/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP010/SRP010084/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP010/SRP010563/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP014/SRP014795/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP017/SRP017044/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP018/SRP018835/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP031/SRP031856/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP040/SRP040659/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP040/SRP040944/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP049/SRP049231/")

# make list with working directory for each experiment as elements
rseq.norm <- list()
for (w in 1:length(wds)){
  wd <- wds[w]
  
  # experiment name format: ERADxxx or SRPxxxxxx
  exp <- unlist(strsplit(wd,"/"))[9]

  # load only the columns with the read counts (default export from subread)
  counts.l <- lapply(rseq.l[[exp]], function(x)read.table(file.path("/media/dianacs/data1/Express_data/RNAseq/SRP/RNAseq_part1/subread_aligned/counts_part1and2/counts/",x), header=TRUE, sep='\t')[,7])   
  names(counts.l) <- gsub(".fastq.gz_subread.sam_counts.txt","",rseq.l[[exp]])

  # make a data frame holding counts (columns) for all sample of one experiment and the transcript IDs
  counts.m <- do.call("cbind", counts.l)
  counts.df <- as.data.frame(cbind(ids,counts.m))
  rownames(counts.df) <- counts.df[,1]
  counts.df[,1] <- NULL
  counts.df.n <- apply(counts.df,2, as.numeric)
  rownames(counts.df.n) <- rownames(counts.df)

  # keep only samples with at least 1 Mio. reads
  counts.df.n <- counts.df.n[,colSums(counts.df.n) > 1e+06]

  # filter out loci with less than 2 counts per million in at least 3 samples
  isexpr <- rowSums(cpm(counts.df.n)>2) >= 3
  counts.ex <- counts.df.n[isexpr,]

  # read in the annotation file (manually curated and assembeled meta-data)
  targets <- readTargets(file.path(paste(wd,"annotation.txt",sep="")))
  
  # proceed with samples having annotation info
  srr <- sapply(targets$File,function(x) unlist(strsplit(x,"\\."))[1])
  targets <- targets[which(srr %in% colnames(counts.ex)),]

  # make the design matrix required for edgeR
  f = factor(targets$Name, levels = unique(targets$Name))
  design = model.matrix(~0 + f)
  colnames(design) = levels(f)
  rownames(design) <- colnames(counts.ex)


  x <- DGEList(counts=counts.ex, group=f)
  
  # TMM normalization (within experiment)
  x <- calcNormFactors(x)
  
  # voom transformation
  v <- voom(x,design)

  # store the TMM normalizated and voom transformed matrix as list element 
  rseq.norm[[w]] <- v
  names(rseq.norm)[w] <- exp
    
  # save the TMM normalized within experiment data in the corresponding experiment folder
  write.table(v$E,file.path(paste(wd,"TMMvoom",exp,".txt",sep="")),sep="\t")

  # plot library size, BVC dist. and voom for each experiment
  pdf(file=file.path(paste(wd,"libSize_MDS_voom",exp,".pdf",sep="")),width=11, height=8.5,pointsize=12, paper='special')
  par(mar=c(8,8,4,4))
  par(mfrow=c(2,2))
  
  barplot(colSums(counts.df.n),ylab="Library size",las=2,col="gold",main=exp,cex.names=0.75)
  plotMDS(x, main="BVC dist.",top=500, col=as.numeric(f),cex=0.75)
  vv <- voom(x,design,plot=TRUE)
  dev.off()
}


save(rseq.norm, file="/media/dianacs/data1/Express_data/RNAseq/SRP/rseq_repEx_Norm")


### experiments without replicates ###

# specify the working directory for all experiments without replicates [contains the .sra files (not used here), the count files and one annotation file for all samples (manually curated and assembeled meta-data)]
wd <-"/media/dianacs/data1/Express_data/RNAseq/SRP/RNAseq_part1/subread_aligned/counts_part1and2/single_doubleSample_exp/"

# load the transcript IDs (are the same and in the same orders in each file; default export from subread)
ids <- as.matrix(read.table(file.path(wd,"ERR185151.fastq.gz_subread.sam_counts.txt"),header=TRUE, sep='\t')[,1])

# load only the columns with the read counts (default export from subread)
fls <- list.files(wd,pattern="counts.txt")
countsNoR.l <- lapply(fls, function(x)read.table(file.path(wd,x), header=TRUE, sep='\t')[,7])   
names(countsNoR.l) <- gsub(".fastq.gz_subread.sam_counts.txt","",fls)


# make a data frame holding counts (columns) for samples from all experiments without replicates and the transcript IDs
counts.m <- do.call("cbind", countsNoR.l)
counts.df <- as.data.frame(cbind(ids,counts.m))
rownames(counts.df) <- counts.df[,1]
counts.df[,1] <- NULL
counts.df.n <- apply(counts.df,2, as.numeric)
rownames(counts.df.n) <- rownames(counts.df)

# keep only samples with at least 1 Mio. reads
counts.df.n <- counts.df.n[,colSums(counts.df.n) > 1e+06]

# filter out loci with less than 2 counts per million in at least 3 samples
isexpr <- rowSums(cpm(counts.df.n)>2) >= 3
counts.ex <- counts.df.n[isexpr,]

# TMM normalization (within experiment)
x <- calcNormFactors(counts.ex)

libSize <- colSums(counts.ex)
normalizedCounts <- log2(10e+06 *(counts.ex / (libSize * x))+1)

# plot TMM normalized expression values distribution
plotDensities(normalizedCounts)

# save the TMM normalized within experiment data (no replicates) 
write.table(normalizedCounts,file.path(paste(wd,"TMM_noReps.txt",sep="")),sep="\t")



### Quantile normalization ###

# move the TMM normalized within experiment files into one folder (for both experiments with or without replicates) and specify this as working directory
wd <- "/media/dianacs/data1/Express_data/RNAseq/SRP/RNAseq_part1/subread_aligned/counts_part1and2/RNAseq_quantNormBetween/input_NormWithin/" 

# read in all TMM normalized files and store as list
fls <- list.files(wd)
norm.rseq = sapply(fls, function(x)read.table(file.path(paste(wd,x,sep="")), header=TRUE, sep='\t'))


# count and plot the actual number of samples per experiment (remained after QC filtering)
ex.no <-cbind(as.numeric(summary(norm.rseq)[,"Length"]),as.numeric(summary(norm.rseq)[,"Length"]))
colnames(ex.no) <- c("NormSampleNo","same")
rownames(ex.no) <- rownames(summary(norm.rseq))

dev.off()
pdf(file=paste(file.path(wd),"NoNormRNAseqExp.pdf",sep=""),width=11, height=8.5,pointsize=12, paper='special')

par(mar=c(8,8,4,4))
xLabLocs <- barplot(t(ex.no[,1]), space=0.4, xaxt='n', ann=FALSE,ylab='Number of Normalized (TMM, voom) RNA-seq samples per Experiment',col='gold',border=NA,width = 1)
axis(1, cex.axis=0.75, las=2, at=xLabLocs, labels=rownames(ex.no))
text(x= xLabLocs, y= ex.no[,1]+1, labels=as.character(ex.no[,1]), xpd=TRUE,cex=0.75,col="magenta")
dev.off()

# make a data frame holding the TMM normalized data for all experiments
# use only common transcript IDs across all experiments (! QC locus filter)
for (l in names(norm.rseq)){
  norm.rseq[[l]]$geneID <- row.names(norm.rseq[[l]])
}

all.IDs <-list()
for (l in names(norm.rseq)){
  all.IDs[[l]] <- row.names(norm.rseq[[l]])
}

comm.IDs <-Reduce(intersect,all.IDs[4:length(all.IDs)])

commIDs.arrays <-list()
for (l in names(norm.rseq)){
  commIDs.arrays[[l]] <- norm.rseq[[l]][which(row.names(norm.rseq[[l]]) %in% comm.IDs),]
}

merge.norm.rseq <- Reduce(function(x,y) merge(x,y, all=T,by.x='geneID',by.y='geneID'),commIDs.arrays, accumulate=F)

row.names(merge.norm.rseq) <- merge.norm.rseq$geneID
row.names(merge.norm.rseq)[1:5]
merge.norm.rseq$geneID <- NULL

merge.norm.rseq <- as.matrix(apply(merge.norm.rseq,2,as.numeric))
merge.norm.rseq[is.na(merge.norm.rseq)] <- 0
row.names(merge.norm.rseq) <- row.names(merge.norm.rseq)

# normalize between experiments (quant. norm)
merge.norm.rseq.qnorm <- normalizeQuantiles(merge.norm.rseq)
merge.norm.rseq.qnorm[1:5,1:5]

# plot expression levels distribution before and after quantile normalization
plot(density(tst.qnorm),col="blue",xlab="mRNA levels (log2)")
lines(density(tst),col="magenta")

# save the TMM normalized within experiment followed by Quantile normalized between experiment data (all experiments with or without replicates)
write.table(tst.qnorm,file.path("/media/dianacs/data1/Express_data/RNAseq/SRP/RNAseq_part1/subread_aligned/counts_part1and2/RNAseq_quantNormBetween/QuantileNorm_RNAseq.txt"),sep="\t")

### average replicates ###

# read in the normalized data (TMM and Quantile Norm.)
rseq.qnorm <- read.table(file="/media/dianacs/data1/Express_data/RNAseq/SRP/RNAseq_part1/subread_aligned/counts_part1and2/RNAseq_quantNormBetween/QuantileNorm_RNAseq.txt",header=TRUE,row.names=1,sep="\t")

# specify the working directory for each experiment with or without replicates [contains the .sra files (not used here) and the annotation file (manually curated and assembeled meta-data)]
wds <- c("/media/dianacs/data1/Express_data/RNAseq/SRP/EERAD124/EERAD124/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP002/SRP002284/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP003/SRP003630/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP005/SRP005483/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP009/SRP009466/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP010/SRP010062/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP010/SRP010084/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP010/SRP010563/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP014/SRP014795/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP017/SRP017044/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP018/SRP018835/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP031/SRP031856/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP040/SRP040659/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP040/SRP040944/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/SRP049/SRP049231/",
         "/media/dianacs/data1/Express_data/RNAseq/SRP/RNAseq_part1/subread_aligned/counts_part1and2/single_doubleSample_exp/")

# make a list with annotation for each experiment as elements
anno.l <- list()
for (w in 1:length(wds)){
  wd <- wds[w]
  anno <- read.table(file.path(paste(wd,"annotation.txt",sep="")), header=TRUE, sep='\t')[,c("Name","File")]
  anno$SampleID <- gsub(".fastq.gz_subread.sam_counts.txt","",anno$File)
  anno.l[[w]] <- anno
  names(anno.l)[w] <- unlist(strsplit(wd,"/"))[9]
}

# store all annotations in a data frame 
anno.all <- do.call("rbind",anno.l)

# for samples with replicates assign to each sample name the replication info
anno.all$Name <- paste(sapply(rownames(anno.all),function(x) unlist(strsplit(x,"\\."))[1]),anno.all$Name,sep="_")

# for samples without replicates assign unique sample name
anno.all[280:360,"Name"] <- paste(anno.all[280:360,"SampleID"],anno.all[280:360,"Name"],sep="_")
anno.all$Name <- gsub("subread_aligned_","",anno.all$Name)
rownames(anno.all)[280:360] <- anno.all[280:360,"SampleID"]

# keep only those sample annotation which match the sampleID is the normalized matrix
commSam <- intersect(anno.all$SampleID,colnames(rseq.qnorm))
anno.all <- anno.all[which(anno.all$SampleID %in% commSam),]  
rseq.qnorm <- rseq.qnorm[,which(colnames(rseq.qnorm) %in% commSam),]

# infer replication for each experiment
rep <- anno.all$Name
names(rep) <- rep
names(rep)[1:3]
uninames <- unique(rep)
names(uninames) <- uninames

# make list for each experiment with the corresponding sample replication
reps <- list()
for (n in names(uninames)){
  reps[n] <- paste(anno.all[which(anno.all$Name == n),"SampleID"],collapse=";")#grep(uninames[[n]],colnames(tst.qnorm), value=TRUE
}


LOCUStoPNtable <- do.call(rbind,reps)
colnames(LOCUStoPNtable) <- "probeNames"

LOCUStoPN <- sapply(LOCUStoPNtable[,"probeNames"], function(x) strsplit(x, ";", fixed = TRUE)); names(LOCUStoPN) <- rownames(LOCUStoPNtable)

# average expression values for replicates
locusData <- list()
for (sample in row.names(rseq.qnorm)) {
	sub <- rseq.qnorm[sample,]
	locusData[[sample]] <- sapply(LOCUStoPN, function(x) ifelse(length(x)>1, mean(as.numeric(sub[x])), as.numeric(sub[x]))) 
}


array.av <- do.call("rbind", locusData)


# save normalized and averaged expression values for all experiments
write.table(array.av,"/media/dianacs/data1/Express_data/RNAseq/SRP/RNAseq_part1/subread_aligned/counts_part1and2/RNAseq_quantNormBetween/ChlamyRNAseq_NormBewteen_AvrgRep.txt",sep="\t")



####################
### microarrays  ###
####################

library(limma)
library(corrplot)

### Between experiment Quantile normalization (limma) ###


workDir <- "/media/dianacs/data1/Express_data/Microarrays/processed/" #"/media/dianacs/data1/Express_data/Chlamy_Arrays/" 

# read in the limma normalized within experiment files for each microarray experiment
filelist <- list.files(path=workDir, pattern = "*meanNormWithin.txt")
allarrays = sapply(filelist, function(x)read.table(file.path(paste(workDir,x,sep="")), header=T, sep='\t'))

# count sample number per experiment
array.count <- 0
for (l in names(allarrays)){
  array.count <- sum(array.count + length(allarrays[[l]]))
}

# use only common transcript IDs across all experiments
for (l in names(allarrays)){
  allarrays[[l]]$geneID <- row.names(allarrays[[l]])
}

all.IDs <-list()
for (l in names(allarrays)){
  all.IDs[[l]] <- row.names(allarrays[[l]])
}

comm.IDs <-Reduce(intersect,all.IDs[4:length(all.IDs)])

# only use experiments interrogating at least 10000 transcripts
ssetarrays <- allarrays[4:length(allarrays)]

commIDs.arrays <-list()
for (l in names(ssetarrays)){
  commIDs.arrays[[l]] <- ssetarrays[[l]][which(row.names(ssetarrays[[l]]) %in% comm.IDs),]
}

# count sample number per experiments interrogating at least 10000 transcripts
array.count <- 0
for (l in names(commIDs.arrays)){
  array.count <- sum(array.count + length(commIDs.arrays[[l]]))
}

# make a data frame hodling within  experiment normalized data for all experiments interrogating at least 10000 transcripts
merge.allarrays <- Reduce(function(x,y) merge(x,y, all=T,by.x='geneID',by.y='geneID'),commIDs.arrays, accumulate=F)
row.names(merge.allarrays) <- merge.allarrays$geneID
merge.allarrays$geneID <- NULL

merge.allarrays <- as.matrix(apply(merge.allarrays,2,as.numeric))
merge.allarrays[is.na(merge.allarrays)] <- 0
row.names(merge.allarrays) <- row.names(merge.allarrays)

# normalize between experiments (quant. norm.)
merge.allarrays.qnorm <- normalizeQuantiles(merge.allarrays)

# plot expression levels distribution before and after quantile normalization
plot(density(tst.qnorm),col="blue",xlab="mRNA levels (log2)")
lines(density(tst),col="magenta")


### average replicates ###

# infer replication for each experiment
rep <- sapply(colnames(tst.qnorm),function(x) paste(unlist(strsplit(x,"_"))[1:(length(unlist(strsplit(x,"_")))-1)],collapse="_"))


names(rep) <- rep
names(rep)[1:3]
uninames <- unique(rep)
names(uninames) <- uninames

# make list for each experiment with the corresponding sample replication
reps <- list()
for (n in names(uninames)){
  reps[n] <- paste(grep(uninames[[n]],colnames(tst.qnorm), value=TRUE),collapse=";")
}


LOCUStoPNtable <- do.call(rbind,reps)
colnames(LOCUStoPNtable) <- "probeNames"


LOCUStoPN <- sapply(LOCUStoPNtable[,"probeNames"], function(x) strsplit(x, ";", fixed = TRUE)); names(LOCUStoPN) <- rownames(LOCUStoPNtable)

# average expression values for replicates
locusData <- list()
for (sample in row.names(tst.qnorm)) {
	sub <- tst.qnorm[sample,]
	locusData[[sample]] <- sapply(LOCUStoPN, function(x) ifelse(length(x)>1, mean(sub[x]), sub[x])) # switch mean with median if you like
}

array.av <- do.call("rbind", locusData)

# save normalized and averaged expression values for all experiments
write.table(array.av,file.path(paste(workDir,"ChlamyArrays_NormBewteen_AvrgRep.txt",sep="")),sep="\t")


##########################
### sample correlation ###
##########################

library(corrplot)

# read in the normalized and repl. averaged data
eset <- as.matrix(read.table(file="/media/dianacs/data1/Express_data/RNAseq/SRP/RNAseq_part1/subread_aligned/counts_part1and2/RNAseq_quantNormBetween/ChlamyRNAseq_NormBewteen_AvrgRep.txt",header=TRUE,row.names=1,sep="\t"))

# calculate the correlation matrix
M <- cor(eset, method="pearson")
colnames(M) <- sapply(colnames(M),function(x) unlist(strsplit(x,"_"))[1])
row.names(M) <- colnames(M)

# specify the color scheme
col2a=colorRampPalette(c("black","darkblue","blue2","blue","aquamarine"),space="rgb")(255)

# set the margins of the figures to be smaller than default
par(mar=c(2,3,3,1)) 

# plot the correlogram
corrplot(M,method="ellipse",hclust.method="average",order="hclust",
	col=col2a,type="full",
	tl.cex=0.5,tl.col="black",
	cl.pos="b")


##############################################
### prepare the input data for bash aracne ###
##############################################

library("genefilter")
library("reshape2")


### zebrafish embryo ###

eset <- as.matrix(read.csv(file="/media/dianacs/data2/ZebrafishEmbryo/subreadCounts/finalSubreadCounts/Archive/normZebrafishEmb.csv",header=TRUE,row.names=1))
dim(eset)

########################
########################


### microarrays ###

# read in the normalized and repl. averaged data
eset <- as.matrix(read.table(file="/media/dianacs/data1/Express_data/Microarrays/processed/ChlamyArrays_NormBewteen_AvrgRep.txt",header=TRUE,row.names=1,sep="\t"))


### RNA-seq ###

# read in the normalized and repl. averaged data
eset <- as.matrix(read.table(file="/media/dianacs/data1/Express_data/RNAseq/SRP/RNAseq_part1/subread_aligned/counts_part1and2/RNAseq_quantNormBetween/ChlamyRNAseq_NormBewteen_AvrgRep.txt",header=TRUE,row.names=1,sep="\t"))

# variance based filtering of loci which do not show variation across samples
m.var <- varFilter(eset, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)  

dim(m.var)

# prepare the input matrix for bash aracne
input.ar <- matrix(nrow=nrow(m.var),ncol=ncol(m.var)+2)
str(rownames(m.var))

input.ar[,1] <- gsub(".v5.5","",rownames(m.var))


input.ar[,2] <- "---"

input.ar[,3:ncol(input.ar)] <- m.var

colnames(input.ar) <- c("IDs","Anno",colnames(m.var))

dim(input.ar)
# remove loci with missing values
input.ar <- as.data.frame(na.omit(input.ar))
dim(input.ar)

# save bash aracne input matrix


#### zebrafish embryo

write.table(input.ar,"/media/dianacs/data2/ZebrafishEmbryo/subreadCounts/finalSubreadCounts/Archive/inAracne_ZebrafishEmbRNAseq.txt",sep="\t",quote=FALSE,row.names=FALSE)


##########

# microarrays
write.table(input.ar,"/media/dianacs/data1/Express_data/Microarrays/processed/inAracne_ChlamyUarrays.txt",sep="\t",quote=FALSE,row.names=FALSE)

# RNAseq
write.table(input.ar,"/media/dianacs/data1/Express_data/RNAseq/SRP/RNAseq_part1/subread_aligned/counts_part1and2/RNAseq_quantNormBetween/inAracne_ChlamyRNAseq.txt",sep="\t",quote=FALSE,row.names=FALSE)

##################################
### Transcription Factors info ###
##################################

# TFs assembled from public DBs have IDs corresponding to different genome releases
# map the TF IDs to the latest chlamy release (JGIv5.5)

# read in the latest chlamy loci annotation (JGIv5.5)
cre.anno <- read.csv("/media/dianacs/data1/chlamyJGI/Creinhardtii_download/Creinhardtii/v5.5/annotation/creanno.csv",header=FALSE)
colnames(cre.anno) <- c("PhytozomeTranscrID","PhytozomeGene","PhytozomeTranscr","PhytozomeProt","PFAM","Panther","KOG","KEGGec","KEGGortho","GO","TAIR10name","TAIR10symb","TAIR10defline")
head(cre.anno,10)
# read in the file with the info for conversion between genome releases
conv <- read.table("/media/dianacs/data1/chlamyJGI/Creinhardtii_download/Creinhardtii/v5.5/annotation/ChlamydomonasTranscriptNameConversionBetweenReleases.Mch12b.txt",header=FALSE,sep="\t")
colnames(conv) <- c("JGIv5_5","JGIv3_1","Genbank","JGIv4","JGIv4_3","Aug5","Aug6","JGIv5_3_1")

# read in the latest chlamy gene names (JGIv5.5)
gene <- read.table("/media/dianacs/data1/chlamyJGI/Creinhardtii_download/Creinhardtii/v5.5/annotation/Creinhardtii_281_v5.5.geneName.txt",header=FALSE,sep="\t")
colnames(gene) <- c("JGIv5_5","genename")

# read in the TFs IDs
tf <- read.table("/media/dianacs/data1/ChlamyTFs/TF.txt",header=TRUE,sep="\t")

# map the TF IDs to the latest chlamy release (JGIv5.5)
tf$JGIv5_5 <- conv[which(tf[,1] %in% conv[,"JGIv5_3_1"]),"JGIv5_5"]
tf$gene_model <- NULL

tf <- merge(tf,gene,by="JGIv5_5",all=T,by.x="JGIv5_5",by.y="JGIv5_5")

# save the annotated TFs
write.csv(tf,"/media/dianacs/data1/ChlamyTFs/TFanno.csv")

write.table(tf$JGIv5_5,"/media/dianacs/data1/ChlamyTFs/TFaracne.txt",sep="\t",quote=FALSE,row.names=FALSE)



################################################
### prepare bash aracne output for Cytoscape ###
################################################


###zebrafish embryo


aracne <- read.table("/media/dianacs/data2/ZebrafishEmbryo/subreadCounts/finalSubreadCounts/Archive/aracneZebrafishEmbryo.txt",header=FALSE,sep="\t",na.strings = " ",fill = TRUE)


dim(aracne)
head(aracne)

##########################
##########################


# microarrays
aracne <- read.csv("/media/dianacs/data1/Express_data/Microarrays/processed/aracneUarrays.csv",header=FALSE,na.strings="")

head(aracne,5)
# RNAseq
aracne <- read.csv("/media/dianacs/data1/Express_data/RNAseq/SRP/RNAseq_part1/subread_aligned/counts_part1and2/RNAseq_quantNormBetween/aracneRNAseq.csv",header=FALSE,na.strings="")

head(aracne)

# remove headers and check for field separators in file: inAracne_ChlamyUarrays_k0.18_t0.17_e0.1.adj (output bash aracne)
aracne.l <- list()
for (l in 1:nrow(aracne)){
  line <- as.vector(aracne[l,])
  line <- line[! is.na(line)]
  linem <- as.data.frame(matrix(nrow=(length(line)-1)/2,ncol=3))
  colnames(linem) <- c("Query","Target","MI")
  linem$Query <- line[1]
  #   line.q <- strsplit(line," ")[1]
  line.tmi <- paste((line.tmi <- strsplit(line," ")[-1]))
  line.dup <- split(line.tmi,ceiling(seq_along(line.tmi)/2))
  linem$Target <- sapply(line.dup,function(x) strsplit(x," ")[1])
  linem$MI <- sapply(line.dup,function(x) strsplit(x," ")[2])
  aracne.l[[l]] <- linem
}

aracne.cytoscape <- as.matrix(do.call("rbind",aracne.l))

head(aracne.cytoscape)
dim(aracne.cytoscape)

#############zebrafish embryo
tst <- aracne.cytoscape[which(aracne.cytoscape[,3] > 0.5),]
dim(tst)
tst


write.table(as.matrix(tst),"/media/dianacs/data2/ZebrafishEmbryo/subreadCounts/finalSubreadCounts/Archive/CytoscapeIn_ZebrafishEmbryoRseq_BASHaracneNetMI05.txt",sep="\t",quote=FALSE,row.names=FALSE)
#############################

######################################################################################################################
### microarray network

# filter network nodes based on differential expression (microarrays based network); union of diffE in all comparisons
wd <- "/media/dianacs/data1/Express_data/Microarrays/processed/diffE/" 


wd <- "/media/dianacs/data1/Express_data/Microarrays/processed/diffE/" 
fls <- list.files(wd,pattern=".csv")
fls <- fls[10:18]
fls


diffE <- sapply(fls, function(x)read.table(file.path(paste(wd,x,sep="")), header=FALSE, sep=',')[1])
summary(diffE)

# total number of diffE genes from all comparisons (duplicate IDs)
sum(as.numeric((summary(diffE)[,"Length"])))

# total number of diffE genes from all comparisons (no duplicate IDs)

all.diffE <- Reduce(union,diffE)[-1]
str(all.diffE)

# make sure the correct network is loaded in the workspace >> microarray
aracne.cytoscape <- as.data.frame(aracne.cytoscape)
dim(aracne.cytoscape)
head(aracne.cytoscape)

#no MI threshold
allMIdiffE.filt <- aracne.cytoscape[which(aracne.cytoscape$Query %in% all.diffE),]
dim(allMIdiffE.filt)

# MI threshold 0.4
MIdiffE.filt <- aracne.cytoscape[which(aracne.cytoscape$MI > 0.4 & aracne.cytoscape$Query %in% all.diffE),]
dim(MIdiffE.filt)
head(MIdiffE.filt)




write.table(as.matrix(allMIdiffE.filt),"/media/dianacs/data1/Express_data/Microarrays/processed/CytoscapeIn_ChlamyUarrays_BASHaracneSigMIdiffE.txt",sep="\t",quote=FALSE,row.names=FALSE)

write.table(as.matrix(MdiffEI.filt),"/media/dianacs/data1/Express_data/Microarrays/processed/CytoscapeIn_ChlamyUarrays_BASHaracneNetMI04diffE.txt",sep="\t",quote=FALSE,row.names=FALSE)


### RNAseq network

# filter network nodes based on differential expression (microarrays based network); union of diffE in all comparisons
wd <- "/media/dianacs/data1/Express_data/RNAseq/SRP/RNAseq_part1/subread_aligned/counts_part1and2/diffE_RNAseq/" 
fls <- list.files(wd,pattern=".csv")
fls
fls <- fls[15:28]
fls


diffE <- sapply(fls, function(x)read.table(file.path(paste(wd,x,sep="")), header=FALSE, sep=',')[1])
summary(diffE)

# total number of diffE genes from all comparisons (duplicate IDs)
sum(as.numeric((summary(diffE)[,"Length"])))

# total number of diffE genes from all comparisons (no duplicate IDs)

all.diffE <- Reduce(union,diffE)[-1]
all.diffE <- gsub(".v5.5","",all.diffE)
str(all.diffE)

# make sure the correct network is loaded in the workspace >> RNAseq
aracne.cytoscape <- as.data.frame(aracne.cytoscape)
dim(aracne.cytoscape)
head(aracne.cytoscape)


#no MI threshold
allMIdiffE.filt <- aracne.cytoscape[which(aracne.cytoscape$Query %in% all.diffE),]
dim(allMIdiffE.filt)


# MI threshold 0.4
MIdiffE.filt <- aracne.cytoscape[which(aracne.cytoscape$MI > 0.4 & aracne.cytoscape$Query %in% all.diffE),]
dim(MIdiffE.filt)
head(MIdiffE.filt)


par(mfrow=c(1,2))
hist(as.numeric(allMIdiffE.filt$MI))
abline(v=0.4,col="magenta")

hist(as.numeric(MIdiffE.filt$MI))

MIdiffE.filt[MIdiffE.filt$Query=="Cre01.g000350.t1.1",]
MIdiffE.filt[MIdiffE.filt$Query=="Cre03.g187450.t1.2",]


write.table(as.matrix(allMIdiffE.filt),"/media/dianacs/data1/Express_data/RNAseq/SRP/RNAseq_part1/subread_aligned/counts_part1and2/RNAseq_quantNormBetween/CytoscapeIn_ChlamyRNAseq_BASHaracneSigMIdiffE.txt",sep="\t",quote=FALSE,row.names=FALSE)


write.table(as.matrix(MIdiffE.filt),"/media/dianacs/data1/Express_data/RNAseq/SRP/RNAseq_part1/subread_aligned/counts_part1and2/RNAseq_quantNormBetween/CytoscapeIn_ChlamyRNAseq_BASHaracneNetMI04diffE.txt",sep="\t",quote=FALSE,row.names=FALSE)




######################################################################################################################


# save the Cytoscape input (microarrays based network)
write.table(aracne.cytoscape,"/media/dianacs/data1/Express_data/Microarrays/processed/CytoscapeIn_ChlamyUarrays_BASHaracneNet.txt",
            sep="\t",quote=FALSE,row.names=FALSE)

# filter network edges based on Mutual Information (microarrays based network)
basha <- read.table("/media/dianacs/data1/Express_data/Microarrays/processed/CytoscapeIn_ChlamyUarrays_BASHaracneNet.txt",sep="\t",header=TRUE)            
basha.s <- basha[which(basha$MI > 0.4),]

# save the Cytoscape input with MI threshold (microarrays based network)
write.table(basha.s,"/media/dianacs/data1/Express_data/Microarrays/processed/CytoscapeIn_ChlamyUarrays_BASHaracneNetMI04.txt",
            sep="\t",quote=FALSE,row.names=FALSE)

          
# save the Cytoscape input (RNA-seq based network)
write.table(aracne.cytoscape,"/media/dianacs/data1/Express_data/RNAseq/SRP/RNAseq_part1/subread_aligned/counts_part1and2/RNAseq_quantNormBetween/CytoscapeIn_ChlamyRNAseq_BASHaracneNet.txt",
            sep="\t",quote=FALSE,row.names=FALSE)
      
# filter network edges based on Mutual Information (RNA-seq based network)
basha <- read.table("/media/dianacs/data1/Express_data/RNAseq/SRP/RNAseq_part1/subread_aligned/counts_part1and2/RNAseq_quantNormBetween/CytoscapeIn_ChlamyRNAseq_BASHaracneNet.txt",sep="\t",header=TRUE)            
basha.s <- basha[which(basha$MI > 0.4),]

# save the Cytoscape input with MI threshold (RNA-seq based network)
write.table(basha.s,"/media/dianacs/data1/Express_data/RNAseq/SRP/RNAseq_part1/subread_aligned/counts_part1and2/RNAseq_quantNormBetween/CytoscapeIn_ChlamyRNAseq_BASHaracneNetMI04.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


#####################################
### igraph and graph_tool(python) ###
#####################################


mthresh <- read.table("/media/dianacs/data1/Express_data/Microarrays/processed/CytoscapeIn_ChlamyUarrays_BASHaracneSigMIdiffE.txt",header=TRUE)
dim(mthresh)
mthresh$edge <- paste(mthresh[,1],mthresh[,2],sep=" ")#,mthresh[,3]
head(mthresh)
write.table(mthresh[,"edge"],"/media/dianacs/data1/Express_data/net/microarray_diffE_igraph.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)


rthresh <- read.table("/media/dianacs/data1/Express_data/RNAseq/SRP/RNAseq_part1/subread_aligned/counts_part1and2/RNAseq_quantNormBetween/CytoscapeIn_ChlamyRNAseq_BASHaracneSigMIdiffE.txt",header=TRUE)
dim(rthresh)
rthresh$edge <- paste(rthresh[,1],rthresh[,2],sep=" ")#,rthresh[,3]
head(rthresh)
write.table(rthresh[,"edge"],"/media/dianacs/data1/Express_data/net/rseq_diffE_igraph.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

# mrthesh <- as.matrix(Reduce(function(x,y) merge(x,y, all=T,by.x='edge',by.y='edge'),list(mthresh,rthresh), accumulate=F))
# dim(mrthesh)
# head(mrthesh)

###################################
# igraph                          #
# install.packages("igraph")      #
# install.packages("d3Network")   #
###################################

rm(list=ls())
ls()

library(igraph)
library(d3Network)
library(rnetcarto)
library(varhandle)
library(lattice)


###zebrafish embryo

m.g <- read.graph("/media/dianacs/data2/ZebrafishEmbryo/subreadCounts/finalSubreadCounts/Archive/igraph_ZebrafishEmbryoRseq_BASHaracneNetMI05.txt",format="ncol",directed=FALSE)
summary(m.g)
##########


#linux
#m.g <- read.graph("/media/dianacs/data1/Express_data/net/microarray_diffE_igraph.txt",format="ncol",directed=FALSE)

m.g <- read.graph("D://Express_data//net//microarray_diffE_igraph.txt",format="ncol",directed=FALSE)


m.g

# remove duplicate edges
is.simple(m.g)
m.g.s <- simplify(m.g, remove.multiple = TRUE, remove.loops = TRUE,edge.attr.comb = getIgraphOpt("edge.attr.comb"))
m.g.s
is.simple(m.g.s)

#linux
# r.g <- read.graph("/media/dianacs/data1/Express_data/net/rseq_diffE_igraph.txt",format="ncol",directed=FALSE)

r.g <- read.graph("D://Express_data//net//rseq_diffE_igraph.txt",format="ncol",directed=FALSE)

r.g

is.simple(r.g)
r.g.s <- simplify(r.g, remove.multiple = TRUE, remove.loops = TRUE,edge.attr.comb = getIgraphOpt("edge.attr.comb"))
r.g.s
is.simple(r.g.s)

# graph union
# mr.g <- m.g %u% r.g

mr.g.s <- graph.union(m.g.s,r.g.s,byname=TRUE)
mr.g.s

V(mr.g.s)$label <- seq_along(V(mr.g.s))

V(mr.g.s)[1:5]

summary(mr.g.s)

vcount(mr.g.s)

is.simple(mr.g.s)

is.connected(mr.g.s)

no.clusters(mr.g.s)

table(clusters(mr.g.s)$csize)

graph.density(mr.g.s)



write.graph(mr.g.s, "/media/dianacs/data1/Express_data/net/mergedRseqMArr_diffE_igraph.txt", format="ncol")

write.graph(r.g.s, "/media/dianacs/data1/Express_data/net/Rseq_diffE_graphml", format="graphml")

write.graph(mr.g.s, "/media/dianacs/data1/Express_data/net/mergedRseqMArr_diffE_graphml", format="graphml")


#### community detection
summary(mr.g.s)


global.ndeg <- degree(mr.g.s, v = V(mr.g.s), mode = "total", normalized = FALSE)
str(global.ndeg)
range(global.ndeg)

colnames(global.ndeg) <- "global.ndeg"

hist(global.ndeg, col="gray")
abline(v=c(10,50),col="magenta")


mr.g.s <- set_vertex_attr(mr.g.s, "nodeDegree", index = V(mr.g.s), global.ndeg)
list.vertex.attributes(mr.g.s)
get.vertex.attribute(mr.g.s,'nodeDegree',1:5)

set.seed(333)
lo <- cluster_louvain(mr.g.s)
str(lo)

lo$modularity
sizes(lo)
names(sizes(lo))

modularity(lo)

mr.g.s <- set_vertex_attr(mr.g.s, "louvain", index = V(mr.g.s), lo$membership)
list.vertex.attributes(mr.g.s)
get.vertex.attribute(mr.g.s,'louvain',1:5)
mr.g.s



#######################################################################################################################
#######################################################################################################################

#select GeneIDs for a certain community

cur <- read.table("clipboard", sep = "\t", header=FALSE)
cur

cur.g <- induced.subgraph(mr.g.s,vids=unlist(neighborhood(mr.g.s,order=1,nodes=V(mr.g.s)$name %in% cur[,1])))
summary(cur.g)

ShortPth <- get.all.shortest.paths(cur.g, "Cre10.g434350.t1.1","Cre09.g390023.t1.1")    
# ShortPth <- get.shortest.paths(cur.g, which(V(cur.g)$name == "Cre10.g434650.t1.1"), which(V(cur.g)$name == "Cre09.g390023.t1.1"))    

as.matrix(unlist(ShortPth[[1]]))
E(cur.g)$color <- "grey"
E(cur.g)$width <- 1
V(cur.g)$color <- "grey"
V(cur.g)$size <- 5

for (p in ShortPth$res) { 
  E(cur.g, path=p)$color <- "orange"
  E(cur.g, path=p)$width <- 5
  V(cur.g)[p]$color <- "orange"
  V(cur.g)[p]$size <- 15
}

plot(cur.g,vertex.frame.color=NA, vertex.label=NA,layout=layout.fruchterman.reingold)

tkplot(cur.g,vertex.frame.color=NA, vertex.label=NA,layout=layout.fruchterman.reingold)




write.graph(cur.g, "D://Express_data//net//SelCTR_mergedRseqMArr_diffE_igraph.txt", format="ncol")



set.seed(333)
lo.cur <- cluster_louvain(cur.g)
str(lo.cur)

modularity(lo)

cur.g <- set_vertex_attr(cur.g, "louvain", index = V(cur.g), lo.cur$membership)
list.vertex.attributes(cur.g)
get.vertex.attribute(cur.g,'louvain',1:5)

ids.a <- names(V(induced.subgraph(cur.g,vids=unlist(neighborhood(cur.g,order=0,nodes=V(cur.g)$louvain == 1)))))
str(ids.a)
ids.b <- names(V(induced.subgraph(cur.g,vids=unlist(neighborhood(cur.g,order=0,nodes=V(cur.g)$louvain == 2)))))
ids.c <- names(V(induced.subgraph(cur.g,vids=unlist(neighborhood(cur.g,order=0,nodes=V(cur.g)$louvain == 3)))))
ids.d <- names(V(induced.subgraph(cur.g,vids=unlist(neighborhood(cur.g,order=0,nodes=V(cur.g)$louvain == 4)))))


#######################################################################################################################
#######################################################################################################################

lo_cl <- as.matrix(membership(lo))
lo_cl <-cbind(rownames(lo_cl),lo_cl)
rownames(lo_cl)<-NULL
colnames(lo_cl) <- c("GeneID","LouvainID_CuRel")
head(lo_cl)

range(unique(lo_cl[,1]))
table(lo_cl[,1])


write.graph(mr.g.s, "/media/dianacs/data1/Express_data/net/chlamy_netoxi_mergedRseqMArr_diffE", format="ncol")
write.table(lo_cl, "/media/dianacs/data1/Express_data/net/noaLouvainCl_chlamy_netoxi_mergedRseqMArr_diffE.txt",quote=FALSE)
write.table(ndeg, "/media/dianacs/data1/Express_data/net/noaNdeg_chlamy_netoxi_mergedRseqMArr_diffE.txt",quote=FALSE)





###

global.ndeg <- degree(mr.g.s, v = V(mr.g.s), mode = "total", normalized = FALSE)

comm.carto <- list()
for (comm in names(sizes(lo))){
#    cat(paste("Community: ",comm, sep=''))
  comm.g <- induced.subgraph(mr.g.s,which(V(mr.g.s)$louvain==comm))
#    cat(paste("Community: ",comm,"\n",summary(comm.g), sep=''))
  ndeg <- degree(comm.g, v = V(comm.g), mode = "total", normalized = FALSE) 
  
  v.carto <- list()
  
  for(v in names(V(comm.g))){  

#     cat(paste("Community: ",comm," Vertex: ",v,sep=''))
#     v.comm <- comm
# 
#     v=names(V(comm.g))[1]
#     v
#     ndeg[v]
#     global.ndeg[v]
    
    zi <- unname((ndeg[v] - mean(ndeg))/sd(ndeg))

#     cat(paste("Community: ",comm, "Vertex: ",v," within_conn: ",zi,sep=''))
    parti <- unname(1-((ndeg[[v]]/global.ndeg[[v]])^2))
    #     cat(paste("Community: ",comm, "Vertex: ",v," modif.participation: ",pi,sep=''))
    

    role <- ifelse(zi<2.5 & parti<0.05,"ultraperipheral",
		  ifelse(zi<2.5 & parti>=0.05 & parti<0.625,"peripheral",
		    ifelse(zi<2.5 & parti>=0.625 & parti<0.8,"connector",
		      ifelse(zi<2.5 & parti>=0.8,"kinless",
			ifelse(zi>2.5 & parti<0.3,"provincialHub",
			  ifelse(zi>2.5 & parti>=0.3 & parti<0.75,"connectorHub",
			    ifelse(zi>2.5 & parti>=0.75,"kinlessHub","na")))))))    	


    v.carto[[v]] <- c(v,comm,zi,parti,role)
    
  }
  
#   pi.ori.fi <- 1-(pi.ori.fi+pi.ori)
  comm.carto[[comm]] <- do.call("rbind",v.carto) 
#   comm.carto[[comm]][[v,4]] <- pi.ori.fi 
}


comm.carto <- as.data.frame(do.call("rbind",comm.carto))

colnames(comm.carto) <- c("geneID","community","within_conn","participation","role")

comm.carto[1:5,]
table(comm.carto$role)

write.table(comm.carto, "/media/dianacs/data1/Express_data/net/commcartoLouvain_chlamy_netoxi_mergedRseqMArr_diffE.txt",sep="\t",quote=FALSE,row.names=FALSE)

cart <- read.table("/media/dianacs/data1/Express_data/net/commcartoLouvain_chlamy_netoxi_mergedRseqMArr_diffE.txt",sep="\t",header=TRUE)
head(cart)
str(cart)
table(cart$role)
table(cart$community)

dev.off()

myPanel <- function(x,y, ...) {
    panel.xyplot(x,y, ...)
#     panel.abline(h=2.5, ...)
#     panel.abline(v=c(0,0.625,0.8), h=rep(),...)
}

       
xyplot(within_conn~participation, group=role, data=cart,pch=18,
       panel=panel.superpose,panel.groups=myPanel,
       auto.key=list(space="right"),
       ylab="Within-module degree",
       xlab="Participation coefficient")
       
levels(cart$role)
con.hub <- nrow(cart[cart$role=="connectHub",])*100/nrow(cart)
con.hub

u.perif <- nrow(cart[cart$role=="ultraperif",])*100/nrow(cart)
u.perif

con <- nrow(cart[cart$role=="connect",])*100/nrow(cart)
con

perif <- nrow(cart[cart$role=="perif",])*100/nrow(cart)
perif


### expression profiles Rseq


rseq <- read.table("/media/dianacs/data1/Express_data/RNAseq/SRP/RNAseq_part1/subread_aligned/counts_part1and2/RNAseq_quantNormBetween/ChlamyRNAseq_NormBewteen_AvrgRep.txt")

rseq[1:3,1:3]
row.names(rseq) <- gsub(".v5.5","",row.names(rseq))
row.names(rseq)[1:3]

colnames(rseq)
cn <- colnames(rseq)
cn[1:3]

cn.n <- sapply(cn,function(x) unlist(strsplit(x,"_"))[1])
cn.n[1:5]

colnames(rseq) <- cn.n
head(rseq,10)

Cre10.g434350.t1.1
Cre12.g536000.t1.2
Cre09.g388850.t1.1
Cre13.g570600.t1.2
Cre12.g505350.t1.2
Cre09.g410050.t1.2
Cre02.g145100.t1.1

ids <- read.table("clipboard", sep = "\t", header=F)
ids
ids.rseq <- as.data.frame(rseq[which(row.names(rseq) %in% ids[,1]),])
head(ids.rseq,10)
class(ids.rseq)
rownames(ids.rseq)[1:4]

library(gplots)
col1=colorRampPalette(c("green","green4","gray","violet","purple"))

mini <-min(ids.rseq)
maxi <-max(ids.rseq)

mini
maxi



br=seq(mini+2,maxi-2,0.5)
br

heatmap.2(as.matrix(ids.rseq),cexRow=0.65,xaxt="n",
          trace="none",col=col1, margins = c(8, 8),
          breaks=br)

###

am <- get.adjacency(mr.g.s,sparse=FALSE)
str(am)
class(am)
am[1:3,1:3]

atBeginning <- Sys.time()
cat(paste("netcarto ","\n", sep=''))
then <- Sys.time()
cart <- netcarto(am)
write.table(cart[[1]], "/media/dianacs/data1/Express_data/net/netcarto_chlamy_netoxi_mergedRseqMArr_diffE.txt",quote=FALSE,sep="\t",row.names=FALSE)
print(Sys.time()-then)



cart <- read.table("/media/dianacs/data1/Express_data/net/netcarto_chlamy_netoxi_mergedRseqMArr_diffE.txt",sep="\t",header=TRUE)
dim(cart)
table(cart$role)
table(cart$module)



plot(cart[,"participation"],cart[,"connectivity"],pch=18,col="gray",xlim=c(0,1),xlab="Participation coefficient, P",ylab="Within-module degree, z")
abline(h=2.5,col="green")
abline(v=c(0,0.625,0.8),col="green")

n.mod <- unique(cart$module)
n.mod

table(cart$module)

levels(cart$role)

con.hub <- nrow(cart[cart$role=="Connector Hub",])*100/nrow(cart)
con.hub

u.perif <- nrow(cart[cart$role=="Ultra peripheral",])*100/nrow(cart)
u.perif

con <- nrow(cart[cart$role=="Connector",])*100/nrow(cart)
con

perif <- nrow(cart[cart$role=="Peripheral",])*100/nrow(cart)
perif

##############################




cre.anno[cre.anno[["PhytozomeTranscr"]] == "Cre02.g143250.t1.2",c("PhytozomeTranscr","TAIR10name","TAIR10defline")]

######


# list.vertex.attributes(subg)
# range(V(subg)$'nodeDegree')
# get.vertex.attribute(subg,'nodeDegree',1:5)
# hist(V(subg)$'nodeDegree')
# 
# 
# # fc <- cluster_fast_greedy(subg)
# # # membership(fc)
# # sizes(fc)
# # 
# # 
# # wt <- cluster_walktrap(subg)
# # # membership(wt)
# # sizes(wt)
# 
# 
# lo <- cluster_louvain(subg)
# sizes(lo)
# typeof(lo)
# 
# sizes(lo) > 25
# 
# as.matrix(membership(lo)[1:10])
# 
# subg <- set_vertex_attr(subg, "louvain", index = V(subg), lo$membership)
# list.vertex.attributes(subg)
# get.vertex.attribute(subg,'louvain',1:5)
# subg
# 
# 
# lo_cl <- as.matrix(membership(lo))
# head(lo_cl)
# colnames(lo_cl) <- "LouvainID"
# 
# 
# ndeg.subg <- as.matrix(degree(subg, v = V(subg), mode = "total", normalized = FALSE))
# head(ndeg.subg)
# colnames(ndeg.subg) <- "ndeg"
# range(ndeg.subg[,"ndeg"])



#######



ndeg <- degree(mr.g.s, v = V(mr.g.s), mode = "total", normalized = FALSE)
str(ndeg)
degree_distribution(mr.g.s, cumulative = FALSE)
hist(ndeg)

###########################################
###########################################
# subset graph based on vertex name/ ndeg

ids <- read.table("clipboard", sep = "\t", header=TRUE,row.names=1)
str(ids)
rownames(ids)


#subg <- delete.vertices(mr.g.s, ndeg < 250)


which(V(mr.g.s)$name %in% rownames(ids))

subg <- induced.subgraph(mr.g.s,vids=unlist(neighborhood(mr.g.s,order=1,nodes=V(mr.g.s)$name %in% rownames(ids))))
                        
summary(mr.g.s)
summary(subg)

lo.subg <- cluster_louvain(subg)

sizes(lo.subg)

modularity(lo.subg)

lo.subg$membership[1:5]

subg <- set_vertex_attr(subg, "louvain", index = V(subg), lo.subg$membership)


get.vertex.attribute(subg,'louvain',1:5)


lo.subg_cl <- as.matrix(membership(lo.subg))
lo.subg_cl <-cbind(rownames(lo.subg_cl),lo.subg_cl)
rownames(lo.subg_cl)<-NULL
colnames(lo.subg_cl) <- c("GeneID","LouvainID_CuRel")

head(lo.subg_cl)



range(unique(lo.subg_cl[,2]))
table(lo.subg_cl[,2])


write.graph(subg, "D://Express_data//net//chlamy_netoxi_mergedRseqMArr_diffE_CuRel", format="ncol")
write.table(lo.subg_cl, "D://Express_data//net//noaLouvainCl_chlamy_netoxi_mergedRseqMArr_diffE_CuRel.txt",quote=FALSE,row.names=FALSE)


#########################
degree_distribution(subg, cumulative = FALSE)
ndeg.subg <- degree(subg, v = V(subg), mode = "total", normalized = FALSE)
hist(ndeg.subg)

cluster_louvain(subg)

fc <- cluster_fast_greedy(subg)
plot_dendrogram(fc,xlab=NA)




# https://rpubs.com/kaz_yos/igraph-individuals
# Check the transitivity of a graph (probability that the adjacent vertices of a vertex are connected)
(transitivityDat <- transitivity(mr.g.s, type = "localaverage",isolates = "zero"))

# Transitivity of a random graph of the same size
g <- erdos.renyi.game(vcount(mr.g.s), ecount(mr.g.s), type="gnm")
transitivity(g)


# Average path length between any two given nodes
(averagePathLength <- average.path.length(mr.g.s))

# Community structure detection based on edge betweenness
atBeginning <- Sys.time()
cat(paste("Community structure detection based on edge betweenness ","\n", sep=''))
then <- Sys.time()
communityEdgeBetwn <- edge.betweenness.community(mr.g.s)
print(Sys.time()-then)


# plots
# Set the seed to get the same result
set.seed("4559")
## Add community indicating background colors
plot(igraphDat,
     vertex.color = communityEdgeBetwn$membership, vertex.size = log(degree(igraphDat) + 1),
     mark.groups = by(seq_along(communityEdgeBetwn$membership), communityEdgeBetwn$membership, invisible))

## Annotate
title("Stanford Facebook data", sub = "http://snap.stanford.edu/data/egonets-Facebook.html")
text(x = -1, y = -1, labels = sprintf("Average path length: %.2f\nTransitivity: %.2f",
                         averagePathLength, transitivityDat))



adjacent.triangles (mr.g.s, vids = V(mr.g.s))

# V(mr.g.s)$name

# V(mr.g.s)$label <- seq_along(V(mr.g.s))

# plots cont.

plot.igraph(mr.g.s,vertex.label=NA, vertex.size = 0,vertex.label=NA, edge.arrow.size=0,vertex.shape="none")


E(mr.g.s)$weight <- runif(ecount(mr.g.s))
head(E(mr.g.s)$weight)

E(mr.g.s)$color <- "grey"
E(mr.g.s)[weight > 0.95]$color <- "magenta"

dev.off()
plot(mr.g.s, vertex.size = log(degree(mr.g.s) + 1, vertex.label=NA, layout=layout.kamada.kawai,edge.width=2+3*E(mr.g.s)$weight))

plot(degree.distribution(mr.g.s),pch=18,col="grey")



### handle network node attribute files ###

mapm.lev<- read.table("/media/dianacs/data1/Express_data/net/mapman_RaracneTF1N.txt",sep="\t",header=TRUE)
mapm.lev <- mapm[mapm[,2] == "1",]

mapm.l <- list()
for (l in 1:nrow(mapm.lev)){
  ids <- strsplit(as.character(mapm.lev[l,5]),",")
  mapm.anno <- data.frame(ids)
  mapm.anno$pathw <- rep(as.character(mapm.lev[l,1]),length(ids))
  colnames(mapm.anno) <- c("ID","pathw")
  mapm.l[[l]] <- mapm.anno
}


mapm.noa <- do.call("rbind",mapm.l)
mapm.noa$ID <- gsub(" ","",mapm.noa$ID)

write.table(mapm.noa,"/media/dianacs/data1/Express_data/net/mapmanL1.noa.txt",sep="\t",quote=FALSE)

#####################
## copper,silver diffE

rm(list = ls())
noa <- read.table("/media/dianacs/data1/Express_data/net/anno.txt",sep="\t",header=TRUE)
dim(noa)
head(noa)

table(noa$ID,noa$diffE)[1:3,-1]

cs <- table(noa$ID,noa$diffE)[,-1]
str(cs)
head(cs)

cs <- cbind(cs,chartr("0,1","A,B",cs[,"Copper"]),chartr("0,1","A,B",cs[,"Silver"]))
head(cs)

cs <- as.data.frame(cs)
head(cs)

cs$assign <- paste(cs[,3],cs[,4],sep="")

cs[which(cs$assign == "BA"),5] <- "Copper" 
cs[which(cs$assign == "BB"),5] <- "CopperSilver" 
cs[which(cs$assign == "AA"),5] <- "other" 
cs[which(cs$assign == "AB"),5] <- "Silver" 

head(cs)


write.table(cs,"/media/dianacs/data1/Express_data/net/CopperSilver_diffEanno.txt",sep="\t",quote=FALSE,col.names=TRUE)


##########


#####################
## transporters

rm(list = ls())
noa <- read.table("/media/dianacs/data1/Express_data/net/Transp.txt",sep="\t",header=TRUE, na.strings="N/A")
dim(noa)
head(noa)
str(noa)

noa[["Family.ID"]] <- as.character(noa[["Family.ID"]])
noa[["Family.Name"]] <- as.character(noa[["Family.Name"]])
noa[["Transporter.Type"]] <- as.character(noa[["Transporter.Type"]])


noa <- na.omit(noa)
dim(noa)

row.names(noa)<-noa$GeneID
row.names(noa)[1:5]
noa$GeneID <- NULL

head(noa)
dim(noa)

locusData <- list()
for (r in row.names(noa)) {
	row <- unlist(strsplit(r,", ",fixed=TRUE))
	rowcol <- list()
	for (i in 1:length(row)){
	  rowcol[[i]] <- cbind(row[i],noa[r,])
	}
	locusData[[r]] <- do.call("rbind",rowcol)	
}

transp <- do.call("rbind", locusData)
dim(transp)
row.names(transp) <- NULL
colnames(transp) <- c("geneID","TpFam","TpFamName","TpType")
head(transp)



write.table(transp,"/media/dianacs/data1/Express_data/net/noaTp.txt",sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)


### is TF or Tp

tf <- read.table("/media/dianacs/data1/Express_data/net/noaTF.txt",header=TRUE)
tf <- tf[,c("geneID","isTF")]
colnames(tf) <- c("geneID","is")

tp <- read.table("/media/dianacs/data1/Express_data/net/noaTp.txt",sep="\t",header=TRUE)
tp <- tp[,c("geneID","isTp")]
colnames(tp) <- c("geneID","is")

head(tf)
head(tp)


tftp <- rbind(tf,tp)
head(tftp

table(tftp$is)
write.table(tftp,"/media/dianacs/data1/Express_data/net/noaIS_TFTP.txt",sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)


# GO analysis
source("http://bioconductor.org/biocLite.R")
biocLite("AnnotationHub")
biocLite("GO.db")


library("topGO")

geneID2GO <- readMappings("D://chlamyJGI//Creinhardtii_download//Creinhardtii//v5.5//annotation//gene_association.txt")
str(head(geneID2GO))

geneNames <- names(geneID2GO)
head(geneNames)

ids <- read.table("clipboard", sep = "\t", header=TRUE)
str(ids[,1])


geneList <- factor(as.integer(geneNames %in% ids[,1]))
names(geneList) <- geneNames
str(geneList)
levels(geneList)
geneList[1:6]

GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata

a <- genes(GOdata) ## obtain the list of genes
head(a)

numGenes(GOdata)


sg <- sigGenes(GOdata)
str(sg)
numSigGenes(GOdata)
graph(GOdata)

ug <- usedGO(GOdata)
head(ug)


num.ann.genes <- countGenesInTerm(GOdata, ug) ## the number of annotated genes
table(num.ann.genes)

ann.genes <- genesInTerm(GOdata, ug) ## get the annotations

#ann.score <- scoresInTerm(GOdata, ug)

ann.score <- scoresInTerm(GOdata, ug, use.names = TRUE)
str(ann.score)


gosig <- termStat(GOdata, ug)
gosig[gosig$Significant > 0,]


resultFisher <- runTest(GOdata,algorithm="classic",statistic="fisher")

allRes <- GenTable(GOdata, classicFisher = resultFisher,
                  ranksOf = "classicFisher", topNodes = 10)


allRes

showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all')

printGraph(GOdata, resultFisher, firstSigNodes = 5, fn.prefix = "sampleFile", useInfo = "all", pdfSW = TRUE)

########################
## PATHWAY ENRICHMENT ##
## Fisher test        ## 
########################
require(XLConnect)

library(corrgram)
library(ellipse)
library(gplots)
library('plotrix')
library(reshape)
library(wordcloud)
library(RColorBrewer)

### geneIDs should be stored as a list 

### extract in R separate geneIDs for:
#     - communities
#     - node role
#     - Co/Ag manual separation 
### store as list elements


### quick&dirty testing input



# ids.a <- read.table("clipboard", sep = "\t", header=TRUE)
# ids.b <- read.table("clipboard", sep = "\t", header=TRUE)
# ids.c <- read.table("clipboard", sep = "\t", header=TRUE)
# ids <- list(as.character(ids.a[,1]),as.character(ids.b[,1]),as.character(ids.c[,1]))
# names(ids) <- c("Clust1","Clust2","Clust3")
# summary(ids)
 
### Read in the pathways  stored in different sheets of an excel file and store as list (sheet=list elem.)


sp.wb <- loadWorkbook(file="D://chlamyJGI//Pathways//pathw.xlsx")
sp.l = readWorksheet(sp.wb, sheet = getSheets(sp.wb))
summary(sp.l)
names(sp.l)
sp.l[[1]]


###########################################################################
## prepare 1 ct=ontingency table for Fisher test for each pathway tested ##
## store as R list (as many cont. tables as pathways tested)             ##
###########################################################################

### N is the vertex number of the starting network

N=9178

u <- ids
names(u)
summary(u)
length(u[[1]])


class(u)

m <- numeric()
m
m <- unlist(lapply(u,function(i) length(i)))
str(m)

m <-numeric() 
for (i in names(u)){ 
  m <- c(m,length(u[[i]]))
}
m

ct.l <- list()# do a ct for each protein(row) and store in a list
for (p in names(sp.l)){
  
#   p=names(sp.l)[1]
#   p
  n <- nrow(sp.l[[p]])
#   n
  
  k.t <- numeric()
  for (i in names(u)){
    k.t <- c(k.t,length(intersect(u[[i]], sp.l[[p]]$GeneID)))
  }
#   k.t
  ct <- as.data.frame(matrix(nrow=length(u),ncol=4))#adjust "nrow" to fit the number of elements in list "u"
  ct[,1] <-k.t
  ct[,2] <-m
  ct[,3] <-n
  ct[,4] <-N
  colnames(ct) <- c('k','m','n','N')
  row.names(ct)<- names(u)
  ct<- as.matrix(apply(ct,2, as.numeric))
  k<-ct[,'k']
  nk<- apply(ct,1, function(x)  x[3]-x[1])
  mk<- apply(ct,1, function(x)  x[2]-x[1])
  nmk<-apply(ct,1, function(x)  x[4]-x[2]-x[3]+x[1])
  ct.fisher <-cbind(k,nk,mk,nmk)
  ct.l[[p]] <- ct.fisher
}


summary(ct.l)#the list with 1 contingency table / pathway tested
ct.l[[3]]

###############################################################
## Fisher test for each ct=cont. table in the pe.l list      ##
## Store the results in one final table:                     ##
##           rows=mutant names                               ##
##           cols=p values (Fisher test) for each pathway    ##
###############################################################

pval.l <- list()
for (i in names(ct.l)){
  fisher.pval.pathway <- apply(ct.l[[i]],1, function(x) fisher.test(matrix(x,nr=2),alternative="greater")$p.value)
  path.enrich <-as.data.frame(fisher.pval.pathway)
  row.names(path.enrich) <- names(u)
  pval.l[[i]] <- path.enrich
}

summary(pval.l)

pval.all<-do.call('cbind',pval.l)
colnames(pval.all) <- names(ct.l)
pval.all


