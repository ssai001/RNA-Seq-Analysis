#Samples and environment settings
library(systemPipeRdata)
genWorkenvir(workflow="rnaseq")
setwd("rnaseq")
library(systemPipeR)
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR") #arranging all FASTQ files into "targets.txt"
targets <- read.delim(targetspath, comment.char = "#")[,1:4] #defines all FASTQ files and sample comparisons of the analysis workflow
targets

#Read preprocessing - Read quality filtering and trimming
args <- systemArgs(sysma="param/trim.param", mytargets="targets.txt")
preprocessReads(args=args, Fct="trimLRPatterns(Rpattern='GCCCGGGTAA', subject=fq)",
                batchsize=100000, overwrite=TRUE, compress=TRUE)
writeTargetsout(x=args, file="targets_trim.txt", overwrite=TRUE)

# FASTQ Quality Report
args <- systemArgs(sysma="param/tophat.param", mytargets="targets.txt")
fqlist <- seeFastq(fastq=infile1(args), batchsize=100000, klength=8)
pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
seeFastqPlot(fqlist)
dev.off()

#Read mapping with Bowtie2/Tophat2
args <- systemArgs(sysma="param/tophat.param", mytargets="targets.txt")
sysargs(args)[1] # Command-line parameters for first FASTQ file
system("bowtie2-build ./data/tair10.fasta ./data/tair10.fasta") #This is the reference genome that the target files are aligned to

#Read and alignment stats
read_statsDF <- alignStats(args=args) 
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
read.table(system.file("extdata", "alignStats.xls", package="systemPipeR"), header=TRUE)[1:4,]

#Read counting with summarizeOverlaps in parallel mode using multiple cores
library("GenomicFeatures"); library(BiocParallel)
txdb <- makeTxDbFromGFF(file="data/tair10.gff", format="gff", dataSource="TAIR", organism="Arabidopsis thaliana")
saveDb(txdb, file="./data/tair10.sqlite")
txdb <- loadDb("./data/tair10.sqlite")
(align <- readGAlignments(outpaths(args)[1])) # Demonstrates how to read bam file into R
eByg <- exonsBy(txdb, by=c("gene"))
bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
counteByg <- summarizeOverlaps(eByg, bfl, mode="Union", 
                               ignore.strand=TRUE, 
                               # preprocess.reads=invertStrand,
                               inter.feature=FALSE, 
                               singleEnd=TRUE)
countDFeByg <- assays(counteByg)$counts
rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=eByg))
write.table(countDFeByg, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")

# Sample-wise correlation analysis












