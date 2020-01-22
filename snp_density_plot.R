#!/usr/bin/Rscript

############################################################
#
# snp_density_plot.R
#
# Draw SNP density to visualize 1000 genome SNP distribution
#
# Author : Thomas Sunghoon Heo
#
############################################################
library(ggplot2)

gene.range <- read.table("./data/gene_range", sep="\t", header=T,stringsAsFactors=F)

# This will be inserted into for loop

for(gene in gene.range$Gene){
	target.range <- gene.range[gene.range$Gene == gene,]
	target.chrom <- target.range$chrom
	target.start <- target.range$Start
	target.end <- target.range$End
	target.strand <- target.range$Strand

	snps <- read.table(paste0("./data/", gene, ".snps"), sep="\t", header=T,stringsAsFactors=F)
	n <- ncol(snps)
	snps$VAF1 <- rowSums(snps[,7:n] >= 1)
	snps$VAF5 <- rowSums(snps[,7:n] >= 5)
	snps$VAF10 <- rowSums(snps[,7:n] >= 10)

	# Build a dataframe

	make.df <- data.frame(matrix(NA, nrow=nrow(snps), ncol=7))
	colnames(make.df) <- c("pos", "rsid", "is_vaf1", "is_vaf5", "is_vaf10", "is_exon","height")
	make.df$pos <- snps$pos
	make.df$rsid <- snps$rsid
	make.df$is_vaf1 <- snps$VAF1 >= 1
	make.df$is_vaf5 <- snps$VAF5 >= 1
	make.df$is_vaf10 <- snps$VAF10 >= 1
	make.df$is_exon <- snps$isExonic
	make.df$height <- 0.3
	make.df[make.df$is_exon=="NO",]$height <- 0.15

	all.X <- make.df$pos
	all.Y <- rep(3,times=length(all.X))
	all.height <- make.df$height


	vaf1.X <- make.df[make.df$is_vaf1,]$pos
	vaf1.Y <- rep(4,times=length(vaf1.X))
	vaf1.height <- make.df[make.df$is_vaf1,]$height

	vaf5.X <- make.df[make.df$is_vaf5,]$pos
	vaf5.Y <- rep(5,times=length(vaf5.X))
	vaf5.height <- make.df[make.df$is_vaf5,]$height

	vaf10.X <- make.df[make.df$is_vaf10,]$pos
	vaf10.Y <- rep(6,times=length(vaf10.X))
	vaf10.height <- make.df[make.df$is_vaf10,]$height

	tiff(paste0(gene, ".tiff"))
	plot(NULL,xlim=c(target.start - 1000, target.end + 1000), ylim = c(0,7),xlab="Genomic coordinate", ylab="", yaxt="n", main=gene,axes=F)
	segments(all.X, all.Y-all.height,all.X, all.Y + all.height)
	segments(vaf1.X, vaf1.Y-vaf1.height,vaf1.X, vaf1.Y + vaf1.height)
	segments(vaf5.X, vaf5.Y-vaf5.height,vaf5.X, vaf5.Y + vaf5.height)
	segments(vaf10.X, vaf10.Y-vaf10.height,vaf10.X, vaf10.Y + vaf10.height)

	exon.info <- read.table(paste0("~/synbioserver/cnv_snp_analysis/data/", gene, ".liftover.converted.txt"), sep="\t", stringsAsFactors=F, header=F)
	exon.info <- exon.info[,1:3]
	colnames(exon.info) <- c("chrom","start","end")

	ybottom <- 0.5
	ytop <- 1.5

	for(i in 1:nrow(exon.info)){
		X <- exon.info[i,2:3]
		rect(X$start, ybottom, X$end, ytop, col="black")
	}

	#segments(target.start, 1, target.end, 1)

	arrow.starts <- seq(target.start, target.end - 1000, by= 1000)
	arrow.ends <- seq(target.start + 1000, target.end, by= 1000)
	arrow.ystarts <- rep(1, length(arrow.starts))
	arrow.yends <- rep(1, length(arrow.ends))

	if(target.strand == 1){
		arrows(arrow.starts,arrow.ystarts, arrow.ends, arrow.yends,length=0.1)
	} else {
		arrows(arrow.ends,arrow.ystarts, arrow.starts, arrow.yends,length=0.1)
	}

	axis(side=1,
	at=seq(target.start-1000,target.end+1000,by=1000),
	labels=seq(target.start-1000,target.end+1000,by=1000)
	)


	axis(2,
	at=c(1,3,4,5,6),
	labels=FALSE,
	tick=FALSE
	)
	text(par("usr")[1], c(1,3,4,5,6), labels=c("Genome","AllSnps", "VAF1", "VAF5", "VAF10"), srt=0, pos=2, tck=-0.02, xpd=TRUE)
	
	dev.off()
}
