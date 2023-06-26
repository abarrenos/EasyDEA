

#' create.output.hierarchy
#'
#' Create the hierarchy of output directories that we need for the analysis
#' 
#' @param	rnaseq.out	the folder where we want to store all our output
#' @param	use.both.reads	whether both reads of a paired-reads sequencing
#'			experiment should be required to map in the analysis
#' 
#' @return	the name of the output folder where results will be stored
#'
#' @usage	out.dir <- create.output.hierarchy(rnaseq.out, use.both.reads)
#' 
#' @examples	out.dir <- create.output.hierarchy('.', TRUE)
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
create.output.hierarchy <- function(rnaseq.out='.', use.both.reads=TRUE)
{
# set the name of the output folder and log file
    if (use.both.reads == T) {
        requireBothEnds <- T
        folder <- paste(rnaseq.out, "both_ends", sep='/')	# whether both ends sshould match or just any end
    } else {
        requireBothEnds <- F
        folder <- paste(rnaseq.out, "any_end", sep='/')
    }

    # create needed directory hierarchy
    # this could go into a separate function for simplicity
    dir.create(rnaseq.out, showWarnings=FALSE)
    dir.create(folder, showWarnings=FALSE)
    dir.create(file.path(folder, "log"), showWarning=FALSE)
    dir.create(file.path(folder, "img"), showWarning=FALSE)
    dir.create(file.path(folder, "annotation"), showWarning=FALSE)
    dir.create(file.path(folder, "edgeR"), showWarning=FALSE)
    dir.create(file.path(folder, "edgeR/img"), showWarning=FALSE)
    dir.create(file.path(folder, "edgeR/go"), showWarning=FALSE)
    dir.create(file.path(folder, "edgeR/pfam"), showWarning=FALSE)
    dir.create(file.path(folder, "edgeR/cluster"), showWarning=FALSE)
    dir.create(file.path(folder, "edgeR/cluster/kmeans"), showWarning=FALSE)
    dir.create(file.path(folder, "edgeR/cluster/pam"), showWarning=FALSE)
    dir.create(file.path(folder, "edgeR/cluster/dbscan"), showWarning=FALSE)
    dir.create(file.path(folder, "DESeq2/raw"), showWarning=FALSE)
    dir.create(file.path(folder, "DESeq2/shrunk"), showWarning=FALSE)
    dir.create(file.path(folder, "DESeq2/signif"), showWarning=FALSE)
    dir.create(file.path(folder, "DESeq2"), showWarning=FALSE)
    dir.create(file.path(folder, "DESeq2/go"), showWarning=FALSE)
    dir.create(file.path(folder, "DESeq2/pfam"), showWarning=FALSE)
    dir.create(file.path(folder, "DESeq2/img"), showWarning=FALSE)
    dir.create(file.path(folder, "DESeq2/cluster"), showWarning=FALSE)
    dir.create(file.path(folder, "DESeq2/cluster/kmeans"), showWarning=FALSE)
    dir.create(file.path(folder, "DESeq2/cluster/pam"), showWarning=FALSE)
    dir.create(file.path(folder, "DESeq2/cluster/dbscan"), showWarning=FALSE)

    return(folder)
}

#' align.fastq
#'
#' Align all fastq files inside a folder against a reference genome
#' 
#' @param path	the path to the folder where the fastq files are stored
#'
#' @param	reference	the name of the reference genome (without extension)
#'
#' @param 	aln.out	the path to a directory where we want to store the alignments
#'
#' @param	save.dir	the path to a directory where we will store summaries
#' 
#' @return	nothing
#'
#' @usage	align.fastq(path, reference, aln.out, save.dir=aln.out)
#' 
#' @examples	align.fastq('./fastq', 'Hsapiens', 'alignments')
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
align.fastq <- function(path, reference, aln.out, save.dir=aln.out) {
    fastq.data <- path
    
    # get the fastq file names
    R1.fastq.files <- list.files(path=fastq.data, pattern='R1', full.names=TRUE)
    R2.fastq.files <- list.files(path=fastq.data, pattern='R2', full.names=TRUE)
    print(R1.fastq.files)
    print(R2.fastq.files)
    
    dbname <- reference

    # we expect a fasta reference file inside a directory named 'ref' in
    # the alignment output directory, named 'reference'.{fna|fa|fas|fasta}
    curWD <- getwd()	# save current location
    setwd(aln.out)	# go to where we want to have the alignments
    reference.fasta <- paste("./ref/", reference, ".fna", sep='')
    if (! file.exists(reference.fasta))
        reference.fasta <- paste("./ref/", reference, ".fa", sep='')
    if (! file.exists(reference.fasta))
        reference.fasta <- paste("./ref/", reference, ".fas", sep='')
    if (! file.exists(reference.fasta))
        reference.fasta <- paste("./ref/", reference, ".fasta", sep='')
    # if all else failed we stop here
    if (! file.exists(reference.fasta))
        cat.err("No suitable reference fasta file found in\n",
                "	", getwd(), "/ref/",
                "(I tried with ", reference, ".fna .fa .fas .fasta)\n", 
                sep='')
        
    # we'll check if the output files exist to avoid repeating
    # work already done
    if (! file.exists(paste(reference, '.0.b.tab', sep=''))) {
        cat.info('\nBUILDING INDEX\n')
        # build the reference index inside the 'aln.out' directory
        buildindex(basename=dbname,reference=reference.fasta)
        dir()
    }
    
    r11 <- basename(R1.fastq.files[1])
    if (! file.exists(paste(aln.out, '/', r11, '.subread.BAM', sep=''))) {
        cat.info('\nALIGNING\n')
        # align the reads
        # IMPORTANT NOTE: WE NEED TWO FILES LISTING ALL THE FASTQ FILES TO ALIGN
        #	R1.fastq.files and R2.fastq.files
        align(index=dbname, readfile1=R1.fastq.files, readfile2=R2.fastq.files)
        
        # align will generate the output in the data directory, we
        # do not want to mix datasets, so we move the alignment results
        # to the corresponding output directory
        system(paste("mv ", fastq.data, "/*.BAM .", sep=""))
        system(paste("mv ", fastq.data, "/*.vcf .", sep=""))
        system(paste("mv ", fastq.data, "/*.summary .", sep=""))
    }
    setwd(curWD)	# go back to where we were
    
    # get the bam file names and inspect them to count the number
    # of reads that map to each genome position
    if ( ! file.exists(paste(save,dir, 'bam_files_stats.txt', sep='/'))) {
        bam.files <- list.files(path=aln.out, pattern='.BAM$', full.names = TRUE)
        print(bam.files)
        props <- propmapped(files=bam.files)
        print(props)
        write.table(props, file=paste(save,dir, 'bam_files_stats.txt', sep='/'), 
    	    row.names=T, col.names=T, sep='\t')
    }

}


#' compute.feature.counts
#'
#' Takes a list of aligned bam files and a reference and computes the times
#' each gene (feature) in the reference is matched by an aligned read. It
#' will try to use a file named 'reference'.gtf or, if one does not exist,
#' a file named 'reference'.gff as a source of annotation with the feature
#' (gene) coordinates. As such, the reference gtf/gff file must exist and
#' match the reference fasta sequence used in the alignment.
#' 
#' @param	bam.files	a list of bam files with reads mapped to the reference
#' @param	reference	the base name of the reference used for aligning and
#'				to be used to search for an annotation source
#' @param	requireBothEnds	whether both ends of a read pair are required to 
#'				match in order to be counted
#' @param	save.dir	the name of the directory where the feature
#'				counts will be saved for future reference
#' 
#' @return	the feature count table
#'
#' @usage	fc <- compute.feature.counts(
#'			bam.files=list.files(path = aln.dir, pattern = '.BAM$', full.names = TRUE),
#'			reference=ref.genome,
#'			requireBothEnds=T,
#'			save.dir='.')
#' 
#' @examples
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
compute.feature.counts <- function(bam.files, reference, requireBothEnds=T, save.dir='.') {

    reference.ann   <- paste("./ref/", reference, ".gtf", sep='')
    is.gtf <- TRUE
    if (! file.exists(reference.ann) ) {
        reference.ann   <- paste("./ref/", reference, ".gff", sep='')
        is.gtf=FALSE
    }
    if (! file.exists(reference.ann) ) {
        cat.err("No annotation file ./ref/.", reference, ".gtf or .gff found\n",
                sep='')
    }
    fc <- featureCounts(bam.files=bam.files, 
	    annot.ext=reference.ann, 
            isGTFAnnotationFile=is.gtf,
	    isPairedEnd=T, 
            requireBothEndsMapped=requireBothEnds, 
            primaryOnly=T, 
            ignoreDup=T, 
            useMetaFeatures=T)
    #  this means: process all BAM files
    #	Use as annotation the external file reference.gtf which is GTF
    #	BAM files contain paired end reads and we will only consider those
    #	where both ends match against the genome
    #    -------------------------------------------------------------
    #          ----->        <-----
    #	We will filter matches so that if a read may match more than one
    #	place in the genome, we will only consider the primary match and
    #	ignore other, duplicate matches
    #	We use meta-features
    #	We match against any feature in the GTF file, not only genes
    #	(this implies we will need to check the annotation carefully later)
    #	We do not remove chimeric fragments, but do actually count them too

    # SAVE FEATURE COUNTS
    # -------------------
    # now the variable fc contains 4 columns: annotation, target, counts and stats
    # but they exist only in the RAM memory, they are not stored somewhere safe,
    # so, we save them and create new variables so as to make it easier for us to 
    # manipulate the data
    write.csv(fc$counts, file=paste(save.dir, 'featureCounts.csv', sep='/'),
	      row.names=T, col.names=T, quote=F)
    write_delim(fc$stat, file=paste(save.dir, 'featureCounts_stat.txt', sep='/'), 
	        delim='      ')
    # save all the contents of 'fc' in an RDS file
    saveRDS(fc, file=paste(save.dir, '/featureCounts.rds', sep=''))
    #	'fc' can later be recovered with: fc <- readRDS(file='featureCounts.rds'
    # and save as well as Rdata file
    save(fc, file=paste(save.dir, '/featureCounts.RData', sep=''))
    #	'fc' can later be recovered with: fc <- load(file='featureCounts.RData')

    return(fc)
}



#' h.cluster.changes
#'
#' Cluster gene expression data using hierarchycal clustering
#' 
#' @param
#' 
#' @return
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
h.cluster.changes <- function(data.table, annotated.data, 
                              distance="euclidean",	# see ?dist()
                              clusters=1,
                              method="complete",	# see ?hclust()
                              estimate=TRUE,
                              interactive=TRUE
                              ) {
                              

    # Normalize the data
    by.row <- 1
    by.col <- 2
    dif <- data.table
    means <- apply(dif, by.col, mean)
    sds <- apply(dif, by.col, sd)		# standard deviations
    nor <- scale(dif,center=means,scale=sds)	# normalized values

    # Calculate distance matrix  
    distan = dist(nor, method=distance)

    # Hierarchical agglomerative clustering  
    cat.info("
    H I E R A R C H I C A L   C L U S T E R I N G
    =============================================
    \n\n")
    
    hc = hclust(distan)
#    if (verbose) {
#        as.png(plot(hc), 'hclust.png')
#        as.png(plot(hc,hang=-1), 'hclust.hang.png')
#        as.png(plot(hc,labels=rownames(data.table),main='Default from hclust')
#               "hclust.labelled.png")
#    }
    if (interactive) {
        print(plot(hc,labels=rownames(data.table),main='Default from hclust'))
        continue.on.enter("Press [ENTER] to continue ")
    } 
       
    # Cluster membership
    if ((clusters <= 1) & (estimate == T)) {
        ans <- continue.on.enter(prompt="How many clusters should we use? ")
        nclust <- as.numeric(ans)
    } else
        nclust <- clusters
    
    member <- cutree(hc, nclust)	# cut to 4 groups
    cat.info("Cluster membership counts\n")
    print(table(member))
    # Characterizing clusters 
    cat('\n')
    cat.info("Means by cluster in the normalized data\n")
    print(aggregate(nor, list(member), mean))
    cat('\n')
    cat.info("Means by cluster in the non-normalized data\n")
    print(aggregate(dif, list(member), mean))

    for (i in 1:nclust) {
        clus.i <- member[ member == i ]
        data.i <- annotated.data[annotated.data$ensembl.gene.id %in% names(clus.i), ]
        # sort by l2FC
        data.i <- data.i[order(data.i$log2FoldChange, decreasing=T), ]
        cat.info("Showing annotation for cluster", i, "(", length(clus.i), " elements)\n")
        
        #View(data.i[ , c("log2FoldChange", "GENENAME", "entrezgene_description", "definition_1006")],
        #     title=paste("Cluster no.", i, "(", length(clus.i), ") elements)"))
        
        # or, using tcltk and gWidget2
        #library(tcltk)
        #library(gWidgets2)
        data.to.show <- data.i[ , c("log2FoldChange", "GENENAME", "entrezgene_description", "definition_1006")]
        # clean up for showing
        data.to.show[ is.na(data.to.show) ] <- 'NA'
        window.name <- paste("Cluster no.", i, "(", length(clus.i), " elements)")
        #window <- gwindow(title=window.name, visible=TRUE)
        #tab <- gtable(data.to.show,
        #       container=window)
        window <- show.data.frame(data.to.show, window.name, visible=F)
        if (interactive) {
          visible(window) <- TRUE	# setting it to FALSE removes window
          keypress()
          visible(window) <- FALSE
        }
    }

    cat.info("
    
    Silhouette Plot for hierarchical clustering with normalized data
    ----------------------------------------------------------------
    Measure similarity of each object to its own cluster (cohesion) compared to
    other clusters (dispersion). Values range from -1 to +1. Large values
    indicate objects well matched to their own cluster and badly to neighboring
    clusters. If many points have low or negative value, the number of clusters
    is too low or too high.
    \n")
    #if (verbose)
    #    as.png(plot(silhouette(cutree(hc, nclust), distan)), "silhouette.png")
    print(plot(silhouette(cutree(hc, nclust), distan)))

    return(hc)
}



#' hcut.cluster changes
#'
#' Cluster gene expression data using hierarchycal clustering with 'hcut'
#' 
#' @param
#' 
#' @return
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
hcut.cluster.changes <- function(data.table, annotated.data,
                                    algorithm="hclust",	# see ?hcut()
                                    distance="euclidean", # see ?hcut()
                                    method="ward.D2",	# see ?hcut()
                                    clusters=1,	# n. of clusters
                                    nstart=1,	# ignored
                                    gap_bootstrap=500,
                                    estimate=FALSE
                                    ) {

    # Normalize the data
    by.row <- 1
    by.col <- 2
    dif <- data.table
    means <- apply(dif, by.col, mean)
    sds <- apply(dif, by.col, sd)
    nor <- scale(dif,center=means,scale=sds)

    cat.info("
    H c u t
    -------
    \n")

    if (estimate) {
        cat("Suggestions for the best number of clusters\n")

        cat("Plot by within-cluster sums of squares\n")
        cat("    Elbow method: look for a knee\n")
        print(fviz_nbclust(nor, hcut, method="wss"))
        continue.on.enter("Press [ENTER] to continue ")

        cat("
        Average Silhouette Method

        The average silhouette approach measures the quality of a clustering. It
        determines how well each observation lies within its cluster.

        A high average silhouette width indicates a good clustering. The average
        silhouette method computes the average silhouette of observations for
        different values of k.
        \n")
        print(fviz_nbclust(nor, hcut, method="silhouette"))
        continue.on.enter("Press [ENTER] to continue ")

        cat("
        Gap Statistic Method

        This approach can be utilized in any type of clustering method (i.e.
        K-means clustering, hierarchical clustering).

        The gap statistic compares the total intracluster variation for different
        values of k with their expected values under null reference distribution of
        the data.

        \n")
        gap_stat <- clusGap(nor, FUN=hcut, nstart=25,
                            K.max=15, B=gap_bootstrap)
        print(fviz_gap_stat(gap_stat))
        continue.on.enter("Press [ENTER] to continue ")
    }

    # Number of clusters
    if (clusters == 1) {
        ans <- continue.on.enter(prompt="How many clusters should we use? ")
        nclust <- as.numeric(ans)
    } else
        nclust <- clusters
    
    # clustering (N groups)
    set.seed(123)
    hc <- hcut(nor, k=nclust, 
               hc_func=algorithm, hc_method=method, hc_metric=distance, is_diss=FALSE)
    print(head(hc))
    member <- hc$cluster
    
    cat.info("Cluster membership counts\n")
    print(table(member))
    # Characterizing clusters 
    cat('\n')
    cat.info("Means by cluster in the normalized data\n")
    print(aggregate(nor, list(member), mean))
    cat('\n')
    cat.info("\nMeans by cluster in the non-normalized data\n")
    print(aggregate(dif, list(member), mean))

    for (i in 1:nclust) {
        clus.i <- member[ member == i ]
        data.i <- annotated.data[annotated.data$ensembl.gene.id %in% names(clus.i), ]
        # sort by l2FC
        data.i <- data.i[order(data.i$log2FoldChange, decreasing=T), ]
        cat("Showing annotation for cluster", i, "(", length(clus.i), ") elements)\n")

        # invoke a spreadsheet-style data viewer
        View(data.i[ , c("log2FoldChange", "GENENAME", "entrezgene_description", "definition_1006")],
             title=paste("Cluster no.", i, "(", length(clus.i), ") elements)"))
 #       keypress()
    }
    
    # plot a PCA of the clusters
    # If there are more than two dimensions (variables) fviz_cluster will perform
    # principal component analysis (PCA) and plot the data points according to
    # the first two principal components that explain the majority of the
    # variance.
    print(fviz_cluster(hc, data=nor, geom="point"))
    continue.on.enter("Press [ENTER] to continue ")
    # to plot gene names instead of ensembl.ids we need to change rownames to genes
    gor <- nor 
    rownames(gor) <- paste(annotated.data[ , "GENENAME"], 1:length(rownames(gor)))
    print(fviz_cluster(hc, data=gor))
    continue.on.enter("Press [ENTER] to continue ")
}



#' k.means.cluster changes
#'
#' Cluster gene expression data using K-means clustering
#' 
#' @param
#' 
#' @return
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
k.means.cluster.changes <- function(data.table, annotated.data,
                                    algorithm="Hartigan-Wong",	# see kmeans()
                                    clusters=1,	# n. of clusters
                                    nstart=1,	# n. of random start sets to choose
                                    estimate=F,
                                    gap_bootstrap=500
                                    ) {

    # Normalize the data
    by.row <- 1
    by.col <- 2
    dif <- data.table
    means <- apply(dif, by.col, mean)
    sds <- apply(dif, by.col, sd)
    nor <- scale(dif,center=means,scale=sds)

    cat.info("
    K - m e a n s
    -------------
    \n")
    
    if (estimate) {
        cat("Suggestions for the best number of clusters\n")

        #
        # Scree Plot
        cat("

        Scree plot using normalized data

        This allows us to evaluate how much variation we account for as we
        consider more clusters and decide what a reasonable number of clusters
        might be.
        We draw here the variances accounted for using up to 20 K-means clusters
        \n") 
        # compute variances by row
        wss <- (nrow(nor)-1)*sum(apply(nor, by.row, var))
        for (i in 2:20) wss[i] <- sum(kmeans(nor, centers=i)$withinss)
        print(plot(1:20, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares") )
        continue.on.enter("Press [ENTER] to continue ")

        # The scree plot will allow us to see the variabilities in clusters, 
        # we expect that if we increase the number of clusters, then the 
        # within-group sum of squares would come down. 
        #

        cat("Plot by within-cluster sums of squares\n")
        cat("    Elbow method: look for a knee\n")
        print(fviz_nbclust(nor, kmeans, method="wss"))
        continue.on.enter("Press [ENTER] to continue ")

        cat("
        Average Silhouette Method

        The average silhouette approach measures the quality of a clustering. It
        determines how well each observation lies within its cluster.

        A high average silhouette width indicates a good clustering. The average
        silhouette method computes the average silhouette of observations for
        different values of k.
        \n")
        print(fviz_nbclust(nor, kmeans, method="silhouette"))
        continue.on.enter("Press [ENTER] to continue ")

        cat("
        Gap Statistic Method

        This approach can be utilized in any type of clustering method (i.e.
        K-means clustering, hierarchical clustering).

        The gap statistic compares the total intracluster variation for different
        values of k with their expected values under null reference distribution of
        the data.

        \n")
        gap_stat <- clusGap(nor, FUN=kmeans, nstart=25,
                            K.max=15, B=gap_bootstrap)
        print(fviz_gap_stat(gap_stat))
        continue.on.enter("Press [ENTER] to continue ")
    }

    # Number of clusters
    if (clusters == 1) {
        ans <- continue.on.enter(prompt="How many clusters should we use? ")
        nclust <- as.numeric(ans)
    } else
        nclust <- clusters
    
    # K-means clustering (N groups)
    kc <- kmeans(nor, centers=nclust, nstart=nstart, algorithm=algorithm)
    print(head(kc))
    member <- kc$cluster
    
    cat("Cluster membership counts\n")
    print(table(member))
    # Characterizing clusters 
    cat("\nMeans by cluster in the normalized data\n")
    print(aggregate(nor, list(member), mean))
    cat("\nMeans by cluster in the non-normalized data\n")
    print(aggregate(dif, list(member), mean))

    for (i in (1:nclust)) {
        clus.i <- member[ member == i ]
        data.i <- annotated.data[annotated.data$ensembl.gene.id %in% names(clus.i), ]
        # sort by l2FC
        data.i <- data.i[order(data.i$log2FoldChange, decreasing=T), ]
        cat("Showing annotation for cluster", i, "(", length(clus.i), ") elements)\n")
        View(data.i[ , c("log2FoldChange", "GENENAME", "entrezgene_description", "definition_1006")],
             title=paste("Cluster no.", i, "(", length(clus.i), ") elements)"))
 #       keypress()
    }
    
    # plot a PCA of the clusters
    # If there are more than two dimensions (variables) fviz_cluster will perform
    # principal component analysis (PCA) and plot the data points according to
    # the first two principal components that explain the majority of the
    # variance.
    print(fviz_cluster(kc, data=nor, geom="point"))
    continue.on.enter("Press [ENTER] to continue ")
    # to plot gene names instead of ensembl.ids we need to change rownames to genes
    gor <- nor 
    rownames(gor) <- paste(annotated.data[ , "GENENAME"], 1:length(rownames(gor)))
    print(fviz_cluster(kc, data=gor))
    continue.on.enter("Press [ENTER] to continue ")
}



#'
#' pam.cluster.changes
#'
#'  Cluster gene expression data using Partition around medoids clustering
#' @param
#' 
#' @return
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
pam.cluster.changes <- function(data.table, annotated.data,
                                    distance="euclidean", 	# see dist()
                                    metric="euclidean",		# see pam()
                                    clusters=1,		# n. of clusters
                                    nstart=1,	# n. of random start sets to choose
                                    estimate=T,
                                    gap_bootstrap=500
                                    ) {
                                    
    if (nstart != 1) medoids="random" else medoids=NULL
        
    # Normalize the data
    by.row <- 1
    by.col <- 2
    dif <- data.table
    means <- apply(dif, by.col, mean)
    sds <- apply(dif, by.col, sd)
    nor <- scale(dif,center=means,scale=sds)

    cat("
    Partition around medoids
    ------------------------
    \n")

    if (estimate) {
        cat("Suggestions for the best number of clusters\n")

        cat("    Plot by within-cluster sums of squares\n\n")
        cat("    Elbow method: look at the knee\n")
        print(fviz_nbclust(nor, pam, method="wss"))
        continue.on.enter("Press [ENTER] to continue ")

        cat("
        Average Silhouette Method

        The average silhouette approach measures the quality of a clustering. It
        determines how well each observation lies within its cluster.

        A high average silhouette width indicates a good clustering. The average
        silhouette method computes the average silhouette of observations for
        different values of k.
        \n")
        print(fviz_nbclust(nor, pam, method="silhouette"))
        continue.on.enter("Press [ENTER] to continue ")

        cat("
        Gap Statistic Method

        This approach can be utilized in any type of clustering method (i.e.
        K-means clustering, hierarchical clustering).

        The gap statistic compares the total intracluster variation for different
        values of k with their expected values under null reference distribution of
        the data.

        \n")
        gap_stat <- clusGap(nor, FUN=pam, nstart=25,
                            K.max=15, B=gap_bootstrap)
        print(fviz_gap_stat(gap_stat))
        continue.on.enter("Press [ENTER] to continue ")
    }


    # Number of clusters
    if (clusters == 1) {
        ans <- continue.on.enter(prompt="How many clusters should we use? ")
        nclust <- as.numeric(ans)
    } else
        nclust <- clusters

    # cluster with partition around medioids for k=N clusters 
    # (data is not dissimilarity but distance)
    # compute distance (we have dissimilarity to a common reference but not
    # between the considered groups).
    cat("Computing clusters (may take some time)...\n")
    eucldist <- dist(nor, method=distance) 
    #cluster
    ### NOTE ### NOTE: may be worth trying to cluster separately with diss=TRUE)
    pam.clus <- pam(eucldist, k=nclust, diss=TRUE)
    print(head(pam.clus))
    member <- pam.clus$clustering

    cat("Cluster membership information\n")
    print(pam.clus$clusinfo)
    cat("\nMeans by cluster in the normalized data\n")
    print(aggregate(nor, list(member), mean))
    cat("\nMeans by cluster in the non-normalized data\n")
    print(aggregate(dif, list(member), mean))
    
    for (i in (1:nclust)) {
        clus.i <- member[ member == i ]
        data.i <- annotated.data[annotated.data$ensembl.gene.id %in% names(clus.i), ]
        # sort by l2FC
        data.i <- data.i[order(data.i$log2FoldChange, decreasing=T), ]
        cat("Showing annotation for cluster", i, "(", length(clus.i), ") elements)\n")
        View(data.i[ , c("log2FoldChange", "GENENAME", "entrezgene_description", "definition_1006")],
             title=paste("Cluster no.", i, "(", length(clus.i), ") elements)"))
 #       keypress()
    }
    
    # plot the partitioning
    print(clusplot(pam.clus, shade = FALSE,labels=F,
	    col.clus="blue",col.p="red",
            span=FALSE,
            main="PAM Cluster Mapping",cex=1.2))
    continue.on.enter("Press [ENTER] to continue ")

    # plot a PCA of the clusters
    # If there are more than two dimensions (variables) fviz_cluster will perform
    # principal component analysis (PCA) and plot the data points according to
    # the first two principal components that explain the majority of the
    # variance.
    print(fviz_cluster(pam.clus, data=nor, geom="point"))
    continue.on.enter("Press [ENTER] to continue ")
    # to plot gene names instead of ensembl.ids we need to change rownames to genes
    gor <- nor 
    rownames(gor) <- paste(annotated.data[ , "GENENAME"], 1:length(rownames(gor)))
    # this seemingly ignores the data argument!!!
    print(fviz_cluster(pam.clus, data=gor))
    continue.on.enter("Press [ENTER] to continue ")
}


#' cluster.changes
#'
#' Cluster gene expression data using any of a variety of methods for
#' clustering
#' 
#' @param
#' 
#' @return
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
cluster.changes <- function(data.table, annotated.data,
                            FUN=hcut,
                            clusters=1,			# n. of clusters
                            algorithm="default",	# see below	
                            distance="default", 	# see below
                            method="default",		# see see below
                            nstart=1,		# n. of ran dom start sets to choose
                            eps=1.0,		# epsilon for DBScan
                            gap_bootstrap=100,
                            normalize=TRUE,
                            estimate=TRUE,	# only if clusters > 1
                            output.folder=NULL	# NULL => no output desired
                                    ) {
# FUN -- one of c(hcut, kmeans, pam)
# algorithm -- for 'hcut' one of c(_"hclust"_, "agnes", "diana")
#              for 'kmeans' one of c(_"Hartigan-Wong"_, "Lloyd", "Forgy", "MacQueen")
#              for 'pam' one of c(_"original"_, "o_1", "o_2", "f_3", "f_4", "f_5", 
#                       "faster")
# distance -- for 'hcut' one of c(_"euclidean"_, "manhattan", "maximum", 
#                       "canberra", "binary", "minkowski", "pearson", "spearman", 
#                       "kendall")
#             for 'kmeans' it is ignored
#             for 'pam' one of c(_"euclidean"_, manhattan")
#             for 'dbscan' it is ignored
# method -- for 'hcut' one of c("ward.D"', "ward.D2", "single", "complete", 
#                       "average" (= UPGMA), "mcquitty" (= WPGMA), 
#                       "median" (= WPGMC), "centroid" (= UPGMC), 
#                       "weighted" (=WPGMA), "flexible", "gaverage")
#             for 'kmeans' ignored
#             for 'pam' ignored
#             for 'dbscan' ignored

    # Normalize the data
    by.row <- 1
    by.col <- 2
    if (normalize == TRUE) {
        means <- apply(data.table, by.col, mean)
        sds <- apply(data.table, by.col, sd)
        nor <- scale(data.table,center=means,scale=sds)
    } else {
        nor <- data.table
    }

    fun.name <- deparse(substitute(FUN))
    cat("
    C L U S T E R I N G    W I T H :   ", fun.name, "
    ---------------------------------------------
    \n")

    # re-activate next line to save computation time
    #if (clusters > 1) estimate <- FALSE		# we already know the number
    if (clusters < 1) estimate <- TRUE
    # we may have clusters == 1 and estimate == FALSE, e.g. to find outliers
    
    if (estimate == TRUE) {
        cat("Suggestions for the best number of clusters using", fun.name, "\n")

        if (normalize == TRUE) {
            cat("Plot by within-cluster sums of squares\n")
            cat("    Elbow method: look for a knee (normalized data)\n")
            print(fviz_nbclust(nor, FUN, method="wss"))
            out.png <- sprintf("%s/DESeq2_%s_wss.png",
                output.folder, fun.name)
            as.png(fviz_nbclust(nor, FUN, method="wss"), out.png)
            continue.on.enter("Press [ENTER] to continue ")
        }

        cat("Plot by within-cluster sums of squares\n")
        cat("    Elbow method: look for a knee (raw data)\n")
        print(fviz_nbclust(data.table, FUN, method="wss"))
        out.png <- sprintf("%s/DESeq2_%s_wss.png",
            output.folder, fun.name)
        as.png(fviz_nbclust(data.table, FUN, method="wss"), out.png)
        continue.on.enter("Press [ENTER] to continue ")

        if (fun.name != "dbscan") {
            cat("
            Average Silhouette Method

            The average silhouette approach measures the quality of a clustering. It
            determines how well each observation lies within its cluster.

            A high average silhouette width indicates a good clustering. The average
            silhouette method computes the average silhouette of observations for
            different values of k.
            \n")
            print(fviz_nbclust(nor, FUN, method="silhouette"))
            out.png <- sprintf("%s/DESeq2_%s_silhouette.png",
                output.folder, fun.name)
            as.png(fviz_nbclust(nor, FUN, method="silhouette"), out.png)

            continue.on.enter("Press [ENTER] to continue ")
        }

        cat("
        Gap Statistic Method

        This approach can be utilized in any type of clustering method (e.g.
        K-means, hierarchical clustering, partition around medoids).

        The gap statistic compares the total intracluster variation for different
        values of k with their expected values under null reference distribution of
        the data.

        We are looking for the highest peak identified.

        \n")
        gap_stat <- clusGap(nor, FUN=FUN, K.max=15, B=gap_bootstrap)
        print(fviz_gap_stat(gap_stat))
        out.png <- sprintf("%s/DESeq2_%s_gap_stat.png",
            output.folder, fun.name)
        as.png(fviz_gap_stat(gap_stat), out.png)
        continue.on.enter("Press [ENTER] to continue ")
    }
    
    # Select the number of clusters
    if ((clusters < 1) & (fun.name != 'dbscan')) {
        ans <- continue.on.enter(prompt="How many clusters should we use? ")
        nclust <- as.numeric(ans)
    } else {
        nclust <- clusters
    }
    
    # Do the clustering (N groups). This is necessarily method-specific
    if (fun.name == "hcut") {
        if (algorithm == "default") algorithm <- 'hclust'
        if (distance == "default") distance <- "euclidean"
        if (method == "default") method <- "complete"
        
        cl <- hcut(nor, k=nclust, 
               hc_func=algorithm, hc_method=method, hc_metric=distance, is_diss=FALSE)
        print(head(cl))
        member <- cl$cluster
    } else if (fun.name == "kmeans") {
        if (algorithm == "default") algorithm <- "Hartigan-Wong"
        if (distance == "default") distance <- "euclidean"
        
        cl <- kmeans(nor, centers=nclust, nstart=nstart, algorithm=algorithm)
        print(head(cl))
        member <- cl$cluster
    } else if (fun.name == "pam") {
        if (algorithm == "default") algorithm <- "faster"
        if (distance == "default") distance <- "euclidean"
        
        # cluster with partition around medioids for k=N clusters 
        # (data is not dissimilarity but distance)
        # compute distance (we have dissimilarity to a common reference but not
        # between the considered groups).
        cat("Computing clusters (may take some time)...\n")
        use_dist_matrix <- FALSE	### NOTE ### fviz_cluster fails, why?
        if (use_dist_matrix == TRUE) {
            ### NOTE ### It may be worth trying to cluster with diss=TRUE)
            # calculate dissimilarity matrix and cluster it
            #distm <- get_dist(nor, method=distance, stand=TRUE) 
            cl <- pam(get_dist(nor, method=distance, stand=TRUE),
                      k=nclust, diss=TRUE, 
                      nstart=nstart, metric=distance, variant=algorithm)
        } else {
            cl <- pam(nor, k=nclust, diss=FALSE, 
                      nstart=nstart, metric=distance, variant=algorithm)
        }

        print(head(cl))
        member <- cl$cluster
        cat("Summary\n")
        print(cl$clusinfo)
    }  else if (fun.name == "dbscan") {
        if (algorithm == "defaut") algorithm <- "hybrid"
        if (distance == "default") distance <- "manhattan"
        cl <- dbscan(nor, eps=eps, MinPts=4, showplot=1)
        # showplot=1 makes it produce a movie plot
        continue.on.enter("Press [ENTER] to continue ")
        print(head(cl))
        member <- cl$cluster
        names(member) <- annotated.data$ensembl.gene.id
        nclust <- max(member)
        plot(cl, nor, main="DBScan")
        out.png <- sprintf("%s/DESeq2_%s_plot.png",
            output.folder, fun.name)
        as.png(plot(cl, nor, main="DBScan"), out.png)
        continue.on.enter("Press [ENTER] to continue ")
    }
        
    # this may take very long and is of little use for now
    if (FALSE) {
        cat("Distance plot\n")
        # plot distances
        distm <- get_dist(nor, distance, stand=TRUE)
        print(fviz_dist(distm,
                  gradient=list(low="blue", mid="white", high="red")))
        out.png <- sprintf("%s/DESeq2_%s_distances.png",
            output.folder, fun.name)
         as.png(fviz_dist(distm,
                    gradient=list(low="blue", mid="white", high="red")),
                    out.png)
        continue.on.enter("Press [ENTER] to continue ")
    }
    
    cat("Cluster membership counts\n")
    print(table(member))
    # Characterizing clusters 
    cat("\nMeans by cluster in the normalized data\n")
    print(aggregate(nor, list(member), mean))
    cat("\nMeans by cluster in the non-normalized data\n")
    print(aggregate(data.table, list(member), mean))

    for (i in (min(member):max(member))) {
        clus.i <- member[ member == i ]
        data.i <- annotated.data[annotated.data$ensembl.gene.id 
                                 %in% names(clus.i), ]
        # sort by l2FC
        data.i <- data.i[order(data.i$log2FoldChange, decreasing=T), ]
        cat("Showing annotation for cluster", i, "(", length(clus.i), ") elements)\n")
        #View(data.i[ , c("log2FoldChange", "GENENAME", "entrezgene_description", "definition_1006")],
        #     title=paste("Cluster no.", i, "(", length(clus.i), ") elements)"))
        d <- data.i[ , c("log2FoldChange", "GENENAME", "entrezgene_description", "definition_1006")]
        d[ is.na(d) ] <- 'NA'
        t <- paste("Cluster no.", i, "(", length(clus.i), ") elements")
        w <- show.data.frame(d, t, TRUE)
        visible(w) <- FALSE
        #visible(w) <- TRUE
#        keypress()
         clus.file <- sprintf("%s/DESeq2_%s_nc=%03d_c=%03d.tab",
             output.folder, fun.name, max(member), i)
         write.table(data.i, file=clus.file,
             row.names=T, col.names=T, sep='\t')
    }
    
    # plot a PCA of the clusters
    # If there are more than two dimensions (variables) fviz_cluster will perform
    # principal component analysis (PCA) and plot the data points according to
    # the first two principal components that explain the majority of the
    # variance.
    print(fviz_cluster(cl, data=nor, geom="point", ellipse.type="convex"))
    out.png <- sprintf("%s/DESeq2_%s_nc=%03d_PCA.png",
            output.folder, fun.name, max(member))
    as.png(fviz_cluster(cl, data=nor, geom="point", ellipse.type="convex"),
           out.png)
    continue.on.enter("Press [ENTER] to continue ")
    # to plot gene names instead of ensembl.ids we need to change rownames to genes
    gor <- nor 
    rownames(gor) <- paste(annotated.data[ , "GENENAME"], 1:length(rownames(gor)))
    rcl <- cl
    if (fun.name == 'hcut') {
        cat("Setting names to gene names\n")
        dimnames(rcl$data)[[1]] <- rownames(gor)
    } else if (fun.name == 'pam') {
        names(cl$cluster <- rownames(gor))
    }
    print(fviz_cluster(rcl, data=gor, 
              ellipse.type="convex", show.clust.cen=FALSE, labelsize=6))
    out.png <- sprintf("%s/DESeq2_%s_nc=%03d_PCA_genes.png",
            output.folder, fun.name, max(member))
    as.png(fviz_cluster(rcl, data=gor, 
              ellipse.type="convex", show.clust.cen=FALSE, labelsize=6),
           out.png)
    continue.on.enter("Press [ENTER] to continue ")
    
    return(cl)
}





#' build.offline.annotation
#'
#' Try different strategies to buid the refernece genome annotation 
#' database online.
#' 
#' @param
#' 
#' @return
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
build.offline.annotation <- function(
                                     folder,
                                     db.dir='EnsDb.Ggallus.v106',
                                     target.organism,
                                     reference.gtf,
                                     release,
                                     ens.version,
                                     user, 
                                     pass
                                     author,
                                     maintainer,
                                     license
                                    ) 
{

    # We'll try to build an EnsDb Package
    # so we can keep it locally for future use
    #EnsDbPackageDir <- paste(folder, 'EnsDb.Ggallus.v106', sep='/')
    #EnsDbPackageDir <- paste(folder, 'EnsDb.Cjaponica.v105', sep='/')
    EnsDbPackageDir <- paste(folder, db.dir, sep='/')

    # source script included with package 'ensembldb'
    scr <- system.file("scripts/generate-EnsDbs.R", package = "ensembldb")
    source(scr)

    # use our local GTF file if that didn't work
    if (dir.exists(EnsDbPackageDir) ) {
        sqlite <- paste(target.organism, ".", release, ".", ens.version, ".sqlite", sep='/')
        DBFile <- paste(folder, '/', sqlite)
    } else {
        # generate SQLite database in place from GTF. Produces e.g.
        #    ./gallus_gallus.GRCg6a.106.sqlite
        gtf.edb <- ensDbFromGtf(gtf=reference.gtf, 
	        organism=target.organism,
                genomeVersion=release,
                version=ens.version,
                destDir=folder)			# lacks entrezid
        #dir()
        
        # file name of the SQLite database created
        #DBFile <- paste(folder, 'gallus_gallus.GRCg6a.106.sqlite', sep='/')
        sqlite <- paste(target.organism, ".", release, ".", ens.version, ".sqlite", sep='/')
        DBFile <- paste(folder, '/', sqlite)

        # we'll select the SQLite database file generated from the GTF to
        # create an EnsDb R package 
        system(paste("mv", sqlite, DBFile))

        makeEnsembldbPackage(DBFile, version=ens.version, 
    		             maintainer=maintainer, 
    		             author=author,
                             destDir=folder, license=license)
    }

    # retry with GFF if that didn't work
    if ( ! dir.exists(EnsDbPackageDir) ) {
        # As an alternative, we may use the GFF3 file (which should contain
        # the same information
        gff.edb <- ensDbFromGff(gff=reference.gff, 
                organism=target.organism,
                genomeVersion=release,
                version=ens.version,
                destDir=folder)
        # but it seems to lack entrezid and transcript_id

        #dir()
        # file name of the SQLite database created
        sqlite <- paste(target.organism, ".", release, ".", ens.version, ".sqlite", sep='/')
        DBFile <- paste(folder, '/', sqlite)

        # we'll select the SQLite database file generated from the GFF to
        # create an EnsDb R package 
        system(paste("mv", sqlite, DBFile))

        makeEnsembldbPackage(DBFile, version=ens.version, 
    		             maintainer=maintainer, 
    		             author=author,
                             destDir=folder, license=license)
    } else {
        sqlite <- paste(target.organism, ".", release, ".", ens.version, ".sqlite", sep='/')
        DBFile <- paste(folder, '/', sqlite)
    }

    # we still have a third option, build it from a MySQL file
    if ( ! file.exists(EnsDbPackageDir) ) {
        # This shouln't be needed because we have already done it above
        # using the GTF/GFF3 file.
        # This is an alternate way to do generate the EnsDb package from
        # MySQL data downloaded from ENSEMBL
        local.mysql.db <- paste("mysql/", 
        	tolower(target.organism), 
                "_core_", ens.version, "_", release, sep='')
        #local.mysql.db <- "mysql/coturnix_japonica_core_104_2"
        #local.mysql.db <- "mysql/gallus_gallus_core_105_6"
        createEnsDbForSpecies(ens_version = ens.version,
                user = user, 
                pass=pass,
                host = "localhost", 
                local_tmp=local.mysql.db, 
                species=target.organism,
                dropDb=FALSE)
        # This should create a .sqlite file like ensDbFromGtf above,
        # named, e.g.
        #	Gallus_gallus.GRCg6a.106.sqlite
        #       ^
        # with similar contents to the one produced from ensDBFromGtf/Gff;
        # if none of these does work, then tweak the process by hand in
        # source("scripts/generate-EnsDBs.R")
        sqlite <- paste(target.organism, ".", release, ".", ens.version, ".sqlite", sep='/')
        DBFile <- paste(folder, '/', sqlite)
        system(paste("mv", sqlite, DBFile))
        makeEnsembldbPackage(DBFile, version=ens.version, 
    		             maintainer=maintainer, 
    		             author=author",
                             destDir=folder, license=license)
                             
    } 
    return(DBFile)
}



#' get.biomart.ensembl.annotation
#'
#' obtain ENSEMBL-related annotation using BiomaRt either from a local
#' cache file previously saved or directly from the network
#' 
#' @param
#' 
#' @return
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
get.biomart.ensembl.annotation <- function(mart.db, folder) {

    if ( ! file.exists( paste(folder, 'biomart.ensembl.tab', sep='/')) ) {
        bm.ensembl.annot <- getBM(attributes=c(
                               "ensembl_gene_id", 
                               "ensembl_transcript_id", 
                               "start_position", "end_position", 
                               "chromosome_name", "gene_biotype", 
                               "description"),
                            mart=mart.db)
        write.table(bm.ensembl.annot, paste(folder, '/biomart.ensembl.tab', sep=''), 
	        row.names=T, col.names=T, sep='\t')
    } else {
        bm.ensembl.annot <- read.table(paste(folder, '/biomart.ensembl.tab', sep=''), 
	        header=T, sep='\t')
    }
    return(bm.ensembl.annot)
}


#' get.biomart.entrez.annotation
#'
#' obtain ENTREZ/NCBI-related annotation from BiomaRt, using either local
#' cached data or directly from the network
#' 
#' @param
#' 
#' @return
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
get.biomart.entrez.annotation <- function(mart.db, folder) {
    if ( ! file.exists(paste(folder, '/biomart.entrez.tab', sep='')) ) {
        bm.entrez.annot <- getBM(attributes=c(
                               "ensembl_gene_id", 
                                "entrezgene_id", "entrezgene_accession", "entrezgene_description"),
                            mart=mart.db)
        write.table(bm.entrez.annot, paste(folder, '/biomart.entrez.tab', sep=''), 
	        row.names=T, col.names=T, sep='\t')
    } else {
        bm.entrez.annot <- read.table(paste(folder, '/biomart.entrez.tab', sep=''), 
	        header=T, sep='\t')
    }
    return(bm.entrez.annot)
}



#' get.biomart.go.annotation
#'
#' obtain GO annotation from BiomaRt
#' 
#' @param
#' 
#' @return
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
get.biomart.go.annotation <- function(mart.db, folder) {
    if ( ! file.exists(paste(folder, '/biomart.go.tab', sep='')) ) {
        bm.go.annot <- getBM(attributes=c(
                               "ensembl_gene_id", 
                                "go_id", "name_1006", "definition_1006", "go_linkage_type", "namespace_1003"), 
                           mart=mart.db)
        write.table(bm.go.annot, paste(folder, '/biomart.go.tab', sep=''), 
	        row.names=T, col.names=T, sep='\t')
    } else {
        bm.go.annot <- read.table(paste(folder, '/biomart.go.tab', sep=''), 
	        header=T, sep='\t')
    }
    return(bm.go.annot)
}



#' get.biomart.goslim.annotation
#'
#' obtain GOslim annotatin from BiomaRt
#' 
#' @param
#' 
#' @return
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
get.biomart.goslim.annotation <- function(mart.db, folder) {
    if ( ! file.exists(paste(folder, '/biomart.goslim.tab', sep='')) ) {
        bm.goslim.annot <- getBM(attributes=c(
                               "ensembl_gene_id", 
                                "goslim_goa_accession", "goslim_goa_description"),
                           mart=mart.db)
        write.table(bm.goslim.annot, paste(folder, '/biomart.goslim.tab', sep=''), 
	        row.names=T, col.names=T, sep='\t')
    } else {
        bm.goslim.annot <- read.table(paste(folder, '/biomart.goslim.tab', sep=''), 
	        header=T, sep='\t')
    }
}


#' get.biomart.family.annotation
#'
#' obtain protein family (PFAM, TIGRFam, etc...) annotation from BiomaRt
#' 
#' @param
#' 
#' @return
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
get.biomart.family.annotation <- function(mart.db, folder) {
    if ( ! file.exists(paste(folder, '/biomart.fam.tab', sep='')) ) {
        bm.fam.annot <- getBM(attributes=c(
                               "ensembl_gene_id",
                               "pfam",
                               "pirsf",
                               "prints",
                               "tigrfam"
                               ), 
                           mart=mart.db)
        write.table(bm.fam.annot, paste(folder, '/biomart.fam.tab', sep=''), 
	        row.names=T, col.names=T, sep='\t')
    } else {
        bm.fam.annot <- read.table(paste(folder, '/biomart.fam.tab', sep=''), 
	        header=T, sep='\t')
    }
    return(bm.fam.annot)
}


#' get.biomart.prosite.annotation
#'
#' obtain PROSITE annotation from BiomaRt
#' 
#' @param
#' 
#' @return
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
get.biomart.prosite.annotation <- function(mart.db, folder) {
    if ( ! file.exists(paste(folder, '/biomart.prosite.tab', sep='')) ) {
        bm.prosite.annot <- getBM(attributes=c(
                               "ensembl_gene_id",
                               "scanprosite",
                               "pfscan" 
                               ), 
                           mart=mart.db)
        write.table(bm.prosite.annot, paste(folder, '/biomart.prosite.tab', sep=''), 
	        row.names=T, col.names=T, sep='\t')
    } else {
        bm.prosite.annot <- read.table(paste(folder, '/biomart.prosite.tab', sep=''), 
	        header=T, sep='\t')
    }
    return(bm.prosite.annot)
}



#' get.biomart.superfamily.annotation
#'
#' obtain SuperFam annotation from BiomaRt
#' 
#' @param
#' 
#' @return
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
get.biomart.superfamily.annotation <- function(mart.db, folder) {
    if ( ! file.exists(paste(folder, '/biomart.sfam.tab', sep='')) ) {
        bm.sfam.annot <- getBM(attributes=c(
                               "ensembl_gene_id",
                               "superfamily"
                               ), 
                           mart=mart.db)
        write.table(bm.sfam.annot, paste(folder, '/biomart.sfam.tab', sep=''), 
	        row.names=T, col.names=T, sep='\t')
    } else {
        bm.sfam.annot <- read.table(paste(folder, '/biomart.sfam.tab', sep=''), 
	        header=T, sep='\t')
    }
    return(bm.sfam.annot)
}


#' get.biomart.extra.annotation
#'
#' obtain additional miscellaneous annotation from BiomaRt
#' 
#' @param
#' 
#' @return
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
get.biomart.extra.annotation <- function(mart.db, folder) {
    if ( ! file.exists(paste(folder, '/biomart.extra.tab', sep='')) ) {
        bm.extra.annot <- getBM(attributes=c(
                               "ensembl_gene_id", 
                               "pdb",
                               #"reactome", 
                               "uniprotswissprot"), 
                           mart=mart.db)
        write.table(bm.extra.annot, paste(folder, '/biomart.extra.tab', sep=''), 
	        row.names=T, col.names=T, sep='\t')
    } else {
        bm.extra.annot <- read.table(paste(folder, '/biomart.extra.tab', sep=''), 
	        header=T, sep='\t')
    }
    return(bm.extra.annot)
}


#' biomart.merge.annotations
#'
#' merge together various annotation datasets into a single, consolidated one
#' 
#' @param
#' 
#' @return
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @noexport
#
biomart.merge.annotations <- function(ann, by="ensembl_gene_id") {
    if (! file.exists(paste(folder, '/biomaRt.annotation.txt', sep=''))) {
        # ann is a list of annotations to merge
        if (length(ann) == 1) return(ann)		# nothing to merge

        # we have at least two
        bm.annot <- merge(ann[[1]], ann[[2]], by=by)

        if (length(ann == 2) return(bm.annot)

        # we have more annotations to merge
        for (i in 3:length(ann))
            bm.annot <- merge(bm.annot, ann[[i]], by=by)
        write.table(bm.annot, file=paste(folder, '/biomaRt.annotation.txt', sep=''), 
                    sep='\t', row.names=T, col.names=T)
    } else {
        bm.annot <- read.table(paste(folder, '/biomaRt.annotation.txt', sep=''), 
	        header=T, sep='\t')}
    }
    return(bm.annot)
}




# # # # edgeR functions

#
# SAVE TOP 'N' GENES
# ------------------
# save the table of the n.genes top expressed genes sorted in various orders
eR.save.top <- function(fit, folder, n.genes=500, sort.by='p', coef=1, p.value=0.01)
{
    for (by in sort.by) {
        # default adjustment is BH
        tt <- topTable(fit, 
                       coef=coef, 
                       number=n.genes, 
                       sort.by=by, 
                       p.value=p.value)
        file <- paste(folder, 
                      '/edgeR/cmp_coef=', coef, 
                      '_top_', n.genes, 
                      '_by', by, 
                      '.tab', 
                      sep='')
        write.table(tt, file, row.names=T, col.names=T, sep='\t')
    }
}


#
# SAVE TOP 'N' GENES ANNOTATED
# ----------------------------
# we can now run again the topTable and we will have all the annotation 
# information linked 
#	n.genes <- 500 ALREADY DEFINED ABOVE
eR.save.top.annotated <- function(fit, folder, n.genes=500, sort.by='p', coef=1, p.value=0.01) {

    for (by in sort.by) {
        tt <- topTable(fit, 
                       coef=coef, 
                       number=n.genes, 
                       sort.by=by, 
                       p.value=p.value)
        file <- paste(folder, 
                      '/edgeR/cmp=', coef, 
                      '_top_', n.genes, 
                      '_by_', by, 
                      '_annotated.tab', 
                      sep='')
        write.table(tt, file, w.names=T, col.names=T, sep='\t')
    }
}



# compare to a threshold and save
# SAVE TOP RESULTS RELATIVE TO THRESHOLD
# --------------------------------------
eR.save.top.treat.ann <- function(fit, folder, n.genes=500, sort.by='p', coef=1, p.value=0.01, threshold=1) {
    for (by in sort.by) {
        #     vvvvvvvv		here we use topTreat instead of topTable
        tt <- topTreat(fit, 
                       coef=coef, 
                       number=n.genes, 
                       sort.by=by, 
                       p.value=p.value)
        file <- paste(folder, 
                      '/edgeR/cmp=', coef, 
                      '_top_', n.genes, 
                      '_by_', by, 
                      '_lfc>=', threshold, 
                      '_annotated.tab', 
                      sep='')
        write.table(tt, file, row.names=T, col.names=T, sep='\t')
    }
}



# save top N genes, annotated with biomart
eR.save.top.ann.thresh <- function(fit, folder, n.genes=500, sort.by='p', coef=1, p.value=0.01, threshold=1) {
    for (by in sort.by) {
        file <- paste(folder, 
                      '/edgeR/cmp=', coef, 
                      '_top_', n.genes, 
                      '_by_', by, 
                      '_lfc>=', threshold, 
                      '_threshold_annotated.tab', sep='')
        tt <- topTable(fit, 
                       coef=coef, 
                       number=n.genes, 
                       sort.by=by, 
                       p.value=p.value)
        write.table(tt, file, row.names=T, col.names=T, sep='\t')
    }
}



# annotation functions for edgeR


eR.fit.annotate <- function(fit, ens.db, biomart) {

     if (str_subrownames(fit, 1, 3) == "ENS"))
         by='GENEID'
     else
     	by='ENTREZID'
        
     ens.ann <- ensembldb::select(ens.db, 
                  column=by, 
                  keytype=by, keys=rownames(fit), 
                  columns= c('SEQNAME', 'SYMBOL', 'DESCRIPTION',
                             'GENENAME', 'GENEID', 'ENTREZID',
                             'TXNAME', 'TXID', 'TXBIOTYPE',
                             'PROTEINID', 'UNIPROTID'
                            ))

#    ann <- ensembldb::select(ens.db, 
#              keytype= 'GENEID', keys=rownames(fit), 
#              columns= c('SEQNAME', 'SYMBOL', 'DESCRIPTION',
#                         'GENENAME', 'GENEID', 'ENTREZID', 
#                          'TXNAME', 'TXID', 'TXBIOTYPE', 
#                          'PROTEINID', 'UNIPROTID'))

    # use only first annotation
    ens.ann.1 <- ens.ann[ ! duplicated(ens.ann$GENEID), ]

    # As long as it works perfectly we continue by connecting the 
    # annotation information to the fit using a new bucket named "genes"
    # (fit is an object of class "MArrayLM" (package "limma"), DGELRT or
    # DGELM (package edgeR) which is a list. I.e. we assign the annotation 
    # to a new element named 'genes' of this list.
    #	This works if they are in the same order
    #fit$genes <- ens.ann.1
    #	This works by matching rownames and gene IDs
    fit$genes <- ens.ann.1[ match(rownames(fit), ens.ann.1$GENEID), ]
    rownames(fit$genes) <- rownames(fit)
    #
    # add biomart annotation
    fit$genes <- merge(fit$genes, biomart, by.x="GENEID", by.y="ensembl_gene_id")
    #head(fit$genes)


    return (fit)
}

eR.save.fit <- function(fit, name) {
    # SAVE THE ANNOTATED FIT
    # ----------------------
    # save all the contents of 'fit' in an RDS file
    saveRDS(fit, file=paste(name,'.rds', sep=''))
    #	'fit' can later be recovered with: 
    #		fit <- readRDS(file=paste(folder, 'edgeR/', name, '.rds', sep=''))
    # and save as well as Rdata file
    save(fit, file=paste(name, '.RData', sep=''))
    #	'fit' can later be recovered with: 
    #		fc <- load(file=paste(folder, '/edgeR/', name, '.RData', sep=''))
}

eR.save.top.fit <- function(fit, file, n.genes=500, sort.by='PValue', p.value=0.01) {    # default adjust.method is BH
    # default p.value is 1 (all genes)
    tt <- topTags(fit, n=n.genes, sort.by=sort.by, p.value=p.value)
    n <- dim(tt)[1]
    # cap at the maximum number of genes requested
    if (n > n.genes) n <- n.genes
    file.name <- paste(file, '_top_', n, '_by_', sort.by, '.tab', sep='')
    write.table(tt$table[1:n, ], file.name, 
	row.names=T, col.names=T, sep='\t')
    # if fit is annotated the annotation will also be saved
}



eR.differential.gene.expression<- function(fc, 
                                           metadata='SampleInfo.txt',
                                           threshold,
                                           ens.db,
                                           folder)
{

    # make convenience names
    counts <- fc$counts
    annotation <- fc$annotation
    stats <- fc$stats
    target <- read.delim(paste(rnaseq.out, metadata, sep='/'))

    # then we proceed to the analysis.. for this, we will need other packages
    #library(edgeR)
    #library(limma)
    #library(RColorBrewer)
    #library(gplots)

    # The next step is to distinguish the genes whose expression is significant
    # from  the ones that have an 'insignificant' expression that could be by
    # chance or irrelevant to the  situation of the cell. So we set up the
    # threshold of expression for each gene to a minimum 10-15 counts. Since
    # we will work with normalized CPM (counts per million) data, we need
    # to check which number of CPM corresponds to  to this amount of counts in
    # order to filter the data

    # convert counts to CPM
    cpm.counts <- cpm(counts)


    # we'll plot the correlation of cpm and counts to see which is the number of
    # cpm that corresponds to 10-15 counts minimum. we see this information
    # graphically we use for example column 1... we could check every file
    # independently but since they are all similar in number of reads there 
    # should be no need

    out.png <- paste(folder, '/edgeR/img/edgeR_CPMcounts.png', sep='')
    as.png( {
            plot(counts[,1:12], cpm.counts[,1:12], xlim=c(0,20), ylim=c(0,5))
            abline(v=10, col=1)
        }, out.png )

    # we will set up the threshold to 0.25 according on the plot we have drawn
    # this command will return a table of trues and falses, then, we want to keep
    # only the rows that exceed the threshold at least in three different
    # cases(experiments .bam files)


    #thres <- cpm.counts > 0.25		# Coturnix
    thres <- cpm.counts > threshold		# 0.5 for Gallus
    keep <- rowSums(thres) >=3
    if(VERBOSE) print(table(keep))


    # then we store in a different variable the genes whose counts exceed the
    # threshold and visualise the content of the new variable to see the amount of
    # remaining genes

    counts.keep <- counts[keep,]
    dim(counts.keep)
    if ( VERBOSE ) {
        # some paranoid manual checks
        # we are using counts instead of cpm
        cpmavg <- data.frame(vd000.0=apply(counts.keep[,13:15],1,mean), 
                             vd000.1=apply(counts.keep[,1:3],1,mean),
                             vd001.0=apply(counts.keep[,10:12],1,mean), 
                             vd010.0=apply(counts.keep[,7:9],1,mean), 
                             vd100.0=apply(counts.keep[,4:6],1,mean) 
                             )

        f.c000.0 <- data.frame(vd000.0_vd000.1=(cpmavg$vd000.0 / cpmavg$vd000.1),
                               vd000.0_vd001.0=(cpmavg$vd000.0 / cpmavg$vd001.0),
                               vd000.0_vd010.0=(cpmavg$vd000.0 / cpmavg$vd010.0),
                               vd000.0_vd100.0=(cpmavg$vd000.0 / cpmavg$vd100.0)
                              ) 
        row.names(f.c000.0) <- row.names(cpmavg)
        l.f.c000.0 <- log2(f.c000.0)
        write.table(f.c000.0, file=paste(out.dir, '/edgeR/hand.fc_000.0.tab', sep=''), sep='\t')
        write.table(l.f.c000.0, file=paste(out.dir, '/edgeR/hand.lfc_000.0.tab', sep=''), sep='\t')
        write.table(cpmavg, file=paste(out.dir, '/edgeR/hand.cpmavg.tab', sep=''), sep='\t')
    }

    # We have manipulated the data discarding whatever is not of high interest
    # and now we need to see the differencial expression and highlight the
    # differences among the cells

    # convert the counts.keep to a DGEList
    dge <- DGEList(counts.keep)


    # do TMM normalization
    dge <- calcNormFactors(dge)
    if (VERBOSE) print(dge$samples)

    # Plot the library size of the different samples.
    out.png <- paste(out.dir, '/edgeR/img/edgeR_sample_lib_size.png', sep='')
    as.png(barplot(dge$samples$lib.size, cex.names= 1, main = "Library Size", col = dge$samples$group, names.arg=dge$samples$group, ylab = "Reads"), out.png, overwrite=TRUE)


    # Now, do some quality control plots, barplots and boxplots
    # we need normalized data counts so we take the logarithm
    # we need to group together all the experiment data that correspond to the
    # same  cell type (target) and asign a different color to each one so as to
    # distinguish them graphically; we know we have 4 groups so we have to specify
    # 4 different colors
    out.png <- paste(folder, '/edgeR/img/edgeR_sample_lib_size.png', sep='')
    as.png(barplot(dge$sample$lib.size), out.png)

    logcpm <- cpm(dge$counts, log=TRUE)
    group.col <- c('red', 'blue', 'green', 'yellow')[target$[ ,design.column]] 

    out.png <- paste(folder, '/edgeR/img/edgeR_log2_cpm.png', sep='')
    as.png( {
            par(mfrow=c(1,1))
            boxplot(logcpm, xlab='', ylab=' Log2 counts per million', 
	            col= group.col, las=2, outline=FALSE)
            abline(h=median(logcpm), col='red')
            title('Boxplots of logCPMs unnormalised')
        }, out.png )

    # then we produce an MDS plot to see any significant difference between the
    # groups
    out.png <- paste(folder, '/edgeR/img/edgeR_mds_plot.png', sep='')
    as.png( {
            par(mfrow= c(1,1))
            plotMDS(dge, col=group.col)
        }, out.png)


    # now is time to see the differences in the expression (variance). We have
    # to apply a funcion that calculates the variance by rows (genes) and then
    # retrieve the n.genes most DE genes

    logcounts <- cpm(dge, log =TRUE)
    var_genes <- apply(logcounts, by.rows, var)
    select_var <- names(sort(var_genes, decreasing=TRUE))[1:n.genes]
    highly_var <- logcounts[select_var,]
    dim(highly_var)

    # we now set up the colours we will use for the heatmap plot 
    # and then we create all the colors in between in the palette
    mypalette <- brewer.pal(11, 'RdYlBu')
    morecolors <- colorRampPalette(mypalette)

    # we plot the heatmap.2 (gplots) without a line (trace), scale by row
    # (difference in color) margins (something about the labels used), also we
    # reverse the colors because by default the red is associated with low
    # expression and we are more used to it meaning "hot"
    out.png <- paste(folder, '/edgeR/img/edgeR_heatmap.png', sep='')
    as.png( {
            #margins <- par("mar")
            #par(mar=c(25, 5, 5, 10))
            heatmap.2(highly_var, 
                      col= rev(morecolors(50)), 
                      trace='none', 
                      ColSideColors=group.col,
                      scale='row', 
                      margins= c(15,5))
            #par(mar=margins)
        }, out.png)

    # This plot is of limited use. We'd better have other names for rows and 
    # columns and plot the n first (most variable) to see them well
    n <- 50
    high_var <- highly_var
    colnames(high_var) <- gsub("_R1*.[Bb][Aa][Mm]", "", colnames(highly_var))
    
    ## NOTE: we should add the annotation here, before doing the next plot
     if (str_subrownames(fit, 1, 3) == "ENS"))
         by='GENEID'
     else
     	by='ENTREZID'
    name <- ensembldb::select(ens.db, keys=rownames(high_var), 
                       column=by, keytype=by, 
                       columns=c('GENENAME'))
    out.png <- paste(folder, '/edgeR/img/edgeR_heatmap.', n,'.png', sep='')
    as.png( {
            #margins <- par("mar")
            #par(mar=c(25, 5, 5, 10))
            heatmap.2(high_var[1:n,1:15], 
                          col=rev(morecolors(50)), 
                          trace='none', 
                          ColSideColors=group.col,
                          scale='row', 
                          margins= c(15,5),
                          labRow=name[1:n, 2])
            #par(mar=margins)
        }, out.png )


    # We can automate the estimation of the dispersion and add it to the dge object
    dge <- estimateCommonDisp(dge)

    # And now we can estimate gene-wise dispersion estimates allowing for a
    # possible trend with average count size. These will allow us to use a GLM
    # instead of a plain LM.
    # This gives us the BCV (Biological Coefficient of Variation) between samples
    dge <- estimateGLMTrendedDisp(dge)
    dge <- estimateTagwiseDisp(dge)

    out.png <- paste(folder, '/edgeR/img/edgeR_BCV_dispersions.png', sep='')
    as.png(plotBCV(dge), out.png)

    return(dge)
}

eR.dge.voom.variation.analysis <- function(dge, design.column, folder)
{
     ############## JR #########################
    # 
    # edgeR user's guide
    # Chapter 3 
    # 
    # Specific experimental designs 3.1 Introduction In this chapter, we outline
    # the principles for setting up the design matrix and forming contrasts for
    # some typical experimental designs.
    # 
    # Throughout this chapter we will assume that the read alignment, normalization
    # and dispersion estimation steps described in the previous chapter have
    # already been completed. We will assume that a DGEList object y has been
    # created containing the read counts, library sizes, normalization factors and
    # dispersion estimates.
    # 
    # 3.2 Two or more groups
    # 
    # 3.2.1 Introduction
    # 
    # The simplest and most common type of experimental design is that in which a
    # number of experimental conditions are compared on the basis of independent
    # biological replicates of each condition. Suppose that there are three
    # experimental conditions to be compared, treatments A, B and C, say. The
    # samples component of the DGEList data object might look like:
    # 
    # > y$samples
    # group lib.size norm.factors
    # Sample1 A 100001 1
    # Sample2 A 100002 1
    # Sample3 B 100003 1
    # Sample4 B 100004 1
    # Sample5 C 100005 1
    # 
    # Note that it is not necessary to have multiple replicates for all the
    # conditions, although it is usually desirable to do so. By default, the
    # conditions will be listed in alphabetical order, regardless of the order that
    # the data were read:
    # 
    # > levels(y$samples$group)
    # [1] "A" "B" "C"
    # 29
    # 
    # 3.2.2 Classic approach
    # 
    # The classic edgeR approach is to make pairwise comparisons between the
    # groups. For example,
    # 
    # > et <- exactTest(y, pair=c("A","B"))
    # > topTags(et)
    # 
    # will find genes differentially expressed (DE) in B vs A. Similarly
    # 
    # > et <- exactTest(y, pair=c("A","C"))
    # 
    # for C vs A, or
    # 
    # > et <- exactTest(y, pair=c("C","B"))
    # 
    # for B vs C.
    # 
    # Alternatively, the conditions to be compared can be specified by number, so
    # that
    # 
    # > et <- exactTest(y, pair=c(3,2))
    # 
    # is equivalent to pair=c("C","B"), given that the second and third levels of
    # group are B and C respectively.
    # 
    # Note that the levels of group are in alphabetical order by default, but can
    # be easily changed.
    # 
    # Suppose for example that C is a control or reference level to which
    # conditions A and B are to be compared. Then one might redefine the group
    # levels, in a new data object, so that C is the first level:
    # 
    # > y2 <- y
    # > y2$samples$group <- relevel(y2$samples$group, ref="C")
    # > levels(y2$samples$group)
    # [1] "C" "A" "B"
    # 
    # Now
    # 
    # > et <- exactTest(y2, pair=c("A","B"))
    # 
    # would still compare B to A, but
    # 
    # > et <- exactTest(y2, pair=c(1,2))
    # 
    # would now compare A to C.
    # 
    # When pair is not specified, the default is to compare the first two group
    # levels, so
    # 
    # > et <- exactTest(y)
    # 
    # compares B to A, whereas
    # 
    # > et <- exactTest(y2)
    # 
    # compares A to C.
    # 
    # 
    # 
    # 3.2.3 GLM approach
    # 
    # The glm approach to multiple groups is similar to the classic approach, but
    # permits more general comparisons to be made. The glm approach requires a
    # design matrix to describe the treatment conditions. We will usually use the
    # model.matrix function to construct the design matrix, although it could be
    # constructed manually. There are always many equivalent ways to define this
    # matrix. Perhaps the simplest way is to define a coefficient for the
    # expression level of each group:
    # 
    # > design <- model.matrix(~0+group, data=y$samples)
    # > colnames(design) <- levels(y$samples$group)
    # > design
    # A B C
    # Sample1 1 0 0
    # Sample2 1 0 0
    # Sample3 0 1 0
    # Sample4 0 1 0
    # Sample5 0 0 1
    # attr(,"assign")
    # [1] 1 1 1
    # attr(,"contrasts")
    # attr(,"contrasts")$group
    # [1] "contr.treatment"
    # 
    # Here, the 0+ in the model formula is an instruction not to include an
    # intercept column and instead to include a column for each group.
    # 
    # One can compare any of the treatment groups using the contrast argument of
    # the glmQLFTest or glmLRT function. For example,
    # 
    # > fit <- glmQLFit(y, design)
    # > qlf <- glmQLFTest(fit, contrast=c(-1,1,0))
    # > topTags(qlf)
    # 
    # will compare B to A. The meaning of the contrast is to make the comparison
    # -1*A + 1*B + 0*C, which is of course is simply B-A.
    # 
    # The contrast vector can be constructed using makeContrasts if that is
    # convenient. The above comparison could have been made by
    # 
    # > BvsA <- makeContrasts(B-A, levels=design)
    # > qlf <- glmQLFTest(fit, contrast=BvsA)
    # 
    # One could make three pairwise comparisons between the groups by
    # 
    # > my.contrasts <- makeContrasts(BvsA=B-A, CvsB=C-B, CvsA=C-A, levels=design)
    # > qlf.BvsA <- glmQLFTest(fit, contrast=my.contrasts[,"BvsA"])
    # > topTags(qlf.BvsA)
    # > qlf.CvsB <- glmQLFTest(fit, contrast=my.contrasts[,"CvsB"])
    # > topTags(qlf.CvsB)
    # > qlf.CvsA <- glmQLFTest(fit, contrast=my.contrasts[,"CvsA"])
    # > topTags(qlf.CvsA)
    # 
    # which would compare B to A, C to B and C to A respectively.
    # 
    # 
    # Any comparison can be made. For example,
    # 
    # > qlf <- glmQLFTest(fit, contrast=c(-0.5,-0.5,1))
    # 
    # would compare C to the average of A and B. Alternatively, this same contrast
    # could have been specified by
    # 
    # > my.contrast <- makeContrasts(C-(A+B)/2, levels=design)
    # > qlf <- glmQLFTest(fit, contrast=my.contrast)
    # 
    # with the same results.
    # 
    ############## NOTE END #########################

    # we know that the variability we see in the expression depends on infection so
    # we have to take this into account and create a model
    # 0+ forces the design to include all groups and not have
    # an intercept (reference) column

    #design <- model.matrix(~0 + target[ , design.column])
    #design
    design <- model.matrix(~ target[ , design.column])
    colnames(design) <- levels(as.factor(target[ , design.column]))
    rownames(design) <- rownames(target)
    if (VERBOSE) print(design)

    # IF WE DO NOT INCLUDE THE "~ 0 + " IN THE FORMULA, MODEL.MATRIX()
    # WILL USE AS REFERENCE THE FIRST ALPHABETICAL ORDER LEVEL !!!
    # 
    # another problem is that we are limited to the comparisons defined by
    # coefficents, if we want more control we need to use makeContrasts
    #
    # the same applies for the following analyses

    # Let us test for differential expression with the DGE data using a GLM:
    #   first, fit genewise GLMs
    gfit <- glmFit(dge, design)
    if (VERBOSE) {
        print(names(gfit))
        print(head(coef(gfit)))
    }
    
    # We can now conduct Likelihood Rato tests and show the top genes
    # for the selected comparison coefficients (reference vs. coeff)
    lrt <- glmLRT(gfit, coef=1)	# coef = 1... length(gfit$coefficients)
    print(topTags(lrt))
    # the problem here us that we are limited to comparisons defined by
    # the fitting coefficients (ref vs. variable-in-coeff).
    # If we want more control we need to use makeContrasts()
    # This is due to the formula used:  ~ table[ , design.column ]
    # If we used ~ 0 + table[ , design.column] we would have diffs for
    # all levels, but no reference.

    # so, for now, let's do a VOOM analysis with the formula employed
    # in the current design

    #voom transform the data 
    v <- voom(dge, design, plot=FALSE)
    out.png <- paste(folder, '/edgeR/img/edgeR_voom.png', sep='')
    as.png( {
            par(mfrow= c(1,1))
            voom(dge, design, plot=TRUE)
        }, out.png )



    # Carry on a summary variation analysis using the VOOM transformed data
    #
    # we fit the results of the voom depending on the model we created with our
    # design, we overwrite the fit with the use of eBayes (stat model), and
    # depending on that we decide abut which tests fit the best and summarise the
    # results: topTable gives us helpful information about everything

    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    results <- decideTests(fit)
    if (VERBOSE) print(summary(results))
    #topTable(fit, coef= 1, sort.by='p')

    eR.save.top(fit, folder, n.genes, sort.by=c('p', 'B', 'logFC', 'AveExpr'))

    # the problem here is that we are limited to the comparisons defined by
    # coefficents, if we want more control we need to use makeContrasts
    #
    # but for now we will leave it here, as we will do it in more detail
    # with DESeq2

    return(fit)

}

eR.fit.annotate.ensembl.biomart.save <- function(fit,
                            ens.db,
                            bm.annot.1,
                            folder,
                            n.genes=1000	# top N genes to save in tables
                            )
{
    #---------------------------------------------------------------------------
    #Now is time to connect all the results we have with the existing
    #information  we know from the literature, so we will retrieve infos from
    #ENSEMBL and connect with the genes  we have kept. We are interested only in
    #genes and transcriptomes (non-characterized genes) from the organism
    #of interest

    # we created two different databases one from the gtf file (small archive)
    # and one from the data we extracted from ENSEMBL and stored it in a sqlite 
    # file

    # we'll use ens.db
    if (str_sub(rownames(fit, 1,3) == 'ENS')
        by <- 'GENEID'
    else
        by <- 'ENTREZID'

    # once we have prepared the database, we now want to extract the annotation 
    # for the genes in 'fit' to add useful information to our genes
    ens.ann <- ensembldb::select(ens.db, 
                      column=by, keytype=by, keys=rownames(fit), 
                      columns= c('SEQNAME', 'SYMBOL', 'DESCRIPTION', 
                                 'GENENAME', 'GENEID', 'ENTREZID',
                                 'TXID', 'TXBIOTYPE',
                                 'PROTEINID', 'UNIPROTID'
                                ))

    #ann <- ensembldb::select(ens.db, 
    #              keytype= 'GENEID', keys=rownames(fit), 
    #              columns= c('SEQNAME', 'SYMBOL', 'DESCRIPTION',
    #                         'GENENAME', 'GENEID', 'ENTREZID', 
    #                          'TXNAME', 'TXBIOTYPE', 
    #                          'PROTEINID', 'UNIPROTID'))



    # SAVE ANNOTATION
    # ---------------
    # this is all the annotation for all the genes in 'fit'
    write.table(ens.ann, file=paste(folder, '/ensembl.annotation.txt', sep=''), 
	    sep='\t', row.names=T, col.names=T)

    # check if the amount of genes we have is the same as the number of 
    # the annotations that we have extracted
    if ( ! table(ens.ann$GENEID==rownames(fit)) ) {
        cat("Annotation does not match genes\n")
        cat("Using only one entry per gene\n")
        ens.ann.1 <- ens.ann[ ! duplicated(ens.ann$GENEID), ]
    } else {
        ens.ann.1 <- ens.ann
    }
    write.table(ens.ann.1, file=paste(folder, '/ensembl.annotation.1st.txt', sep=''), 
	    sep='\t', row.names=T, col.names=T)


    # As long as it works perfectly we continue by connecting the 
    # annotation information to the fit using a new bucket named "genes"
    # (fit is an object of class "MArrayLM" (package "limma") which is a
    # list. I.e. we assign the annotation to a new element named 'genes'  
    # of this list.
    #
    #	This checks that genes and annotation go in the same order)
    fit$genes <- ens.ann.1[ match(rownames(fit), ens.ann.1$GENEID), ]

    
    # and now we will also add to 'fit$genes' the biomaRt annotation
    fit$genes <- merge(fit$genes, bm.annot.1, by.x="GENEID", by.y="ensembl_gene_id")
    if (VERBOSE) print(head(fit$genes))

    # SAVE THE FULL ANNOTATED FIT
    # ---------------------------
    # save all the contents of 'fit' in an RDS file
    saveRDS(fit, file=paste(folder, '/edgeR/annotatedVOOMfit.rds', sep=''))
    #	'fit' can later be recovered with: 
    #		fit <- readRDS(file=paste(folder, '/annotatedVOOMfit.rds', sep=''))
    # and save as well as Rdata file
    save(fit, file=paste(folder, '/edgeR/annotatedVOOMfit.RData', sep=''))
    #	'fit' can later be recovered with: 
    #		fc <- load(file=paste(folder, '/annotatedVOOMfit.RData', sep=''))


    # save top annotated genes
    eR.save.top.annotated(fit, folder, n.genes, 
                          sort.by=c('p', 'logFC', 'AveExpr'))
    

    return(fit)

}


eR.fit.treat <- function(fit,
                         threshold=1,
                         folder,
                         verbose=T)
{
    #Testing relative to a threshold: 1 means a 2x fold change
    fit.thres <- treat(fit, lfc=threshold)
    if (verbose) {
        res.thres <- decideTests(fit.thres)
        print(summary(res.thres))
    }
    if (VERBOSE)
        print(topTreat(fit.thres, coef=1, sort.by='p'))
    
    # save topTreat data
    eR.save.top.treat.ann(fit.thres, folder, n.genes, 
                          sort.by=c('p', 'logFC', 'AveExpr'))

    # also save fitted to a threshold as tables
    eR.save.top.ann.thresh(fit.thres, folder, n.genes, 
                           sort.by=c('p', 'B', 'logFC', 'AveExpr'))
    
    return(fit.thres)
}

eR.dge.all.comparisons <- function(dge, 
                               design.column,
                               ens.db,
                               biomart,
                               folder)
{

    # ---------------------------------------------------------------------
    # Do all comparisons at once

    # for Coturnix we use 'viral.dose' as basis for comparison
    #design.column <- 'viral.dose'
    # for Gallus we will use 'Src'
    design <- model.matrix(~0 + target[ , design.column])
    colnames(design) <- levels(as.factor(target[ , design.column]))
    grps <- levels(as.factor(target[ , design.column]))
    colnames(design) <- grps
    n.grps <- length(grps)

    qlfit <- glmQLFit(dge, design, robust=TRUE, abundance.trend=TRUE)
    png.file <- paste(folder, '/edgeR/img/edgeR_QLFit_', design.column, '.png', sep='')
    as.png(plotQLDisp(qlfit), png.file)

    # if we wanted to apply a log-fold-change theshold, we could do
    # the calculations using qlfit here before testing for contrasts 
    # using e.g.
    # qlf.cmp <- glmTreat(glfit, cef=1..ncol(glfit$design), lfc=threshold)
    # or
    # qlf.cmp <- glmTreat(glfit, contrast=cmp, lfc=threshold)
    #

    eR.data <- list()
    for (i in grps) {
        for (j in grps) {
            if (i == j) next	# ignore self-comparisons
            formula <- paste(i, '-', j)
            cat("\nComputing DGE:", i, '-', j, '\n')
            cmp <- makeContrasts(formula, levels=design)
            # glmQLFTest is similar to glmLRT except it uses Bayes quasi-likelihood
            # the P-values are always >= those produced by glmLRT
            qlf.cmp <- glmQLFTest(qlfit, contrast=cmp)
            # or if a threshold has been defined
            #qlf.cmp <- glmTreat(glfit, contrast=cmp, lfc=threshold)

            # get summary of up/down regulated genes
            if (VERBOSE) {
                print(summary(decideTests(qlf.cmp)))
            }
            png.file <- paste(folder, "/edgeR/img/edgeR_QLF_MD_", 
                              design.column, '_', i, '-', j, '.png', sep='')
            as.png( { 
                    plotMD(qlf.cmp)
                    abline(h=c(-1, 1), col="darkgreen")
                    }, png.file)
                    
            if (INTERACTIVE) {
                plotMD(qlf.cmp)
                abline(h=c(-1, 1), col="darkgreen")
                continue.on.enter("Press [RETURN] to continue: ")
            }
            
            # ANNOTATE the results without saving them
            ### NOTE consider using eR.fit.annotate.ensembl.biomart
            qlf.cmp <- eR.fit.annotate(qlf.cmp, ens.db, biomart)

            name <- paste(folder, '/edgeR/fit_', design.column, '_', i, '_-_', j, '_annot', sep='')
            # qlf.cmp is a list of tables, if we want to save it,
            # we'll need to save the whole object
            # will add .rds and .RData to the files created
            #
            eR.save.fit(qlf.cmp, name)
            # defaults: n.genes=500, sort.by='PValue', p.valu=0.01
            # we'll save all (<=100.000) significant genes

            name <- paste(folder, '/edgeR/comp_', design.column, '_', i, '_-_', j, '_annot', sep='')
	    # qlf.cmp$table is a table with logFC, logCPM, F and PValue
            # that is what weill be saved when using topTags and write.table
            # will add "_top_" n "_by_" sort.by
            # defaults: n.genes=500, sort.by='PValue', p.value=0.01
            eR.save.top.fit(qlf.cmp, file=name, n.genes=100000)

            # we can use limma to test for over-representation of gene
            # ontology (GO) terms or KEGG pathways with goana() or kegga()
            # using the entrez.gene.ids of DE genes, to an FDR of 0.05 (default)
	    # we can use species.KEGG="gga" or "cjo"
	    # ( see https://www.kegg.jp/kegg/catalog/org_list.html )
	    #
            # eR.go <- goana(qlf.cmp, species="Cj")
	    # eR.go <- goana(qlf.cmp, species="Gg")
	    # eR.kegg <- kegga(qlf.cmp, species.KEGG="cjo")
	    # eR.kegg <- kegga(qlf.cmp, species.KEGG="gga")
	    # topGO(go, sort="up", number=n.genes)
	    # topKEGG(keg, sort="up", number=n.genes)

            qlf.result <- list(
	                      eR.cmp=qlf.cmp
                              #, eR.go=eR.go
			      #, eR.kegg=eR.kegg
			      )
            eR.data[[formula]] <- qlf.result 
            #print(topTags(qlf.cmp))
        }
    }
    
    return(eR.data)
}



##### DESEQ2



ds2.get.annotation.ens.org <- function(fit, ens.db, org.db)
{
    if (sub_str(rownames(fit), 1, 3) == 'ENS')

    # annotate with ensembl ens.db
    ens.ann <- ensembldb::select(ens.db, 
                      column='GENID', keytype= 'GENEID', keys=rownames(fit), 
                      columns= c('SEQNAME', 'SYMBOL', 'DESCRIPTION', # no longer available
                                 'GENENAME', 'GENEID', 'ENTREZID', # empty
                                 'TXID', 'TXBIOTYPE', # these make the call fail
                                 'PROTEINID', 'UNIPROTID' # no longer available
                                ))

    ann <- ensembldb::select(ens.db, 
                  keytype= 'GENEID', keys=rownames(fit), 
                  columns= c('SEQNAME', 'SYMBOL', 'DESCRIPTION',
                             'GENENAME', 'GENEID', 'ENTREZID', 
                              'TXNAME', 'TXBIOTYPE', 
                              'PROTEINID', 'UNIPROTID'))

    ens.ann.1 <- ens.ann[ ! duplicated(ens.ann$GENEID), ]

    # Add annotation to dds to keep everything in one place
    #dds$ens.annot.1 <- ens.ann.1
    #dds$bm.annot.1 <- b,.annot.1

    #
    # Extract annotation using Org object if avaiable
    #
    if ( ! is.null(org.db) ) {

        # and now we should be able to use the Org package if we successfully built it
        # at the beginning.
        geneSymbols <- mapIds(org.db, 
                              keys=as.character(dfwt_vs_0.1$genes$ENTREZID), 
                              column="ENTREZID", 
                              keytype="ENTREZID", 
                              multiVals="first")


        # retrieve go ids
        go.ann <- AnnotationDbi::select(org.Ggallus.eg.db, 
	        keys=as.character(dfwt_vs_0.1$genes$ENTREZID), 
                columns=c("ENTREZID", "GO", "GOALL", "ONTOLOGY","ONTOLOGYALL"), 
                keytype="GID", 
                multiVals="CharacterList")

        # retrieve corresponding descriptions
        GOdescription <- AnnotationDbi::select(GO.db, keys=go.ann$GO, 
                         columns= c("GOID", "TERM", "DEFINITION", "ONTOLOGY"), 
                         keytype= "GOID")

    }
    
    return(list(ens.ann, ens.ann.1, geneSymbols, go.ann, GOdescription))

}



ds2.dds.compare.annotate.save <- function(dds, 
	column, 
        x, 
        y, 
        filterFun=ihw, 
        alpha=0.01, 
        ensembl.ann, 
        biomart.ann,
        outDir='.' ) {

    out.file.base <- paste("raw", column, x, "vs", y, sep='_')
    # obtain comparison results
    cmp <- results(dds, 
                   contrast=c(column, x, y),
                   filterFun=filterFun,
                   alpha=alpha)
    # convert to data frame and add row names (ensembl.gene.id) 
    # as an additional column
    cmp.df.a <- data.frame(ensembl.gene.id=rownames(cmp), cmp)
    # ensembl.ann matches all ENSMEBL-ID in the cmp results 
    # (but lacks description)
    cmp.df.a <- cbind(cmp.df.a, 
                      ensembl.ann[ match(cmp.df.a$ensembl.gene.id, ensembl.ann$GENEID), ])
    # bm.annot fails to annotate some entries (why?)
    cmp.df.a <- cbind(cmp.df.a, 
                      biomart.ann[ match(cmp.df.a$ensembl.gene.id, biomart.ann$ensembl_gene_id), ])

    # and now save
    # save summary
    sink(paste(outDir, '/DESeq2/', out.file.base, '_summary.txt', sep=''), split=T)
    summary(cmp)
    sink()
    # unnanotated results object
    saveRDS(cmp, paste(outDir, "/DESeq2/", out.file.base, ".rds", sep=""))
    # annotated results as data frame (table)
    write.table(cmp.df.a,
	    paste(outDir, "/DESeq2/", out.file.base, "_annotated.tab", sep=""),
	    row.names=T, col.names=T, sep='\t')
    # histogram plot
    out.png <- paste(folder, '/DESeq2/img/DESeq2_', out.file.base, '_hist.png', sep='')
    as.png( {
        margins <- par("mar")
        par(mar=c(5, 5, 5, 5))
        hist(cmp.df.a$pvalue, main=out.file.base, breaks=1/alpha, xlab="p-value")
        par(mar=margins)
    }, out.png )

    return (cmp.df.a)
}




ds2.dds.plot.and.save <- function(dds, 
                          column, x, y, 
                          filterFun=ihw, alpha=0.01,
        		  ensembl.ann, biomart.ann,
                          outDir='.',
                          save=TRUE ) {
    # base output file name
    out.base <- paste(column, ':_', x, '_x_', y, sep='')

    # Do the comparison and get the raw results (with baseMean, log2FC, p, padjusted)
    raw <- results(dds, 
                   contrast=c(column, x, y),
                   filterFun=filterFun,
                   alpha=alpha)
    # raw contains the data after comparing x and y as a DESEq2 result object
    # convert to a data frame and add row names (ensembl.gene.id) 
    # as an additional column
    raw.df <- data.frame(ensembl.gene.id=rownames(raw), raw)
    # ensembl.ann matches all ENSMEBL-ID in the raw results 
    # (but lacks description)
    raw.df <- cbind(
                raw.df, 
                ensembl.ann[ match(raw.df$ensembl.gene.id, ensembl.ann$GENEID), ]
                )
    # bm.annot fails to annotate some entries (why?)
    raw.df <- cbind(
                raw.df, 
                biomart.ann[ 
                  match(raw.df$ensembl.gene.id, biomart.ann$ensembl_gene_id),
                   ]
                )

    # and now save
    # save summary
    sink(paste(outDir, '/DESeq2/raw/raw_', out.base, '_summary.txt', sep=''), split=T)
    summary(raw)
    sink()
    # unnanotated results object
    saveRDS(raw, paste(outDir, "/DESeq2/raw/raw_", out.base, ".rds", sep=""))
    # annotated results as data frame (table)
    write.table(raw.df,
	    paste(outDir, "/DESeq2/raw/raw_", out.base, "_annotated.tab", sep=""),
	    row.names=T, col.names=T, sep='\t')
    # histogram plot
    out.png <- paste(folder, '/DESeq2/img/DESeq2_raw', out.base, '_hist.png', sep='')
    as.png( {
        margins <- par("mar")
        par(mar=c(5, 5, 5, 5))
        hist(raw.df$pvalue, main=out.base, breaks=1/alpha, xlab="p-value")
        par(mar=margins)
    }, out.png )


    # now we'll shrink the data to improve visualization and ranking
    shrunk.lfc <- lfcShrink(dds, contrast=c(column, x, y), type="ashr")
    if (save) {
        ofile <- paste(outDir, "/DESeq2/img/DESeq2_", out.base, "_raw+shrunk_MA.png", sep='')
        cat("    plotting", ofile, '\n')
    } else ofile <- NULL
    as.png( {
        dim <- par("mfrow")
        par(mfrow=c(1,2))
        # plotMA shows the log2 fold changes attributable to a given variable
        # over the mean of normalized counts for all the samples in the dataset
        # Points above alpha are colored, outliers are shown as directed triangles    
        plotMA(raw, alpha=alpha)
        # it is useful to look at the shrunk l2fc values, which removes the noise
        # associated with l2fc changes from low-count genes
        plotMA(shrunk.lfc, alpha=alpha)
        # after plotMA, one may identify interesting points interactively using
        # identify() and clicking on them:
        # idx <- identify(res$baseMean, res$log2FoldChange)
        # rownames(res[idx, ])
        #
        # Alternatively, looking at the plot and deciding which coordinates are
        # of interest, one can use, e.g. 
        # res_wt_vs_PC[ (res_wt_vs_PC$log2FoldChange) > 4) 
        #               & (res_wt_vs_PC$baseMean > 1000), ]
        # or
        # res_wt_vs_PC[  (abs(res_wt_vs_PC$log2FoldChange) > 4) 
        #              & (res_wt_vs_PC$baseMean > 1000), ]
        par(mfrow=dim)
    }, ofile )

    # prepare for ggplot:
    #	convert to data frame
    #	add rownames as two new columns GENEID and ensembl_gene_id
    #   add annotation from EnsDb using GENEID
    #   add biomaRt annotation using ensembl_gene_id
    #   rename some columns for plotting
    shrunk.lfc$ensembl_gene_id <- rownames(shrunk.lfc)
    ann.shrunk <- as.data.frame(shrunk.lfc) %>%
    rownames_to_column("GENEID") %>% 
    left_join(ensembl.ann, "GENEID") %>% 
    left_join(biomart.ann, "ensembl_gene_id") %>% 
    rename(logFC=log2FoldChange, FDR=padj)   

    # ggplot does not work inside a function, so this code seems useless
    if (FALSE) {
        if (save) {
            ofile <- paste(outDir, "/DESeq2/img/DESeq2_", out.base, "_l2FCshrunk.png", sep='')
            cat("    plotting", ofile, '\n')
        } else ofile <- NULL
        as.png( {
            ggplot(ann.shrunk, 
              aes(x = log2(baseMean), y=logFC),
              environment=environment()) + # this is supposed to make it work in a local env
                geom_point(aes(colour=FDR < alpha), shape=20, size=0.5) +
                geom_text(data=~top_n(.x, 10, wt=-FDR), aes(label=SYMBOL)) +
                labs(x="mean of normalised counts", y="log fold change")
        }, ofile )
    }

    if (save) {
        ofile <- paste(outDir, "/DESeq2/shrunk/shrunk_", out.base, ".rds", sep='')
        cat("    saving", ofile, '\n')
        saveRDS(shrunk.lfc, file=ofile)
        ofile <- paste(outDir, "/DESeq2/shrunk/shrunk_", out.base, "_annotated.tab", sep='')
        cat("    saving", ofile, '\n')
        write.table(ann.shrunk, 
                    file=ofile,
                    row.names=T, col.names=T, sep='\t')
    }

    # find statistically significant changes
    signif <- raw[ raw$padj < alpha, ]
    signif$abs_lfc <- abs(signif$log2FoldChange)
    if (save) {
        ofile <- paste(folder, '/DESeq2/signif/signif_', out.base, "_α<", alpha, ".tab", sep='')
        cat("    saving", ofile, '\n')
        write.table(signif, 
                file=ofile, 
                sep='\t', row.names=T, col.names=T)
    }

    # order the data, firstly by the decreasing absolute 
    # value of log2FoldChange and secondly by the increasing pvalue....
    srt <- signif[ order(signif$abs_lfc,
                        signif$padj,
                        decreasing=c(T, F)), ]
    if (save) {
        ofile <- paste(folder, "/DESeq2/signif/signif_sorted_", out.base, ".tab", sep='')
        cat("    saving", ofile,'\n')
        write.table(srt, 
                file=ofile, 
	        sep='\t', row.names=T, col.names=T)
    }
    # annotate the sorted significant data
    srt.df <- data.frame(ensembl.gene.id=rownames(srt), srt)
    srt.df <- cbind(srt.df, 
                ensembl.ann[ match(srt.df$ensembl.gene.id, ensembl.ann$GENEID), ])
    srt.df <- cbind(srt.df, 
                biomart.ann[ match(srt.df$ensembl.gene.id, biomart.ann$ensembl_gene_id), ])
    if (save) {
        ofile <- paste(folder, "/DESeq2/signif/signif_sorted_", out.base, "_annotated.tab", sep='')
        cat("    saving", ofile,'\n')
        write.table(srt.df, 
                file=ofile, 
	        sep='\t', row.names=T, col.names=T)
    }


    return(list(result=raw, 
                shrunk=shrunk.lfc, 
                shrunk.annot=ann.shrunk, 
                signif=srt, 
                signif.annot=srt.df))
}



ds2.interactively.print.n.most.significant.de.genes <- function(ds.data, n=10) 
{
    continue.on.enter("You may want to maximize your terminal before continuing ")

    options(width=200)
    n <- 10
    for (i in names(ds.data)) {
        cat("\nMost significant', n, 'genes for", i, '\n')
        print(head(ds.data[[i]]$signif.annot[ , c("log2FoldChange", "entrezgene_description")] ), n)
        continue.on.enter("Press [ENTER] to continue ")
    }
    options(width=80)
    continue.on.enter("Done, you can restore your terminal now ")
}

ds2.interactively.plot.top.up.down.regulated.gene <- function(dds,
                                                              ds.data, 
                                                              design.column)
{
    # Plot counts of the gene with maximal l2FC (up or down)
    threshold <- 0
    for (n in names(ds.data)) {
        res <- ds.data[[n]]$signif.annot
        up <- res[ res$log2FoldChange > threshold, ]	# up regulated
        down <- res[ res$log2FoldChange < -threshold, ]	# down regulated

        most.up <- rownames(up)[which.max(up$log2FoldChange)]	# max up
        cat("\ncounts of most overexpressed gene in", n, ":\n",
            most.up, 
            up$log2FoldChange[which.max(up$log2FoldChange)], 
            up$entrezgene_description[which.max(up$log2FoldChange)], 
            "\n")
        #plotCounts(dds, gene=rownames(res)[which.max(res$log2FoldChange)], intgroup="PFU")
        print(plotCounts(dds, gene=most.up, intgroup=design.column))
        k <- continue.on.key()
        if (k == "q") break

        most.down <- rownames(down)[which.min(down$log2FoldChange)]	# min down
        cat("\ncounts of most underexpressed gene in", n, ":\n",
            most.down, 
            down$log2FoldChange[which.min(down$log2FoldChange)], 
            down$entrezgene_description[which.min(down$log2FoldChange)],
            "\n")
        #plotCounts(dds, gene=rownames(res)[which.min(res$log2FoldChange)], intgroup="PFU")
        print(plotCounts(dds, gene=most.down, intgroup=design.column))
        continue.on.key()
        if (k == "q") break
    }
    cat("\n")

}


ds2.analyze.go.representation <- function(ds.data, bm.go.annot, folder) 
{ 
    # Analyze GO representation
    # -------------------------
    #
    # we use the biomaRt database from above
    #	mart.db might have been already assigned above, but not
    #	necessarily
    #ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
    #mart.db <- useMart("ensembl", mart.name)
    #listAttributes(mart.db)

    for (n in names(ds.data)) {
        cat("Processing GOs for", n, '\n')
        sig <- as.data.frame(ds.data[[n]]$signif)
        sig.a <- as.data.frame(ds.data[[n]]$signif.annot)
        sig$ensembl_gene_id <- rownames(sig)
        ourFilterType <- "ensembl_gene_id"
        filterValues <- sig$ensembl_gene_id
        # this will give us ALL GO annotations
        # we should use the saved searches from above to save bandwidth.
    #    gos <- getBM(attributes=c("ensembl_gene_id", 
    #                              "go_id", "name_1006", 
    #                              "definition_1006" ), 
    #                 mart=mart.db,
    #                 filters=ourFilterType,
    #                 values=filterValues)
        gos <- bm.go.annot[ bm.go.annot$ensembl_gene_id %in% filterValues,
                            c("ensembl_gene_id", 
                              "go_id", "name_1006", 
                              "definition_1006" ) ]
        sig.a <- merge(sig, gos, by="ensembl_gene_id")
        table(sig.a$go_id)	# this gives counts, but we'd like to multiply those
                                # counts by log2FC
        #sig.a[ , c(1, 3, 10)]
        # aggregate log2FC by go_id and sum the values
        print(
        aggregate(sig.a$log2FoldChange, 
                  by=list(Category=sig.a$go_id),
                  FUN=sum)
        )
        # or with formula interface
        gosums <- aggregate(log2FoldChange ~ go_id, sig.a, sum)
        gosums <- gosums[ order(gosums$log2FoldChange), ]
        # annotate gosums
        gosums.a <- cbind(gosums, bm.annot.1[ match(gosums$go_id, bm.annot.1$go_id), ])

        out.dir<- paste(folder, '/DESeq2/go', sep='')
        dir.create(out.dir, showWarning=FALSE)
        out.file <- paste(out.dir, '/GO_', n, "_sum.tab", sep='')
        out.file <- paste(folder, '/DESeq2/GO_', n, "_sum.tab", sep='')
        write.table(gosums.a, out.file, sep='\t')

        goavgs <- aggregate(log2FoldChange ~ go_id, sig.a, mean)
        goavgs <- gosums[ order(gosums$log2FoldChange), ]
        # annotate gosums
        goavgs.a <- cbind(goavgs, bm.annot.1[ match(goavgs$go_id, bm.annot.1$go_id), ])
        out.file <- paste(out.dir, '/GO_', n, "_average.tab", sep='')
        write.table(goavgs.a, out.file, sep='\t')
    }

}


annotate.go.ancestry <- function (annot.data, l2fc.threshold, outDir) {
    res <- annot.data

    goBPanc <- as.list(GOBPANCESTOR)
    # remove GO terms that do not have any ancestor
    goBPanc <- goBPanc[ ! is.na(goBPanc) ]

    goCCanc <- as.list(GOCCANCESTOR)
    # remove GO terms that do not have any ancestor
    goCCanc <- goCCanc[ ! is.na(goCCanc) ]

    goMFanc <- as.list(GOMFANCESTOR)
    # remove GO terms that do not have any ancestor
    goMFanc <- goMFanc[ ! is.na(goMFanc) ]

    res <- res[ abs(res$log2FoldChange) > l2fc.threshold, ]
    go.l2fc <- res[ , c("ensembl_gene_id", "log2FoldChange", "GENENAME", "go_id") ]

    for (i in 1:nrow(res)) {
        if ((i %% 100) == 0) cat(".")
        if (is.null(res[i, "go_id"])) next
        if (is.na(res[i, "go_id"])) next
        if (res[i, "go_id"] == '') next
        anc <- goBPanc[[ res[i, "go_id"] ]]
        if (! is.null(anc) & ! length(anc) == 0) {
            if (! is.na(anc) & anc != '') {
            # add  these GoIds to the table of gene/l2FC/GoIds 
                for (anc.goid in anc) {
                    go.l2fc <- rbind(go.l2fc, data.frame(
                          ensembl_gene_id=res[i, "ensembl.gene.id"],
                          log2FoldChange=res[i, "log2FoldChange"],
                          GENENAME=res[i, "GENENAME"],
                          go_id=anc.goid))
                }
            }
        }
        anc <- goCCanc[[ res[i, "go_id"] ]]
        if (! is.null(anc) ) {
            if (! is.na(anc) & anc != '') {
            # add  these GoIds to the table of gene/l2FC/GoIds 
                for (anc.goid in anc) {
                    go.l2fc <- rbind(go.l2fc, data.frame(
                          ensembl_gene_id=res[i, "ensembl.gene.id"],
                          log2FoldChange=res[i, "log2FoldChange"],
                          GENENAME=res[i, "GENENAME"],
                          go_id=anc.goid))
                }
            }
        }
        anc <- goMFanc[[ res[i, "go_id"] ]]
        if (! is.null(anc) ) {
            if (! is.na(anc) & anc != '') {
            # add  these GoIds to the table of gene/l2FC/GoIds 
                for (anc.goid in anc) {
                    go.l2fc <- rbind(go.l2fc, data.frame(
                          ensembl_gene_id=res[i, "ensembl.gene.id"],
                          log2FoldChange=res[i, "log2FoldChange"],
                          GENENAME=res[i, "GENENAME"],
                          go_id=anc.goid))
                }
            }
        }
    }
    # aggregate data by GOID (will result in c(Category, x) columns
    go.l2fc.sum <- aggregate(go.l2fc$log2FoldChange, 
               by=list(Category=go.l2fc$go_id), FUN=sum)
    go.l2fc.avg <- aggregate(go.l2fc$log2FoldChange, 
               by=list(Category=go.l2fc$go_id), FUN=mean)

    # rename Category, x to go_id, l2fc
    colnames(go.l2fc.sum) <- c('go_id', 'sum.l2fc')
    colnames(go.l2fc.avg) <- c('go_id', 'avg.l2fc')


    # annotate
    godesc <- AnnotationDbi::select(GO.db, keys=go.l2fc$go_id, 
                  columns= c("GOID", "TERM", "DEFINITION", "ONTOLOGY"), 
                  keytype= "GOID")
    go.l2fc <- cbind(go.l2fc, godesc)
    godesc <- AnnotationDbi::select(GO.db, keys=go.l2fc.sum$go_id, 
                  columns= c("GOID", "TERM", "DEFINITION", "ONTOLOGY"), 
                  keytype= "GOID")
    go.l2fc.sum <- cbind(go.l2fc.sum, godesc)
    godesc <- AnnotationDbi::select(GO.db, keys=go.l2fc.avg$go_id, 
                  columns= c("GOID", "TERM", "DEFINITION", "ONTOLOGY"), 
                  keytype= "GOID")
    go.l2fc.avg <- cbind(go.l2fc.avg, godesc)


    # sort by l2FC and save
    go.l2fc <- go.l2fc[ order(go.l2fc$log2FoldChange, decreasing=T), ]
    go.l2fc.sum <- go.l2fc.sum[ order(go.l2fc.sum$sum.l2fc, decreasing=T), ]
    go.l2fc.avg <- go.l2fc.avg[ order(go.l2fc.avg$avg.l2fc, decreasing=T), ]
    # sort by abs(log2FC)
    go.l2fc.sum.abs <- go.l2fc.sum[ order(abs(go.l2fc.sum$sum.l2fc), decreasing=T), ]
    go.l2fc.avg.abs <- go.l2fc.avg[ order(abs(go.l2fc.avg$avg.l2fc), decreasing=T), ]

    # save
    out.file <- paste(outDir, '/GOANC_', n, ".tab", sep='')
    write.table(go.l2fc, out.file, sep='\t')
    out.file <- paste(outDir, '/GOANC_', n, "_sum.tab", sep='')
    write.table(go.l2fc.sum, out.file, sep='\t')
    out.file <- paste(outDir, '/GOANC_', n, "_average.tab", sep='')
    write.table(go.l2fc.avg, out.file, sep='\t')
    out.file <- paste(outDir, '/GOANC_', n, "_sum_abs.tab", sep='')
    write.table(go.l2fc.sum.abs, out.file, sep='\t')
    out.file <- paste(outDir, '/GOANC_', n, "_average_abs.tab", sep='')
    write.table(go.l2fc.avg.abs, out.file, sep='\t')

}



ds2.analyze.go.ancestry <- function(ds.data, l2fc.threshold=0, folder)
{
    # Reannotate GO with ancestry
    # library(GO.db)

    l2fc.threshold <- 0

    out.dir <- paste(folder, '/Deseq2/goanc')
    dir.create(out.dir, showWarning=FALSE)

    for (n in names(ds.data)) {
        cat("\nTracing GO ancestry for", n, '\n')
        res <- ds.data[[n]]$signif.annot

        annotate.go.ancestry(res, l2fc.threshold, out.dir)

        #ans <- continue.on.enter("Press RETURN to continue: ")
        #if (ans == "q") break
    }
}



ds2.analyze.pfam.representation <- function(ds.data, bm.fam.annot, folder)
{
   # Analyze PFAM representation
    # ---------------------------
    #
    # for PFAM, we can use
    #
    # Get PFAM database and use AC -> DE mapping
    #library(PFAM.db)
    db <- PFAMDE
    mk <- PFAMDE[mappedkeys(PFAMDE)]
    xx <- as.list(mk)
    pfam.table <- toTable(PFAMDE)

    out.dir <- paste(folder, '/DESeq2/pfam', sep='')
    dir.create(out.dir, showWarning=FALSE)

    #for (i in pfam.ids) print(xx[[i]])

    # Get PFAM families and sort them by their representation in the dataset
    options(width=200)
    for (n in names(ds.data)) {
        sig <- as.data.frame(ds.data[[n]]$signif) 
        names(sig)
        sig$ensembl_gene_id <- rownames(sig)
        sig.a <- cbind(sig, 
               bm.fam.annot[ match(sig$ensembl_gene_id, bm.fam.annot$ensembl_gene_id), ])
        pfam.ids <- sig.a$pfam[ ! is.na(sig.a$pfam) ]
        aggregate(sig.a$log2FoldChange, by=list(Category=sig.a$pfam), FUN=sum)

        sig.a <- merge(sig, bm.fam.annot, by="ensembl_gene_id")
        head(sig.a, 10)

        sig.sum <- aggregate(sig.a$log2FoldChange, by=list(Category=sig.a$pfam), FUN=sum)
        sig.sum <- sig.sum[ order(sig.sum$x, decreasing=T), ]
        sig.sum <- sig.sum[ ! is.na(sig.sum$Category), ]
        sig.sum <- sig.sum[ sig.sum$Category != '', ]
        names(sig.sum) <- c("pfam", "sum.l2FC")
        # annotate
        for (i in 1:length(sig.sum$pfam)) { 
            pf <- sig.sum$pfam[i] ; 
            if (! is.null(xx[[pf]])) { 
                sig.sum$pfam.de[i] <- xx[[pf]] 
            } else { 
                sig.sum$pfam.de[i] <- '' 
            }
        }
        cat("\nMost represented PFAM families in", n, "\n")
        print(head(sig.sum, 20)) ; print(tail(sig.sum, 20))
        out.file <- paste(folder, '/DESeq2/PFAM_', n, "_sum.tab", sep='')
        write.table(sig.sum, out.file, sep='\t')

        sig.avg <- aggregate(sig.a$log2FoldChange, by=list(Category=sig.a$pfam), FUN=mean)
        sig.avg <- sig.avg[ order(sig.avg$x, decreasing=T), ]
        sig.avg <- sig.avg[ ! is.na(sig.avg$Category), ]
        sig.avg <- sig.avg[ sig.avg$Category != '', ]
        names(sig.avg) < - c("pfam", "avg.l2FC")
        # annotate
        for (i in 1:length(sig.avg$pfam)) { 
            pf <- sig.avg$pfam[i] ; 
            if (! is.null(xx[[pf]])) { 
                sig.avg$pfam.de[i] <- xx[[pf]] 
            } else { 
                sig.avg$pfam.de[i] <- '' 
            }
        }
        out.file <- paste(folder, '/DESeq2/PFAM_', n, "_average.tab", sep='')
        write.table(sig.avg, out.file, sep='\t')

        # save also the files sorted by abs(l2fc)
        sig.sum.abs <- sig.sum[ order(abs(sig.sum$sum.l2FC), decreasing=T), ]
        sig.avg.abs <- sig.avg[ order(abs(sig.avg$avg.l2FC), decreasing=T), ]
        out.file <- paste(folder, '/DESeq2/PFAM_', n, "_sum_abs.tab", sep='')
        write.table(sig.sum.abs, out.file, sep='\t')
        out.file <- paste(folder, '/DESeq2/PFAM_', n, "_average_abs.tab", sep='')
        write.table(sig.avg.abs, out.file, sep='\t')

        if (INTERACTIVE) {
            #continue.on.key()
            ans <- continue.on.enter(prompt="Press return to continue ")
            if (ans == "q") break
        }
    }
    cat('\n')
    options(width=80)
}




GO_fgsea <- function (ann.shrunk.lfc, 
		      max.size=250,
                      out.dir=paste(folder, 'go_fgsea', sep='/'), 
                      out.name='GO_fgsea',
                      use.description=TRUE,
                      top.n=20,
                      top.biblio=5,
                      verbose=FALSE) {
    # Do GSEA on GO terms using fgsea

    # Rank all genes on their fold change.
    #	Here we exclude genes for which we have no EntrezID and
    #	use shrunk LFC values
#    gseaDat <- filter(ann.shrunk.lfc, !is.na(ENTREZID))
    gseaDat <- filter(ann.shrunk.lfc, !is.na(GENEID))

    ranks <- gseaDat$lfc
    #names(ranks) <- gseaDat$ENTREZID
    names(ranks) <- gseaDat$GENEID
    head(ranks)
    #uranks <- ranks[!duplicated(sort(names(ranks)))]

    # plot all the ranked fold changes
    out.png <- paste(out.dir, '/', out.name, '.barplot.png', sep='')
    as.png(barplot(sort(ranks, decreasing=T)), out.png)
    
    # load pathways
    #pathways <- ann.shrunk.lfc$ENTREZID
    #pathways <- ann.shrunk.lfc$go_id
    #names(pathways) <- gseaDat$ENTREZID
    #names(pathways) <- gseaDat$GENEID
    #upathways <- pathways[!duplicated(sort(names(pathway)))]
    
    # create a list of go_id terms where each term contains a vector
    # of emsembl_gene_id in that term
    #   first recover the annotation (in case we are re-run and do not
    #   have it
    if ( ! exists(substitute(bm.go.annot)) ) {
        bm.go.annot <- read.table(paste(folder, '/biomart.go.tab', sep=''), 
	            header=T, sep='\t')
    }
    #if ( ! exists(substitute(bm.goslim.annot)) ) {
    #    bm.goslim.annot <- read.table(paste(folder, '/biomart.goslim.tab', sep=''), 
    #	            header=T, sep='\t')
    #}

    if (use.description == FALSE) {
        # split() will divide the gene-ids by their go-id
        pathways <- split(bm.go.annot$ensembl_gene_id, bm.go.annot$go_id)
    } else {
        # Do the same but with large gene ontology names
        # split() will divide the gene-ids by their go-name
        pathways.go <- split(bm.go.annot$ensembl_gene_id, bm.go.annot$name_1006)
    }
    # do analysis
    # the resulting table contains enrichment scores and p-values
    out.file <- sprintf("%s/go_fgsea_10-%d.RData", out.dir, max.size)
    if (file.exists(out.file)) {
        load(file=out.file)
    } else {
        fgseaRes <- fgsea(pathways=pathways.go, 
                          stats=ranks, 
		          minSize=10, 
		          maxSize=max.size, 
                          nPermSimple=100000 
                          )
        save(fgseaRes, file=out.file)
    }
    if (verbose == TRUE)
        head(fgseaRes[order(padj, -abs(NES)), ], n=10)

    if (verbose == TRUE) {
        # plot enrichment score
        sorted.fgsea.res <- fgseaRes[order(padj, -abs(NES)), ]
        sfr.names <- sorted.fgsea.res$pathway
        for (i in 1:top.n) {
            if (use.description == FALSE)
                descr <- bm.go.annot[bm.go.annot$go_id == sfr.names[i], "name_1006"]
            else
                descr <- sfr.names[i]
            print(
                plotEnrichment(pathways.go[[ sfr.names[i] ]], ranks) +
                    labs(title=descr)
                )
            ans <- readline("Press RETURN to continue: ")
            if (ans == "q") break
        }
    }
    # gsea table plot of top.n gene families
    #   top_n() is now deprecated
    topUp <- fgseaRes %>%
        filter(ES > 0) %>%
        top_n(top.n, wt=-padj)

    topDown <- fgseaRes %>%
        filter(ES < 0) %>%
        top_n(top.n, wt=-padj)

    topPathways <- bind_rows(topUp, topDown) %>%
        arrange(-ES)
    # last resort when "pos" is used    
    #topPathways <- sorted.fgsea.res[1:top.n, ]

    # do the plots and save descriptions
    out.file <- paste(out.dir, "/topUp.", top.n, '.txt', sep='')
    for (i in 1:top.n) {
        name <- topUp[i]$pathway
        if (use.description == FALSE)
            descr <- unique(bm.go.annot[bm.go.annot$go_id == name, "name_1006"])
        else
            descr <- name
        cat(i, descr, file=out.file, '\n', sep='\t', append=TRUE)
        #print(name)
        #out.png <- paste(out.dir, '/', out.name, '.topUp', i, '.enrichment.png', sep='')
        out.png <- sprintf("%s/%s.topUp.%03d.enrichment.png", 
        	out.dir, out.name, i)
        out.plt <-
            plotEnrichment(pathways.go[[ name ]] , ranks) +
            labs(title=descr)
        ggsave(out.png, out.plt)
    }
    out.file <- paste(out.dir, "/topDn.", top.n, '.txt', sep='')
    for (i in 1:top.n) {
        name <- topDown[i]$pathway
        if (use.description == FALSE)
            descr <- unique(bm.go.annot[bm.go.annot$go_id == name, "name_1006"])
        else
            descr <- name
        cat(i, descr, file=out.file, '\n', sep='\t', append=TRUE)
        #print(name)
        #out.png <- paste(out.dir, '/', out.name, '.topDown', i, '.enrichment.png', sep='')
        out.png <- sprintf("%s/%s.topDn.%03d.enrichment.png", 
        	out.dir, out.name, i)
        out.plt <-
            plotEnrichment(pathways.go[[ name ]] , ranks) +
            labs(title=descr)
        ggsave(out.png, out.plt)
    }

    out.png <- paste(out.dir, '/', out.name, '.GSEAtable.png', sep='')
    as.png(
        plotGseaTable(pathways.go[topPathways$pathway], 
                      ranks, 
                      fgseaRes, 
                      gseaParam = 0.5)
    , out.png, width=1024, height=100*top.n, overwrite=TRUE)
    
    
    # and now do a plot of the interest in citations during
    # the last ten years for the top 5 sets
    out.png <- paste(out.dir, '/', out.name, '.EUPMC.png', sep='')
    cur.year <- as.integer(format(Sys.Date(), "%Y"))
    terms <- topDown[1:top.biblio]$pathway
    as.png(
        pmcplot(terms, (cur.year-10):cur.year, proportion=FALSE),
        out.png)

    
    # finally, return fgseaRes
    return(fgseaRes)
}



GO_KEGG_clusterProfiler <- function(ann.shrunk.lfc, 
		      max.size=250,
                      out.dir=paste(folder, 'go_cProf', sep='/'), 
                      out.name='GO_cProf',
                      use.description=TRUE,
                      OrgDb = org.Cjaponica.eg.db,
                      kegg_organism = "cjo",	# (cjo = coturnix japonica)
                                             # (gga = gallus gallus)
                      top.n=10,
                      top.biblio=5,
                      verbose=FALSE) {
        
    gseaDat <- filter(ann.shrunk.lfc, !is.na(ENTREZID))
    ranks <- gseaDat$lfc
    names(ranks) <- gseaDat$ENTREZID
    ranks<-na.omit(ranks)

    s.ranks <- sort(ranks, decreasing=T)
    
    ### G O annotation
    #
    gse.out.file <- paste(out.dir, '/', out.name, '.topGO.Rdata', sep='')
    if (VERBOSE) cat("Doing GO GSEA with Cluster profiler\n")
    if ( ! file.exists(gse.out.file) ) {
        gse <- gseGO(geneList=s.ranks, 
                     ont ="ALL", 
                     keyType = "ENTREZID", 
                     minGSSize = 3, 
                     maxGSSize = max.size, 
                     pvalueCutoff = 0.05, 
                     verbose = TRUE, 
                     OrgDb = OrgDb, 
                     pAdjustMethod = "fdr")
                     #pAdjustMethod = "none")
        if (dim(gse)[1] == 0) {
            # p.adjusting may have failed to produce any result:
            # Try without correction issuing a warning
            gse <- gseGO(geneList=s.ranks, 
                     ont ="ALL", 
                     keyType = "ENTREZID", 
                     minGSSize = 3, 
                     maxGSSize = max.size, 
                     pvalueCutoff = 0.05, 
                     verbose = TRUE, 
                     OrgDb = OrgDb, 
                     #pAdjustMethod = "fdr")
                     pAdjustMethod = "none")
            out.name <- paste(out.name, '.raw_p', sep='')
        }
        if (dim(gse)[1] > 1) {
            saveRDS(gse, file=gse.out.file)
            gse.tab.file <- paste(out.dir, '/', out.name, '.topGO.tab', sep='')
            write.table(gse, file=gse.tab.file, row.names=T, col.names=T, sep='\t')
        }
    } else {
        if (file.exists(gse.out.file)) {
            # default is adjusted p
            gse <- readRDS(gse.out.file)
        } else {
            # try raw p
            out.name <- paste(out.name, '.raw_p', sep='')
            gse.out.file <- paste(out.dir, '/', out.name, '.topGO.Rdata', sep='')
            if (file.exists(gse.out.file)) {
                gse <- readRDS(gse.out.file)
            } else {
                # use an empty table so next check can be done
                gse <- table(c())
            }
        }
    }
    # if dim(gse)[1] == 0 then we won't save anything
    #	hopefully, if re-run again it might work next time by chance
    # else
    if (dim(gse)[1] > 0) {
	if (VERBOSE) cat("Plotting GO GSEA with Cluster profiler\n")
        # DO GENERIC PLOTS
        # **Dotplot**: For each group shows if it is up or down 
        # regulated and to which extent. The circle size is proportional to the
        # size (the number of genes contained) of the group, and the color to the
        # p-value.
        out.png <- paste(out.dir, '/', out.name, '.GOdotplot.png', sep='')
        as.png(
            dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign)
            , out.png)

        # **Enrichment map**: organizes terms in a network with edges 
        # connecting overlapping gene sets (i.e. shows which gene sets contain common
        # genes). Dot diameter represents the size of the gene set (the number of genes
        # it contains) and color represents the adjusted probability.
        out.png <- paste(out.dir, '/', out.name, '.GOemmaplot.png', sep='')
        as.png(
            emapplot(pairwise_termsim(gse), showCategory = top.n)
            , out.png)

        # **Ridgeplot**: density plots grouped by gene set depicting
        # the frequency of fold change values per gene within each set. Helps 
        # interpret up/down-regulated pathways.
        out.png <- paste(out.dir, '/', out.name, '.GOridgeplot.png', sep='')
        as.png(
            ridgeplot(gse) + labs(x = "enrichment distribution")
            , out.png)

        # **PubMed trend** of enriched terms
        # Plots the number/proportion of publications trend based on 
        # the query result from PubMed Central.
        out.png <- paste(out.dir, '/', out.name, '.GOpmcplot.png', sep='')
        cur.year <- as.integer(format(Sys.Date(), "%Y"))
        terms <- gse$Description[1:top.biblio]
        as.png(
            pmcplot(terms, (cur.year-10):(cur.year-1), proportion=FALSE)
            , out.png)
        
        # DO DETAILED PLOTS
        for (i in 1:dim(gse)[1] ) {
            # **GSEA plot**: __Plot of the Running Enrichment Score__ (green
            # line) for a gene set as the analysis walks down the ranked gene list,
            # including the location of the maximum enrichment score (the red line).
            # The black lines in the Running Enrichment Score show where the members of the
            # gene set appear in the ranked list of genes, indicating the leading edge
            # subset.
            # 
            # __Ranked list metric__ shows the value of the ranking
            # metric (log2 fold change) as you move down the list of ranked genes. The
            # ranking metric measures a gene’s correlation with a phenotype.

            # Use the `Gene Set` param for the index in the title, and as the value for geneSetId
            out.png <- paste(out.dir, '/', 
                    out.name, '.GOgseaplot.', i, '.', gse$ID[i], '.png', sep='')
            as.png(
                gseaplot(gse, by = "all", title = gse$Description[i], geneSetID = i)
                , out.png)

            # **Category Netplot**: shows the
            # linkage between genes and biological concepts as a network (helpful to see
            # which genes are involved in enriched pathways and genes that may belong to
            # multiple annotation categories).
            out.png <- paste(out.dir, '/', 
                    out.name, '.GOcnetplot.', i, '.', gse$ID[i], '.png', sep='')
            # categorySize can be either 'pvalue' or 'geneNum'
            as.png(
                cnetplot(gse, categorySize="pvalue", foldChange=s.ranks, 
                         showCategory=i)
                , out.png)
        }
        #out.file <- paste(out.dir, '/', out.name, '.topGO.tab', sep='')
        out.file <- gse.out.file
        write.table(gse, file=gse.out.file, row.names=T, col.names=T, sep='\t')
    }
    
    ### K E G G annotation
    # let's try with KEGG (from the ENTREZID which is the same a ncbi-genid)
    
    if (VERBOSE) cat("Doing KEGG GSEA with ClusterProfiler\n")
    kse.out.file <- paste(out.dir, '/', out.name, '.topKEGG.tab', sep='')
    if (! file.exists(kse.out.file)) {
        kse <- gseKEGG(geneList     = s.ranks,
                   organism     = kegg_organism,
                   #nPerm        = 10000,
                   minGSSize    = 3,
                   maxGSSize    = max.size,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "fdr",
                   keyType       = "ncbi-geneid",
                   nPermSimple = 100000)

        if (dim(kse)[1] == 0) {
            # Try without correction issuing a warning
            kse <- gseKEGG(geneList     = s.ranks,
                       organism     = kegg_organism,
                       #nPerm        = 10000,
                       minGSSize    = 3,
                       maxGSSize    = max.size,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "none",
                       keyType       = "ncbi-geneid",
                       nPermSimple = 100000)
            out.name <- paste(out.name, 'raw_p', sep='')
        }
        if (dim(kse)[1] > 1) {
            saveRDS(kse, file=kse.out.file)
            kse.tab.file <- paste(out.dir, '/', out.name, '.topGO.tab', sep='')
            # save also as text for user access
            write.table(kse, file=kse.tab.file, row.names=T, col.names=T, sep='\t')
        }
    } else {
        if (file.exists(kse.out.file)) {
            # default is adjusted p
            kse <- readRDS(kse.out.file)
        } else {
            # try raw p
            out.name <- paste(out.name, '.raw_p', sep='')
            kse.out.file <- paste(out.dir, '/', out.name, '.topGO.Rdata', sep='')
            if (file.exists(kse.out.file)) {
                kse <- readRDS(kse.out.file)
            } else {
                # use an empty table so next check can be done
                kse <- table(c())
            }
        }
    }
    
    # if dim(kse)[1] == 0 then we won't save anything
    #	hopefully, if re-run again it might work next time by chance
    # else
    if (dim(kse)[1] > 0) {
        if (VERBOSE) cat("Plotting KEGG GSEA with ClusterProfiler\n")
        # <p>A <strong>Dotplot</strong>: For each group shows if it is up or down 
        # regulated and to which extent. The circle size is proportional to the
        # size (the number of genes contained) of the group, and the color to the
        # p-value.</p>
        out.png <- paste(out.dir, '/', out.name, '.KEGGdotplot.png', sep='')
        as.png( 
            dotplot(kse, showCategory = 10, title = "Enriched Pathways" , 
                    split=".sign") + facet_grid(.~.sign)
            , out.png)

        # <p>B <strong>Enrichment map</strong>: organizes terms in a network with edges 
        # connecting overlapping gene sets (i.e. shows which gene sets contain common
        # genes). Dot diameter represents the size of the gene set (the number of genes
        # it contains) and color represents the adjusted probability.</p>
        out.png <- paste(out.dir, '/', out.name, '.KEGGemmaplot.png', sep='')
        as.png(
            emapplot(pairwise_termsim(kse))
            , out.png)


        # <p>D <strong>Category Netplot</strong>: shows the
        # linkage between genes and biological concepts as a network (helpful to see
        # which genes are involved in enriched pathways and genes that may belong to
        # multiple annotation categories).</p>
        # categorySize can be either 'pvalue' or 'geneNum'
        out.png <- paste(out.dir, '/', out.name, '.KEGGcnetplot.png', sep='')
        as.png(
            cnetplot(kse, categorySize="pvalue", foldChange=s.ranks)
            , out.png)

        # <p>C <strong>Ridgeplot</strong>: density plots grouped by gene set depicting
        # the frequency of fold change values per gene within each set. Helps 
        # interpret up/down-regulated pathways.</p>
        out.png <- paste(out.dir, '/', out.name, '.KEGGridgeplot.png', sep='')
        as.png(
            ridgeplot(kse) + labs(x = "enrichment distribution")
            , out.png)

        # PubMed trend of enriched terms
        # 
        # Plots the number/proportion of publications trend based on the query result
        # from PubMed Central.
        out.png <- paste(out.dir, '/', out.name, '.KEGGpmcplot.png', sep='')
        cur.year <- as.integer(format(Sys.Date(), "%Y"))
	    as.png(
            pmcplot(kse$description[1:5], (cur.year-10):(cur.year-1), proportion=FALSE)
            , out.png)

        cur.dir <- getwd()
        for (i in 1:dim(kse)[1]) {
            # for each of the pathways in kse
            # Use the `Gene Set` param for the index in the title, and as the value for geneSetId
            
            # GSEA plot 
            # Plot of the Running Enrichment Score (green
            # line) for a gene set as the analysis walks down the ranked 
            # gene list, including the location of the maximum enrichment 
            # score (the red line). The black lines in the Running Enrichment 
            # Score show where the members of the gene set appear in the 
            # ranked list of genes, indicating the leading edge subset.
            #
            # The Ranked list metric shows the value of
            #  the ranking metric (log2 fold change) as you move down the 
            # list of ranked genes. The ranking metric measures a gene’s 
            # correlation with a phenotype.
            out.png <- paste(out.dir, '/', 
                out.name, '.KEGGgseaplot.', i, '.', kse$ID[i], '.png', sep='')
            as.png(
              gseaplot(kse, by = "all", title = kse$Description[i], geneSetID = i)
              , out.png)
            
            # **Pathview**
            # This will create a PNG and a __different__ PDF of the enriched 
            # KEGG pathway in the current working directory.
            setwd(out.dir)	# change to the appropriate directory
            # Produce the native KEGG plot (PNG)
            dme <- pathview(gene.data=s.ranks, pathway.id=kse$ID[i], species = kegg_organism)

            # Produce a different plot (PDF) (different from the previous one)
            dme <- pathview(gene.data=s.ranks, pathway.id=kse$ID[i], species = kegg_organism, kegg.native = F)
            setwd(cur.dir) # and back
        }
    }
}

