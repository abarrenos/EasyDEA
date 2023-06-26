#!/usr/bin/env Rscript
#
# We start from data that has already been cleaned by performing a quality
# check with FastQC and subsequent edge trimming.
#
# The next step is to align the reads and calculate the counts of reads
# that map to each gene. This has been done as well previously, obtaining
# alignments in BAM format: the .bam files (binary files containing the
# reads aligned to the reference genome sequence) and the .bam.bai files
# containing the indexes for each one). This step is mostly a matter of
# time and has already been done:
# 
#   Paired-end Illumina short-reads were aligned against Coturnix japonica
#   genome (v2.0 primary assembly) using RNA-STAR (1) (--outReadsUnmapped
#   Fastx; --alignIntronMax 10000; -- alignMatesGapMax 10000). PCR and optical
#   duplicates were marked using the MarkDuplicates function of Picard-Tools 
#   (GATK)(2) (TAGGING_POLICY=ALL). Alignment results, saved as BAM files, 
#   were sorted and indexed using samtools (3).
# 
# We take the alignments produced by the Bioinformatics for Proteomics and
# Genomics Service of CNB directly, which are stored in the folder
# 'Alignments_Coturnix', and the reference genome used by them which is in
# the folder 'refGenomes/Coturnix_Japonica'. The reference files correspond
# to the 2.0 primary assembly. This is important for we will need their
# indexes for the next step.
#
# Now we need to calculate the counts of the reads that map to gene regions,
# summing up the reads that match each gene. For this, we use each .bam file
# and process it with the function featureCounts from R package Rsubread, and
# the reference genome data, which is in file Cjaponica.gtf. This last file
# contains the information about the features annotated in the genome of
# Coturnix japonica, including the coordinates of each feature in the reference
# genome. The function featureCounts will use these coordinates to know to
# which feature each read maps.
#


# we need this function here so we can include additional files
# until we make this into a package
getScriptPath <- function()
{
     
    # this works if we were called with 'source()'
    src.path <- getSrcFilename(function() {}, full.names=T)
    if (! is.null(src.path) && src.path != '')
        return(src.path)

    # this works for Rscript, it may match more than one path
    # if Rscript was used with several --file= arguments, in
    # which case we will return *only* the first one
    cmd.args <- commandArgs()
    m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
    script.path <- regmatches(cmd.args, m)
    if (length(script.path) >= 1) return(script.path[1])
    if (length(script.path) == 1) return(script.path)	# may return multiple matches

    # this works for 'R -f', it will return only the first -f argument
    for (i in 1:length(cmd.args) ) {
        print(i); print(cmd.args[i])
        if (cmd.args[i] == '-f') return(cmd.args[i+1])
    }
    
    # if we arrive here, we didn't match anything, turn to last resort
    return (sys.frame(1)$ofile)
}

# get my location
mydir <- dirname(getScriptPath())
# inside my location there should be a 'lib' directory with the needed
# auxiliary scripts: we'll source them all
sourceDir(paste(mydir, "lib", sep='/'))


use.package(optparse)
use.package(ensembldb)		# needs to be first to avoid S4 inconsistencies
use.package(Rsubread)		# for read mapping
use.package(AnnotationHub)		# to seek annotation packages
use.package(AnnotationForge)	# to build our annotation package
use.package(GO.db)
use.package(PFAM.db)

use.package("biomaRt")		# an alternate approach to retrieve annotation

use.package(tibble)			# general tools
use.package(tidyr)
use.package(stringr)
use.package(dplyr)
use.package(readr)
use.package(keypress)

use.package(edgeR)			# RNAseq with edgeR
use.package(limma)
use.package(RColorBrewer)
use.package(gplots)

use.package(DESeq2)			# RNAseq with DESeq2
use.package(IHW)			# for p-value adjustment with IHW
use.package(ggplot2)

use.package(cluster) 
use.package(factoextra)
use.package(fpc)
use.package("NbClust")

use.package(tcltk)
use.package(gWidgets2)


options <- get.options()

# for convenience, we will assign options to specific names
#	we could get a similar efect if instead of a list, options
#	were a data.frame and then it would suffice to use attch()
#	but this makes it evident which variables correspond to 
# 	options
ALIGN <-                  options$ALIGN
BOTH <-					  options$BOTH
USE.ONLINE.ANNOTATION <-  options$USE.ONLINE.ANNOTATION
USE.EDGER <-              options$USE.EDGER
USE.DESEQ2 <-             options$USE.DESEQ2
reference <-              options$reference		# reference genom sequence base name
release <-                options$release		# release name
target.organism <-        options$target.organism
ens.version <-            options$ens.version	# version of ENSEMBL to use
mart.name <-              options$mart.name      # mart from BiomaRt to use
org.package <-            options$org.package
n.genes <-                options$n.genes		# number of top genes to revise
fastq.data <-             options$fastq.data     # directory with fastq files
alignment.dir <-          options$alignment.dir
rnaseq.out <-             options$rnaseq.out
my.name <-                options$my.name		# used to identify maintainer of
my.email <-               options$my.email		# any created annotation package
my.user <-                options$my.user        # used to access MySQL
my.password <-            options$my.password
metadata <-               options$metadata
cpm.threshold <-          options$cpm.threshold
significance.threshold <- options$significance.threshold
design.column <-          options$design.column
config.file <-            options$config.file
INTERACTIVE <-            options$INTERACTIVE
VERBOSE <-                options$VERBOSE

# convenience variables
by.rows=1
by.columns=2

##############################################################################
#
# DO THE WORK
#
##############################################################################


short.title('RNAseq')		# Print a visible title

folder <- create.output.hierarchy(rnaseq.out, use.both.reads=BOTH)

# next are for reproducibility

# make a copy of the metadata into the rnaseq folder
system(paste("cp", metadata, rnaseq.out))
system(paste("cp", metadata, folder))
# this could give problem if we are run on our installation folder for
# we would get an error and might confuse the user
#system(paste("cp -R", mydir, folder))

# save the options so we keep a record of how the analysis was generated
save_options(options, paste(folder, 'RNASEQ_OPTIONS.R', sep='/'))

# keep a log file for tracking and reporting
logfile <- paste(folder, "log", "RNAseq.log", sep='/')
log <- openLogFile(logfile)	# opens with sink(), may be closed
						    # specifically with closeLogFile(log)
                            
# when a script ends or when stop() is called, there are a number
# of housekeeping tasks to do. on.exit() allows us to add additional
# tasks so that they, too, are executed at the end of the script;
# adding this here we do not need to keep track of all the sink()s
on.exit(sink.titanic, add=T)





##############################################################################
#
#  PREPARE ANNOTATION SO IT IS READY WHEN WE GET THE RESULTS
#
##############################################################################
#---------------------------------------------------------------------------

# Now is time to prepare to connect all the results we get with the existing
# information  we know from the literature, so we will retrieve infos from
# ENSEMBL and connect with the genes  we have kept. We are interested only in
# genes and transcriptomes (non-characterized genes) from the organism
# Gallus galllus. there is no database available for this organism so we
# have to create it by ourselves

# we create two different databases one from the gtf file (small archive)
# and one from the data we extracted from ENSEMBL and stored it in a sqlite 
# file

# we need to load the package that enables us to connect a file as an
# external database ensembldb. After that we create a variable for the DB
# storing it in the memory simple process: we extract the files from ENSEMBL
# and store the data to MySQL account and then with the script
# 'generate-EnsDBs.R' we create a sqlite file that has all the information
# needed to continue

short.title('annotation (ensembl)')

# NOTE: to be used to store/search annotation from now on
annotation.dir <- file.path(folder, "annotation")

if (USE.ONLINE.ANNOTATION == TRUE) {
    # get access to annotation hub
    ah <- AnnotationHub()
    # get ENSEMBL annotations for our target
    qr <- query(ah, c("EnsDb", target.organism))
    if (length(qr > 0) {
        # there is at least one, choose the most recent (last) one
        #	NOTE: this might not be desired sometimes
        last.ref <- names(qr)[length(names(qr))]
        net.edb <- qr[[last.ref]]
    }

    # with this we do not need to build it from GTF/GFF/mysql,
    ens.db <- net.db

} else {    
    DBFile <- build.offline.annotation( 
                              annotation.dir,
                              db.dir='EnsDb.Ggallus.v106',
                              target.organism=target.organism,
                              reference.gtf=reference.gtf,
                              relese=release,
                              ens.version=ens.version,
                              user=my.user,
                              pass=my.password,
                              author=my.name,
                              maintainer=paste(my.name, my.email), 
                              license="Artistic-2.0")

    ens.db <- EnsDb(DBFile)
}

if (VERBOSE) {
    # check it
    print(ens.db)
    print(head(keys(ens.db, 'GENEID')))
    columns(ens.db)
}


###################################################
###   A L T E R N A T E   A N N O T A T I O N   ###
###################################################

# ---------------------------------------------------------------
# M A K E     O R G     P A C K A G E
# ---------------------------------------------------------------

short.title('annotation (Org.db)')

org.db <- NULL

if (USE.ONLINE.ANNOTATION == TRUE) {
    # use AnnotationHub to seek a suitable package
    #ah <- AnnotationHub()	# we have already open the connection
    # check if an Org.db is available for our organism)
    qo <- query(ah, "Orgdb")
    if ( last.ref %in% names(qo) ) {
        org.db <- query(ah, "Orgdb")[[ last.ref ]]
    } 
}

# if that didn't work, try offline
if (is.null(org.db)) {
     # then try to build off-line annotation
     if ( ! require(org.package, character.only=T)) {
        #-----------------------------------------------
        # Create Org Package
        #-----------------------------------------------
        #
        # The datasets were first downloaded by hand from NCBI
        # and SwissProt into org.Gg.eg.db.
        # Then the following command had to be used:
        #
        makeOrgPackageFromNCBI(
                author  = "J. R. Valverde <jrvalverde@cnb.csic.es>", 
                maintainer = "J. R. Valverde <jrvalverde@cnb.csic.es>", 
                tax_id  = "9031", # from NCBI Taxonomy browser (ncbi:txid9031)
                genus   = "Gallus", 
                species = "gallus", 
                version = "6a", 
                #tax_id = "93934", 
                #genus = "Coturnix", 
                #species = "japonica", 
	        #version = "2.0", 
                outputDir = "./org", 
                NCBIFilesDir = "./ncbi"#, 
                #rebuildCache=FALSE
                )
        # We specify a firectory to save locally the files used (and 
        # retrieved from NCBI) to create the Org database.
        # Normally, if the files are older than one day, they will be
        # downloaded again. Since they are very large, the download may
        # take too long and be interrupted frequently. This implies that
        # if we are unable to download everything in a single day, we are
        # doomed. The 'rebuildCache' option should avoid the re-downloading,
        # but it is described as an option for "internal use only and for
        # testing", and it may imply that no files are downloaded at all
        # and only local files in the cache directory are used, so we will
        # not use it unless strictly necessary.
        # GET
        # ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/*
        # ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/ARCHIVE/gene2unigene
        # ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz
        # if that fails try using https://... instead  <<< PREFERRED!!!
        # if that fails try curl

        # We can download the files by hand first and then run this command
        # so it uses the already downloaded files.
        
        # install the package
        install.packages(org.package, repos=NULL)
        if (require(org.package, character.only=T)) {
            org.db <- eval(parse(text=org.package))
            columns(org.db)
        } else {
            org.db <- NULL
            cat.warn('\nCould not generate org.db\n')
        }
        ### NOTE
        # This fails in Gallus gallus due to download failures from NCBI
        # because of excesses in download times.
        # We need to use the latest GitHub version installed with
        # library(devtools)
        # install_github("Bioconductor/AnnotationHub")
    }
}


# at this point ens.db should contain the ENSEMBL data and 
# org.db the Org type data.




# ---------------------------------------------------------------
# O B T A I N   B I O M A R T   A N N O T A T I O N
# ---------------------------------------------------------------
short.title('annotation (biomaRt)')

if ( file.exists(paste(annotation.dir, 'biomaRt.annotation.1st.txt', sep='/')) ) {
   # we have aready retrieved and saved the annotation, use it
    bm.annot.1 <- read.table(
    		file=paste(annotation.dir, '/biomaRt.annotation.1st.txt', sep=''), 
	        sep='\t', 
                header=T)
} else {
    # get all the annotation

    marts <- listMarts()
    if (VERBOSE) head(marts)
    datasets <- listDatasets(useMart("ensembl"))
    mart.db <- useMart("ensembl", mart.name)

    if (VERBOSE) {
        head(datasets)
        ## set up connection to ensembl database
        ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
    
        # list the available datasets (species) for the record
        listDatasets(ensembl) %>%  filter(str_detect(description, release))
        attributes <- listAttributes(mart.db)
        head(attributes)
        filters <- listFilters(mart.db)
        head(filters)

        # check the available "filters" - things you can filter for
        listFilters(mart.db) %>% 
            filter(str_detect(name, "ensembl"))

        # check the available "attributes" - things you can retreive
        attr <- listAttributes(mart.db)[,1:2]	# the first 200 are general
        listAttributes(mart.db)[,1:2] %>% 
            head(20)
    }
    
    # we cannot get all the annotation at once because it times out
    #full.annot <- getBM(attributes=
    #                       c("ensembl_gene_id", "ensembl_transcript_id", 
    #			   "start_position", "end_position", 
    #                          "chromosome_name", "gene_biotype", 
    #                          "description", 
    #                          "entrezgene_id", "entrezgene_accession", "entrezgene_description", 
    #                          "go_id", "name_1006", "definition_1006", "go_linkage_type", "namespace_1003", 
    #                          "goslim_goa_accession", "goslim_goa_description", 
    #                          "pdb", 
    #                          "reactome", "uniprotswissprot"), 
    #                       mart=mart.db)
    #
    # so we will retrieve the data in pieces, including ensembl_gene_id in
    # each piece so we can use it as key for merging the annotation and
    # saving it in a file to avoid future downloads
    # each of these calls will cache locally the annotation ro speed up
    # subsequent accesses
    bm.ensembl.annot <- get.biomart.ensembl.annotation(mart.db, annotation.dir)

    bm.entrez.annot <- get.biomart.entrez.annotation(mart.db, annotation.dir)

    bm.go.annot <- get biomart.go.annotation(mart.db, annotation.dir)

    bm.goslim.annot <- get.biomart.goslim.annotation(mart.db, annotation.dir)

    bm.fam.annot <- get.biomart.family.annotation(mart.db, annotation.dir)

    bm.prosite.annot <- get.biomart.prosite.annotation(mart.db, annotation.dir)

    bm.sfam.annot <- get.biomart.superfamily.annotation(mart.db, folder)

    bm.extra.annot <- get.biomart.extra.annotation(mart.db, folder)

    # Now that we have all the pieces, merge them all together
    # into a single annotation variable
    #	THIS TAKES TOO LONG AND TOO MUCH MEMORY, COMMENTED FOR NOW
    #bm.annot <- biomart.merge.annotations(list(
    #                  bm.ensembl.annot,
    #                  bm.entrez.annot,
    #                  bm.go.annot,
    #                  bm.goslim.annot,
    #                  bm.fam.annot,
    #                  bm.prosite.annot,
    #                  bm.sfam.annot,
    #                  bm.extra.annot
    #                ),
    #                by="ensembl_gene_id", folder)
    #
    #write.table(bm.annot, file=paste(folder, '/biomaRt.annotation.txt', sep=''), 
    #        sep='\t', row.names=T, col.names=T)


    # Now that we have the annotation we can select unique entries
    # for our dataset

    # One possible way to do it would be to filter the queries above
    #	to retrieve the annotation matching ensembl_ids
    #
    # We can set a field to use to filter the output data
    # Set the filter type and values
    #ourFilterType <- "ensembl_gene_id"
    # and the values to select from that field
    #filterValues <- rownames(fit)
    #
    # and then obtain the specified annotation from records that match the values
    # specified in the filter field
    #fit.bm.extra.annot <- getBM(attributes=c(
    #                       "ensembl_gene_id", 
    #                       "pdb",
    #                       "reactome", 
    #                       "uniprotswissprot"), 
    #                   mart=mart.db,
    #                   filters=ourFilterType,
    #                   values=filterValues)
    #                   
    # deduplicate selecting the first annotation
    #fit.bm.extra.annot.1 <- fit.bm.extra.annot[ ! duplicated(fit.bm.extra.annot$ensembl_gene_id), ]
    # then we would repeat this for each annotation subset and merge all of them
    # at the end...

    # or we could deduplicate everything first and match aftwerards
    # this has the advantage that we keep all the annotation at hand and
    # can reuse it for any gene dataset instead of annotating each
    # specific dataset separately
    bm.ensembl.annot.1 <- bm.ensembl.annot[ ! duplicated(bm.ensembl.annot$ensembl_gene_id), ]
    bm.entrez.annot.1 <- bm.entrez.annot[ ! duplicated(bm.entrez.annot$ensembl_gene_id), ]
    bm.go.annot.1 <- bm.go.annot[ ! duplicated(bm.go.annot$ensembl_gene_id), ]
    bm.goslim.annot.1 <- bm.goslim.annot[ ! duplicated(bm.goslim.annot$ensembl_gene_id), ]
    bm.extra.annot.1 <- bm.extra.annot[ ! duplicated(bm.extra.annot$ensembl_gene_id), ]

    # this should now be manageable (many entries should have been removed)
    bm.annot.1 <- merge(bm.ensembl.annot.1, bm.entrez.annot.1, by="ensembl_gene_id")
    bm.annot.1 <- merge(bm.annot.1, bm.go.annot.1, by="ensembl_gene_id")
    bm.annot.1 <- merge(bm.annot.1, bm.goslim.annot.1, by="ensembl_gene_id")
    bm.annot.1 <- merge(bm.annot.1, bm.extra.annot.1,  by="ensembl_gene_id")

    write.table(bm.annot.1, 
                file=paste(annotation.dir, '/biomaRt.annotation.1st.txt', sep=''), 
	        sep='\t', row.names=T, col.names=T)
    # we save it to avoid repeating this in the future

    # or even do it all at once? if we had enough power for building bm.annot:
    ##bm.annot.1 <- bm.annot[ ! duplicated(bm.annot$ensembl_gene_id), ]
    ##write.table(bm.annot.1, 
    ##            file=paste(folder, '/biomart.annotation.1st.txt', sep=''), 
    ##            sep='\t', row.names=T, col.names=T)
}

# at this point we have all the different annotation subsets we would
# like to use.


# And now we are ready with our reference info at hand...


#################################################################################
#
# Get the alignments and feature counts
#
#################################################################################

short.title("aligning")

if ( ALIGN == TRUE) {
    out <- alignment.dir
    align.fastq(path=path, reference=reference, aln.out=out, save.dir=folder) 
} else {
    out <- alignment.dir
    #out <- 'Alignments_Coturnix'
}

bam.files <- list.files(path = out, pattern = '.BAM$', full.names = TRUE)

# get the feature counts
#	we use EnsEMBL annotation from release 6a
#	this will tae the position counts and match them against the 
#	genes annotated in the GFF3 file for the reference, obtaining
#	a list of genes and the number of reads that mapped to them
#	as an indicator of their expression level

short.title("features")

cat('\nCOMPUTING FEATURE COUNTS\n')

if (file.exists(paste(folder, '/featureCounts.rds', sep=''))) {
    fc <- readRDS(file=paste(folder, '/featureCounts.rds', sep=''))
} else {
    # compute and save fc in cache
    fc <- compute.feature.counts(files=bam.files, 
                           reference=reference, 
                           requireBothEnds=requireBothEnds, 
                           save.dir=folder) 
}

#################################################################
#
#################################################################
                                                                #
#################################################################
#
#################################################################
                                                                #
#################################################################

# We are ready; we have
#	annotation (ensembldb, org.db, biomar)
#	feature counts table

if (USE.EDGER) {

    # ---------------------------------------------------------------
    # A N A L Y S I S     W I T H     E D G E R
    # ---------------------------------------------------------------

    short.title('edgeR')

    dge <- eR.differential.gene,expression(fc, 
    			metadata=metadata,
                threshold=cpm.threshold, 
                ens.db=ens.db, 
                folder=folder)
    
    fit <- eR.dge.voom.variation.analysis(dge, design.column, folder)

    fit <- eR.fit.annotate.ensembl.biomart.save(fit,
    					   				   		ens.db, 
                                           		bm.annot.1, 
                                           		folder, 
                                           		n.genes)

    # now create a volcano plot for only the top 1/2 genes using gene
    # annotation
    out.png <- paste(folder, '/edgeR/img/edgeR_volcanoplot.png', sep='')
    as.png(volcanoplot(fit, highlight=n.genes/2, coef=1, names=fit$genes$SYMBOL),
        out.png)

    # testing relative to a threshold
    fit.thres <- eR.fit.treat(fit, threshold=signif.threshold, folder)
    
    # save all the contents of 'fit' in an RDS file
    saveRDS(fit, file=paste(folder, '/edgeR/annotatedVOOMfit+.rds', sep=''))
    saveRDS(fit.thres, 
            file=paste(folder, '/edgeR/annotatedVOOMfit+.gt.',
                       significance.threshold, '.rds', sep=''))
    #	'fit' can later be recovered with: 
    #		fit <- readRDS(file=paste(folder, '/annotatedVOOMfit+.rds', sep=''))
    # and save as well as Rdata file
    save(fit, file=paste(folder, '/edgeR/annotatedVOOMfit+.RData', sep=''))
    save(fit, 
         file=paste(folder, '/edgeR/annotatedVOOMfit+.gt.',
                    significance.threshold, '.RData', sep=''))
    #	'fit' can later be recovered with: 
    #		fc <- load(file=paste(folder, '/annotatedVOOMfit+.RData', sep=''))


    # -----------------------------------------------------------------
    # Getting beyond here is easier with DESeq2
    # -----------------------------------------------------------------
    # 
    # # redefine the design so we can make any kind of comparison by
    # # using a null reference ("0 + ")
    # 
    # this carries out all comparisons and annotates each fit obtained
    # saving all of them in a list
    eR.data <- eR.dge.all.comparisons(dge, design.column, ens.db, bm.annot.1, folder)

    # now eR.data is a list where each element is a comparison A-B (A minus B),
    # i.e. each A-B is a list of tables, one of them named "table" and containing
    # logFC, logCPM, F and PValue
    #
    # We would like to have FDR-corrected p-values as well, which can be got by
    # defining a threshold. But that implies we know of a meaningful one, which
    # we don't yet.

}    # end if (USE.EDGER)



#################################################################
#
#################################################################
                                                                #
#################################################################
#
#################################################################
                                                                #
#################################################################

# ---------------------------------------------------------------
# A N A L Y S I S     W I T H     D E S E Q 2
# ---------------------------------------------------------------

if (USE.DESEQ2) {
    # get counts
    if (exists(get.name.of(fc))
        countData <- fc$counts
    else {
        ### NOTE @@
        # this should not be needed
        countDataFile <- paste(folder, "/featureCounts.csv", sep='/')

        countData <- read.csv(countDataFile,
		        row.names=1) %>%
		        as.matrix()
    }
    countData <- countData[rowSums(countData)>1, ]

    if (VERBOSE) print(head(countData))

    # get metadata
    colData <- read.delim(paste(rnaseq.out, metadata, sep='/'))

    name <- paste(folder, "DESeq2/dds.DESeq2", sep='/')
    if ( ! file.exists(paste(name, "rds", sep='.')) ) {

        dds <- DESeqDataSetFromMatrix(countData, 
                                      colData,  
                                      design=paste('~', design.column, sep='') 
                                      tidy=F)

        dds <- DESeq(dds)

        # SAVE DESEQ ANALYSIS
        # -------------------
        name <- paste(folder, "/DESeq2/dds.DESeq2", sep='')
        save(dds, file=paste(name, "RData", sep='.'))
        #	'dds' can later be recovered with: 
        #		dds <- load(file=paste(folder, "/dds.DESeq2.RData", sep=''))
        saveRDS(dds,  file=paste(name, "rds", sep='.'))
        # read as follows:    
    } else {
        dds <- readRDS(file=paste(name, "rds", sep='.'))
    }

    cat("DESeq2 analysis produced the following fit")
    print(resultsNames(dds))

    annot <- ds2.get.annotation.ens.org(dds, ens.db, org.db)
    ens.ann <- annot$ens.ann		### NOTE we could use attach here
    ens.ann.1 <- annot$ens.ann.1	# but this is more explicit
    geneSymbols <- annot$geneSymbols
    go.ann <- annot$go.ann
    GOdescription <- annot$GOdescription
    
    ### NOTE
    # change by function
    # annotate with ensembl ens.db
    ens.ann <- ensembldb::select(ens.db, 
                      column='GENEID', keytype= 'GENEID', keys=rownames(dds), 
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


    # SAVE SELECTED COMPARISONS
    # -------------------------
    # We will get all comaprisons
    for (a in levels(as.factor(target[ , design.column])) ) {
        for (b in levels(as.factor(target[ , design.column])) ) {
            print(paste(a, b))
            if (a == b) next
            ds2.dds.compare.annotate.save( 
                            dds,
                            column=design.column,
                            x=a, y=b,
                            filterFun=ihw, alpha=0.01,
                            ensembl.ann=ens.ann.1,
                            biomart.ann=bm.annot.1,
                            outDir=folder
                            )

        }
    }



    # SAVE GENES WITH SIGNIFICANT CHANGES
    # -----------------------------------

    ds.data <- list()
    #x <- 0
    contr <- design.column
    grps <- levels(colData[ , contr ])
    for (a in grps) {
        for (b in grps) {
            print(paste(a, b))
            if (a == b) next
            res <- ds2.dds.plot.and.save( 
                            dds, contr,
                            x=a, y=b,
                            filterFun=ihw, alpha=0.01,
                            ensembl.ann=ens.ann.1,
                            biomart.ann=bm.annot.1,
                            outDir=folder
                            )
            print(names(res))
	        name <- paste(contr, "_", a, "_", b, sep='')
            #x <- x + 1
            #ds.data[[x]] <- res
            ds.data[[name]] <- res
            #names(ds.data)[x] <- name
            #stop()
        }
    }
    print(names(ds.data))

    if (INTERACTIVE) {
        ds2.interactively.print.n.most.significant.de.genes(ds.data, n=10)

        ds2.interactively.plot.top.up.down.regulated.genes(dds, 
                                                           ds.data, 
                                                           design.column)
    }


    # -------------------------------
    # Do Gene Set Enrichment Analysis
    # -------------------------------


    # we do not know what is the best upper limit, so we'll try several
    ms <- 500 # default value in GSEA, should work as well as the others
    #for (ms in c(50, 100, 250, 500)) {
    for (ms in c(500)) {					# use this to save computation time
        for (cmp.name in names(ds.data)) {
            # for each comparison name
            cmp.data <- ds.data[[ cmp.name ]]
            cat('Doing GSEA on', cmp.name, '\n')
            ann.shrunk.lfc <- cmp.data[[ "shrunk.annot" ]]

    #        #ms <- 500 
    #        #out.dir <- paste(folder, "DESeq2/GO_fgsea", cmp.name, sep='/')
    #        out.dir <- sprintf("%s/DESeq2/GO_fgsea/max.size=%03d/%s",
    #                folder, ms, cmp.name)
    #        dir.create(out.dir, showWarning=FALSE, recursive=TRUE)
    #        gogsea <- GO_fgsea(ann.shrunk.lfc, 
    #                   max.size=ms,
    #                   out.dir=out.dir, 
    #                   out.name='fgsea',
    #                   use.description=TRUE)
    #        # we should add fgsea results to the cmp.data list
    #        ds.data[[cmp.name]][['go.fgsea']] <- gogsea

            #ms <- 500 
            out.dir <- sprintf("%s/DESeq2/GO+KEGG_cProf/max.size=%03d/%s",
                    folder, ms, cmp.name)
            dir.create(out.dir, showWarning=FALSE, recursive=TRUE)
            GO_KEGG_clusterProfiler(
                          ann.shrunk.lfc, 
		          max.size=ms,
                          out.dir=out.dir, 
                          out.name='cProf',
                          use.description=TRUE,
                          verbose=FALSE)
        }
    }


    # -------------------------
    # Analyze GO representation
    # -------------------------
    ds2.analyze.go.representation(ds.data, bm.go.annot, folder)

    ds2.analyze.go.ancestry(ds.data, l2fc.threshold=0, folder)

    ds2.analyze.pfam.representation(ds.data, bm.fam.annot, folder)


    # ------------------------------
    # DO PCA AND CLUSTERING ANALYSES
    # ------------------------------

    ## CREATE THE COLUMNS FOR THE PCA ANALYSIS 
    # we need to retrieve from the dataframes only the log2Fold change and the
    # gene names
    # we will take the data from the unsorted tables signif_* 
    # we save the rownames as a distinct column and we pass it as a column 
    # we delete column1 (gene-names) before clustering because it is not needed

    # now that we have the data we need to keep only the common rows to all of them
    # we use the function intersect by pairs and then all together
    # the common genes are rows from common2
    # finally we create a dataframe with everything we have


    # This is here in case we decide to loop over several columns later 
    contrasts.column <- design.column
    references <- levels(colData[ , contrasts.column ])

    for (ref in references) {
        # Find genes that change w.r.t. the reference strain

        # create a convenience text variable to simplify/unify filenames below
        ccol_ref <- paste(contrasts.column, ref, sep='_')

        # find all the genes common to all samples
        name <- paste(ccol_ref, levels(as.factor(target[ , contrasts.column]))[1], sep='_')
        common <- rownames(ds.data[[name]]$signif)
        for (i in levels(as.factor(target[ , contrasts.column]))) {
            if (i == ref) next
            name <- paste(ccol_ref, i, sep='_')
            print(name)
            common <- intersect(common, rownames(ds.data[[name]]$signif))
        }
        length(common)	# 681 in Coturnix, 1670 in Gallus

        data.table <- data.frame(genes=common)
        # compare all against all other samples
        for (i in levels(as.factor(target[ , contrasts.column]))) {
            if (i == ref) next
            name <- paste(ccol_ref, i, sep='_')
            print(name)
            data.table[name] <- ds.data[[name]]$signif[common, "log2FoldChange"]
        }
        rownames(data.table) <- common

        data.annot <- ds.data[[name]]$signif.annot[common, ]

        par(mfrow=c(1,1))

        set.seed(1963)

        # First have a general look at the methods to get a feeling for the
        # best number of clusters


        dif <- data.table[ , -1]
        by.row <- 1
        by.col <- 2
        means <- apply(dif, by.col, mean)
        sds <- apply(dif, by.col, sd)
        nor <- scale(dif,center=means,scale=sds)
        if (INTERACTIVE) {
            # Do a scatterplot matrix
            car::scatterplotMatrix(dif)
        }
        out.png <- paste(folder, "/DESeq2/img/", 
                   'scatterplot_matrix_', ccol_ref, ".png", 
                   sep='')
        as.png( {
                print(car::scatterplotMatrix(dif))
                }, out.png)


        # Try to guess the optimum number of K-means clusters
        # NBClust
        out.log <- paste(folder, "/DESeq2/cluster/NBClust_", ccol_ref, ".log", sep='')
        openlog(out.log)
        #library("NbClust")
        # predict best number of clusters for hierarchical and k-means clustering
        nbc <- NbClust(dif, diss=NULL, 
                distance="euclidean", method="complete", 
                min.nc=3, max.nc=10, 
                index="all", alphaBeale=0.1)
        # 4 for the C. japonica wt vs others data
        # 5 for G. gallus wt vs infected.data
        sink()


        # Let the user see the various clusters and make a decision
        #
        out.log <- paste(folder, "/DESeq2/cluster/kmeans/DESeq2_kmeans_all_", 
                         ccol_ref, ".log", sep='')
        openlog(out.log)
        for (i in 2:10) {
            cat("Clustering with K-means (", i, " clusters)\n")
            cl <- kmeans(nor, i)				# NOTE: nor
            print(table(cl$cluster))
            print(fviz_cluster(cl, geom = "point", data=nor))
            out.png <- sprintf(
                    "%s/DESeq2/cluster/kmeans/DESeq2_kmeans_%s_nc=%03d.png", 
        	    folder, ccol_ref, i)
            as.png(fviz_cluster(cl, geom = "point", data=nor,
                        main=paste("K-means clsutering nc =", i) ),
                    out.png)
            if (INTERACTIVE)
                continue.on.enter("Press [ENTER] to continue ")
        }
        sink()

        cat("The log and plots for K-means clustering have been saved to\n", 
            "\t", folder, "/DESeq2/cluster/kmeans/\n", 
            "please, inspect them and select the best number of clusters\n",
            sep='')

        #kmeans.nc <- readline("How many clusters should I use for k-means? ")
        #kmeans.nc <- as.numeric(kmeans.nc)
        # for gg
        kmeans.nc <- 5


        out.log <- paste(folder, "/DESeq2/cluster/pam/DESeq2_pam_all_", 
                         ccol_ref, ".log", sep='')
        openlog(out.log)
        for (i in 2:10) {
            cat("Clustering with PAM (", i, " clusters)\n")
            cl <- pam(nor, i)					# NOTE:nor
            print(table(cl$cluster))
            print(fviz_cluster(cl, geom = "point"))
            out.png <- sprintf("%s/DESeq2/cluster/pam/DESeq2_pam_%s_nc=%03d.png", 
                folder, ccol_ref, i)
            as.png(fviz_cluster(cl, geom = "point",
                        main=paste("Partitioning Around Medoids clsutering nc =", i) ), 
                    out.png)
            #continue.on.enter("Press [ENTER] to continue ")
        }
        sink()

        cat("The log and plots for PAM clustering have been saved to\n", 
            "\t", folder, "/DESeq2/cluster/pam/\n", 
            "please, inspect them and select the best number of clusters\n",
            sep='')

        #pam.nc <- readline("How many clusters should I use for PAM? ")
        #pam.nc <- as.numeric(pam.nc)
        # for gg
        pam.nc <- 7

        out.log <- paste(folder, "/DESeq2/cluster/dbscan/DESeq2_dbscan_all_", 
                         ccol_ref, ".log", sep='')
        openlog(out.log)
        for (e in c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2.0)) {
            cat("\nClustering with DBScan ( eps =", e, " )\n")
            cl <- dbscan(dif, eps=e, MinPts=10, showplot=1)	# NOTE: dif
            print(table(cl$cluster))
            print(fviz_cluster(cl, data=dif, 
                         show.clust.cent=FALSE, labelsize=4,
                         ellipse=TRUE, ellipse.type="convex"))
            out.png <- sprintf("%s/DESeq2/cluster/dbscan/DESeq2_dbscan_%s_eps=%03.2f.png", 
                         folder, ccol_ref, e)
            as.png(fviz_cluster(cl, data=dif, 
                          show.clust.cent=FALSE, 
                          geom="point", #labelsize=4,
                          ellipse=TRUE, ellipse.type="convex",
                          main=paste("Density based clustering eps =", e) ),
                    out.png)
            #continue.on.enter("Press [ENTER] to continue ")
        }
        sink()

        cat("The log and plots for DBscan clustering have been saved to\n", 
            "\t", folder, "/DESeq2/cluster/dbscan/\n", 
            "please, inspect them and select the best epsilon value\n",
            sep='')

        #dbscan.eps <- readline("Which epsilon should I use for DBscan? ")
        #dbscan.eps <- as.numeric(dbscan.eps)
        # for gg
        dbscan.eps <- 1.5


        hclust.nc <- 5	# we'll set it by hand for now

        # Now, proceed in detail, method by method, with a deeper and more 
        # detailed analysis

        # possible distances are c("euclidean", "maximum", "manhattan", "canberra",
        #	"binary", "minkowski", "pearson", "spearman", "kendall")
        # possible methods are c("ward.D", "ward.D2", "single", "complete",
        #	"average" (UPGMA), "mcquitty" (WPGMA), "median" (WPGMC), "centroid" (UPGMC))
        out.log <- paste(folder, "/DESeq2/cluster/hcluster/DESeq2_hcluster_", contrasts.column, "_", ref, ".txt", sep='')
        sink(out.log, split=T)
        h.cluster.changes(
            #data.table=pca.table[ , -1],
            #annotated.data=ds.data[[name]]$signif.annot[common, ],
            clusters=hclust.nc,
            data.table=dif,
            annotated.data=data.annot
        )
        sink()


        boot <- 100

        hclust_folder <- paste(folder, "/DESeq2/cluster/hclust_", ccol_ref, sep='')
        dir.create(hclust_folder, showWarnings=FALSE)

        hcut_folder <- paste(folder, "/DESeq2/cluster/hcut_", ccol_ref, sep='')
        dir.create(hcut_folder, showWarnings=FALSE)

        out.log <- paste(hcut_folder, "/DESeq2_hcut_", 
                         ccol_ref, ".txt", sep='')
        sink(out.log, split=T)
        hc.cl <- cluster.changes(
            #data.table=pca.table[ , -1],
            #annotated.data=sa.data[[name]]$signif.annot[common, ],
            data.table=dif,
            annotated.data=data.annot,
            FUN=hcut,
            clusters=hclust.nc,
            estimate=TRUE,
            gap_bootstrap=boot,
            output.folder=hcut_folder
            )
        sink()


        kmeans_folder <- paste(folder, "/DESeq2/cluster/kmeans_", ccol_ref, sep='')
        dir.create(kmeans_folder, showWarnings=FALSE)

        out.log <- paste(kmeans_folder, "/DESeq2_kmeans_", 
                         ccol_ref, ".txt", sep='')
        sink(out.log, split=T)
        km.cl <- cluster.changes(
            #data.table=pca.table[ , -1],
            #annotated.data=ds.data[[name]]$signif.annot[common, ],
            data.table=dif,
            annotated.data=data.annot,
            FUN=kmeans,
            algorithm="Hartigan-Wong",
            clusters=kmeans.nc,
            nstart=kmeans.nc*10,
            estimate=T,
            gap_bootstrap=boot,
            output.folder=kmeans_folder
            )
        sink()


        pam_folder <- paste(folder, "/DESeq2/cluster/pam_", ccol_ref, sep='')
        dir.create(pam_folder, showWarnings=FALSE)

        out.log <- paste(pam_folder, "/DESeq2_pam_", 
                         ccol_ref, ".txt", sep='')
        sink(out.log, split=T)
        pam.cl <- cluster.changes(
            #data.table=pca.table[ , -1],
            #annotated.data=ds.data[[name]]$signif.annot[common, ],
            data.table=dif,
            annotated.data=data.annot,
            FUN=pam,
            distance="euclidean",	# see dist()
            clusters=pam.nc,	# n. of clusters
            nstart=pam.nc*10,	# n. of random start sets to choose
            estimate=T,
            gap_bootstrap=boot,
            output.folder=pam_folder
            )
        sink()


        dbscan_folder <- paste(folder, "/DESeq2/cluster/dbscan_", ccol_ref, sep='')
        dir.create(dbscan_folder, showWarnings=FALSE)

        out.log <- paste(dbscan_folder, "/DESeq2_dbscan_", 
                         ccol_ref, ".txt", sep='')
        sink(out.log, split=T)
        dbs.cl <- cluster.changes(
            #data.table=pca.table[ , -1],
            #annotated.data=ds.data[[name]]$signif.annot[common, ],
            data.table=dif,
            annotated.data=data.annot,
            FUN=dbscan,
            distance="euclidean", # see dist()
            eps=dbscan.eps,
            normalize=F,
            estimate=T,
            gap_bootstrap=boot,
            output.folder=dbscan_folder
            )
        sink()



        # Now cluster by experiment
        # -------------------------

    ## ALREADY DONE ABOVE
    ##    # Find genes that change w.r.t. the reference strain
    ##
    ##     name <- paste(contrasts.column, ref, levels(as.factor(target[ , contrasts.column]))[1], sep='_')
    ##     common <- rownames(ds.data[[name]]$signif)
    ##     for (i in levels(as.factor(target[ , contrasts.column]))) {
    ##         if (i == ref) next
    ##         name <- paste(contrasts.column, ref, i, sep='_')
    ##         print(name)
    ##         common <- intersect(common, rownames(ds.data[[name]]$signif))
    ##     }
    ##     length(common)	# 681 in Coturnix, 1670 in Gallus
    ## 
    ##     data.table <- data.frame(genes=common)
    ##     #for (i in c("wt", "PC", "P", "1.0", "0.1")) {
    ##     for (i in levels(as.factor(target[ , contrasts.column]))) {
    ##     # we do not want to include PC this time
    ##     #for (i in c("P", "1.0", "0.1")) {
    ##         if (i == ref) next
    ##         name <- paste(contrasts.column, ref, i, sep='_')
    ##         print(name)
    ##         data.table[name] <- ds.data[[name]]$signif[common, "log2FoldChange"]
    ##     }
    ##     rownames(data.table) <- common
    ## 

        fid <- t(dif)
        means <- apply(fid, by.col, mean)
        sds <- apply(fid, by.col, sd)
        ron <- scale(fid,center=means,scale=sds)

        # here we have a small number of rows and can set a maximum number of clusters
        maxclust <- nrow(fid) - 1

        # hclust
        distan = dist(ron, method="euclidean")
        hcl <- hclust(distan)
        plot(hcl,labels=rownames(fid),main='Default from hclust')
        out.png <- sprintf(
                "%s/DESeq2_hcluster_grps_%s.png",
	        hclust_folder, ccol_ref)
        as.png(plot(hcl,labels=rownames(fid),main='Default from hclust'), out.png)

        # kmeans
        fviz_nbclust(ron, kmeans, method="silhouette", k.max=maxclust)
        out.png <- sprintf(
                "%s/DESeq2_kmeans_grps_silhouette_%s.png",
	        kmeans_folder, ccol_ref, ccol_ref)
        as.png(fviz_nbclust(ron, kmeans, method="silhouette", k.max=maxclust), out.png)
        kcl <- kmeans(ron, centers=3, nstart=100)
        fviz_cluster(kcl, data=ron)
        out.png <- sprintf(
                "%s/DESeq2_kmeans_grps_%s_nc=%03d.png",
	        kmwans_folder, ccol_ref, ccol_ref, 3)
        as.png(fviz_cluster(kcl, data=ron), cout.png)

        # pam
        fviz_nbclust(ron, pam, method="silhouette", k.max=maxclust)
        out.png <- sprintf(
                "%s/DESeq2_pam_grps_silhouette_%s.png",
	        pam_folder, ccol_ref)
        as.png(fviz_nbclust(ron, pam, method="silhouette", k.max=maxclust), out.png)
        pcl <- pam(ron, k=3, diss=F)
        fviz_cluster(pcl, data=ron)
        out.png <- sprintf(
                "%s/DESeq2_pam_grps_%s_nc=%03d.png",
	        folder, ccol_ref, 3)
        as.png(fviz_cluster(pcl, data=ron), out.png)

        # dbscan
        # this results in the same two PCs but one cluster
        fviz_nbclust(ron, dbscan, method="silhouette", k.max=maxclust)
        out.png <- sprintf(
                "%s/DESeq2_dbscan_grps_silhouette_%s_%s.png",
	        dbscan_folder, ccol_ref)
        as.png(fviz_nbclust(ron, dbscan, method="silhouette", k.max=maxclust), out.png)
        ### JR ### eps should be tuned for each experiment
        dcl <- dbscan(ron, eps=0.2, MinPts=2, showplot=1)
        fviz_cluster(dcl, data=ron)
        out.png <- sprintf(
                "%s/DESeq2_dbscan_grps_%s_eps=%03.2f.png",
	        dbscan_folder, ccol_ref, 0.2)
        as.png(fviz_cluster(dcl, data=ron), out.png)

        # this fails
        tryCatch( {
            NbClust(ron, diss=NULL, 
                    distance="euclidean", method="complete", 
                    min.nc=3, max.nc=10, 
                    index="all", alphaBeale=0.1)
        } )
    }

}

while (sink.number() > 0) sink()


