##==========================##
##   ANNOTATION FUNCTIONS   ##
##==========================##

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
                                     pass,
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
    		             author=author,
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

        if (length(ann == 2)) return(bm.annot)

        # we have more annotations to merge
        for (i in 3:length(ann))
            bm.annot <- merge(bm.annot, ann[[i]], by=by)
        write.table(bm.annot, file=paste(folder, '/biomaRt.annotation.txt', sep=''), 
                    sep='\t', row.names=T, col.names=T)
    } else {
        bm.annot <- read.table(paste(folder, '/biomaRt.annotation.txt', sep=''), 
	        header=T, sep='\t')
    }
    return(bm.annot)
}

