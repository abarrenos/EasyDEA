
# 1) ANNOTATE RESULTS

annotateDESeqResults <- function(results,
				ensembl.db,
				biomart.db,
				org.db = NULL,
				annotation.dir,		# For saving annotation dataset
				save = TRUE,
				out.file.base,
				out.dir ){
	
	cat("\n	Annotating genes ...\n\n")
	
	## Get gene annotation from ENSEMBL databsase and BioMart database

	if ( ! file.exists(paste(annotation.dir, 'ensembl.annotation.1st.txt', sep = '/'))) {
		ensembl.ann <- annotateGenesFromENSEMBL(ensembl.db = ensembl.db,
							gene.ids = rownames(results),
							save = save,
							out.dir = annotation.dir)
	} else {
		ensembl.ann <- read.table(paste(annotation.dir, 'ensembl.annotation.1st.txt', sep = '/'))
	}
	
	if ( ! is.null(org.db) ){
		## Get gene annotation from ENSEMBL databsase and BioMart database
		if ( ! file.exists(paste(annotation.dir, 'orgDb.GO.annotation.1st.txt', sep = '/'))) {
			org.ann <- annotateGenesFromORG(org.db = org.db,
								gene.ids = rownames(results),
								save = save,
								out.dir = annotation.dir)
		} else {
			org.ann <- read.table(paste(annotation.dir, 'orgDb.GO.annotation.1st.txt', sep = '/'))
		}
	}
	
	# Remove duplicated records preserving only first entry

	ensembl.ann <- ensembl.ann[ ! duplicated(ensembl.ann$GENEID), ]
	biomart.ann <- biomart.db[ ! duplicated(biomart.db$ensembl_gene_id), ]
	if ( ! is.null(org.db) ) org.ann <- org.ann[ ! duplicated(org.ann$ENTREZID), ]


	## Add annotation to DESeq results data.frame to keep everything in one place

	# Add row names (ensembl.gene.id) as an additional column
	res.ann <- results
	res.ann$gene.id <- rownames(res.ann)
	
	### Identify keys for ENSEMBL ###
	if (all(str_sub(res.ann$gene.id, 1, 3) == "ENS")){
        by='GENEID'
    }else{
    	by='GENENAME'}

	# ENSEMBL GENENAMES match the GENE.IDs from our results dataframe
	#res.ann <- cbind(res.ann, ensembl.ann[ match(res.ann$ensembl.gene.id, ensembl.ann$GENENAME), ])
	for (i in colnames(ensembl.ann)){
		res.ann[i] <- ensembl.ann[match(res.ann$gene.id, ensembl.ann[,by], nomatch = NA), i]
	}	
	
	# BIOMART ENTREZ ACCESSION match the GENE.IDs from our results dataframe
	#res.ann <- cbind(res.ann, biomart.ann[ match(res.ann$ensembl.gene.id, biomart.ann$entrezgene_accession), ])
	for (i in colnames(biomart.ann)){
		res.ann[i] <- biomart.ann[match(res.ann$GENEID, biomart.ann$ensembl_gene_id, nomatch = NA), i]
	}

	if ( ! is.null(org.db) ){

		### Identify keys for NCBI (OrgDb Package) ###
		if (all(str_sub(res.ann$gene.id, 1, 3) == "ENS")){
        	by = 'ENSEMBL'
    	} else if (! any(is.na(as.numeric(res.ann$gene.id)))){
    		by = 'ENTREZID'
		} else {
			by = 'SYMBOL'}

		# NCBI SYMBOL match the GENE.IDs from our results dataframe
		for (i in colnames(org.ann)){
			res.ann[i] <- org.ann[match(res.ann$gene.id, org.ann[,by], nomatch = NA), i]
		}
	}

	## Save annotated results as data frame (table)
	if (save){
		write.table(data.frame(res.ann), sep = "\t",
					file = paste(out.dir, '/DESeq2/', out.file.base, '_ann_results.tab', sep = ''),
					row.names = TRUE, col.names = TRUE)

		saveRDS(res.ann, paste(out.dir, "/DESeq2/", out.file.base, "_ann_results.rds", sep=""))
	}
	return(res.ann)
}

# 2) COMPARE AND ANNOTATE RESULTS FROM DESEQ2 OBJECT
ddsCompareAnnotate <- function(	dds, 
								contrast,			# c(factor, numerator, denominator) 
								filterFun = ihw,	# IHW increases statistical power
								alpha = 0.01,
								annotate = TRUE, 
								ensembl.db,
								biomart.db,
								org.db = NULL,
								annotation.dir,		# For saving annotation tables
								save = TRUE,
								out.file.base,
								out.dir ) {

	results <- results(dds,
					contrast = contrast, 
					pAdjustMethod="BH",
					filterFun = filterFun,
					alpha = alpha)
	
	print(head(results))
	print(mcols(results)$description)	# Description of the Results fields

    ## Provide default output base name
    if (missing(out.file.base)) {
		out.file.base <- paste(contrast[1], contrast[2], 'vs', contrast[3], sep = '_')
	}

	if (save){
	
		## Save comparison results
		write.table(data.frame(results),
				file = paste(out.dir, '/DESeq2/raw/', out.file.base, '_raw_results.tab', sep = ''),
				row.names = TRUE, col.names = TRUE)

		saveRDS(results, paste(out.dir, "/DESeq2/raw/", out.file.base, "_raw_results.rds", sep=""))
	
		## Save summary
		sink(paste(out.dir, '/DESeq2/', out.file.base, '_summary.txt', sep=''), split=T)
		summary(results)
		sink()
	}
	
	if (annotate){
		results <- annotateDESeqResults(results = results,
										ensembl.db = ensembl.db,
										biomart.db = biomart.db,
										org.db = org.db,
										annotation.dir = annotation.dir,
										save = save,
										out.file.base = out.file.base,
										out.dir = out.dir )
	}
	return(results)
}



# 3) IDENTIFY TOP N DIFFERENTIALLY EXPRESSED GENES
top_n <- function(ds, n, col = "padj", decreasing = F) {
	
	df <- as.data.frame(ds)
	ordered.idx <- order(abs(df[,col]), decreasing = T)
	ordered.df <- df[ordered.idx, ]
	
	if (dim(ordered.df)[1] < n) n <- dim(ordered.df)[1]
	top.n <- ordered.df[1:n, ]
	
	return(top.n)
}


# 4) PERFORM DESEQ2 COMPARISON, ANNOTATE AND PLOT THE RESULTS

analysePlotDESeq <- function(	#results,
								dds, 
								contrast,			# c(factor, numerator, denominator) 
								filterFun = ihw,	# IHW increases statistical power
								alpha = 0.01,
								shrnk.type,
								annotate = TRUE,
								ensembl.db,
								biomart.db,
								org.db = NULL,
								annotation.dir,		# For saving annotation data
								save = TRUE,
								overwrite = FALSE,
								out.file.base,
								out.dir ) {

    ## Provide default output base name
    if (missing(out.file.base)) {
		out.file.base <- paste(contrast[1], contrast[2], 'vs', contrast[3], sep = '_')
	}

    ## Check if final comparison results already exist. If so 
	if ( ! file.exists(paste(out.dir, "/DESeq2/", out.file.base, "_final_result_list.rds", sep = "")) || (overwrite == TRUE)){

		## Check if DESeqResults object is provided
		
		results <- ddsCompareAnnotate(	
							dds = dds, 
							contrast = contrast,
							filterFun = ihw,
							alpha = alpha,
							annotate = annotate, 
							ensembl.db = ensembl.db,
							biomart.db = biomart.db,
							org.db = org.db,
							annotation.dir = annotation.dir,	
							save = save,
							out.file.base = out.file.base,
							out.dir = out.dir)

		## Histogram plot p-value distribution
		if (save) {
			out.png <- paste(out.dir, '/DESeq2/img/DESeq2_', out.file.base, '_padj_hist.png', sep='')
		} else out.png <- NULL

		as.png( {
				margins <- par("mar")
				par(mar = c(5, 5, 5, 5))
				hist(results$padj, main=paste("Adj. p-value distribution", "\n", out.file.base, sep = ''),
				breaks = 1/alpha, xlab = "Adj. p-value")
				par(mar = margins)
				},
		out.png, overwrite = TRUE)
		
		
		## Shrinkage of data to improve visualization and ranking
		cat("\n\tShrinking Log 2 Fold Change ...\n\n")    
		coef <- paste(contrast[1], contrast[2], 'vs', contrast[3], sep = '_')
		shrunk.lfc <- lfcShrink(dds = dds, coef = coef, type = shrnk.type)
		
		if (save) {
		    out.png <- paste(out.dir, "/DESeq2/img/DESeq2_", out.file.base, "_shrunk_", shrnk.type, "_MA.png", sep='') 
		} else out.png <- NULL
		
		as.png( {

			dim <- par("mfrow")
			par(mfrow = c(2,1))
		   
			# plotMA shows the log2 fold changes attributable to a given variable
			# over the mean of normalized counts for all the samples in the dataset
			# Points above alpha are colored, outliers are shown as directed triangles    
			DESeq2::plotMA(results, alpha = alpha, main = "Log2 Fold Change")
			
			# it is useful to look at the shrunk l2fc values, which removes the noise
			# associated with l2fc changes from low-count genes
			DESeq2::plotMA(shrunk.lfc, alpha = alpha, main = paste("Shrunken LFC (", shrnk.type, ")", sep = ""))
			
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
			par(mfrow = dim)
			},
			
		out.png, overwrite = TRUE)
		
		## Save shrunken results
		if (save) {
		
		    out.file <- paste(out.dir, "/DESeq2/shrunk/", out.file.base, "_shrunk_", shrnk.type, "_raw_results.rds", sep='')
		    cat("\n\tSaving", out.file, '\n\n')
		    saveRDS(shrunk.lfc, file=out.file)
		
			out.file <- paste(out.dir, "/DESeq2/shrunk/", out.file.base, "_shrunk_", shrnk.type, "_raw_results.tab", sep='')
			cat("\tSaving", out.file, '\n')
			write.table(shrunk.lfc, 
				        file=out.file,
				        row.names=T, col.names=T, sep='\t')
		}
		
		## Annotate shrunken data for visualization in MA plot	
		if (annotate) {
			
			ann.shrunk <- annotateDESeqResults(results = shrunk.lfc,
									ensembl.db = ensembl.db,
									biomart.db = biomart.db,
									org.db = org.db,
									annotation.dir = annotation.dir,
									save = FALSE)	
			if (save) {

			## Save shrunken annotated results		
				out.file <- paste(out.dir, "/DESeq2/shrunk/", out.file.base, "_shrunk_", shrnk.type, "_ann_results.rds", sep='')
				cat("\n\tSaving", out.file, '\n\n')
				saveRDS(ann.shrunk, file=out.file)
			
				out.file <- paste(out.dir, "/DESeq2/shrunk/", out.file.base, "_shrunk_", shrnk.type, "_ann_results.tab", sep='')
				cat("\tSaving", out.file, '\n')
				write.table(ann.shrunk, 
					        file=out.file,
					        row.names=T, col.names=T, sep='\t')
			}
		# Else, unnannotated shrunk.lfc will be represented in 
		} else ann.shrunk <- shrunk.lfc
		
		if (FALSE) {
		    if (save) {
		        out.file <- paste(out.dir, "/DESeq2/img/DESeq2_", out.file.base, "_shrunk_", shrnk.type, "_l2FC.png", sep='')
		        cat("\n\tPlotting", out.file, '\n\n')
		    } else out.file <- NULL
		    
		    as.png( {
					ggplot2::ggplot(data.frame(ann.shrunk), 
					ggplot2::aes(x = log2(baseMean), y=log2FoldChange),
					environment=environment()) + 	# this is supposed to make it work in a local env
					
					ggplot2::geom_point(ggplot2::aes(colour = padj < alpha), shape = 10, size = 1) +
					
					#ggplot2::geom_text(data=~top_n(.x, 10, wt=-padj),
					ggplot2::geom_text(data = top_n(ann.shrunk, 10, "padj"),
					ggplot2::aes(label = rownames(ann.shrunk))) +
					ggplot2::labs(x="Log2 Mean of Normalised Counts", y="Log2 Fold Change")
		    
		    }, out.file , overwrite = TRUE)
		}
		
		# Find statistically significant changes. Ignore genes with NAN p-value.
		signif <- results[ (results$padj < alpha) & ! is.na(results$padj) , ]
		signif$abs_lfc <- abs(signif$log2FoldChange)
		
		if (save) {
		    out.file <- paste(out.dir, '/DESeq2/signif/signif_', out.file.base, "_α<", alpha, ".tab", sep='')
		    cat("\n\tSaving", out.file, '\n')
		    write.table(data.frame(signif), file = out.file, sep = '\t', row.names = T, col.names = T)
		}
		
		# Order the data, firstly by the decreasing absolute 
		# value of log2FoldChange and secondly by the increasing pvalue....
		srt.signif <- signif[ order( signif$abs_lfc,
								signif$padj,
								decreasing=c(T, F)), ]
		if (save) {
		    out.file <- paste(out.dir, "/DESeq2/signif/signif_sorted_", out.file.base, ".tab", sep='')
		    cat("\tSaving", out.file,'\n\n')
		    write.table(srt.signif, file = out.file, sep = '\t', row.names = T, col.names = T)
		}
		
		## If shrunken data is unnanotated, return null.
		if ( ! annotate) ann.shrunk <- NULL
		
		# Generate output list containing results    
		res.list <- list(	result.ann = results, 
		            		shrunk = shrunk.lfc, 
		            		shrunk.ann = ann.shrunk,
							signif = srt.signif[,1:7],
		            		signif.ann = srt.signif)
		
		# Save results as an object to avoid redundant analyses.
		cat("\n\tSaving final results ...\n\n")		
		saveRDS(res.list, file=paste(out.dir, "/DESeq2/", out.file.base, "_final_result_list.rds", sep=''))
		save(res.list, file=paste(out.dir, "/DESeq2/", out.file.base, "_final_result_list.RData", sep=''))
		
		cat("\n =========================== \n")
		cat(  "| A N A L Y S I S   D O N E |\n")
		cat(  " =========================== \n")
        
	} else {
    
		cat("\n\t Loading from file ", basename(out.file.base), "_final_result_list.rds", "\n\n", sep = "")
    	res.list <- readRDS(paste(out.dir, "/DESeq2/", out.file.base, "_final_result_list.rds", sep = ""))
		#load(paste(out.dir, "/DESeq2/", out.file.base, "_final_result_list.RData", sep = ""))
		cat("\n\t DONE! \n\n")
        }
        return(res.list)
}
