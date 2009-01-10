
#############################################

# The FunNet library (C) by Corneliu Henegar <corneliu@henegar.info>

#############################################

.packageName <- "FunNet"

.funnet.version <- "1.00-6"

try(Sys.setlocale("LC_ALL", "en_US.utf8"), silent = TRUE)

.First.lib <- function(lib, pkg, ...)
{
  verbose <- .Options$Hverbose
  if(!length(verbose) || verbose)

cat(paste("\nThis is FunNet package ",.funnet.version," built and maintained by Corneliu Henegar.\n",
	"Using Gene Ontology and KEGG annotations updated on ",annot.date,".\n\n",
	"FunNet(wd='', org=c('HS','MM','RN','SC'), two.lists=TRUE, up.frame=NULL, down.frame=NULL,\n",
			"\t genes.frame=NULL, restrict=FALSE, ref.list=NULL, logged=TRUE,\n",
			"\t discriminant=FALSE, go.bp=TRUE, go.cc=TRUE, go.mf=TRUE, kegg=TRUE,\n",
			"\t annot.method=c('specificity','terminological','decorrelated'),\n",
			"\t annot.details=TRUE, direct=FALSE, enriched=TRUE, fdr=5, build.annot.net=TRUE,\n",
			"\t coexp.matrix=NULL, coexp.method=c('spearman','pearson','kendall','euclid'),\n",
			"\t estimate.th=FALSE, hard.th=NA, soft.th=NA, topological = FALSE, keep.sign=FALSE,\n",
			"\t level=NA, annot.clust.method=c('umilds','ucknn'), annot.prox.measure=c('dynamical',\n",
			"\t 'unilat.pond.norm.mean','unilat.norm.sum','norm.sum','pond.norm.mean'),\n",
			"\t test.recovery=FALSE, test.robust=FALSE, replace.annot=NA,\n", 
			"\t build.gene.net=FALSE, gene.clust.method='hclust', gene.net.details=FALSE,\n",
			"\t gene.clusters=NA, alpha=0.05, RV=0.90, sigma=NA, keep.rdata=FALSE, zip=TRUE)\n\n",sep=''))
  invisible()
}


#############################################################################################
#
# 1. Function FunNet() -> Assures the global control of the analytic process
#
#############################################################################################

			
FunNet <- function(wd="", org="HS", two.lists=TRUE, up.frame=NULL, down.frame=NULL,
	  		genes.frame=NULL, restrict=FALSE, ref.list=NULL, logged=TRUE,
	  		discriminant=FALSE, go.bp=TRUE, go.cc=TRUE, go.mf=TRUE, kegg=TRUE,
			annot.method="specificity", annot.details=TRUE, 
	  		direct=FALSE, enriched=TRUE, fdr=NA, build.annot.net=TRUE,
	  		coexp.matrix=NULL, coexp.method="spearman", estimate.th=FALSE, 
	  		hard.th=NA, soft.th=NA, topological = FALSE, keep.sign=FALSE, level=NA, 
	  		annot.clust.method="umilds", annot.prox.measure="unilat.pond.norm.mean",
	  		test.recovery=FALSE, test.robust=FALSE, replace.annot=NA, 
	  		build.gene.net=FALSE, gene.clust.method="hclust", gene.net.details=FALSE,
			gene.clusters=NA, alpha=0.05, RV=0.90, sigma=NA, keep.rdata=FALSE, zip=TRUE){


	# "two.lists" means 2 lists of genes (i.e. UP & DOWN) to be analyzed instead of just one 
	# "restrict" means restrict the enrichment calculus to a list of genes (in oposition with the whole genome by default)
	# "annot.clust.method" can take values "ucknn" or "umilds"
	# "coexp.method" can take values "spearman", "pearson", "kendall" or "euclid" (which uses a rescaled euclidean distance to compute co-expression)
	# "annot.method" can take values "specificity", "terminological", or "decorrelated"
	# "keep.sign" indicate considering separately co- and inverse expression relations among genes

	# keep the evidence of the parameters
	
	parameter.list <- list(analysis.date=date(), package.version=.funnet.version, org=org, annot.date=annot.date, annot.method=annot.method, 
			annot.clust.method=annot.clust.method, annot.prox.measure=annot.prox.measure, direct=direct, enriched=enriched, fdr=fdr, 
			two.lists=two.lists, restrict=restrict, go.bp=go.bp, go.cc=go.cc, go.mf=go.mf, kegg=kegg, discriminant=discriminant, 
			logged=logged, annot.details=annot.details, estimate.th=estimate.th, hard.th=hard.th, soft.th=soft.th, 
			coexp.method=coexp.method, topological = topological, keep.sign=keep.sign, build.annot.net=build.annot.net, level=level,
			test.recovery=test.recovery, test.robust=test.robust, replace.annot=replace.annot, build.gene.net=build.gene.net, 
			gene.clust.method=gene.clust.method, gene.clusters=gene.clusters, gene.net.details=gene.net.details, alpha=alpha, 
			RV=RV, sigma=sigma, keep.rdata=keep.rdata, zip=zip)
	
	# check parameters for errors
	
	.check.parameters(parameter.list, coexp.matrix, up.frame, down.frame, ref.list, genes.frame)
	
	if(discriminant){
		two.lists <- TRUE
		restrict <- TRUE
	}
	
	if(!is.null(up.frame) & !is.null(down.frame) & is.null(genes.frame)){
		two.lists <- TRUE
	}
	
	if(!is.null(genes.frame) & (is.null(up.frame) | is.null(down.frame))){
		two.lists <- FALSE
		discriminant <- FALSE
	}

	
	cat(paste("\n\tFunNet started at: ",date(),sep=""))
	cat(paste("\n\t\tUsing annotations updated on: ",annot.date,sep=""))
	
	if(wd != ""){setwd(wd)} # set working directory
		
	# create results directories
	results.dir <- paste("Results_",format(Sys.time(), "%Y_%b_%d_%H-%M-%S"),sep="")
	dir.create(paste(getwd(),"/",results.dir,sep=""))
	dir.create (paste(getwd(),"/",results.dir,"/html",sep=""))
	dir.create (paste(getwd(),"/",results.dir,"/images",sep=""))
	
	try(write.table(as.matrix(print(parameter.list)),col.names=F,file=paste(getwd(),"/",results.dir,"/","parameters_list.txt",sep=""),sep="\t"))
	
# Load data files	

	wd <- getwd()	# working directory
	#locus.name <- NULL
	#locus.symbol <- NULL

	# EntrezGene ID's	
	if(org == "HS"){locus.name <- HS.locus.name[,1:2]}
	if(org == "MM"){locus.name <- MM.locus.name[,1:2]}
	if(org == "RN"){locus.name <- RN.locus.name[,1:2]}
	if(org == "SC"){locus.name <- SC.locus.name[,1:2]}
	if(org == "HS"){locus.symbol <- HS.locus.name[,c(1,3)]}
	if(org == "MM"){locus.symbol <- MM.locus.name[,c(1,3)]}
	if(org == "RN"){locus.symbol <- RN.locus.name[,c(1,3)]}
	if(org == "SC"){locus.symbol <- SC.locus.name[,c(1,3)]}
	
	rownames(locus.name) <- locus.name[,1]
	rownames(locus.symbol) <- locus.symbol[,1]
	
	# prepare expression data


	if(two.lists){	# when 2 lists of genes are analysed in oposition (UP & DOWN)

		#if(is.null(up.frame)){try(up.frame <- read.table(paste(wd,"/","up.txt",sep=""),sep="\t"))}	# genes UP
		#if(is.null(down.frame)){try(down.frame <- read.table(paste(wd,"/","down.txt",sep=""),sep="\t"))}	# genes DOWN

		up.down <- .filter.genes(up.frame=up.frame,down.frame=down.frame,two.lists=TRUE,locus.name=locus.name,logged=logged)
		up.frame <- up.down$up.frame
		down.frame <- up.down$down.frame
		rm(up.down)
		
		if(discriminant){	# specifing a reference list of genes for the calculation of annotations enrichment
			ref.list <- c(as.character(up.frame[,1]),as.character(down.frame[,1]))
		}else if(restrict & !discriminant){
			ref.list <- .filter.genes(restrict=TRUE,ref.list=ref.list,locus.name=locus.name)
		}else if(!restrict & !discriminant){
			ref.list <- NULL
		}
		

	}else{
		#if(is.null(genes.frame)){genes.frame <- read.table(paste(wd,"/","genes.txt",sep=""),sep="\t")}	# a unique data frame of genes

		genes.frame <- .filter.genes(genes.frame=genes.frame,two.lists=FALSE,locus.name=locus.name)
		
		if(restrict){
			ref.list <- .filter.genes(restrict=TRUE,ref.list=ref.list,locus.name=locus.name)
		}else{
			ref.list <- NULL
		}

	}
	
	# save start-up environment in the working directory as well as in the results directory
	cat(paste("\n\tSaving start-up environment... ",format(Sys.time(), "%X"),sep=""))
	save(up.frame,down.frame,ref.list,genes.frame,parameter.list,coexp.matrix,locus.name,
			file=paste(getwd(),"/",results.dir,"/","start-up_environment.RData",sep=""),compress=T)
	save(up.frame,down.frame,ref.list,genes.frame,parameter.list,coexp.matrix,locus.name,
			file=paste(getwd(),"/","start-up_environment.RData",sep=""),compress=T)
	
	
	if(estimate.th){
	
		datas <- NULL
		
		if(!is.null(up.frame) & !is.null(down.frame)){datas <- rbind(up.frame,down.frame)}
		if(!is.null(genes.frame) & (is.null(up.frame) | is.null(down.frame))){datas <- genes.frame}
		
		rownames(datas) <- datas[,1]
		datas <- datas[,2:ncol(datas)]
		cat("\n\tHard thresholding...\n")
		try(hard.th <- .PickHardThreshold(datExpr1=t(datas),coexp.method=coexp.method))
				
		try(write.table(hard.th$tablou,file = paste(getwd(),"/",results.dir,"/",coexp.method,"_hard_threshold.txt",sep=""),append=FALSE, col.names=TRUE,,row.names=F,sep="\t"))
		try(write(paste("\n\nHard threshold estimate: ",hard.th$estimate,"\n",sep=""),file = paste(getwd(),"/",results.dir,"/",coexp.method,"_hard_threshold.txt",sep=""),append=TRUE))
		try(write(paste("Do not trust this automated estimation without checking it!\n",
				"Do not hesitate to select another threshold depending on the associated connectivity values.\n",
				"Then please restart FunNet interaction analysis with your selected threshold.\n",sep=""),file = paste(getwd(),"/",results.dir,"/",coexp.method,"_hard_threshold.txt",sep=""),append=TRUE))


		cat("\n\tSoft thresholding...\n")
		try(soft.th <- .PickSoftThreshold(datExpr1=t(datas),coexp.method=coexp.method))
				
		try(write.table(soft.th$tablou,file = paste(getwd(),"/",results.dir,"/",coexp.method,"_soft_threshold.txt",sep=""),append=FALSE, col.names=TRUE,,row.names=F,sep="\t"))
		try(write(paste("\n\nSoft threshold estimate: ",soft.th$estimate,"\n",sep=""),file = paste(getwd(),"/",results.dir,"/",coexp.method,"_soft_threshold.txt",sep=""),append=TRUE))
		try(write(paste("Do not trust this automated estimation without checking it!\n",
				"Do not hesitate to select another threshold depending on the associated connectivity values.\n",
				"Then please restart FunNet interaction analysis with your selected threshold.\n",sep=""),file = paste(getwd(),"/",results.dir,"/",coexp.method,"_soft_threshold.txt",sep=""),append=TRUE))


		print("Estimation of the co-expression threshold finished!") 
		print("Please restart FunNet with your chosen threshold.")
		
		if(!keep.rdata){
			try(unlink(paste(getwd(),"/",results.dir,"/",list.files(path=paste(getwd(),"/",results.dir,"/",sep=""),pattern="[:print:]*.RData"),
				sep=""),recursive=TRUE))
		}

		if(zip){
			try(unlink(paste(getwd(),"/",results.dir,"/html",sep=""),recursive=TRUE))
			try(unlink(paste(getwd(),"/",results.dir,"/images",sep=""),recursive=TRUE))
			try(system(command=paste("zip -r9q ",results.dir,".zip ", "./",results.dir,"/*",sep="")))
			try(unlink(paste(getwd(),"/",results.dir,sep=""),recursive=TRUE))
		}
		
		options(show.error.messages=FALSE)
		stop()
	}

	if(is.null(coexp.matrix) & (build.annot.net | build.gene.net)){
		cat(paste("\n\tComputing co-expression matrix... ",format(Sys.time(), "%X"),sep=""))
		
		datas <- NULL
		
		if(!is.null(up.frame) & !is.null(down.frame)){datas <- rbind(up.frame,down.frame)}
		if(!is.null(genes.frame) & (is.null(up.frame) | is.null(down.frame))){datas <- genes.frame}

		rownames(datas) <- datas[,1]
		datas <- datas[,2:ncol(datas)]
		if(coexp.method %in% c("spearman","pearson","kendall")){
			coexp.matrix <- rcorr(t(datas),type=coexp.method)$r
		}else if(coexp.method == "euclid"){
			coexp.matrix <- 1 - (as.matrix(dist(datas))/max(dist(datas),na.rm=TRUE))
		}
		try(save(coexp.matrix,file=paste(getwd(),"/",results.dir,"/","brut_coexp_matrix.RData",sep=""),compress=T))
	}

	
	
	
	# compute the gene adjacency matrix
	
	if(build.annot.net | build.gene.net){
	
		sign.matrix <- coexp.matrix/abs(coexp.matrix)
		coexp.matrix <- abs(coexp.matrix)

		if(annot.clust.method %in% c("umilds","spectral")){

			if(!is.na(hard.th)){
				coexp.matrix[coexp.matrix >= hard.th] <- 1
				coexp.matrix[coexp.matrix < hard.th] <- 0
			}else if(!is.na(soft.th)){
				coexp.matrix <- coexp.matrix^soft.th
			}

			if(topological){
				coexp.matrix <- 1 - .TOMdist(adjmat1=coexp.matrix)
			}else if(keep.sign){
				coexp.matrix <- coexp.matrix * sign.matrix
			}

		}else if(annot.clust.method == "ucknn" & !is.na(hard.th)){
			coexp.matrix[coexp.matrix < hard.th] <- 0
		}

		try(save(coexp.matrix,sign.matrix, parameter.list, file=paste(getwd(),"/",results.dir,"/","gene_adj_matrix.RData",sep="")))
	}
	


# KEGG Annotations
	if(kegg){
		terms.name <- KEGG.terms.name	# KEGG structure
		rownames(terms.name) <- terms.name[,1]
		# KEGG annotations file
		if(org == "HS"){file.annot <- HS.KEGG.file.annot}
		if(org == "MM"){file.annot <- MM.KEGG.file.annot}
		if(org == "RN"){file.annot <- RN.KEGG.file.annot}
		if(org == "SC"){file.annot <- SC.KEGG.file.annot}
		taxoname <- "KEGG"
	
	
		if(two.lists){
			.main.loop(file.annot=file.annot,taxoname=taxoname,annot.method=annot.method,terms.name=terms.name,direct=direct,fdr=fdr,
				go=FALSE,results.dir=results.dir,alpha=alpha,locus.name=locus.name,annot.clust.method=annot.clust.method,
				up.frame=up.frame,down.frame=down.frame,restrict=restrict,ref.list=ref.list,annot.details=annot.details,level=NA,
				build.annot.net=build.annot.net,test.recovery=test.recovery,test.robust=test.robust,replace.annot=replace.annot,
				locus.symbol=locus.symbol,annot.prox.measure=annot.prox.measure,coexp.matrix=coexp.matrix,parameter.list=parameter.list,
				org=org,gene.net.details=gene.net.details,RV=RV,sigma=sigma)
		}else{
			.main.loop(file.annot=file.annot,taxoname=taxoname,annot.method=annot.method,terms.name=terms.name,direct=direct,fdr=fdr,
				go=FALSE,results.dir=results.dir,annot.clust.method=annot.clust.method,alpha=alpha,locus.name=locus.name,
				genes.frame=genes.frame,restrict=restrict,ref.list=ref.list,annot.details=annot.details,level=NA,build.annot.net=build.annot.net,
				test.recovery=test.recovery,test.robust=test.robust,replace.annot=replace.annot,locus.symbol=locus.symbol,
				annot.prox.measure=annot.prox.measure,coexp.matrix=coexp.matrix,parameter.list=parameter.list,org=org,
				gene.net.details=gene.net.details,RV=RV,sigma=sigma)
		}
	}

# Gene Ontology Annotations

	go.name <- c("GO Biological Process","GO Cellular Component","GO Molecular Function") # names of GO branches
	terms.name <- GO.terms.name # GO structure 
	rownames(terms.name) <- terms.name[,1]


	if(go.bp){
	# GO annotations file BP
		if(org == "HS"){file.annot <- HS.GO.DIR.BP.file.annot}
		if(org == "MM"){file.annot <- MM.GO.DIR.BP.file.annot}
		if(org == "RN"){file.annot <- RN.GO.DIR.BP.file.annot}
		if(org == "SC"){file.annot <- SC.GO.DIR.BP.file.annot}

		taxoname <- go.name[1]	# name of the GO branch considered for gene annotations

		if(two.lists == TRUE){
			.main.loop(file.annot=file.annot,taxoname=taxoname,annot.method=annot.method,terms.name=terms.name,direct=direct,fdr=fdr,
				go=TRUE,results.dir=results.dir,alpha=alpha,locus.name=locus.name,annot.clust.method=annot.clust.method,
				up.frame=up.frame,down.frame=down.frame,restrict=restrict,ref.list=ref.list,annot.details=annot.details,level=level,
				build.annot.net=build.annot.net,test.recovery=test.recovery,test.robust=test.robust,replace.annot=replace.annot,
				locus.symbol=locus.symbol,annot.prox.measure=annot.prox.measure,coexp.matrix=coexp.matrix,parameter.list=parameter.list,
				org=org,gene.net.details=gene.net.details,RV=RV,sigma=sigma)
		}else{
			.main.loop(file.annot=file.annot,taxoname=taxoname,annot.method=annot.method,terms.name=terms.name,direct=direct,fdr=fdr,
				go=TRUE,results.dir=results.dir,annot.clust.method=annot.clust.method,alpha=alpha,locus.name=locus.name,
				genes.frame=genes.frame,restrict=restrict,ref.list=ref.list,annot.details=annot.details,level=level,build.annot.net=build.annot.net,
				test.recovery=test.recovery,test.robust=test.robust,replace.annot=replace.annot,locus.symbol=locus.symbol,
				annot.prox.measure=annot.prox.measure,coexp.matrix=coexp.matrix,parameter.list=parameter.list,org=org,
				gene.net.details=gene.net.details,RV=RV,sigma=sigma)
		}		
	}

	if(go.cc){
	# GO annotations file CC
		if(org == "HS"){file.annot <- HS.GO.DIR.CC.file.annot}
		if(org == "MM"){file.annot <- MM.GO.DIR.CC.file.annot}
		if(org == "RN"){file.annot <- RN.GO.DIR.CC.file.annot}
		if(org == "SC"){file.annot <- SC.GO.DIR.CC.file.annot}

		taxoname <- go.name[2]	# name of the GO branch considered for gene annotations

		if(two.lists == TRUE){
			.main.loop(file.annot=file.annot,taxoname=taxoname,annot.method=annot.method,terms.name=terms.name,direct=direct,fdr=fdr,
				go=TRUE,results.dir=results.dir,alpha=alpha,locus.name=locus.name,annot.clust.method=annot.clust.method,
				up.frame=up.frame,down.frame=down.frame,restrict=restrict,ref.list=ref.list,annot.details=annot.details,level=level,
				build.annot.net=build.annot.net,test.recovery=test.recovery,test.robust=test.robust,replace.annot=replace.annot,
				locus.symbol=locus.symbol,annot.prox.measure=annot.prox.measure,coexp.matrix=coexp.matrix,parameter.list=parameter.list,
				org=org,gene.net.details=gene.net.details,RV=RV,sigma=sigma)
		}else{
			.main.loop(file.annot=file.annot,taxoname=taxoname,annot.method=annot.method,terms.name=terms.name,direct=direct,fdr=fdr,
				go=TRUE,results.dir=results.dir,annot.clust.method=annot.clust.method,alpha=alpha,locus.name=locus.name,
				genes.frame=genes.frame,restrict=restrict,ref.list=ref.list,annot.details=annot.details,level=level,build.annot.net=build.annot.net,
				test.recovery=test.recovery,test.robust=test.robust,replace.annot=replace.annot,locus.symbol=locus.symbol,
				annot.prox.measure=annot.prox.measure,coexp.matrix=coexp.matrix,parameter.list=parameter.list,org=org,
				gene.net.details=gene.net.details,RV=RV,sigma=sigma)
		}		
	}

	if(go.mf){
	# GO annotations file MF
		if(org == "HS"){file.annot <- HS.GO.DIR.MF.file.annot}
		if(org == "MM"){file.annot <- MM.GO.DIR.MF.file.annot}
		if(org == "RN"){file.annot <- RN.GO.DIR.MF.file.annot}
		if(org == "SC"){file.annot <- SC.GO.DIR.MF.file.annot}

		taxoname <- go.name[3]	# name of the GO branch considered for gene annotations

		if(two.lists == TRUE){
			.main.loop(file.annot=file.annot,taxoname=taxoname,annot.method=annot.method,terms.name=terms.name,direct=direct,fdr=fdr,
				go=TRUE,results.dir=results.dir,alpha=alpha,locus.name=locus.name,annot.clust.method=annot.clust.method,
				up.frame=up.frame,down.frame=down.frame,restrict=restrict,ref.list=ref.list,annot.details=annot.details,level=level,
				build.annot.net=build.annot.net,test.recovery=test.recovery,test.robust=test.robust,replace.annot=replace.annot,
				locus.symbol=locus.symbol,annot.prox.measure=annot.prox.measure,coexp.matrix=coexp.matrix,parameter.list=parameter.list,
				org=org,gene.net.details=gene.net.details,RV=RV,sigma=sigma)
		}else{
			.main.loop(file.annot=file.annot,taxoname=taxoname,annot.method=annot.method,terms.name=terms.name,direct=direct,fdr=fdr,
				go=TRUE,results.dir=results.dir,annot.clust.method=annot.clust.method,alpha=alpha,locus.name=locus.name,
				genes.frame=genes.frame,restrict=restrict,ref.list=ref.list,annot.details=annot.details,level=level,build.annot.net=build.annot.net,
				test.recovery=test.recovery,test.robust=test.robust,replace.annot=replace.annot,locus.symbol=locus.symbol,
				annot.prox.measure=annot.prox.measure,coexp.matrix=coexp.matrix,parameter.list=parameter.list,org=org,
				gene.net.details=gene.net.details,RV=RV,sigma=sigma)
		}		
	}
	
		
	# build conventional co-expression net if asked to

	if(build.gene.net & !is.null(coexp.matrix)){

		try(clusters <- .build.coexp.net(coexp.matrix=coexp.matrix, locus.name=locus.name, locus.symbol=locus.symbol, 
					gene.clust.method=gene.clust.method, gene.clusters=gene.clusters))

		try(save(clusters,parameter.list, file=paste(getwd(),"/",results.dir,"/","co-expression_clusters.RData",sep=""),compress=T))

		net.matrix <- coexp.matrix
		rownames(net.matrix) <- locus.symbol[rownames(net.matrix),2]
		colnames(net.matrix) <- locus.symbol[colnames(net.matrix),2]

		try(.cyto.sym(net.matrix=net.matrix,file.net=paste(getwd(),"/",results.dir,"/","co-expression_net.txt",sep=""),
			diagonal=FALSE,thresh=NULL))
		rm(net.matrix)

		try(centrality <- .genes.centrality(adj.matrix=coexp.matrix,clusters=clusters,taxoname=taxoname,locus.symbol=locus.symbol,
			results.dir=results.dir,coexp=TRUE))
		if(two.lists){
			up.down <- rbind(matrix(1,nrow(up.frame),1),matrix(0,nrow(down.frame),1))
			rownames(up.down) <- c(as.character(up.frame[,1]),as.character(down.frame[,1]))

			try(write.table(cbind(rownames(clusters$gene.connect),as.vector(locus.symbol[rownames(clusters$gene.connect),2]),
				as.vector(locus.name[rownames(clusters$gene.connect),2]),up.down,clusters$gene.connect,
				centrality[rownames(clusters$gene.connect),]),file=paste(getwd(),"/",results.dir,"/","co-expression_net_info.txt",sep=""),
				sep="\t",col.names=c("geneid","symbol","name","up(1)_down(0)",colnames(clusters$gene.connect),colnames(centrality)),
 				row.names=F))
			try(rm(clusters,up.down,centrality))
		}else{
			try(write.table(cbind(rownames(clusters$gene.connect),as.vector(locus.symbol[rownames(clusters$gene.connect),2]),
				as.vector(locus.name[rownames(clusters$gene.connect),2]),clusters$gene.connect,
				centrality[rownames(clusters$gene.connect),]),file=paste(getwd(),"/",results.dir,"/","co-expression_net_info.txt",sep=""),
				sep="\t",col.names=c("geneid","symbol","name",colnames(clusters$gene.connect),colnames(centrality)),
				row.names=F))
			try(rm(clusters,centrality))
		}


		cat(paste("\n\tCo-expression net building finished... ",date(),sep=""))
		rm()

	}
	
	if(!keep.rdata){
		try(unlink(paste(getwd(),"/",results.dir,"/",list.files(path=paste(getwd(),"/",results.dir,"/",sep=""),pattern="[:print:]*.RData"),
			sep=""),recursive=TRUE))
	}
	
	if(zip){
		try(system(command=paste("zip -r9q ",results.dir,".zip ", "./",results.dir,"/*",sep="")))
		try(unlink(paste(getwd(),"/",results.dir,sep=""),recursive=TRUE))
	}

	cat(paste("\n\tEnd  of treatment at: ",date(),"\n",sep=""))	
	rm()


}


#############################################################################################
#
# 2. Function .check.parameters() -> Verifies the correctness of the FunNet parameters
#
#############################################################################################


.check.parameters <- function(parameter.list,coexp.matrix,up.frame,down.frame,ref.list,genes.frame){

	cat(paste("\n\tChecking parameters for inconsistencies... \n",sep=""))	

	if(!(parameter.list$org %in% c("HS","MM","RN","SC"))){
		cat(paste("\n\t\tParameter org = '",parameter.list$org,"' is incorrect!\n",sep=""))
		cat(paste("\n\t\tPossible values: 'HS' | 'MM' | 'RN' | 'SC'\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("org")
	}

	if(!(parameter.list$annot.clust.method %in% c("umilds","ucknn","spectral"))){
		cat(paste("\n\t\tParameter annot.clust.method = '",parameter.list$annot.clust.method,"' is incorrect!\n",sep=""))
		cat(paste("\n\t\tPossible values: 'umilds' | 'ucknn'\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("annot.clust.method")
	}

	if(!(parameter.list$annot.prox.measure %in% c("dynamical","unilat.pond.norm.mean","unilat.norm.sum","norm.sum","pond.norm.mean"))){
		cat(paste("\n\t\tParameter annot.prox.measure = '",parameter.list$annot.prox.measure,"' is incorrect!\n",sep=""))
		cat(paste("\n\t\tPossible values: 'norm.sum' | 'unilat.norm.sum' | 'pond.norm.mean' | 'unilat.pond.norm.mean'\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("annot.prox.measure")
	}

	if(!(parameter.list$annot.method %in% c("specificity","terminological","decorrelated"))){
		cat(paste("\n\t\tParameter annot.method = '",parameter.list$annot.method,"' is incorrect!\n",sep=""))
		cat(paste("\n\t\tPossible values: 'specificity' | 'terminological' | 'decorrelated'\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("annot.method")
	}

	if(!is.logical(parameter.list$direct)){
		cat(paste("\n\t\tParameter direct = ",parameter.list$direct," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should have a logical value.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("direct")
	}
	
	if(!is.logical(parameter.list$enriched)){
		cat(paste("\n\t\tParameter enriched = ",parameter.list$enriched," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should have a logical value.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("enriched")
	}

	if(!is.logical(parameter.list$go.bp)){
		cat(paste("\n\t\tParameter go.bp = ",parameter.list$go.bp," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should have a logical value.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("go.bp")
	}

	if(!is.logical(parameter.list$go.cc)){
		cat(paste("\n\t\tParameter go.cc = ",parameter.list$go.cc," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should have a logical value.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("go.cc")
	}

	if(!is.logical(parameter.list$go.mf)){
		cat(paste("\n\t\tParameter go.mf = ",parameter.list$go.mf," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should have a logical value.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("go.mf")
	}

	if(!is.logical(parameter.list$kegg)){
		cat(paste("\n\t\tParameter kegg = ",parameter.list$kegg," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should have a logical value.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("kegg")
	}
	
	if(!is.na(parameter.list$fdr) && (parameter.list$fdr < 0 | parameter.list$fdr > 30)){
		cat(paste("\n\t\tParameter fdr = ",parameter.list$fdr," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should be NA or take numerical values in between 0 and 30.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("fdr")
	}
	
	if(!is.na(parameter.list$hard.th) && (parameter.list$hard.th < 0 | parameter.list$hard.th > 1)){
		cat(paste("\n\t\tParameter hard.th = ",parameter.list$hard.th," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should be NA or take numerical values in between 0 and 1.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("hard.th")
	}
	
	if(!is.na(parameter.list$soft.th) && (parameter.list$soft.th < 1 | parameter.list$soft.th > 20)){
		cat(paste("\n\t\tParameter soft.th = ",parameter.list$soft.th," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should be NA or take numerical values in between 1 and 20.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("soft.th")
	}
	
	if(!is.logical(parameter.list$two.lists)){
		cat(paste("\n\t\tParameter two.lists = ",parameter.list$two.lists," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should have a logical value.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("two.lists")
	}
	
	if(!is.logical(parameter.list$restrict)){
		cat(paste("\n\t\tParameter restrict = ",parameter.list$restrict," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should have a logical value.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("restrict")
	}

	if(!is.logical(parameter.list$annot.details)){
		cat(paste("\n\t\tParameter annot.details = ",parameter.list$annot.details," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should have a logical value.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("annot.details")
	}
	
	if(!is.logical(parameter.list$estimate.th)){
		cat(paste("\n\t\tParameter estimate.th = ",parameter.list$estimate.th," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should have a logical value.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("estimate.th")
	} 

	if(!(parameter.list$coexp.method %in% c("euclid","spearman","pearson","kendall"))){
		cat(paste("\n\t\tParameter coexp.method = '",parameter.list$coexp.method,"' is incorrect!\n",sep=""))
		cat(paste("\n\t\tPossible values: 'euclid' | 'spearman' | 'pearson' | 'kendall'\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("coexp.method")
	}

	if(!is.logical(parameter.list$discriminant)){
		cat(paste("\n\t\tParameter discriminant = ",parameter.list$discriminant," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should have a logical value.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("discriminant")
	}
	
	if(parameter.list$alpha < 0 || parameter.list$alpha > 1){
		cat(paste("\n\t\tParameter alpha = ",parameter.list$alpha," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should take numerical values in between 0 and 1.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("alpha")
	}

	if(!is.logical(parameter.list$build.annot.net)){
		cat(paste("\n\t\tParameter build.annot.net = ",parameter.list$build.annot.net," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should have a logical value.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("build.annot.net")
	}
	
	if(!is.na(parameter.list$level) && (parameter.list$level < 1 | parameter.list$level > 16)){
		cat(paste("\n\t\tParameter level = ",parameter.list$level," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should be NA or take numerical values in between 1 and 16.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("level")
	}
	
	if(!is.logical(parameter.list$topological)){
		cat(paste("\n\t\tParameter topological = ",parameter.list$topological," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should have a logical value.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("topological")
	}
	
	if(!is.logical(parameter.list$test.recovery)){
		cat(paste("\n\t\tParameter test.recovery = ",parameter.list$test.recovery," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should have a logical value.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("test.recovery")
	}
	
	if(!is.logical(parameter.list$test.robust)){
		cat(paste("\n\t\tParameter test.robust = ",parameter.list$test.robust," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should have a logical value.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("test.robust")
	}
	
	if(!is.na(parameter.list$replace.annot) && (parameter.list$replace.annot < 0 | parameter.list$replace.annot >= 100)){
		cat(paste("\n\t\tParameter replace.annot = ",parameter.list$replace.annot," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should be NA or take numerical values in between 0 and 100.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("replace.annot")
	}
	
	if(!is.logical(parameter.list$keep.sign)){
		cat(paste("\n\t\tParameter keep.sign = ",parameter.list$keep.sign," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should have a logical value.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("keep.sign")
	}
		
	if(!is.logical(parameter.list$build.gene.net)){
		cat(paste("\n\t\tParameter build.gene.net = ",parameter.list$build.gene.net," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should have a logical value.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("build.gene.net")
	}
	
	if(parameter.list$gene.clust.method != "hclust"){
		cat(paste("\n\t\tParameter gene.clust.method = '",parameter.list$gene.clust.method,"' is incorrect!\n",sep=""))
		cat(paste("\n\t\tPossible values: 'hclust'\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("gene.clust.method")
	}
	
	if(!is.na(parameter.list$gene.clusters) && !is.numeric(parameter.list$gene.clusters)){
		cat(paste("\n\t\tParameter gene.clusters = ",parameter.list$gene.clusters," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should be NA or take integer values smaller than the total number of analyzed genes.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("gene.clusters")
	}
	
	if(!is.logical(parameter.list$logged)){
		cat(paste("\n\t\tParameter logged = ",parameter.list$logged," is incorrect!",sep=""))
		cat(paste("\n\t\tThis parameter should have a logical value.",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("logged")
	}
	
	if(!is.logical(parameter.list$gene.net.details)){
		cat(paste("\n\t\tParameter gene.net.details = ",parameter.list$gene.net.details," is incorrect!",sep=""))
		cat(paste("\n\t\tThis parameter should have a logical value.",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("gene.net.details")
	}
	
	if(parameter.list$RV < 0 || parameter.list$RV > 1){
		cat(paste("\n\t\tParameter RV = ",parameter.list$RV," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should take numerical values in between 0 and 1.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("RV")
	}
	
	if(!is.na(parameter.list$sigma) && !is.numeric(parameter.list$sigma)){
		cat(paste("\n\t\tParameter sigma = ",parameter.list$sigma," is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should take numerical values.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("sigma")
	}
	
	if(!is.null(coexp.matrix) && !(is.matrix(coexp.matrix) & ncol(coexp.matrix) == nrow(coexp.matrix)
		& max(coexp.matrix) <= 1 & min(coexp.matrix) >= -1)){
		cat(paste("\n\t\tParameter coexp.matrix is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should be NULL or a square matrix containing genes co-expression similarity.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("coexp.matrix")
	}
		
	if(!is.null(up.frame) && !is.data.frame(up.frame)){
		cat(paste("\n\t\tParameter up.frame is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should be NULL or a data.frame with the expression levels of the up-regulated genes.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("up.frame")
	}
		
	if(!is.null(down.frame) && !is.data.frame(down.frame)){
		cat(paste("\n\t\tParameter down.frame is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should be NULL or a data.frame with the expression levels of the down-regulated genes.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("down.frame")
	}
			
	if(!is.null(ref.list) && !is.data.frame(ref.list)){
		cat(paste("\n\t\tParameter ref.list is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should be NULL or a data.frame with the reference list of the analyzed genes.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("ref.list")
	}
	
	if(!is.null(genes.frame) && !is.data.frame(genes.frame)){
		cat(paste("\n\t\tParameter genes.frame is incorrect!\n",sep=""))
		cat(paste("\n\t\tThis parameter should be NULL or a data.frame with the expression levels of the analyzed genes.\n",sep=""))
		cat(paste("\n\t\tPlease check again and restart FunNet.\n",sep=""))
		stop("genes.frame")
	}

}	




#############################################################################################
#
# 3. Function .main.loop() -> Controls in detail the analytic algorithm
#
#############################################################################################

.main.loop <- function(file.annot, taxoname, annot.method, terms.name, direct=FALSE, fdr=FALSE,
			go=FALSE, results.dir, alpha, locus.name, locus.symbol, annot.clust.method,
			up.frame=NULL, down.frame=NULL, genes.frame=NULL, restrict=NULL, ref.list=NULL, 
			annot.prox.measure=NULL, annot.details=FALSE, level=NA, 
			build.annot.net=FALSE, test.recovery=FALSE, test.robust=FALSE, 
			replace.annot=NA, coexp.matrix=NULL, parameter.list, org="HS",
			gene.net.details=FALSE, RV, sigma){

annot.matrix <- NULL
up.annot.matrix <- NULL
down.annot.matrix <- NULL
genes.annot.matrix <- NULL


# annotate genes

if(!is.null(up.frame) & !is.null(down.frame) & is.null(genes.frame)){
	
	try(up.annot <- .annotate.for.net(file.annot=file.annot,fdr=fdr,go=go,direct=direct,annot.method=annot.method,taxoname=taxoname,
			terms.name=terms.name,restrict=restrict,genes.frame=up.frame,ref.list=ref.list,nom="UP"))
	try(down.annot <- .annotate.for.net(file.annot=file.annot,fdr=fdr,go=go,direct=direct,annot.method=annot.method,taxoname=taxoname,
			terms.name=terms.name,restrict=restrict,genes.frame=down.frame,ref.list=ref.list,nom="DOWN"))
	try(save(up.annot, down.annot, terms.name, parameter.list, file=paste(getwd(),"/",results.dir,"/",taxoname,"_annot_data.RData",sep="")),
			silent=FALSE)
		

	if(annot.details & exists("down.annot") & exists("up.annot")){
	
		try(.resmat(annot.matrix=up.annot$annot.matrix,print.data=up.annot$print.data,go=go,terms.name=terms.name,taxoname=taxoname,
			nom="UP",locus.name=locus.name,results.dir=results.dir,bgcolor="red",org=org))
		try(.resmat(annot.matrix=down.annot$annot.matrix,print.data=down.annot$print.data,go=go,terms.name=terms.name,taxoname=taxoname,
			nom="DOWN",locus.name=locus.name,results.dir=results.dir,bgcolor="green",org=org))
			
		try(.profil.plot.two(up.annot=up.annot, down.annot=down.annot, terms.name=terms.name, taxoname=taxoname, 
			dev ="pdf", extra=paste(getwd(),"/",results.dir,"/images/",sep="")))
		try(.profil.plot.two(up.annot=up.annot, down.annot=down.annot, terms.name=terms.name, taxoname=taxoname, 
			dev ="png", extra=paste(getwd(),"/",results.dir,"/images/",sep="")))
	}

	if(build.annot.net & exists("down.annot") & exists("up.annot")){		
		
		
		if(go & !is.na(level) & annot.method != "decorrelated"){
			if(level > min(length(up.annot$annot.matrix),length(down.annot$annot.matrix))){level <- 1}
			up.annot.matrix <- up.annot$annot.matrix[[level]]
			down.annot.matrix <- down.annot$annot.matrix[[level]]
		}else if(go & is.na(level) & annot.method != "decorrelated"){
			up.annot.matrix <- up.annot$annot.matrix[[1]]
			down.annot.matrix <- down.annot$annot.matrix[[1]]
		}else{
			up.annot.matrix <- up.annot$annot.matrix
			down.annot.matrix <- down.annot$annot.matrix
		}
		
		# bind together UP and DOWN annotating categories
		
		if(!is.null(up.annot.matrix) & !is.null(down.annot.matrix)){
			up.sup <- matrix(0,nrow(up.annot.matrix),ncol(down.annot.matrix))
			colnames(up.sup) <- colnames(down.annot.matrix)
			rownames(up.sup) <- rownames(up.annot.matrix)

			down.sup <- matrix(0,nrow(down.annot.matrix),ncol(up.annot.matrix))
			colnames(down.sup) <- colnames(up.annot.matrix)
			rownames(down.sup) <- rownames(down.annot.matrix)
			up.annot <- cbind(up.annot.matrix,up.sup)
			down.annot <- cbind(down.sup,down.annot.matrix)
			
			doubles <- rownames(up.annot[rownames(up.annot) %in% rownames(down.annot),])	# search for double annotations
					
			annot.matrix <- rbind(up.annot[!(rownames(up.annot) %in% rownames(down.annot)),],down.annot[!(rownames(down.annot) %in% rownames(up.annot)),])
			rownames(annot.matrix) <- c(rownames(up.annot)[!(rownames(up.annot) %in% rownames(down.annot))],rownames(down.annot)[!(rownames(down.annot) %in% rownames(up.annot))])

			if(length(doubles) > 0){
				for(i in 1:length(doubles)){
					x <- up.annot[doubles[i],] + up.annot[doubles[i],]
					x[x>1] <- 1
					annot.matrix <- rbind(annot.matrix,x)
					rownames(annot.matrix)[nrow(annot.matrix)] <- doubles[i]
				}
			}
		
		}else if(!is.null(up.annot.matrix) & is.null(down.annot.matrix)){
			annot.matrix <- up.annot.matrix
		}else if(is.null(up.annot.matrix) & !is.null(down.annot.matrix)){
			annot.matrix <- down.annot.matrix
		}
		
		if(!is.null(annot.matrix)){annot.matrix <- t(annot.matrix)}
	
	}


}else if(!is.null(genes.frame)){


	try(genes.annot <- .annotate.for.net(file.annot=file.annot,fdr=fdr,go=go,direct=direct,annot.method=annot.method,taxoname=taxoname,
			terms.name=terms.name,restrict=restrict,genes.frame=genes.frame,ref.list=ref.list,nom="GENES"))
	try(save(genes.annot, terms.name, parameter.list, file=paste(getwd(),"/",results.dir,"/",taxoname,"_annot_data.RData",sep="")),silent=FALSE)
		

	if(annot.details & exists("genes.annot")){
	
		try(.resmat(annot.matrix=genes.annot$annot.matrix,print.data=genes.annot$print.data,go=go,terms.name=terms.name,taxoname=taxoname,
			nom="GENES",locus.name=locus.name,results.dir=results.dir,bgcolor="blue",org=org))
		try(.profil.plot.one(genes.annot=genes.annot,terms.name=terms.name, taxoname=taxoname, 
			dev ="pdf", extra=paste(getwd(),"/",results.dir,"/images/",sep="")))
		try(.profil.plot.one(genes.annot=genes.annot, terms.name=terms.name, taxoname=taxoname, 
			dev ="png", extra=paste(getwd(),"/",results.dir,"/images/",sep="")))
	}

	if(build.annot.net & exists("genes.annot")){
		
		if(go & !is.na(level) & annot.method != "decorrelated"){
			if(level > length(genes.annot$annot.matrix)){level <- 1}
			genes.annot.matrix <- genes.annot$annot.matrix[[level]]
		}else if(go & is.na(level) & annot.method != "decorrelated"){
			genes.annot.matrix <- genes.annot$annot.matrix[[1]]
		}else{
			genes.annot.matrix <- genes.annot$annot.matrix
		}
		
		annot.matrix <- t(genes.annot.matrix)
	}
}


# cluster annotations

if(build.annot.net & !is.null(annot.matrix)){

	clusters <- NULL
	
	# dynamic spectral
	if(annot.clust.method == "umilds"){
		if(test.robust & !is.na(replace.annot)){
		
			try(robust <- .annot.robustness(annot.matrix=annot.matrix, adj.matrix=coexp.matrix, terms.name=terms.name, 
					RV=RV, method="sum", annot.prox.measure="dynamical", taxoname=taxoname, spectral=TRUE, 
					replace.annot=replace.annot, sigma=sigma))
			try(save(robust, parameter.list, file=paste(getwd(),"/",results.dir,"/",taxoname,"_annot_robustness_test.RData",sep="")))
			
		}else if(test.recovery & !is.na(replace.annot)){
		
			try(recovery <- .annot.recovery(annot.matrix=annot.matrix, adj.matrix=coexp.matrix, RV=RV, method="sum", 
					replace.annot=replace.annot))
			try(save(recovery, parameter.list, file=paste(getwd(),"/",results.dir,"/",taxoname,"_annot_recovery_test.RData",sep="")))
			
		}else{
		
			try(proximity <- .dynamic.proximity(annot.matrix=annot.matrix, adj.matrix=coexp.matrix, RV=RV, 
					method="sum", annotated.only=TRUE, sigma=sigma))
			try(save(proximity,file=paste(getwd(),"/",results.dir,"/",taxoname,"_proximity_matrix.RData",sep=""),compress=T))
			
			try(clusters <- .spectral.cluster(annot.proximity=proximity$annot.proximity, terms.name=terms.name, 
					adj.matrix=coexp.matrix, marker.matrix=proximity$marker.matrix[[length(proximity$marker.matrix)]],					
					annot.matrix=t(annot.matrix), iter.max=20000))

			try(save(clusters,parameter.list, file=paste(getwd(),"/",results.dir,"/",taxoname,"_spectral_clusters.RData",sep=""),compress=T))
		}		
	
	# spectral with other proximity measures
	}else if(annot.clust.method == "spectral" & annot.prox.measure != "dynamical"){
		
		if(test.robust & !is.na(replace.annot)){

			try(robust <- .annot.robustness(annot.matrix=annot.matrix, adj.matrix=coexp.matrix, terms.name=terms.name, 
					RV=RV, method="sum", annot.prox.measure=annot.prox.measure,taxoname=taxoname,spectral=TRUE, 
					replace.annot=replace.annot, sigma=sigma))
			try(save(robust,parameter.list, file=paste(getwd(),"/",results.dir,"/",taxoname,"_annot_robustness_test.RData",sep="")))
		}else{			
			try(annot.distance <- .annotation.distance(annot.matrix=t(annot.matrix), coexp.matrix=coexp.matrix, measure=annot.prox.measure))
			try(proximity <- .presumptive.proximity(annot.distance=annot.distance, sigma=sigma))
			
			try(save(proximity,file=paste(getwd(),"/",results.dir,"/",taxoname,"_proximity_matrix.RData",sep=""),compress=T))
			
			try(clusters <- .spectral.cluster(annot.proximity=proximity, terms.name=terms.name, 
					adj.matrix=coexp.matrix, marker.matrix=NULL, annot.matrix=t(annot.matrix), 
					iter.max=20000))
			try(save(clusters,parameter.list, file=paste(getwd(),"/",results.dir,"/",taxoname,"_spectral_clusters.RData",sep=""),compress=T))
		}
	
	# knn with asymetric proximity measures
	}else if(annot.clust.method == "ucknn"){

		annot.matrix <- t(annot.matrix)

		if(test.robust & !is.na(replace.annot)){

			try(robust <- .annot.robustness(annot.matrix=t(annot.matrix), adj.matrix=coexp.matrix, terms.name=terms.name, 
					RV=RV, method="sum", annot.prox.measure=annot.prox.measure,taxoname=taxoname,spectral=FALSE, 
					replace.annot=replace.annot))
			try(save(robust,parameter.list, file=paste(getwd(),"/",results.dir,"/",taxoname,"_annot_robustness_test.RData",sep="")))
		}else{			

			try(clusters <- .uc.knn(annot.matrix=annot.matrix, coexp.matrix=coexp.matrix, alpha = 0.05, taxoname=taxoname, 
					annot.prox.measure = annot.prox.measure, normal = TRUE,terms.name=terms.name))
			try(save(clusters, parameter.list, file=paste(getwd(),"/",results.dir,"/",taxoname,"_knn_clusters.RData",sep=""),compress=T))
		}

	}
	
	
# save clustering results

	if(!is.null(clusters)){
		if(exists("down.annot") & exists("up.annot")){

			try(.cytoscape.two(clusters, file.net=paste(getwd(),"/",results.dir,"/",taxoname,"_annotations_net.txt",sep=""),
				file.net.info=paste(getwd(),"/",results.dir,"/",taxoname,"_annotations_net_info.txt",sep=""),
				file.param=paste(getwd(),"/",results.dir,"/",taxoname,"_annotations_net_proximity.txt",sep=""),
				up.annot.matrix=up.annot.matrix,down.annot.matrix=down.annot.matrix,annot.clust.method=annot.clust.method))
			
			try(.central.plot.two(up.annot=up.annot.matrix, down.annot=down.annot.matrix, clusters=clusters, taxoname=taxoname, 
				dev ="png", extra=paste(getwd(),"/",results.dir,"/images/",sep=""), 
				annot.clust.method=annot.clust.method))
			try(.central.plot.two(up.annot=up.annot.matrix, down.annot=down.annot.matrix, clusters=clusters, taxoname=taxoname, 
				dev ="pdf", extra=paste(getwd(),"/",results.dir,"/images/",sep=""), 
				annot.clust.method=annot.clust.method))
			
			if(gene.net.details){
				up.down <- rbind(matrix(1,nrow(up.frame),1),matrix(0,nrow(down.frame),1))
				rownames(up.down) <- c(up.frame[,1],down.frame[,1])

				try(centrality <- .genes.centrality(adj.matrix=coexp.matrix,clusters=clusters,taxoname=taxoname,locus.symbol=locus.symbol,
					results.dir=results.dir))

				try(write.table(cbind(rownames(clusters$gene.connect),as.vector(locus.symbol[rownames(clusters$gene.connect),2]),
					as.vector(locus.name[rownames(clusters$gene.connect),2]),up.down,clusters$gene.connect,
					centrality[rownames(clusters$gene.connect),]),
					file=paste(getwd(),"/",results.dir,"/",taxoname,"_genes_net_info.txt",sep=""),
					sep="\t",col.names=c("geneid","symbol","name","up(1)_down(0)",colnames(clusters$gene.connect),colnames(centrality)),
					row.names=F))
					
				try(rm(up.down,centrality))
			}
			
			try(rm(clusters))

		}else if(exists("genes.annot")){

			try(.cytoscape.one(clusters, file.net=paste(getwd(),"/",results.dir,"/",taxoname,"_annotations_net.txt",sep=""),
				file.net.info=paste(getwd(),"/",results.dir,"/",taxoname,"_annotations_net_info.txt",sep=""),
				file.param=paste(getwd(),"/",results.dir,"/",taxoname,"_annotations_net_proximity.txt",sep=""),
				annot.clust.method=annot.clust.method))

			try(.central.plot.one(genes.annot=genes.annot.matrix, clusters=clusters, taxoname=taxoname, 
				dev ="png", extra=paste(getwd(),"/",results.dir,"/images/",sep=""), 
				annot.clust.method=annot.clust.method))
			try(.central.plot.one(genes.annot=genes.annot.matrix, clusters=clusters, taxoname=taxoname, 
				dev ="pdf", extra=paste(getwd(),"/",results.dir,"/images/",sep=""), 
				annot.clust.method=annot.clust.method))

			if(gene.net.details){
				try(centrality <- .genes.centrality(adj.matrix=coexp.matrix,clusters=clusters,taxoname=taxoname,locus.symbol=locus.symbol,
					results.dir=results.dir))

				try(write.table(cbind(rownames(clusters$gene.connect),as.vector(locus.symbol[rownames(clusters$gene.connect),2]),
					as.vector(locus.name[rownames(clusters$gene.connect),2]),clusters$gene.connect,
					centrality[rownames(clusters$gene.connect),]),
					file=paste(getwd(),"/",results.dir,"/",taxoname,"_genes_net_info.txt",sep=""),
					sep="\t",col.names=c("geneid","symbol","name",colnames(clusters$gene.connect),colnames(centrality)),
					row.names=F))
				try(rm(centrality))
			}
			
			try(rm(clusters))
					
		}	
	}
	
	cat(paste("\n\t",taxoname," themes proximity net building finished... ",date(),sep=""))
	rm()

}


}




#############################################################################################
#
# 4. Function .build.coexp.net() -> Builds a "conventional" gene co-expression net if requested
#
#############################################################################################

.build.coexp.net <- function(coexp.matrix, locus.name, locus.symbol, gene.clust.method="hclust", gene.clusters=NA){
	
	cat(paste("\n\tBuilding conventional gene co-expression net... ",format(Sys.time(),"%X"),sep=""))	
	
	dist.genes <- as.dist(1 - coexp.matrix)
	
	best.partition <- NULL
	sil.part <- NULL
	sil.cluster <- NULL
	sil <- NULL
	cluster.length <- NULL
	best.index <- NULL
	gene.cluster <- NULL
	gene.connect <- NULL
	
	if(gene.clust.method == "hclust"){
		
		hc <- hclust(dist.genes, method = "ward")
		.sil.hc <- as.vector(NULL)	# vector of silhouettes

		# calculate cluster silhouettes

		for(i in 2:(ncol(coexp.matrix)-1)){ 
			# cut the tree 
			memb <- cutree(hc, k = i)
			sil <- silhouette(memb, dist.genes)
			.sil.hc <- c(.sil.hc, mean(summary(sil)$clus.avg.width))
		}

		names(.sil.hc) <- 2:(ncol(coexp.matrix)-1)	# choose the max silhouette

		# find the best partition 

		best.index <- as.numeric(names(.sil.hc[.sil.hc == max(.sil.hc)]))

		cat(paste("\n\t\tBest index: ",best.index,sep=""))
		
		if(!is.na(gene.clusters)){
			best.partition <- cutree(hc, k = gene.clusters)
		}else{
			best.partition <- cutree(hc, k = best.index)
		}

		sil <- silhouette(best.partition, dist.genes)
		sil.part <- mean(summary(sil)$clus.avg.width)
		sil.cluster <- summary(sil)$clus.avg.width
		sil <- .sil.hc
	
	}
	
	if(!is.null(best.partition)){
		gene.cluster <- list(NULL)

		for(i in 1:max(best.partition)){

			gene.cluster[[i]] <- as.character(names(best.partition[best.partition == i]))
			cluster.length <- c(cluster.length, length(gene.cluster[[i]]))

		}
		
		gene.connect <- matrix(0,nrow(coexp.matrix),length(gene.cluster)+2)
		rownames(gene.connect) <- rownames(coexp.matrix)
		colnames(gene.connect) <- c("co-expression_module",as.character(1:length(gene.cluster)),"total_net")

		for(i in 1:length(gene.cluster)){
			gene.connect[as.character(gene.cluster[[i]]),"co-expression_module"] <- i

		}

		for(i in 1:nrow(gene.connect)){

			for(j in 1:length(gene.cluster)){
				x <- coexp.matrix[rownames(gene.connect)[i],as.character(gene.cluster[[j]])]
				gene.connect[i,as.character(j)] <- sum(x)

			}

			gene.connect[i,"total_net"] <- sum(gene.connect[i,as.character(1:length(gene.cluster))])

		}

		for(i in 1:length(gene.cluster)){
			lngth <- sqrt(sum(gene.connect[,as.character(i)] * gene.connect[,as.character(i)]))
			gene.connect[,as.character(i)] <- gene.connect[,as.character(i)]/lngth
			rm(lngth)
		}

		# total_net
		lngth <- sqrt(sum(gene.connect[,"total_net"] * gene.connect[,"total_net"]))
		gene.connect[,"total_net"] <- gene.connect[,"total_net"]/lngth
		rm(lngth)

		
		clusters <- list(dist.genes=dist.genes,coexp.matrix=coexp.matrix,sil.part=sil.part,sil.cluster=sil.cluster,
				sil=sil,cluster.length=cluster.length,best.index=best.index,gene.cluster=gene.cluster,
				best.partition=best.partition,gene.connect=gene.connect,locus.name=locus.name,
				locus.symbol=locus.symbol)
	}else{
		clusters <- NULL
		print("No clusters were created...")
	}
	
	return(clusters)

}


#############################################################################################
#
# 5. Function .filter.genes() -> Filtering provided expression profiles
#
#############################################################################################


.filter.genes <- function(up.frame=NULL,down.frame=NULL,genes.frame=NULL,two.lists=TRUE,restrict=FALSE,ref.list=NULL,locus.name,logged=FALSE){
	
	if(restrict == TRUE && !is.null(ref.list)){
		
		cat(paste("\n\tFiltering reference list started... ",format(Sys.time(),"%X"),sep=""))
		
		genes.locus <- NULL
		if(!is.data.frame(ref.list)){ref.list <- as.data.frame(ref.list)}
		
		for(i in 1:nrow(ref.list)){
			if(regexpr(";",as.character(ref.list[i,1])) != -1){
				x <- strsplit(gsub(" ","",as.character(ref.list[i,1])),";")[[1]]
							
				for(j in 1:length(x)){
					genes.locus <- c(genes.locus,x[j])				
				}
			}else{
				genes.locus <- c(genes.locus,as.character(ref.list[i,1]))				
			}		
		}
		
		ref.list <- unique(as.character(genes.locus))
		ref.list <- ref.list[ref.list != ""]
		names(ref.list) <- ref.list
		ref.list <- ref.list[ref.list %in% as.character(locus.name[,1])]
		
		return(ref.list)
	}



if(two.lists == FALSE && !is.null(genes.frame)){	

	cat(paste("\n\tFiltering genes list started... ",format(Sys.time(), "%X"),sep=""))

	
	genes.matrix <- NULL
	genes.locus <- NULL
	genes.colnames <- colnames(genes.frame)
	
	for(i in 1:nrow(genes.frame)){
		if(regexpr(";",as.character(genes.frame[i,1])) != -1){
			x <- strsplit(gsub(" ","",as.character(genes.frame[i,1])),";")[[1]]
						
			for(j in 1:length(x)){
				genes.locus <- c(genes.locus,x[j])				
				genes.matrix <- rbind(genes.matrix,as.double(genes.frame[i,2:ncol(genes.frame)]))	
			}
		}else{
			genes.locus <- c(genes.locus,as.character(genes.frame[i,1]))				
			genes.matrix <- rbind(genes.matrix,as.double(genes.frame[i,2:ncol(genes.frame)]))
		}		
	}
	

	genes.frame <- data.frame(genes.locus,genes.matrix)
	
	genes.locus <- unique(as.character(genes.locus))
	genes.locus <- genes.locus[genes.locus %in% as.character(locus.name[,1])]
	
	genes.frame.new <- NULL
	
	for(i in 1:length(genes.locus)){
		
		datas <- genes.frame[as.character(genes.frame[,1]) == genes.locus[i],2:ncol(genes.frame)]
		
		if(logged & !is.vector(datas)){
			datas <- apply(datas,2,function(x) 2^as.double(x))
		}else if(logged & is.vector(datas)){
			datas <- 2^as.double(datas)
		}
			
		if(!is.vector(datas)){
			datas <- apply(datas,2,function(x) mean(as.double(x),na.rm=TRUE))
			
		}else{
			datas <- as.double(datas)
		}			
		
		if(logged){
			datas <- log(datas,2)
		}
		
		datas[is.nan(datas)] <- NA
		
		genes.frame.new <- rbind(genes.frame.new,datas)
		
	}
	
	genes.frame.new <- data.frame(genes.locus,genes.frame.new)
	
	genes.frame <- NULL
	
	for(i in 2:ncol(genes.frame.new)){
		
		genes.frame <- cbind(genes.frame,as.numeric(genes.frame.new[,i]))
	
	}
	
	genes.frame <- data.frame(genes.frame.new[,1],genes.frame)
	colnames(genes.frame) <- genes.colnames
	
	return(genes.frame)



}else if(two.lists == TRUE && !is.null(up.frame) && !is.null(down.frame)){
	
	cat(paste("\n\tFiltering genes UP/DOWN started... ",format(Sys.time(), "%X"),sep=""))
	
	up.matrix <- NULL
	up.locus <- NULL
	up.colnames <- colnames(up.frame)
	
	for(i in 1:nrow(up.frame)){
		if(regexpr(";",as.character(up.frame[i,1])) != -1){
			x <- strsplit(gsub(" ","",as.character(up.frame[i,1])),";")[[1]]
						
			for(j in 1:length(x)){
				up.locus <- c(up.locus,x[j])				
				up.matrix <- rbind(up.matrix,as.double(up.frame[i,2:ncol(up.frame)]))	
			}
		}else{
			up.locus <- c(up.locus,as.character(up.frame[i,1]))				
			up.matrix <- rbind(up.matrix,as.double(up.frame[i,2:ncol(up.frame)]))
		}		
	}

	up.frame <- data.frame(as.character(up.locus),up.matrix)

	down.matrix <- NULL
	down.locus <- NULL
	down.colnames <- colnames(down.frame)
		
	for(i in 1:nrow(down.frame)){
		if(regexpr(";",as.character(down.frame[i,1])) != -1){
			x <- strsplit(gsub(" ","",as.character(down.frame[i,1])),";")[[1]]
							
			for(j in 1:length(x)){
				down.locus <- c(down.locus,x[j])				
				down.matrix <- rbind(down.matrix,as.double(down.frame[i,2:ncol(down.frame)]))	
			}
		}else{
			down.locus <- c(down.locus,as.character(down.frame[i,1]))				
			down.matrix <- rbind(down.matrix,as.double(down.frame[i,2:ncol(down.frame)]))
		}		
	}
		
	down.frame <- data.frame(as.character(down.locus),down.matrix)
	
	up.locus <- unique(as.character(up.locus))
	names(up.locus) <- up.locus
	down.locus <- unique(as.character(down.locus))
	names(down.locus) <- down.locus
	up.locus <- up.locus[!(up.locus %in% down.locus)]
	down.locus <- down.locus[!(down.locus %in% up.locus)]

	up.locus <- up.locus[up.locus %in% as.character(locus.name[,1])]
	down.locus <- down.locus[down.locus %in% as.character(locus.name[,1])]
	
	up.frame.new <- NULL
		
	for(i in 1:length(up.locus)){
			
		datas <- up.frame[as.character(up.frame[,1]) == up.locus[i],2:ncol(up.frame)]
		
		if(logged & !is.vector(datas)){
			datas <- apply(datas,2,function(x) 2^as.double(x))
		}else if(logged & is.vector(datas)){
			datas <- 2^as.double(datas)
		}
			
		if(!is.vector(datas)){
			datas <- apply(datas,2,function(x) mean(as.double(x),na.rm=TRUE))
			
		}else{
			datas <- as.double(datas)
		}			
		
		if(logged){
			datas <- log(datas,2)
		}
		
		datas[is.nan(datas)] <- NA
		
		up.frame.new <- rbind(up.frame.new,datas)
			
	}
	
	up.frame.new <- data.frame(up.locus,up.frame.new)
		
	up.frame <- NULL
		
	for(i in 2:ncol(up.frame.new)){
			
		up.frame <- cbind(up.frame,as.numeric(up.frame.new[,i]))
		
	}
		
	up.frame <- data.frame(up.frame.new[,1],up.frame)
	colnames(up.frame) <- up.colnames
	
	down.frame.new <- NULL
		
	for(i in 1:length(down.locus)){
			
		datas <- down.frame[as.character(down.frame[,1]) == down.locus[i],2:ncol(down.frame)]
		
		if(logged & !is.vector(datas)){
			datas <- apply(datas,2,function(x) 2^as.double(x))
		}else if(logged & is.vector(datas)){
			datas <- 2^as.double(datas)
		}
			
		if(!is.vector(datas)){
			datas <- apply(datas,2,function(x) mean(as.double(x),na.rm=TRUE))
			
		}else{
			datas <- as.double(datas)
		}			
		
		if(logged){
			datas <- log(datas,2)
		}
		
		datas[is.nan(datas)] <- NA
			
		down.frame.new <- rbind(down.frame.new,datas)
			
	}
	
	down.frame.new <- data.frame(down.locus,down.frame.new)
	
	down.frame <- NULL
		
	for(i in 2:ncol(down.frame.new)){
			
		down.frame <- cbind(down.frame,as.numeric(down.frame.new[,i]))
		
	}
		
	down.frame <- data.frame(down.frame.new[,1],down.frame)
	colnames(down.frame) <- down.colnames
	
	return(list(up.frame=up.frame,down.frame=down.frame))
}


	rm()

}




#############################################################################################
#
# 6. Function .annotate.for.net() -> Routine for annotating genes and identifying relevant biological themes
#
#############################################################################################

.annotate.for.net <- function(file.annot,fdr=FALSE,go=TRUE,direct=FALSE,taxoname,annot.method="specificity",terms.name,restrict=FALSE,
				genes.frame=NULL,ref.list=NULL,nom="",enriched=TRUE){
	
	# annot.method can take values "specificity", "terminological", or "decorrelated"
	
	cat(paste("\n\t",taxoname," annotation of genes ",nom," started... ",format(Sys.time(), "%X"),sep=""))


	if(restrict & !is.null(ref.list)){

		file.annot <- file.annot[as.character(file.annot[,1]) %in% as.character(ref.list),]
	
	}
	
	file.annot <- file.annot[as.character(file.annot[,2]) %in% as.character(terms.name[,1]),]
	file.annot <- file.annot[!(file.annot[,2] %in% c("GO:0008150","GO:0000004","GO:0007582","GO:0005575","GO:0008372","GO:0003674","GO:0005554")),]
	
	ref.annot <- unique(as.character(file.annot[,1])) # list of all annotated genes within the considered taxonomical system among the reference list
	pop.total <- length(ref.annot)	# the number of annotated ref genes

	#exp.total <- 0	
		
	annot.matrix <- NULL	# annot matrix for the analyzed genes
	annot.rank <- NULL	# annot rank matrix (GO annotation level estimator)
	print.data <- NULL	# matrix with the annotation data for printing purposes


############################################################################################################

if(go & annot.method=="specificity"){	# GO annotation on separate ontological levels routine whitout keeping genes on superior levels

	go.tree <- GO.terms.hierarchy
	#go.tree <- go.tree[go.tree[,1] %in% as.vector(terms.name[,1]),]
	#go.tree <- go.tree[go.tree[,2] %in% as.vector(terms.name[,1]),]
			
	lowest.level <- as.character(go.tree[!(as.character(go.tree[,1]) %in% as.character(go.tree[,2])),1])	# GO terms which are on the lowest level of GO available
	
	if(nrow(file.annot) > 0 && !is.null(genes.frame)){
		
		genes <- as.character(genes.frame[,1])
		annot <- file.annot[as.character(file.annot[,1]) %in% genes,]	# the list of couples gene/annotation present in the analysed list of genes
		
		
		genes.annot <- unique(as.character(annot[,1])) # list of all annotated genes within the considered taxonomical system among the analyzed list
		exp.total <- length(genes.annot)	# the number of annotated analyzed genes
		
		terms.annot <- unique(as.character(annot[,2]))
		
		
		if(length(terms.annot)>0){	# if there are annotating terms available
		
		
		continue <- TRUE
		rang <- 1	# the GO specificity rank
		
		x <- terms.annot
		
		
		while(continue){
			
			cat(paste("\n\t\tGO specificity level: ",rang,sep=""))
			y <- NULL	# matrix of couples term/pvalue
			
			for(i in 1:length(x)){
				
				pval.i <- .pvalues(data.frame(exp.nr = length(unique(annot[annot[,2] == x[i],1])), exp.total = exp.total, pop.hits = length(unique(file.annot[file.annot[,2] == x[i],1])), pop.total = pop.total))
				
				y <- rbind(y,c(x[i],pval.i))
			}
			
			z <- NULL	# vector of significant annotations
			
			if(enriched){
				if(!is.na(fdr)){z <- as.character(y[.fdr.adjust(pvalues=as.numeric(y[,2]))<=(fdr/100),1])}	# vector with GO terms which are significantly enriched as corected by fdr
				if(is.na(fdr)){z <- as.character(y[as.numeric(y[,2])<=0.05,1])}
			}else{
				z <- as.character(y[,1])
			}

			if(length(z)>0){	# we insert significant annotations in annotation matrix and transfer genes of unsignificant annotations to the superior level of GO
			
			
				if(is.null(annot.matrix)){	# we build the annotation matrix
					annot.matrix <- list(NULL)
				}

				if(is.null(print.data)){	# we build the print matrix
					print.data <- list(NULL)
				}

				annot.matrix[[rang]] <- matrix(0,length(z),length(genes.annot))
				rownames(annot.matrix[[rang]]) <- z	# annotations on the rows
				colnames(annot.matrix[[rang]]) <- genes.annot	# genes on the columns

				print.data[[rang]] <- matrix(0,length(z),5)
				rownames(print.data[[rang]]) <- z	# annotations on the rows
				colnames(print.data[[rang]]) <- c("exp.nr","exp.total","pop.hits","pop.total","pval")	# data on the columns


				for(i in 1:length(z)){
			
					a <- rownames(annot.matrix[[rang]])[i]
					b <- as.character(annot[as.character(annot[,2]) == a,1])
					annot.matrix[[rang]][i,b] <- 1		
					print.data[[rang]][i,"exp.nr"] <- length(unique(as.vector(annot[as.character(annot[,2]) == z[i],1])))
					print.data[[rang]][i,"exp.total"] <- exp.total
					print.data[[rang]][i,"pop.hits"] <- length(unique(as.vector(file.annot[as.character(file.annot[,2]) == z[i],1])))
					print.data[[rang]][i,"pop.total"] <- pop.total
					print.data[[rang]][i,"pval"] <- .pvalues(data.frame(exp.nr = length(unique(as.vector(annot[as.character(annot[,2]) == z[i],1]))), exp.total = exp.total, 
										pop.hits = length(unique(as.vector(file.annot[as.character(file.annot[,2]) == z[i],1]))), pop.total = pop.total))
	
					rm(a,b)
				}

				w <- x[x %in% lowest.level] # we select current annotations (significant or not) which are on the lowest level
			
				if(length(w)>0){	# if there are such annotations we transfer their genes to the superior level of GO
			
					for(i in 1:length(w)){
					
						a <- rbind(NULL,annot[as.character(annot[,2]) == w[i],])
						b <- rbind(NULL,go.tree[as.character(go.tree[,1]) == w[i],])
# 						d <- rbind(NULL,file.annot[file.annot[,2] == w[i],])
					
						for(j in 1:nrow(b)){
							a[,2] <- b[j,2]
							annot <- rbind(annot,a)
# 							d[,2] <- b[j,2]
# 							file.annot <- rbind(file.annot,d)
						}
						rm(a,b)
					}
			
				}
				rm(w)
				
				w <- unique(as.character(file.annot[as.character(file.annot[,2]) %in% lowest.level,2])) # we select annotations which are on the current lowest level
			
				if(length(w)>0){	# if there are such annotations we transfer their genes to the superior level of GO
			
					for(i in 1:length(w)){
					
						b <- rbind(NULL,go.tree[as.character(go.tree[,1]) == w[i],])
						d <- rbind(NULL,file.annot[as.character(file.annot[,2]) == w[i],])
					
						for(j in 1:nrow(b)){
							d[,2] <- b[j,2]
							file.annot <- rbind(file.annot,d)
						}
						rm(b,d)
					}
			
				}
				rm(w)
				
			# we eliminate the lowest level terms from the GO tree and the annotation matrices
			
				annot <- annot[!(as.character(annot[,2]) %in% lowest.level),]
					
				file.annot <- file.annot[!(as.character(file.annot[,2]) %in% lowest.level),]

				go.tree <- go.tree[!(as.character(go.tree[,1]) %in% lowest.level),]
					

			# we recompose a new x and a new lowest.level
				x <- NULL
				
				terms.annot <- unique(as.character(annot[,2]))
				
				
				if(length(unique(go.tree[,2]))>2){
					lowest.level <- as.character(go.tree[!(as.character(go.tree[,1]) %in% as.character(go.tree[,2])),1])
					x <- terms.annot
					rang <- rang + 1
				}else{
					continue <- FALSE
				}
				
			
			}else{	# we simply transfer the genes annotated by the unsignificantly enriched terms to the superior level of GO
				annot.matrix[[rang]] <- NULL
				print.data[[rang]] <- NULL
				w <- x[x %in% lowest.level] # we select current annotations (significant or not) which are on the lowest level
			
				if(length(w)>0){	# if there are such annotations we transfer their genes to the superior level of GO
			
					for(i in 1:length(w)){
					
						a <- rbind(NULL,annot[as.character(annot[,2]) == w[i],])
						b <- rbind(NULL,go.tree[as.character(go.tree[,1]) == w[i],])
# 						d <- rbind(NULL,file.annot[file.annot[,2] == w[i],])
					
						for(j in 1:nrow(b)){
							a[,2] <- b[j,2]
							annot <- rbind(annot,a)
# 							d[,2] <- b[j,2]
# 							file.annot <- rbind(file.annot,d)
						}
						rm(a,b)
					}
			
				}
				rm(w)
				
				w <- unique(as.character(file.annot[as.character(file.annot[,2]) %in% lowest.level,2])) # we select annotations which are on the current lowest level
			
				if(length(w)>0){	# if there are such annotations we transfer their genes to the superior level of GO
			
					for(i in 1:length(w)){
					
						b <- rbind(NULL,go.tree[as.character(go.tree[,1]) == w[i],])
						d <- rbind(NULL,file.annot[as.character(file.annot[,2]) == w[i],])
					
						for(j in 1:nrow(b)){
							d[,2] <- b[j,2]
							file.annot <- rbind(file.annot,d)
						}
						rm(b,d)
					}
			
				}
				rm(w)
			
				annot <- annot[!(as.character(annot[,2]) %in% lowest.level),]
				file.annot <- file.annot[!(as.character(file.annot[,2]) %in% lowest.level),]
				go.tree <- go.tree[!(as.character(go.tree[,1]) %in% lowest.level),]
	
				# we recompose a new x
				x <- NULL
							
				terms.annot <- unique(as.character(annot[,2]))
				
							
				if(length(unique(go.tree[,2]))>2){
					lowest.level <- as.character(go.tree[!(as.character(go.tree[,1]) %in% as.character(go.tree[,2])),1])
					x <- terms.annot
					rang <- rang + 1
				}else{
					continue <- FALSE
				}
			
			
			}	# end if(length(z)>0)
			
		
		
		}	# end while()
		
		}	# end if(length(terms.annot)>0)
	
	
	
	}


#####################################################################################################################

}else if(go & annot.method=="terminological"){ # respect GO terminological levels

	go.tree <- GO.terms.hierarchy
	#go.tree <- go.tree[go.tree[,1] %in% as.vector(terms.name[,1]),]
	#go.tree <- go.tree[go.tree[,2] %in% as.vector(terms.name[,1]),]
	
	terms.level <- as.vector(matrix(0,nrow(terms.name),1))
	names(terms.level) <- terms.name[,1]
	
	continue <- TRUE
	n <- 1
	
	while(continue == TRUE){
		
		terms.level[go.tree[!(go.tree[,2] %in% go.tree[,1]),2]] <- n
		go.tree <- go.tree[go.tree[,2] %in% go.tree[,1],]
		n <- n + 1
		if(nrow(go.tree) == 0){continue <- FALSE}
		
	}
	
	go.tree <- GO.terms.hierarchy
			
	
	if(nrow(file.annot) > 0 && !is.null(genes.frame)){
		
		genes <- as.character(genes.frame[,1])
		annot <- file.annot[file.annot[,1] %in% genes,]	# the list of couples gene/annotation present in the analysed list of genes
		
		
		genes.annot <- unique(as.character(annot[,1])) # list of all annotated genes within the considered taxonomical system among the analyzed list
		exp.total <- length(genes.annot)	# the number of annotated analyzed genes
 		
		terms.annot <- unique(as.character(annot[,2]))
		
		rang <- max(terms.level)	# the GO terminological level
		
		
			
		for(i in 1:(max(terms.level) - 1)){
			
			cat(paste("\n\t\tGO terminological level: ",rang,sep=""))
			x <- terms.annot[terms.annot %in% names(terms.level)[terms.level == rang]]
 			
 			#print(length(x))
			
			if(length(x) > 0){	# if there are annotating terms available


			
			
				y <- NULL	# matrix of couples term/pvalue
			
				for(j in 1:length(x)){
	 				#print(data.frame(exp.nr = length(unique(annot[annot[,2] == x[j],1])), exp.total = exp.total, pop.hits = length(unique(file.annot[file.annot[,2] == x[j],1])), pop.total = pop.total))

					pval.j <- .pvalues(data.frame(exp.nr = length(unique(annot[annot[,2] == x[j],1])), exp.total = exp.total, pop.hits = length(unique(file.annot[file.annot[,2] == x[j],1])), pop.total = pop.total))

					y <- rbind(y,c(x[j],pval.j))
				}
			
				z <- NULL	# vector of significant annotations
				
				if(enriched){
					if(!is.na(fdr)){z <- as.character(y[.fdr.adjust(pvalues=as.numeric(y[,2]))<=(fdr/100),1])}	# vector with GO terms which are significantly enriched as corected by fdr
					if(is.na(fdr)){z <- as.character(y[as.numeric(y[,2])<=0.05,1])}
				}else{
					z <- y[,1]
				}
				
				#print(length(z))

				if(length(z)>0){	# we insert significant annotations in annotation matrix and transfer genes of unsignificant annotations to the superior level of GO


					if(is.null(annot.matrix)){	# we build the annotation matrix
						annot.matrix <- list(NULL)
					}

					if(is.null(print.data)){	# we build the print matrix
						print.data <- list(NULL)
					}
					
					annot.matrix[[i]] <- matrix(0,length(z),length(genes.annot))
					rownames(annot.matrix[[i]]) <- z	# annotations on the rows
					colnames(annot.matrix[[i]]) <- genes.annot	# genes on the columns

					print.data[[i]] <- matrix(0,length(z),5)
					rownames(print.data[[i]]) <- z	# annotations on the rows
					colnames(print.data[[i]]) <- c("exp.nr","exp.total","pop.hits","pop.total","pval")	# data on the columns



					for(j in 1:length(z)){

						a <- rownames(annot.matrix[[i]])[j]
						b <- annot[annot[,2] == a,1]
						annot.matrix[[i]][j,b] <- 1		
						rm(a,b)

						print.data[[i]][j,"exp.nr"] <- length(unique(annot[annot[,2] == z[j],1]))
						print.data[[i]][j,"exp.total"] <- exp.total
						print.data[[i]][j,"pop.hits"] <- length(unique(file.annot[file.annot[,2] == z[j],1]))
						print.data[[i]][j,"pop.total"] <- pop.total
						print.data[[i]][j,"pval"] <- .pvalues(data.frame(exp.nr = length(unique(annot[annot[,2] == z[j],1])), exp.total = exp.total, 
										pop.hits = length(unique(file.annot[file.annot[,2] == z[j],1])), pop.total = pop.total))
						
					}
	 				#print("ok1")


					# we transfer the genes to the superior level of GO

					for(k in 1:length(x)){

						a <- rbind(NULL,annot[annot[,2] == x[k],])
						b <- rbind(NULL,go.tree[go.tree[,1] == x[k],])

						for(j in 1:nrow(b)){
							a[,2] <- b[j,2]
							annot <- rbind(annot,a)
						}
						rm(a,b)
					}

					
					
	 				#print("ok2")

					w <- unique(as.character(file.annot[file.annot[,2] %in% names(terms.level)[terms.level == rang],2])) # we select annotations which are on the current lowest level
	 				#print(length(w))

					if(length(w)>0){	# if there are such annotations we transfer their genes to the superior level of GO

						for(k in 1:length(w)){

							b <- rbind(NULL,go.tree[go.tree[,1] == w[k],])
							d <- rbind(NULL,file.annot[file.annot[,2] == w[k],])

							for(j in 1:nrow(b)){
								d[,2] <- b[j,2]
								file.annot <- rbind(file.annot,d)
							}
							rm(b,d)
						}

					}
					rm(w)
	 				#print("ok3")

	# we eliminate the lowest level terms from the GO tree and the annotation matrices

	# 				
					annot <- annot[!(annot[,2] %in% names(terms.level)[terms.level == rang]),]

					file.annot <- file.annot[!(file.annot[,2] %in% names(terms.level)[terms.level == rang]),]

					go.tree <- go.tree[!(go.tree[,1] %in% names(terms.level)[terms.level == rang]),]


	# we recompose a new x and a new lowest.level
					x <- NULL

					terms.annot <- unique(as.character(annot[,2]))

					rang <- rang - 1

	#print(x)

				}else{	# we simply transfer the genes annotated by the unsignificantly enriched terms to the superior level of GO
					annot.matrix[[i]] <- NULL
					print.data[[i]] <- NULL

					# we select current annotations (significant or not) which are on the lowest level

					# we transfer their genes to the superior level of GO
					#print(x)
					for(k in 1:length(x)){

						a <- rbind(NULL,annot[annot[,2] == x[k],])
						#print(a)
						b <- rbind(NULL,go.tree[go.tree[,1] == x[k],])
						#print(b)

						for(j in 1:nrow(b)){
							a[,2] <- b[j,2]
							annot <- rbind(annot,a)
						}
						rm(a,b)
					}

	 				#print("ok4")

					w <- unique(as.character(file.annot[file.annot[,2] %in% names(terms.level)[terms.level == rang],2])) # we select annotations which are on the current lowest level
	 				#print(length(w))

					if(length(w)>0){	# if there are such annotations we transfer their genes to the superior level of GO

						for(k in 1:length(w)){

							b <- rbind(NULL,go.tree[go.tree[,1] == w[k],])
							d <- rbind(NULL,file.annot[file.annot[,2] == w[k],])

							for(j in 1:nrow(b)){
								d[,2] <- b[j,2]
								file.annot <- rbind(file.annot,d)
							}
							rm(b,d)
						}

					}
					rm(w)
	 				#print("ok5")

					annot <- annot[!(annot[,2] %in% names(terms.level)[terms.level == rang]),]
					file.annot <- file.annot[!(file.annot[,2] %in% names(terms.level)[terms.level == rang]),]
					go.tree <- go.tree[!(go.tree[,1] %in% names(terms.level)[terms.level == rang]),]

	# we recompose a new x
					x <- NULL

					terms.annot <- unique(as.character(annot[,2]))


					rang <- rang - 1

#print(x)
			
				}	# end if(length(z)>0)
			
		}else{ # if(length(x) >0)
			
			annot.matrix[[i]] <- NULL
			print.data[[i]] <- NULL
							
			# we select current annotations (significant or not) which are on the lowest level
						
			# we transfer their genes to the superior level of GO
						
							
			w <- unique(as.character(file.annot[file.annot[,2] %in% names(terms.level)[terms.level == rang],2])) # we select annotations which are on the current lowest level
 			#print(length(w))
						
			if(length(w)>0){	# if there are such annotations we transfer their genes to the superior level of GO
						
				for(k in 1:length(w)){
								
					b <- rbind(NULL,go.tree[go.tree[,1] == w[k],])
					d <- rbind(NULL,file.annot[file.annot[,2] == w[k],])
								
					for(j in 1:nrow(b)){
						d[,2] <- b[j,2]
						file.annot <- rbind(file.annot,d)
					}
					rm(b,d)
				}
						
			}
			rm(w)
			#print("ok6")
						
			annot <- annot[!(annot[,2] %in% names(terms.level)[terms.level == rang]),]
			file.annot <- file.annot[!(file.annot[,2] %in% names(terms.level)[terms.level == rang]),]
			go.tree <- go.tree[!(go.tree[,1] %in% names(terms.level)[terms.level == rang]),]
				
# we recompose a new x
			x <- NULL
										
			terms.annot <- unique(as.character(annot[,2]))
							
			rang <- rang - 1
		
		
		
		}# end if(length(x) >0)
		
		}	# end for()
		
	
	} # end if(nrow...
#####################################################################################################################

}else if(go & annot.method=="decorrelated"){	# GO decorrelated annotation routine

	go.tree <- GO.terms.hierarchy
	#go.tree <- go.tree[go.tree[,1] %in% as.vector(terms.name[,1]),]
	#go.tree <- go.tree[go.tree[,2] %in% as.vector(terms.name[,1]),]
		
			
	lowest.level <- go.tree[!(go.tree[,1] %in% go.tree[,2]),1]	# GO terms which are on the lowest (most specific) annotation level of GO available
	lowest.level <- lowest.level[!(lowest.level %in% c("GO:0008150","GO:0005575","GO:0003674"))]
	
	if(nrow(file.annot) > 0 && !is.null(genes.frame)){
		
		genes <- as.character(genes.frame[,1])
		annot <- file.annot[file.annot[,1] %in% genes,]	# the list of couples gene/annotation present in the analysed list of genes
		
		
		genes.annot <- unique(as.character(annot[,1])) # list of all annotated genes within the considered taxonomical system among the analyzed list
		exp.total <- length(genes.annot)	# the number of annotated analyzed genes
		
		terms.annot <- unique(as.character(annot[,2]))
		
		
		if(length(terms.annot)>0){	# if there are annotating terms available
		
		
		continue <- TRUE
		rang <- 1	# the GO annotation specificity rank
		
		x <- NULL
		direct.stock <- NULL	# a vector to store significant direct annotations which are not on the lowest level of GO
		
		if(direct == FALSE){x <- terms.annot[terms.annot %in% lowest.level]}	# vector with annotating GO terms which are not subsuming others
		if(direct == TRUE){x <- terms.annot}	# consider all direct annotations first
		
		
		
		while(continue){
		
		if(length(x) > 0){ # test that x is not an empty vector
		
			y <- NULL	# matrix of couples term/pvalue
			
			for(i in 1:length(x)){
				
				pval.i <- .pvalues(data.frame(exp.nr = length(unique(annot[annot[,2] == x[i],1])), exp.total = exp.total, pop.hits = length(unique(file.annot[file.annot[,2] == x[i],1])), pop.total = pop.total))
				
				y <- rbind(y,c(x[i],as.character(pval.i)))
			}
			
			z <- NULL	# vector of significant annotations
			
			if(!is.na(fdr)){z <- as.character(y[.fdr.adjust(pvalues=as.numeric(y[,2]))<=(fdr/100),1])}	# vector with GO terms which are significantly enriched as corected by fdr
			if(is.na(fdr)){z <- as.character(y[as.numeric(y[,2])<=0.05,1])}
		
		
		#print(length(z[!is.na(z)]))
		#print(rang)
		if(length(z[!is.na(z)])>0){	# we insert significant annotations in annotation matrix and transfer genes of unsignificant annotations to the superior level of GO
			
			
			if(is.null(annot.matrix)){	# we build the annotation matrix
				
				annot.matrix <- matrix(0,length(z),length(genes.annot))
				rownames(annot.matrix) <- z	# annotations on the rows
				colnames(annot.matrix) <- genes.annot	# genes on the columns

				print.data <- matrix(0,length(z),5)	# we build the print matrix
				rownames(print.data) <- z	# annotations on the rows
				colnames(print.data) <- c("exp.nr","exp.total","pop.hits","pop.total","pval")	# data on the columns


				for(i in 1:length(z)){
			
					a <- rownames(annot.matrix)[i]
					b <- annot[annot[,2] == a,1]
					annot.matrix[i,b] <- 1	
					
					print.data[i,"exp.nr"] <- length(unique(annot[annot[,2] == z[i],1]))
					print.data[i,"exp.total"] <- exp.total
					print.data[i,"pop.hits"] <- length(unique(file.annot[file.annot[,2] == z[i],1]))
					print.data[i,"pop.total"] <- pop.total
					print.data[i,"pval"] <- .pvalues(data.frame(exp.nr = print.data[i,"exp.nr"], exp.total = exp.total, 
										pop.hits = print.data[i,"pop.hits"], pop.total = pop.total))
		
				}
				
				annot.rank <- as.vector(matrix(rang,1,length(z)))
				names(annot.rank) <- z
				
				if(rang == 1 && length(z[!(z %in% lowest.level)])>0){
					direct.stock <- z[!(z %in% lowest.level)]
				}	# we store significant direct annotations which are not on the lowest level

			
			}else{	# we complete the annotation matrix
				
				if(length(direct.stock)>0 & length(z[z %in% direct.stock]) > 0){
					
					#print(direct.stock)
					#print("OK")
					
					m <- z[!(z %in% direct.stock)]
					n <- z[z %in% direct.stock]
					
					# we treat first annotations which were not stored in direct.stock
					if(length(m)>0){
					
						annot.matrix1 <- matrix(0,length(m),length(genes.annot))
						rownames(annot.matrix1) <- m	# annotations on the rows
						colnames(annot.matrix1) <- genes.annot	# genes on the columns
						
						print.data1 <- matrix(0,length(m),5)	# we build the print matrix
						rownames(print.data1) <- m	# annotations on the rows
						colnames(print.data1) <- c("exp.nr","exp.total","pop.hits","pop.total","pval")	# data on the columns
																
						for(i in 1:length(m)){
																	
							a <- rownames(annot.matrix1)[i]
							b <- annot[annot[,2] == a,1]
							annot.matrix1[i,b] <- 1	
							
							print.data1[i,"exp.nr"] <- length(unique(annot[annot[,2] == m[i],1]))
							print.data1[i,"exp.total"] <- exp.total
							print.data1[i,"pop.hits"] <- length(unique(file.annot[file.annot[,2] == m[i],1]))
							print.data1[i,"pop.total"] <- pop.total
							print.data1[i,"pval"] <- .pvalues(data.frame(exp.nr = length(unique(annot[annot[,2] == m[i],1])), exp.total = exp.total, 
											pop.hits = length(unique(file.annot[file.annot[,2] == m[i],1])), pop.total = pop.total))
																
						}
														
						annot.matrix <- rbind(annot.matrix,annot.matrix1)	# we add the new significantly enriched annotations to the annotation matrix
						print.data <- rbind(print.data,print.data1)
						
						annot.rank1 <- as.vector(matrix(rang,1,length(m)))
						names(annot.rank1) <- m
						annot.rank <- c(annot.rank,annot.rank1)
					}
					
					# we treat here annotations which are in direct.stock and thus already in annot.matrix by adding new genes
						
					if(length(n)>0){	
						
						for(i in 1:length(n)){
																						
							b <- annot[annot[,2] == n[i],1]
							annot.matrix[n[i],b] <- 1	
							# we update available data on these annotations
							print.data[n[i],"exp.nr"] <- length(unique(annot[annot[,2] == n[i],1]))
							print.data[n[i],"pop.hits"] <- length(unique(file.annot[file.annot[,2] == n[i],1]))
							print.data[n[i],"pval"] <- .pvalues(data.frame(exp.nr = length(unique(annot[annot[,2] == n[i],1])), exp.total = exp.total, 
											pop.hits = length(unique(file.annot[file.annot[,2] == n[i],1])), pop.total = pop.total))
																					
						}
					}
				
				}else{
					#print(z)
					annot.matrix1 <- matrix(0,length(z),length(genes.annot))
					rownames(annot.matrix1) <- z	# annotations on the rows
					colnames(annot.matrix1) <- genes.annot	# genes on the columns

					print.data1 <- matrix(0,length(z),5)	# we build the print matrix
					rownames(print.data1) <- z	# annotations on the rows
					colnames(print.data1) <- c("exp.nr","exp.total","pop.hits","pop.total","pval")	# data on the columns

					for(i in 1:length(z)){
												
						a <- rownames(annot.matrix1)[i]
						b <- annot[annot[,2] == a,1]
						annot.matrix1[i,b] <- 1	

						print.data1[i,"exp.nr"] <- length(unique(annot[annot[,2] == z[i],1]))
						print.data1[i,"exp.total"] <- exp.total
						print.data1[i,"pop.hits"] <- length(unique(file.annot[file.annot[,2] == z[i],1]))
						print.data1[i,"pop.total"] <- pop.total
						print.data1[i,"pval"] <- .pvalues(data.frame(exp.nr = length(unique(annot[annot[,2] == z[i],1])), exp.total = exp.total, 
										pop.hits = length(unique(file.annot[file.annot[,2] == z[i],1])), pop.total = pop.total))

											
					}
									
					annot.matrix <- rbind(annot.matrix,annot.matrix1)	# we add the new significantly enriched annotations to the annotation matrix
					print.data <- rbind(print.data,print.data1)

					annot.rank1 <- as.vector(matrix(rang,1,length(z)))
					names(annot.rank1) <- z
					annot.rank <- c(annot.rank,annot.rank1)
				
				}
				
			
			}
			
			# print("OK")
			
			w <- x[!(x %in% z)]	# select unsignificant annotations
			
			if(direct == FALSE){w <- w[w %in% lowest.level]}	# select unsignificant annotations which are on the lowest level
			
			if(length(w)>0){	# if there are unsignificant annotations left transfer their genes to the superior level of GO
			
				for(i in 1:length(w)){
				
					a <- rbind(NULL,annot[annot[,2] == w[i],])
					b <- rbind(NULL,go.tree[go.tree[,1] == w[i],])
					d <- rbind(NULL,file.annot[file.annot[,2] == w[i],])
				
					for(j in 1:nrow(b)){
						a[,2] <- b[j,2]
						annot <- rbind(annot,a)
						d[,2] <- b[j,2]
						file.annot <- rbind(file.annot,d)
					}
			
				}
			
			}
		
			# we eliminate all the analysed terms from the GO tree
			
			if(direct == FALSE){
			
				a <- x[x %in% lowest.level]	# we eliminate significant and unsignificant annotations situated on the lowest level of GO; 
				b <- z[z %in% lowest.level]	# significant annotations situated on the lowest level
				d <- annot[annot[,2] %in% b,1]	# vector of genes annotated by significant annotations situated on the lowest level
				
				annot <- annot[!(annot[,1] %in% d),]
				annot <- annot[!(annot[,2] %in% a),]
				
				file.annot <- file.annot[!(file.annot[,1] %in% d),]
				file.annot <- file.annot[!(file.annot[,2] %in% a),]
				
				go.tree <- go.tree[!(go.tree[,1] %in% lowest.level),]
				
				
				
			}else{
				b <- annot[annot[,2] %in% z,1]	# vector of genes annotated by significant annotations
				annot <- annot[!(annot[,1] %in% b),]	# we eliminate all genes already annotated by significant annotations
				annot <- annot[!(annot[,2] %in% x),]	# we eliminate all analysed annotations significant and unsignificant
				file.annot <- file.annot[!(file.annot[,1] %in% b),]
				file.annot <- file.annot[!(file.annot[,2] %in% x),]			
				go.tree <- go.tree[!(go.tree[,1] %in% lowest.level),]
				direct.stock <- direct.stock[!(direct.stock %in% lowest.level)]

			}
			
			# we recompose a new x
			x <- NULL
			
			terms.annot <- unique(as.character(annot[,2]))
			
			
			if(length(unique(go.tree[,2][!is.na(go.tree[,2])]))>2){
				#print(unique(go.tree[,2]))
				lowest.level <- go.tree[!(go.tree[,1] %in% go.tree[,2]),1]
				x <- terms.annot[terms.annot %in% lowest.level]
				rang <- rang + 1
			}else{
				continue <- FALSE
			}
			
			#print(x)
			
		}else{	# we simply transfer the genes annotated by the unsignificantly enriched terms to the superior level of GO
			
			w <- x
			
			if(direct == FALSE){w <- w[w %in% lowest.level]}
			if(direct == TRUE){direct.stock <- direct.stock[!(direct.stock %in% lowest.level)]}
			
			for(i in 1:length(w)){
				
				a <- rbind(NULL,annot[annot[,2] == w[i],])
				b <- rbind(NULL,go.tree[go.tree[,1] == w[i],])
				d <- rbind(NULL,file.annot[file.annot[,2] == w[i],])
				
				for(j in 1:nrow(b)){
					a[,2] <- b[j,2]
					annot <- rbind(annot,a)
					d[,2] <- b[j,2]
					file.annot <- rbind(file.annot,d)
				}
			
			}
		
			annot <- annot[!(annot[,2] %in% w),]
			file.annot <- file.annot[!(file.annot[,2] %in% w),]
			go.tree <- go.tree[!(go.tree[,1] %in% lowest.level),]

			# we recompose a new x
			x <- NULL
						
			terms.annot <- unique(as.character(annot[,2]))
			
						
			if(length(unique(go.tree[,2][!is.na(go.tree[,2])]))>2){
				#print(unique(go.tree[,2]))
				lowest.level <- go.tree[!(go.tree[,1] %in% go.tree[,2]),1]
				x <- terms.annot[terms.annot %in% lowest.level]
				rang <- rang + 1
			}else{
				continue <- FALSE
			}
		
			#print(x)
		
		}	# end if(length(z)>0)
		
		
		}else{
			continue <- FALSE
		}# end if(length(x)>0)
		}	# end while()
		
		}	# end if(length(terms.annot)>0)
	
	
	
	}
	
####################################################################################################################


}else if(!go){	# if KEGG annotations are used (or other non ontological annotations)
	
	genes <- as.character(genes.frame[,1])
	annot <- file.annot[as.character(file.annot[,1]) %in% genes,]	# the list of couples gene/annotation present in the analysed list of genes
			
			
	genes.annot <- unique(as.character(annot[,1]))	# list of all annotated genes within the considered taxonomical system among the analyzed list
	exp.total <- length(genes.annot)	# the number of annotated analyzed genes
			
	terms.annot <- unique(as.character(annot[,2]))
					
	if(length(terms.annot)>0){
		
		y <- NULL	# matrix of couples term/pvalue
	
		for(i in 1:length(terms.annot)){

			pval.i <- .pvalues(data.frame(exp.nr = length(unique(as.character(annot[as.character(annot[,2]) == terms.annot[i],1]))), 
						exp.total = exp.total, 
						pop.hits = length(unique(as.character(file.annot[as.character(file.annot[,2]) == terms.annot[i],1]))), 
						pop.total = pop.total))

			y <- rbind(y,c(terms.annot[i],pval.i))
		}
				
		z <- NULL

		if(enriched){
			if(!is.na(fdr)){z <- as.character(y[.fdr.adjust(pvalues=as.numeric(y[,2]))<=(fdr/100),1])}	# vector with KEGG terms which are significantly enriched as corected by fdr
			if(is.na(fdr)){z <- as.character(y[as.numeric(y[,2])<=0.05,1])}
		}else{
			z <- as.character(y[,1])
		}
	
	
	if(length(z)>0){	# we insert significant annotations in annotation matrix and transfer genes of unsignificant annotations to the superior level of GO
					
					
						
		annot.matrix <- matrix(0,length(z),length(genes.annot))
		rownames(annot.matrix) <- z	# annotations on the rows
		colnames(annot.matrix) <- genes.annot	# genes on the columns

		print.data <- matrix(0,length(z),5)	# we build the print matrix
		rownames(print.data) <- z	# annotations on the rows
		colnames(print.data) <- c("exp.nr","exp.total","pop.hits","pop.total","pval")	# data on the columns
	
		for(i in 1:length(z)){
					
			a <- rownames(annot.matrix)[i]
			b <- as.character(annot[as.character(annot[,2]) == a,1])
			annot.matrix[i,b] <- 1	
			
			print.data[i,"exp.nr"] <- length(unique(as.character(annot[as.character(annot[,2]) == z[i],1])))
			print.data[i,"exp.total"] <- exp.total
			print.data[i,"pop.hits"] <- length(unique(as.character(file.annot[as.character(file.annot[,2]) == z[i],1])))
			print.data[i,"pop.total"] <- pop.total
			print.data[i,"pval"] <- .pvalues(data.frame(exp.nr = length(unique(as.character(annot[as.character(annot[,2]) == z[i],1]))), 
								exp.total = exp.total, 
								pop.hits = length(unique(as.character(file.annot[as.character(file.annot[,2]) == z[i],1]))), 
								pop.total = pop.total))

				
		}
						
				
	}
	}	# end if(length(terms.annot)>0)

}	# end if(go == TRUE)

	return(list(annot.matrix=annot.matrix,print.data=print.data,annot.rank=annot.rank))
	rm()

}


#############################################################################################
#
# 7. Function .dynamic.proximity() -> Computes biological themes proximity within a convergent dynamic system build from annotations and a transcriptional co-expression matrix
#
#############################################################################################


.dynamic.proximity <- function(annot.matrix=NULL, adj.matrix=NULL, RV=0.90, sigma=NA, method="sum", annotated.only=TRUE){

	
	if(!is.null(annot.matrix) & !is.null(adj.matrix)){	# genes on rows, categories (future markers) on columns
		
		# we add non annotated transcripts to the annotation matrix
		
		n <- length(rownames(adj.matrix)[!(rownames(adj.matrix) %in% rownames(annot.matrix))])
		n <- matrix(0, n,ncol(annot.matrix))
		rownames(n) <- rownames(adj.matrix)[!(rownames(adj.matrix) %in% rownames(annot.matrix))]
		colnames(n) <- colnames(annot.matrix)

		annot.matrix <- rbind(annot.matrix,n)
		annot.matrix <- annot.matrix[rownames(adj.matrix),]

		# initiate the dynamic system
		steps <- 1
		continue <- TRUE
		
		marker.matrix <- as.list(NULL)	
		marker.matrix[[1]] <- annot.matrix
		lngth <- apply(marker.matrix[[1]],2, function(x) sqrt(sum(x*x)))
		marker.matrix[[1]] <- t(t(marker.matrix[[1]])/lngth)
		rm(lngth)
		gc(); gc();
				
		while(continue){
			
			steps <- steps + 1
			
			marker.matrix[[steps]] <- marker.matrix[[steps-1]] # genes on lines colors on columns
			marker.matrix[[steps]][,] <- 0
			
			for(k in 1:ncol(adj.matrix)){
				
				if(length(grep(pattern="[[:digit:]]00", x=k)) > 0){print(paste("Iteration ",steps," : ",k,sep=""))}
				if(method == "sum"){
					marker.matrix[[steps]][k,] <- apply(adj.matrix[,k] * marker.matrix[[steps-1]],2,function(x) sum(x,na.rm=T))
				}else if(method == "max"){
					marker.matrix[[steps]][k,] <- apply(adj.matrix[,k] * marker.matrix[[steps-1]],2,function(x) max(x,na.rm=T))
				}else if(method == "prod"){
					marker.matrix[[steps]][k,] <- apply(adj.matrix[,k] * marker.matrix[[steps-1]],2,function(x) prod(x[x!=0],na.rm=T))
				}else if(method == "mean"){
					marker.matrix[[steps]][k,] <- apply(adj.matrix[,k] * marker.matrix[[steps-1]],2,function(x) mean(x,na.rm=T))
				}
			}
			
			lngth <- apply(marker.matrix[[steps]],2, function(x) sqrt(sum(x*x)))
			marker.matrix[[steps]] <- t(t(marker.matrix[[steps]])/lngth)
			rm(lngth)
						
			gc(); gc();
			
			# don't stop on stop(.) in case complex numbers comparison by coinertia() will produce an error 
			options(error = expression(NULL))
			options(warn = -1)
			count.error <- 0
			
			# perform a coinertia analysis in between the actual and the last configurations of the dynamical system to test for convergence
			xt <- NULL
			xtt <- NULL
			#try(xt <- dudi.pca(t(marker.matrix[[steps - 1]]),scale=TRUE,scan=FALSE,nf=ncol(marker.matrix[[steps - 1]])))
			#try(xtt <- dudi.pca(t(marker.matrix[[steps]]),scale=TRUE,scan=FALSE,nf=ncol(marker.matrix[[steps]])))
			try(xt <- dudi.pco(dist(t(marker.matrix[[steps - 1]])),scannf=FALSE,full=TRUE))
			try(xtt <- dudi.pco(dist(t(marker.matrix[[steps]])),scannf=FALSE,full=TRUE))
			
			if(exists("xt") & exists("xtt") & !is.null("xt") & !is.null("xtt")){
			
				try(coin <- coinertia(xt,xtt,scan=FALSE,nf=ncol(marker.matrix[[steps]])))
				
				try(print(try(paste("Coinertia RV: ",coin$RV,sep=""))))
				#try(plot(coin))

				# perform also a Mantel test on distance matrixes 
				try(mantel <- mantel.rtest(dist(t(marker.matrix[[steps - 1]])),dist(t(marker.matrix[[steps]])), nrepet=100))
				#try(plot(mantel))
								
				if(exists("mantel")){
					try(print(try(paste("Mantel's test RV: ",mantel$obs,sep=""))))
					try(print(try(paste("Mantel's test p-value: ",mantel$pvalue,sep=""))))
					
					if(!is.na(mantel$obs) & !is.na(coin$RV)){
						if(coin$RV >= RV & mantel$obs >= RV){
							print("Convergence was attaint!")
							continue <- FALSE
						}
					}else{
						if(coin$RV >= RV){
							print("Convergence was attaint!")
							continue <- FALSE
						}
					}
				}else{
					if(coin$RV >= RV){
						print("Convergence was attaint!")
						continue <- FALSE
					}
				}
				rm(coin)
			
			}else{
				count.error <- count.error + 1
				if(count.error >= 3){
					continue <- FALSE
					print("Coinertia computation error! Stopping iterations and computing proximity.")
				}
			}
			
			# revert to the default way of error handling 
			options(error = NULL)	
			options(warn = 0)
			
			#an old way around...
			#p.val <- wilcox.test(as.numeric(marker.matrix[[steps]]),as.numeric(marker.matrix[[steps-1]]))$p.value
			#print(p.val)
			#if(p.val >= pval){continue <- FALSE}
			
		}
		
		dist.annotation <- NULL
		annot.proximity <- NULL
		
		if(!is.null(marker.matrix[[steps]])){
			mark.matrix <- marker.matrix[[steps]]
			mark.matrix[is.nan(mark.matrix)] <- 0
			mark.matrix[is.na(mark.matrix)] <- 0

			if(annotated.only){
				marked.transcripts <- apply(mark.matrix,1,function(x) sum(x,na.rm=TRUE))
				marked.transcripts <- names(marked.transcripts[marked.transcripts > 0])
				dist.annotation <- as.matrix(dist(t(mark.matrix[marked.transcripts,]), method="euclid"))
			}else{
				dist.annotation <- as.matrix(dist(t(mark.matrix), method="euclid"))
			}
			
								
			# annot.proximity <- exp(-(dist.annotation))
			if(is.na(sigma)){
				# the old way...
				# sigma <- median(dist.annotation[lower.tri(dist.annotation,diag=FALSE)])
				# annot.proximity <- exp(-(dist.annotation^2)/(2*(sigma^2)))
				
				# the new way with local sigma estimation:
				sigma.vect <- apply(dist.annotation,1,function(x) median(x[x >= median(x)]))
				
				annot.proximity <- (dist.annotation^2)*(1/sigma.vect)
				annot.proximity <- t(annot.proximity)*(1/sigma.vect)
				annot.proximity <- exp(-annot.proximity)
				diag(annot.proximity) <- 0
				
			}else{
				annot.proximity <- exp(-(dist.annotation^2)/(sigma^2))
				diag(annot.proximity) <- 0
			}

		}
		
		return(list(annot.proximity=annot.proximity, dist.annotation=dist.annotation, marker.matrix=marker.matrix, steps=steps, 
				annot.matrix=annot.matrix, adj.matrix=adj.matrix, method=method, RV=RV, sigma=sigma))
		
	}
	
	
}



#############################################################################################
#
# 8. Function .spectral.cluster() -> Spectral clustering of biological themes based on their dynamic diffusion profiles in the co-expression network
#
#############################################################################################


.spectral.cluster <- function(annot.proximity=NULL, marker.matrix=NULL, adj.matrix=NULL, terms.name=NULL, annot.matrix=NULL, iter.max=20000, k=NULL, 
				nstart=10){
	
	
	if(!is.null(annot.proximity) & !is.null(annot.matrix) & !is.null(terms.name) & !is.null(adj.matrix)){
		
				
		if(ncol(annot.proximity) > 2){
			
			diag.matrix <- matrix(0,nrow(annot.proximity),nrow(annot.proximity))
			dimnames(diag.matrix) <- dimnames(annot.proximity)
			diag(diag.matrix) <- apply(annot.proximity,1,sum)

			e <- eigen(solve(diag.matrix))
			V <- e$vectors
			rownames(V) <- rownames(diag.matrix)
			sqrt.diag <- V %*% diag(e$values) %*% t(V)
			
			l.matrix <- sqrt.diag %*% annot.proximity %*% sqrt.diag
			rm(e, V, sqrt.diag)
			
			# l.matrix <- (solve(diag.matrix)^0.5) %*% annot.proximity %*% (solve(diag.matrix)^0.5) # old icml/ecml way...
			# l.matrix <- (solve(diag.matrix)^0.5) %*% (diag.matrix - annot.proximity) %*% (solve(diag.matrix)^0.5)

			x.matrix <- eigen(l.matrix)$vectors	# eigenvectors matrix (eigen on columns) containing ALL eigenvectors
			rownames(x.matrix) <- rownames(l.matrix)
			#plot(eigen(l.matrix,only.values = TRUE)$values,pch=21,col="blue",bg="blue",xlab="SCREE plot", ylab="Normalized eigenvalue")
			print(eigen(l.matrix,only.values = TRUE)$values)
			print("First eigenvector:")
			print(x.matrix[,1])
			eigen.centrality <- x.matrix[,1]


			# start clustering procedure
			cluster.tree <- vector("list", ncol(x.matrix)-2)
			names(cluster.tree) <- as.character(2:(ncol(x.matrix)-1))

			#bic.vector <- as.vector(matrix(NA, ncol(x.matrix)-2,1))
			#names(bic.vector) <- as.character(2:(ncol(x.matrix)-1))

			for(i in 2:(ncol(x.matrix)-1)){

				k.eigen <- x.matrix[,1:i]	# first k eigenvectors from x.matrix
				lngth <- apply(k.eigen,1, function(x) sqrt(sum(x*x)))
				y.matrix <- k.eigen/lngth	# eigenvectors matrix normalized on the rows
				y.matrix[is.nan(y.matrix)] <- 0
				y.matrix[is.na(y.matrix)] <- 0
				rm(lngth)

				tmp <- kmeans(y.matrix, centers=i, iter.max = iter.max, nstart = nstart, algorithm="Hart")
				cluster.tree[[i-1]] <- tmp$cluster
				#bic.vector[k-1] <- nrow(x.matrix)*log(sum(tmp$withinss)/nrow(x.matrix)) + k*log(nrow(x.matrix))
				rm(tmp)
			}
			#print(bic.vector)

			# calculate cluster silhouettes  
			sil.k <- as.vector(NULL)	# vector of silhouettes
			dmatrix <- 1-annot.proximity	# dissimilarity matrix computed from the afinity matrix

			for(i in 1:length(cluster.tree)){ 
				# cut the tree 
				memb <- cluster.tree[[i]]
				sil <- silhouette(memb, dmatrix=dmatrix)
				sil.k <- c(sil.k, mean(summary(sil)$clus.avg.width))
			}

			names(sil.k) <- 2:(ncol(x.matrix)-1)	# choose the max silhouette
			print(sil.k)

			# find the best partition 

			if(!is.null(k)){
				best.index <- k
				print(paste("Best index forced to: ",k,sep=""))
			}else{
				if(length(sil.k) > 6){
					best.index <- as.numeric(names(sil.k[sil.k == max(sil.k[1:6])]))
				}else{
					best.index <- as.numeric(names(sil.k[sil.k == max(sil.k)]))
				}
				print(paste("Best index: ",best.index,sep=""))		
			}


			if(best.index == 1){
				best.partition <- as.vector(matrix(1,nrow(x.matrix),1))
				names(best.partition) <- rownames(x.matrix)

				sil.part <- NULL
				sil.cluster <- NULL		
			}else{

				k.eigen <- x.matrix[,1:best.index]
				lngth <- apply(k.eigen,1, function(x) sqrt(sum(x*x)))
				y.matrix <- k.eigen/lngth
				
				y.matrix[is.nan(y.matrix)] <- 0
				y.matrix[is.na(y.matrix)] <- 0
				
				tmp <- kmeans(y.matrix, centers=best.index[1], iter.max = iter.max, nstart = nstart, algorithm="Hart")
				print(tmp)
				rm(tmp,k.eigen,lngth,y.matrix)

				best.partition <- cluster.tree[[best.index-1]]

				sil <- silhouette(best.partition, dmatrix=dmatrix)
				sil.part <- mean(summary(sil)$clus.avg.width)
				sil.cluster <- summary(sil)$clus.avg.width
			}


			cluster.length <- as.vector(NULL)


			id.cluster <- list(NULL)	
			term.cluster <- list(NULL)
			gene.cluster <- list(NULL)

			for(i in 1:best.index){

				id.cluster[[i]] <- as.character(names(best.partition[best.partition == i]))
				cluster.length <- c(cluster.length, length(id.cluster[[i]]))

			}


			if(length(id.cluster) > 0){

				rownames(terms.name) <- terms.name[,1]

				for(i in 1:length(id.cluster)){
					term.cluster[[i]] <- terms.name[id.cluster[[i]],2]

					gene.vect.matrix <- annot.matrix[id.cluster[[i]],]	# annotations on rows, genes on columns
					gene.vect <- NULL

					if(is.matrix(gene.vect.matrix)){

						for(j in 1:nrow(gene.vect.matrix)){
							gene.vect <- c(gene.vect, names(gene.vect.matrix[j,][gene.vect.matrix[j,] == 1]))			
						}
					}else{
						gene.vect <- c(gene.vect, names(gene.vect.matrix[gene.vect.matrix == 1]))	
					}

					gene.cluster[[i]] <- unique(gene.vect)
				}	

				
				
				# extract connectivity information
				
				gene.connect <- matrix(0,nrow(adj.matrix),length(id.cluster)+3)
				rownames(gene.connect) <- rownames(adj.matrix)
				
				colnames(gene.connect) <- c("annotation_module","assigned_module",as.character(1:length(id.cluster)),"total_net")

				for(i in 1:length(id.cluster)){
					gene.connect[as.character(gene.cluster[[i]]),"annotation_module"] <- i

				}

				gene.connect[,"assigned_module"] <- gene.connect[,"annotation_module"]

				
				if(!is.null(marker.matrix)){	# if the proximity was dynamically inferred...
				
					for(i in 1:nrow(gene.connect)){

						for(j in 1:length(id.cluster)){
							x <- marker.matrix[rownames(gene.connect)[i],as.character(id.cluster[[j]])]
							gene.connect[i,as.character(j)] <- sum(x)

						}

						gene.connect[i,"total_net"] <- sum(gene.connect[i,as.character(1:length(id.cluster))])

					}
					
				}else{
					for(i in 1:nrow(gene.connect)){
						if(rownames(gene.connect)[i] %in% colnames(annot.matrix)){
							for(j in 1:length(id.cluster)){
								x <- annot.matrix[as.character(id.cluster[[j]]),rownames(gene.connect)[i]]
								gene.connect[i,as.character(j)] <- sum(x)

							}

							gene.connect[i,"total_net"] <- sum(gene.connect[i,as.character(1:length(id.cluster))])
						}
					}

				}
				
				
				for(i in 1:nrow(gene.connect)){
					if(gene.connect[i,"assigned_module"] == 0){
						x <- gene.connect[i,as.character(1:length(id.cluster))]
						names(x) <- 1:length(x)
						gene.connect[i,"assigned_module"] <- as.numeric(names(x)[x == max(x)][1])
					}
				}


				for(i in 1:length(id.cluster)){
					lngth <- sqrt(sum(gene.connect[,as.character(i)] * gene.connect[,as.character(i)]))
					gene.connect[,as.character(i)] <- gene.connect[,as.character(i)]/lngth
					rm(lngth)
				}

				# total_net
				lngth <- sqrt(sum(gene.connect[,"total_net"] * gene.connect[,"total_net"]))
				gene.connect[,"total_net"] <- gene.connect[,"total_net"]/lngth
				rm(lngth)

				names(best.partition) <- terms.name[names(best.partition),2]
				print(sort(best.partition))

				clusters <- list(annot.proximity=annot.proximity, sil.part=sil.part, sil=sil.k, sil.cluster=sil.cluster, 
							cluster.length=cluster.length, best.index=best.index, id.cluster=id.cluster, 
							term.cluster=term.cluster, gene.cluster=gene.cluster, best.partition=best.partition, 
							gene.connect=gene.connect, iter.max=iter.max, terms.name=terms.name, 
							eigen.centrality=eigen.centrality)

			}else{
				clusters <- NULL
				print("No clusters were created...")
			}
		
		}else{
			clusters <- NULL
			print("No clusters were created...")
		}

	}else{
		
		clusters <- NULL
		print("No clusters were created...")
	}	
	
	return(clusters)

}



#############################################################################################
#
# 9. Function .genes.centrality() -> Buiding co-expression nets and calculating centrality measures in gene co-expression nets
#
#############################################################################################

.genes.centrality <- function(adj.matrix,clusters=NULL,taxoname,locus.symbol=NULL,results.dir=NULL,coexp=FALSE){
	
	centrality <- NULL
	

	if(!is.null(clusters) & !coexp){

		annot.genes <- rownames(clusters$gene.connect)[as.numeric(clusters$gene.connect[,"annotation_module"])>0]
		connect.genes <- rownames(clusters$gene.connect)[as.numeric(clusters$gene.connect[,"total_net"])>0]

		adj.matrix.annot <- adj.matrix[annot.genes,annot.genes]
		adj.matrix.connect <- adj.matrix[connect.genes,connect.genes]

		rownames(adj.matrix.annot) <- locus.symbol[rownames(adj.matrix.annot),2]
		colnames(adj.matrix.annot) <- locus.symbol[colnames(adj.matrix.annot),2]

		rownames(adj.matrix.connect) <- locus.symbol[rownames(adj.matrix.connect),2]
		colnames(adj.matrix.connect) <- locus.symbol[colnames(adj.matrix.connect),2]

		.cyto.sym(net.matrix=adj.matrix.annot,file.net=paste(getwd(),"/",results.dir,"/",taxoname,"_annotated_genes_net.txt",sep=""),
			diagonal=FALSE,thresh=NULL)
		.cyto.sym(net.matrix=adj.matrix.connect,file.net=paste(getwd(),"/",results.dir,"/",taxoname,"_connected_genes_net.txt",sep=""),
			diagonal=FALSE,thresh=NULL)

		centrality <- matrix(NA,nrow(adj.matrix),8)
		rownames(centrality) <- rownames(adj.matrix)
		colnames(centrality) <- c("modular_degree","scaled_modular_degree","total_degree","scaled_total_degree",
						"modular_betweenness","scaled_modular_betweenness","total_betweenness","scaled_total_betweenness")

		adj.matrix.modular <- adj.matrix
		adj.matrix.modular[!(rownames(adj.matrix.modular) %in% annot.genes),!(colnames(adj.matrix.modular) %in% annot.genes)] <- 0

		centrality[,1] <- degree(adj.matrix.modular, gmode="graph", diag=FALSE, rescale=FALSE)
		centrality[,2] <- degree(adj.matrix.modular, gmode="graph", diag=FALSE, rescale=TRUE)
		centrality[,3] <- degree(adj.matrix, gmode="graph", diag=FALSE, rescale=FALSE)
		centrality[,4] <- degree(adj.matrix, gmode="graph", diag=FALSE, rescale=TRUE)
		centrality[,5] <- betweenness(adj.matrix.modular, gmode="graph", diag=FALSE, cmode="undirected", rescale=FALSE)
		centrality[,6] <- betweenness(adj.matrix.modular, gmode="graph", diag=FALSE, cmode="undirected", rescale=TRUE)
		centrality[,7] <- betweenness(adj.matrix, gmode="graph", diag=FALSE, cmode="undirected", rescale=FALSE)
		centrality[,8] <- betweenness(adj.matrix, gmode="graph", diag=FALSE, cmode="undirected", rescale=TRUE)
	
	}else if(coexp){
		centrality <- matrix(NA,nrow(adj.matrix),4)
		rownames(centrality) <- rownames(adj.matrix)
		colnames(centrality) <- c("total_degree","scaled_total_degree","total_betweenness","scaled_total_betweenness")
		
		centrality[,1] <- degree(adj.matrix, gmode="graph", diag=FALSE, rescale=FALSE)
		centrality[,2] <- degree(adj.matrix, gmode="graph", diag=FALSE, rescale=TRUE)
		centrality[,3] <- betweenness(adj.matrix, gmode="graph", diag=FALSE, cmode="undirected", rescale=FALSE)
		centrality[,4] <- betweenness(adj.matrix, gmode="graph", diag=FALSE, cmode="undirected", rescale=TRUE)
	}
	
	
	return(centrality)
}



#############################################################################################
#
# 10. Function .annot.recovery() -> Experimental: tests the recovery of known functional annotations by the dynamic system
#
#############################################################################################


.annot.recovery <- function(annot.matrix=NULL, adj.matrix=NULL, RV=0.90, method="sum", replace.annot=10){
	
		
	connectivity <- apply(adj.matrix,1,sum)
	connectivity <- connectivity[connectivity > 1]
	
	annotated <- apply(annot.matrix,1,sum)
	annotated <- annotated[annotated > 0]
	
	annot.connect <- annotated[names(annotated) %in% names(connectivity)]
	
	annot.connect[] <- runif(n=length(annot.connect))
	
	if(replace.annot >0){
		select.genes <- names(annot.connect[annot.connect<=(replace.annot/100)])
	}else{
		select.genes <- NULL
	}
	
	try(print(select.genes))
	
	discovered <- 0 # discovered removed true annotations
	initial <- 0 # initial true annotations
	supplementary <- 0 # newly infered annotations
	
	if(replace.annot > 0 && length(select.genes) >0){
		annot.matrix.new <- annot.matrix
		annot.matrix.new[select.genes,] <- 0
	
		marker.matrix <- .dynamic.proximity(annot.matrix=annot.matrix.new, adj.matrix=adj.matrix, RV=RV, method="sum")$marker.matrix
		marker.matrix <- marker.matrix[[length(marker.matrix)]]
	
		marker.matrix[marker.matrix >0] <- 1
		discovered <- sum((marker.matrix[select.genes,] * annot.matrix[select.genes,]),na.rm=TRUE)
		initial <- sum(annot.matrix[select.genes,],na.rm=TRUE)
		supplementary <- sum((marker.matrix[select.genes,] - annot.matrix[select.genes,]),na.rm=TRUE)
		
		print("")
		try(print(paste("Recovered: ",discovered * 100/initial,sep="")))
		print("")
		try(print(paste("Newly infered: ",supplementary*100/initial,sep="")))
		print("")
		
	}else if(replace.annot == 0){
		marker.matrix <- .dynamic.proximity(annot.matrix=annot.matrix.new, adj.matrix=adj.matrix, RV=RV, method="sum")$marker.matrix
		marker.matrix <- marker.matrix[[length(marker.matrix)]]
	
		marker.matrix[marker.matrix >0] <- 1
		initial <- sum(annot.matrix,na.rm=TRUE)
		supplementary <- sum((marker.matrix - annot.matrix),na.rm=TRUE)
		
		print("")
		try(print(paste("Newly infered: ",supplementary*100/initial,sep="")))
		print("")
		
	}else{
		
		print("No endpoint")	
	}
	return(list(discovered=discovered,initial=initial,supplementary=supplementary))

}



#############################################################################################
#
# 11. Function .annot.robustness() -> Experimental: tests the robustness to missing functional annotations
#
#############################################################################################


.annot.robustness <- function(annot.matrix=NULL, adj.matrix=NULL, terms.name=NULL, taxoname, RV=0.90, method="sum", 
			spectral=TRUE, replace.annot=0, annot.prox.measure=NULL, sigma=NA){
	
	
	#if(annot.prox.measure=="dynamical"){diag(adj.matrix) <- 0}
	
	connectivity <- apply(adj.matrix,1,sum)
	connectivity <- connectivity[connectivity > 1]
	
	annotated <- apply(annot.matrix,1,sum)
	annotated <- annotated[annotated > 0]
	
	annot.connect <- annotated[names(annotated) %in% names(connectivity)]
	annot.connect[] <- runif(n = length(annot.connect))
	
	# start replacement of genes annotations and compute the resulting spectral partitions
	
	best.part.list <- NULL
	
	for(i in 0:replace.annot){
		
		no <- i/100
		
		print("")
		print(no * 100)
		print("")

		if(no > 0){
			select.genes <- names(annot.connect[annot.connect <= no])
		}else{
			select.genes <- NULL
		}

	
		if(no > 0 && length(select.genes) >0){
			annot.matrix.new <- annot.matrix
			annot.matrix.new[select.genes,] <- 0
			
			best.part <- NULL

			if(spectral){
				
				#annot.proximity <- NULL
			
				if(annot.prox.measure=="dynamical"){
					try(proximity <- .dynamic.proximity(annot.matrix=annot.matrix.new, adj.matrix=adj.matrix, RV=RV, 
									method="sum", annotated.only=TRUE, sigma=sigma))				
				
					try(best.part <- .spectral.cluster(annot.proximity=proximity$annot.proximity, terms.name=terms.name, 
									adj.matrix=adj.matrix, nstart=10, annot.matrix=t(annot.matrix.new), 
									marker.matrix=proximity$marker.matrix[[length(proximity$marker.matrix)]],					
									iter.max=20000)$best.partition)
				}else{
					try(annot.distance <- .annotation.distance(annot.matrix=t(annot.matrix.new), coexp.matrix=adj.matrix, measure=annot.prox.measure))
					try(proximity <- .presumptive.proximity(annot.distance=annot.distance, sigma=sigma))
					
					try(best.part <- .spectral.cluster(annot.proximity=proximity, terms.name=terms.name, 
									adj.matrix=adj.matrix, nstart=10, annot.matrix=t(annot.matrix.new), 
									marker.matrix=NULL, iter.max=20000)$best.partition)
	
				}
				
			}else{
				
				try(best.part <- .uc.knn(annot.matrix=t(annot.matrix.new), coexp.matrix=adj.matrix, alpha = 0.05, taxoname=taxoname, 
						annot.prox.measure = annot.prox.measure, normal = TRUE, terms.name=terms.name)$best.partition)
				
			}

			if(is.null(best.part.list)){

				best.part.list <- list(NULL)
				best.part.list[[1]] <- best.part
			}else{


				best.part.list[[length(best.part.list)+1]] <- best.part

			}

			rm(best.part)
		

		}else if(no == 0){
			
			best.part <- NULL

			if(spectral){
				
				#annot.proximity <- NULL
			
				if(annot.prox.measure=="dynamical"){
					try(proximity <- .dynamic.proximity(annot.matrix=annot.matrix, adj.matrix=adj.matrix, RV=RV, 
									method="sum", annotated.only=TRUE, sigma=sigma))				

					try(best.part <- .spectral.cluster(annot.proximity=proximity$annot.proximity, terms.name=terms.name, 
									adj.matrix=adj.matrix, nstart=10, annot.matrix=t(annot.matrix), 
									marker.matrix=proximity$marker.matrix[[length(proximity$marker.matrix)]],					
									iter.max=20000)$best.partition)
				}else{
					try(annot.distance <- .annotation.distance(annot.matrix=t(annot.matrix), coexp.matrix=adj.matrix, 
											measure=annot.prox.measure))
					try(proximity <- .presumptive.proximity(annot.distance=annot.distance, sigma=sigma))
					
					try(best.part <- .spectral.cluster(annot.proximity=proximity, terms.name=terms.name, 
									adj.matrix=adj.matrix, nstart=10, annot.matrix=t(annot.matrix), 
									marker.matrix=NULL, iter.max=20000)$best.partition)

				}
				
			}else{
				
				try(best.part <- .uc.knn(annot.matrix=t(annot.matrix), coexp.matrix=adj.matrix, alpha = 0.05, taxoname=taxoname, 
						annot.prox.measure = annot.prox.measure, normal = TRUE, terms.name=terms.name)$best.partition)
				
			}

			if(is.null(best.part.list)){

				best.part.list <- list(NULL)
				best.part.list[[1]] <- best.part
			}else{


				best.part.list[[length(best.part.list)+1]] <- best.part

			}

			rm(best.part)
			
		}

	}
	
	part.sym <- as.vector(matrix(NA,length(best.part.list),1))
	
	
	for(i in 1:length(best.part.list)){
	
		if(!is.null(best.part.list[[i]])){
			try(part.sym[i] <- .part.distance(part1=best.part.list[[1]],part2=best.part.list[[i]]))
		}else{
			part.sym[i] <- NA
		}
		
	}
	
	try(plot(part.sym))
	
	return(list(part.sym=part.sym, best.part.list=best.part.list))
	

}



#############################################################################################
#
# 12. Function .hc.cluster() -> Experimental: hierarchical clustering of annotations diffusion profiles
#
#############################################################################################


# experimental!!!

.hc.cluster <- function(marker.matrix=NULL,steps=NULL,terms.name=NULL,annot.matrix=NULL){
	
	
	if(!is.null(marker.matrix) && !is.null(annot.matrix) && !is.null(terms.name)){
		
		if(is.null(steps)){steps <- length(marker.matrix)}
		
		marker.matrix <- marker.matrix[[steps]]

		color <- marker.matrix
		
		dist.annotation <- dist(t(color), method="euclid")
		hc <- hclust(dist.annotation, method = "complete")
  		.sil.hc <- as.vector(NULL)	# vector of silhouettes

 		# calculate cluster silhouettes
  
  		for(i in 2:(ncol(color)-1)){ 
    			# cut the tree 
    			memb <- cutree(hc, k = i)
    			sil <- silhouette(memb, dist.annotation)
    			.sil.hc <- c(.sil.hc, mean(summary(sil)$clus.avg.width))
  		}
		
		names(.sil.hc) <- 2:(ncol(color)-1)	# choose the max silhouette
		
		# find the best partition 
		
		best.index <- as.numeric(names(.sil.hc[.sil.hc == max(.sil.hc)]))
		
		print(paste("Best index: ",best.index,sep=""))
		
		best.partition <- cutree(hc, k = best.index)
	
		sil <- silhouette(best.partition, dist.annotation)
		sil.part <- mean(summary(sil)$clus.avg.width)
		sil.cluster <- summary(sil)$clus.avg.width
		cluster.length <- as.vector(NULL)
		
		
		id.cluster <- list(NULL)	
		term.cluster <- list(NULL)
		gene.cluster <- list(NULL)
		
		for(i in 1:best.index){
			
			id.cluster[[i]] <- as.character(names(best.partition[best.partition == i]))
			cluster.length <- c(cluster.length, length(id.cluster[[i]]))
		
		}
		
		if(length(id.cluster) > 0){
			
			rownames(terms.name) <- terms.name[,1]
			
			for(i in 1:length(id.cluster)){
				term.cluster[[i]] <- terms.name[id.cluster[[i]],2]
				
				gene.vect.matrix <- annot.matrix[id.cluster[[i]],]	# annotations on rows, genes on columns
				gene.vect <- NULL
				
				
				if(is.matrix(gene.vect.matrix)){
				
					for(j in 1:nrow(gene.vect.matrix)){
						gene.vect <- c(gene.vect, names(gene.vect.matrix[j,][gene.vect.matrix[j,] == 1]))			
					}
				}else{
					gene.vect <- c(gene.vect, names(gene.vect.matrix[gene.vect.matrix == 1]))	
				}
				
				gene.cluster[[i]] <- unique(gene.vect)
			}	
			
			gene.connect <- matrix(0,nrow(marker.matrix),length(id.cluster)+2)
			rownames(gene.connect) <- rownames(marker.matrix)
			colnames(gene.connect) <- c("module",1:length(id.cluster),"total_net")
				
			for(i in 1:length(id.cluster)){
				gene.connect[as.character(gene.cluster[[i]]),1] <- i
				
			}
				
				
			for(i in 1:nrow(gene.connect)){
					
				for(j in 1:length(id.cluster)){
					x <- marker.matrix[rownames(gene.connect)[i],as.character(id.cluster[[j]])]
					gene.connect[i,j+1] <- sum(x)
						
				}
				
				gene.connect[i,length(gene.cluster)+2] <- sum(gene.connect[i,2:(length(id.cluster)+1)])
				
			}
			
			names(best.partition) <- terms.name[names(best.partition),2]

			clusters <- list(dist.annotation=dist.annotation, annot.proximity=(1-dist.annotation), marker.matrix=marker.matrix,
						sil.part=sil.part,sil.cluster=sil.cluster,cluster.length=cluster.length,best.index=best.index,
						id.cluster=id.cluster,term.cluster=term.cluster,gene.cluster=gene.cluster,best.partition=best.partition,
						gene.connect=gene.connect)
			
		}else{
			clusters <- NULL
			print("No clusters were created...")
		}
			

	}else{
		
		clusters <- NULL
		print("No clusters were created...")
	}	
	
	return(clusters)
}



#############################################################################################
#
# 13. Function .sym.annotation() -> Experimental: computes a similarity matrix among transcripts based on co-annotations
#
#############################################################################################


.sym.annotation <- function(annot.matrix=NULL){
	
	adj.matrix <- NULL
	
	if(!is.null(annot.matrix)){ # transcripts on lines annotations on columns
		
		adj.matrix <- matrix(0,nrow(annot.matrix),nrow(annot.matrix))
		colnames(adj.matrix) <- rownames(annot.matrix)
		rownames(adj.matrix) <- rownames(annot.matrix)
		
		for(i in 1:nrow(adj.matrix)){
			
			for(j in 1:nrow(adj.matrix)){
				
				x <- annot.matrix[rownames(adj.matrix)[i],] + annot.matrix[colnames(adj.matrix)[j],]
				
				if(length(x[x==2])>0){adj.matrix[i,j] <- length(x[x==2])}
			
			}
		
		
		}
		diag(adj.matrix) <- 0
		threshold <- mean(as.vector(adj.matrix)) + 2*sd(as.vector(adj.matrix))
		
		adj.matrix[adj.matrix >= threshold] <- 1
		adj.matrix[adj.matrix < threshold] <- 0
	}

	return(adj.matrix)
}



#############################################################################################
#
# 14. Function .part.distance() -> Compute a similarity between two partition vectors involving the same objects based on Jaccard's index
#
#############################################################################################


.part.distance <- function(part1=NULL, part2=NULL){
	
	p.dist <- NULL
	#try(part1 <- part1[as.character(unique(names(part1)))])
	#try(part2 <- part1[as.character(unique(names(part2)))])
	
	if(!is.null(part1) && !is.null(part2) && (length(part1) == length(part2))){
		
		m1 <- matrix(0,length(part1),length(part1))
		m2 <- matrix(0,length(part2),length(part2))
		rownames(m1) <- names(part1)
		colnames(m1) <- names(part1)
		rownames(m2) <- names(part2)
		colnames(m2) <- names(part2)
		
		for(i in 1:nrow(m1)){
			
			for(j in 1:ncol(m1)){
			
				if(part1[rownames(m1)[i]] == part1[colnames(m1)[j]]){m1[i,j] <- 1}
				if(part2[rownames(m2)[i]] == part2[colnames(m2)[j]]){m2[i,j] <- 1}
			}
		
		
		}
		
		m2 <- m2[rownames(m1),colnames(m1)]
		
		#diag(m1) <- NA
		#diag(m2) <- NA
		
		m <- abs(m2-m1)
		
		p.dist <- length(m[m==0]) / (length(m[m==0]) + length(m[m==1]))
		
		print(paste("Similarity: ",p.dist,sep=""))
	
	
	}

	return(p.dist)

}




#############################################################################################
#
# 15. Functions .cytoscape.two() & .cytoscape.one() & cyto.sim() -> Internal routines which write Cytoscape format network files
#
#############################################################################################


# for results obtained for two lists of genes

.cytoscape.two <- function(clusters,file.net,file.net.info,file.param,up.annot.matrix,down.annot.matrix,annot.clust.method,threshold=NULL){
	
	# file.param to store proximity parameters
	# file.annot to store network supplementary info: column 1 = categories names, column two = cluster belonging, column 3 = Up(1) or Down(0)
	# file.net to store the network
	
	terms.name <- clusters$terms.name
	rownames(terms.name) <- terms.name[,1]

	annot.proximity <- clusters$annot.proximity
	colnames(annot.proximity) <- clusters$terms.name[colnames(annot.proximity),2]
	rownames(annot.proximity) <- clusters$terms.name[rownames(annot.proximity),2]
	
	#### normalize the connectivity matrix
	norm.connect <- .normalize.connect(clusters=clusters)
	norm.connect$strength.matrix <- norm.connect$strength.matrix[rownames(annot.proximity),colnames(annot.proximity)]
	norm.connect$intra.matrix <- norm.connect$intra.matrix[rownames(annot.proximity),colnames(annot.proximity)]
	
	.cyto.sym(net.matrix=annot.proximity,file.net=file.net,thresh=threshold,diagonal=FALSE,norm.connect=norm.connect)
	
	up.down <- as.vector(matrix(NA,length(clusters$best.partition),1))
	names(up.down) <- names(clusters$best.partition)
	
	if(!is.null(up.annot.matrix) & !is.null(down.annot.matrix)){
	
		rownames(up.annot.matrix) <- terms.name[rownames(up.annot.matrix),2]
		rownames(down.annot.matrix) <- terms.name[rownames(down.annot.matrix),2]

		for(i in 1:length(up.down)){
			if(names(up.down)[i] %in% rownames(up.annot.matrix)){
				up.down[i] <- 1
			}else if(names(up.down)[i] %in% rownames(down.annot.matrix)){
				up.down[i] <- 0
			}

			if(names(up.down)[i] %in% rownames(down.annot.matrix) & names(up.down)[i] %in% rownames(up.annot.matrix)){
				up.down[i] <- 2
			}
		}
	
	}else if(is.null(up.annot.matrix)){
		up.down[] <- 0
	}else if(is.null(down.annot.matrix)){
		up.down[] <- 1
	}
	
	
	adj.matrix <- clusters$annot.proximity
	
	rownames(adj.matrix) <- terms.name[rownames(adj.matrix),2]
	colnames(adj.matrix) <- terms.name[colnames(adj.matrix),2]

	centrality <- matrix(NA,nrow(adj.matrix),4)
	rownames(centrality) <- rownames(adj.matrix)
	colnames(centrality) <- c("degree","scaled_degree","betweenness","scaled_betweenness")
	
	if(annot.clust.method %in% c("umilds","spectral")){
		centrality[,"degree"] <- degree(adj.matrix, gmode="graph", diag=FALSE, rescale=FALSE)
		centrality[,"scaled_degree"] <- degree(adj.matrix, gmode="graph", diag=FALSE, rescale=TRUE)
		centrality[,"betweenness"] <- betweenness(adj.matrix, gmode="graph", diag=FALSE, cmode="undirected", rescale=FALSE)
		centrality[,"scaled_betweenness"] <- betweenness(adj.matrix, gmode="graph", diag=FALSE, cmode="undirected", rescale=TRUE)
	}else if(annot.clust.method=="ucknn"){
		centrality[,"degree"] <- degree(adj.matrix, gmode="digraph", diag=FALSE, rescale=FALSE)
		centrality[,"scaled_degree"] <- degree(adj.matrix, gmode="digraph", diag=FALSE, rescale=TRUE)
		centrality[,"betweenness"] <- betweenness(adj.matrix, gmode="digraph", diag=FALSE, cmode="directed", rescale=FALSE)
		centrality[,"scaled_betweenness"] <- betweenness(adj.matrix, gmode="digraph", diag=FALSE, cmode="directed", rescale=TRUE)
	}
	
	
	centrality[,"scaled_betweenness"][is.na(centrality[,"scaled_betweenness"])] <- 0
	centrality[,"scaled_betweenness"][centrality[,"scaled_betweenness"] == 0] <- (100 - sum(centrality[,"scaled_betweenness"]))/length(centrality[,"scaled_betweenness"][centrality[,"scaled_betweenness"] == 0])
	
	annot.info <- data.frame(rownames(centrality),clusters$best.partition[rownames(centrality)],up.down[rownames(centrality)],centrality)
	write.table(annot.info,file=file.net.info,col.names=c("name","module","up(1)_down(0)",colnames(centrality)),row.names=F,sep="\t")
	
	
	cat(paste("The median proximity among annotations is: ",median(clusters$annot.proximity),"\n",sep=""),file=file.param,append=TRUE)
	cat(paste("The upper quartile of the proximity among annotations is: ",
			median(clusters$annot.proximity[clusters$annot.proximity >= median(clusters$annot.proximity)]),sep=""),
			file=file.param,append=TRUE)

	
}

# for results obtained for one list of genes

.cytoscape.one <- function(clusters,file.net,file.net.info,file.param,terms.name,annot.clust.method,threshold=NULL){

	# file.param to store proximity parameters
	# file.annot to store network supplementary info: column 1 = categories names, column two = cluster belonging
	# file.net to store the network
	
	terms.name <- clusters$terms.name
	rownames(terms.name) <- terms.name[,1]

	annot.proximity <- clusters$annot.proximity
	colnames(annot.proximity) <- clusters$terms.name[colnames(annot.proximity),2]
	rownames(annot.proximity) <- clusters$terms.name[rownames(annot.proximity),2]
	
	#### normalize the connectivity matrix
	norm.connect <- .normalize.connect(clusters=clusters)
	norm.connect$strength.matrix <- norm.connect$strength.matrix[rownames(annot.proximity),colnames(annot.proximity)]
	norm.connect$intra.matrix <- norm.connect$intra.matrix[rownames(annot.proximity),colnames(annot.proximity)]
		
	.cyto.sym(net.matrix=annot.proximity,file.net=file.net,thresh=threshold,diagonal=FALSE,norm.connect=norm.connect)
	
	
	adj.matrix <- clusters$annot.proximity

	rownames(adj.matrix) <- terms.name[rownames(adj.matrix),2]
	colnames(adj.matrix) <- terms.name[colnames(adj.matrix),2]

	centrality <- matrix(NA,nrow(adj.matrix),4)
	rownames(centrality) <- rownames(adj.matrix)
	colnames(centrality) <- c("degree","scaled_degree","betweenness","scaled_betweenness")

	if(annot.clust.method %in% c("umilds","spectral")){
		centrality[,"degree"] <- degree(adj.matrix, gmode="graph", diag=FALSE, rescale=FALSE)
		centrality[,"scaled_degree"] <- degree(adj.matrix, gmode="graph", diag=FALSE, rescale=TRUE)
		centrality[,"betweenness"] <- betweenness(adj.matrix, gmode="graph", diag=FALSE, cmode="undirected", rescale=FALSE)
		centrality[,"scaled_betweenness"] <- betweenness(adj.matrix, gmode="graph", diag=FALSE, cmode="undirected", rescale=TRUE)
	}else if(annot.clust.method=="ucknn"){
		centrality[,"degree"] <- degree(adj.matrix, gmode="digraph", diag=FALSE, rescale=FALSE)
		centrality[,"scaled_degree"] <- degree(adj.matrix, gmode="digraph", diag=FALSE, rescale=TRUE)
		centrality[,"betweenness"] <- betweenness(adj.matrix, gmode="digraph", diag=FALSE, cmode="directed", rescale=FALSE)
		centrality[,"scaled_betweenness"] <- betweenness(adj.matrix, gmode="digraph", diag=FALSE, cmode="directed", rescale=TRUE)
	}
	
	centrality[,"scaled_betweenness"][is.na(centrality[,"scaled_betweenness"])] <- 0
	centrality[,"scaled_betweenness"][centrality[,"scaled_betweenness"] == 0] <- (100 - sum(centrality[,"scaled_betweenness"]))/length(centrality[,"scaled_betweenness"][centrality[,"scaled_betweenness"] == 0])
			
	annot.info <- data.frame(rownames(centrality),clusters$best.partition[rownames(centrality)],centrality)
	write.table(annot.info,file=file.net.info,col.names=c("name","module",colnames(centrality)),row.names=F,sep="\t")
	
	cat(paste("The median proximity among annotations is: ",median(clusters$annot.proximity),"\n",sep=""),file=file.param,append=TRUE)
	cat(paste("The upper quartile of the proximity among annotations is: ",
			median(clusters$annot.proximity[clusters$annot.proximity >= median(clusters$annot.proximity)]),sep=""),
			file=file.param,append=TRUE)

	
}

# writing a Cytoscape format net file

.cyto.sym <- function(net.matrix,file.net,diagonal=FALSE,thresh=NULL,link.type="pp",norm.connect=NULL){
	
	if(!diagonal){diag(net.matrix) <- 0}
	
	# file header
	if(!(is.null(norm.connect))){
		write(paste("node_1","link_type",
				"node_2","link_strength","normalized_strength","intra_modular_link",sep="\t"),
				file = file.net,append=TRUE)
	}else{
		write(paste("node_1","link_type",
				"node_2","link_strength",sep="\t"),
				file = file.net,append=TRUE)
	}
	
	# file content
	for (i in 1:dim(net.matrix)[1]) {
		for (j in i:dim(net.matrix)[2]) {
			if(is.null(thresh)){
				if(net.matrix[i,j] != 0 & !(is.null(norm.connect))){
					write(paste(dimnames(net.matrix)[[1]][i],link.type,
							dimnames(net.matrix)[[2]][j],net.matrix[i,j],norm.connect$strength.matrix[i,j],norm.connect$intra.matrix[i,j],sep="\t"),
							file = file.net,append=TRUE)
				}else if(net.matrix[i,j] != 0){
					write(paste(dimnames(net.matrix)[[1]][i],link.type,
							dimnames(net.matrix)[[2]][j],net.matrix[i,j],sep="\t"),
							file = file.net,append=TRUE)
				}
			}else if(abs(net.matrix[i,j]) >= thresh){
				if(net.matrix[i,j] != 0 & !(is.null(norm.connect))){
					write(paste(dimnames(net.matrix)[[1]][i],link.type,
							dimnames(net.matrix)[[2]][j],net.matrix[i,j],norm.connect$strength.matrix[i,j],norm.connect$intra.matrix[i,j],sep="\t"),
							file = file.net,append=TRUE)
				}else if(net.matrix[i,j] != 0){
					write(paste(dimnames(net.matrix)[[1]][i],link.type,
							dimnames(net.matrix)[[2]][j],net.matrix[i,j],sep="\t"),
							file = file.net,append=TRUE)
				}
			}
		}
	}



}


# routine for normalizing the connectivity matrix

.normalize.connect <- function(clusters=NULL){

	strength.matrix <- NULL
	intra.matrix <- NULL
	
	
	if(!is.null(clusters)){
		term.cluster <- clusters$term.cluster
		annot.proximity <- clusters$annot.proximity
		terms.name <- clusters$terms.name
		terms.name <- terms.name[rownames(annot.proximity),]
		rownames(annot.proximity) <- terms.name[rownames(annot.proximity),2]
		colnames(annot.proximity) <- terms.name[colnames(annot.proximity),2]
		
		betw <- betweenness(annot.proximity, gmode="graph", diag=FALSE, cmode="undirected", rescale=FALSE)
		names(betw) <- rownames(annot.proximity)
		betw[betw==0] <- 1
		
		annot.proximity <- t(annot.proximity*betw)*betw

		strength.matrix <- matrix(0,nrow(annot.proximity),ncol(annot.proximity))
		dimnames(strength.matrix) <- dimnames(annot.proximity)

		intra.matrix <- matrix(0,nrow(annot.proximity),ncol(annot.proximity))
		dimnames(intra.matrix) <- dimnames(annot.proximity)
		

		for(i in 1:length(term.cluster)){

			prox.matrix <- annot.proximity[term.cluster[[i]],term.cluster[[i]]]
			med <- median(prox.matrix[upper.tri(prox.matrix)])
			upperquart <- sort(prox.matrix[upper.tri(prox.matrix)],decreasing=TRUE)[ceiling(length(prox.matrix[upper.tri(prox.matrix)])/4)]
			upperdec <- sort(prox.matrix[upper.tri(prox.matrix)],decreasing=TRUE)[ceiling(length(prox.matrix[upper.tri(prox.matrix)])/10)]
			upperfive <- sort(prox.matrix[upper.tri(prox.matrix)],decreasing=TRUE)[ceiling(length(prox.matrix[upper.tri(prox.matrix)])/20)]
			
			strength.matrix[rownames(prox.matrix),colnames(prox.matrix)][prox.matrix < med] <- 0
			strength.matrix[rownames(prox.matrix),colnames(prox.matrix)][prox.matrix >= med] <- 1
			strength.matrix[rownames(prox.matrix),colnames(prox.matrix)][prox.matrix >= upperquart] <- 2
			strength.matrix[rownames(prox.matrix),colnames(prox.matrix)][prox.matrix >= upperdec] <- 3
			strength.matrix[rownames(prox.matrix),colnames(prox.matrix)][prox.matrix >= upperfive] <- 4
			
			test.connect <- strength.matrix[rownames(prox.matrix),colnames(prox.matrix)]
			diag(test.connect) <- 0
			test.connect <- apply(test.connect,1,sum)
			if(min(test.connect) == 0){			
				for(j in 1:length(test.connect)){
					if(test.connect[j]==0){
						strength.matrix[rownames(prox.matrix),colnames(prox.matrix)][j,prox.matrix[j,]==max(prox.matrix[j,])] <- 1
						strength.matrix[rownames(prox.matrix),colnames(prox.matrix)][prox.matrix[,j]==max(prox.matrix[,j]),j] <- 1
					}
				}
			}

			intra.matrix[rownames(prox.matrix),colnames(prox.matrix)] <- 1
		}

		med <- median(annot.proximity[intra.matrix == 0])
		upperquart <- sort(annot.proximity[intra.matrix == 0],decreasing=TRUE)[ceiling(length(annot.proximity[intra.matrix == 0])/4)]
		upperdec <- sort(annot.proximity[intra.matrix == 0],decreasing=TRUE)[ceiling(length(annot.proximity[intra.matrix == 0])/10)]
		upperfive <- sort(annot.proximity[intra.matrix == 0],decreasing=TRUE)[ceiling(length(annot.proximity[intra.matrix == 0])/20)]
			
		strength.matrix[intra.matrix == 0 & annot.proximity < med] <- 0
		strength.matrix[intra.matrix == 0 & annot.proximity >= med] <- 1
		strength.matrix[intra.matrix == 0 & annot.proximity >= upperquart] <- 2
		strength.matrix[intra.matrix == 0 & annot.proximity >= upperdec] <- 3
		strength.matrix[intra.matrix == 0 & annot.proximity >= upperfive] <- 4
		
	}
	
	return(list(strength.matrix=strength.matrix,intra.matrix=intra.matrix))	

}



#############################################################################################
#
# 16. Function .resmat() -> Organizing and saving of the results (isolated terms) taken from an annotation matrix
#
#############################################################################################

.resmat <- function(annot.matrix,print.data,terms.name,taxoname,nom,locus.name,results.dir="html",go=FALSE,bgcolor="blue",org="HS"){

	#	"nom" may be "UP" or "DOWN"
	#	"bgcolor" indicates the background color for KEGG map objects
	#	"go" is used to treat GO annotations on a one level annot.matrix (i.e. true matrix not list of matrixes)
	#	"org" indicates the organism "HS", "MM", "SC" for choosing the specific KEGG map

	cat(paste("\n\tSaving ",taxoname," ",nom," annotation results... ",format(Sys.time(), "%X"),sep=""))
	
	if(is.null(.funnet.version)){.funnet.version <- "experimental version"}
	rownames(terms.name) <- terms.name[,1]

if(is.list(annot.matrix)){
	
	for(i in 1:length(annot.matrix)){
		if(!is.null(annot.matrix[[i]])){

			if(nrow(annot.matrix[[i]]) > 1){
				print.data[[i]] <- print.data[[i]][order(print.data[[i]][,"exp.nr"],decreasing=TRUE),]
				annot.matrix[[i]] <- annot.matrix[[i]][rownames(print.data[[i]]),]
			}
			
			
			# list (as a string) of genes for each term
			
			list.genes <- matrix("",nrow(annot.matrix[[i]]),1)
			list.names <- matrix("",nrow(annot.matrix[[i]]),1)

			for(k in 1:nrow(annot.matrix[[i]])){
				ligne <- annot.matrix[[i]][k,]
				a <- names(ligne[ligne == 1])
				
				"" -> b
				"" -> d
				for(j in 1:length(a)){
					paste(b,a[j],'\n',sep="") -> b
					paste(d,"<option value='",locus.name[locus.name[,1] == a[j],2],"'>",locus.name[locus.name[,1] == a[j],2],"</option>",sep="") -> d
				
				}
				b -> list.genes[k,1]
				d -> list.names[k,1]
				rm(a,b,d)
			}
			
			# routine for saving the results of treatment for GO annotations as HTML files
			
				wd <- getwd()
				exp.term.lnk <- ""
				exp.id.lnk <- ""
			
			
				head <- paste("<html>
			
				<head>
				<meta http-equiv='Content-Type' content='text/html; charset=windows-1252'>
				<title>Genes ",nom," - ",taxoname," Annotations</title>
				</head>
			
				<body link='#0000FF' vlink='#0000FF' alink='#0000FF'>
			
				<table id='table0' style='width: 700px; border-collapse: collapse' cellPadding='0' border='0'>
					<tr>
						<td style='font-weight: 700; font-size: 10pt; vertical-align: middle; color: navy; font-style: normal; font-family: Arial, sans-serif; white-space: nowrap; height: 38pt; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
						<span lang='en-us' style='LETTER-SPACING: 2pt'><font size='3'>Genes ",nom," -
						",taxoname," Annotations</font></span></td>
					</tr>
					<tr>
						<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
						&nbsp;</td>
					</tr>
					<tr>
						<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
						<font face='Times New Roman' size='3'>
						<b><span lang='en-us'>Used notations:</span></b></font></td>
					</tr>
					<tr>
						<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
						&nbsp;</td>
					</tr>
					<tr>
						<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
						<ul>
							<li><font face='Times New Roman' size='3'><span lang='en-us'><b>List Hits</b> - the number of genes
							annotated by the considered ",taxoname," category or annotation cluster within
							the analyzed list of target genes</span> </font> </li>
						</ul>
						</td>
					</tr>
					<tr>
						<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
						<ul>
							<li><font face='Times New Roman' size='3'><span lang='en-us'><b>List </b></span><b>Total</b><span lang='en-us'>
							- the number of genes within the analyzed list of target genes
							having at least one ",taxoname," annotation</span> </font> </li>
						</ul>
						</td>
					</tr>
					<tr>
						<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
						<ul>
							<li><font face='Times New Roman' size='3'><b>Population<span lang='en-us'> Hits</span></b><span lang='en-us'>
							- the number of genes, available on the entire microarray, annotated
							by the considered ",taxoname," category or annotation cluster</span>
							</font>
							</li>
						</ul>
						</td>
					</tr>
					<tr>
						<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
						<ul>
							<li><font face='Times New Roman' size='3'><b>Population Total</b> - the number of genes available on the
							entire microarray and having at least one ",taxoname," annotation </font>
							</li>
						</ul>
						</td>
					</tr>
					<tr>
						<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
						<ul>
							<li><font face='Times New Roman' size='3'><span lang='en-us'><b>P-value</b> - the significance p-value of
							the gene enrichment of the considered ",taxoname," category or annotation
							cluster, calculated with a unilateral Fisher exact test</span>
							</font> </li>
						</ul>
						</td>
					</tr>
			
				</table>
				<p>",sep="")

				write(head, file = paste(wd,"/",results.dir,"/html/Genes ",nom," - ",taxoname," Annotations - Specificity Level ",i,".htm",sep=""), append = TRUE, sep = " ")
				
					table.head <- paste("
					<p>
				..........................................................................................................................</p>
					<p><b>Terms</b></p>
					<table border='1' id='table2' style='border-collapse: collapse' bordercolor='#808080'>
						<tr>
							<td align='center' width='45'><b>Rank</b></td>
							<td align='center' width='70'><b>ID</b></td>
							<td align='center' width='200'><b>Term</b></td>
							<td align='center' width='45'>
							<p style='margin-top: 0; margin-bottom: 0'><b>List </b></p>
							<p style='margin-top: 0; margin-bottom: 0'><b>Hits</b></td>
							<td align='center' width='45'>
							<p style='margin-top: 0; margin-bottom: 0'><b>List </b></p>
							<p style='margin-top: 0; margin-bottom: 0'><b>Total</b></td>
							<td align='center' width='70'>
							<p style='margin-top: 0; margin-bottom: 0'><b>Population </b></p>
							<p style='margin-top: 0; margin-bottom: 0'><b>Hits</b></td>
							<td align='center' width='70'>
							<p style='margin-top: 0; margin-bottom: 0'><b>Population </b></p>
							<p style='margin-top: 0; margin-bottom: 0'><b>Total</b></td>
							<td align='center' width='100'><b>P-value</b></td>
							<td align='center' width='100'><b>Genes ID's</b></td>
							<td align='center' width='200'><b>Genes Names</b></td>
					</tr>",sep="")
				
				write(table.head, file = paste(wd,"/",results.dir,"/html/Genes ",nom," - ",taxoname," Annotations - Specificity Level ",i,".htm",sep=""), append = TRUE, sep = " ")
				
				for(k in 1:nrow(print.data[[i]])){
				if(!is.na(terms.name[rownames(print.data[[i]])[k],2])){
					exp.term.lnk <- paste("<a target='_blank' href='http://www.ebi.ac.uk/ego/DisplayGoTerm?id=",rownames(print.data[[i]])[k],"' style='text-decoration: none'>
									<font size='3'>",terms.name[rownames(print.data[[i]])[k],2],"</font></a>",sep="")
					exp.id.lnk <- paste("<a target='_blank' href='http://www.ebi.ac.uk/ego/DisplayGoTerm?id=",rownames(print.data[[i]])[k],"' style='text-decoration: none'>
									<font size='3'>",rownames(print.data[[i]])[k],"</font></a>",sep="")
				
				
				
				
					list.terms <- paste("
						<tr>
						<td align='center' width='45'><font size='3'>",k,"</font></td>
						<td align='center' width='70'>
						",exp.id.lnk,"</td>
						<td align='center' width='200'>
						",exp.term.lnk,"</td>
						<td align='center' width='45'><font size='3'>",print.data[[i]][k,"exp.nr"],"</font></td>
						<td align='center' width='45'><font size='3'>",print.data[[i]][k,"exp.total"],"</font></td>
						<td align='center' width='70'><font size='3'>",print.data[[i]][k,"pop.hits"],"</font></td>
						<td align='center' width='70'><font size='3'>",print.data[[i]][k,"pop.total"],"</font></td>
						<td align='center' width='70'><font size='3'>",format(print.data[[i]][k,"pval"], digits = 3, scientific = TRUE),"</font></td>
						<td align='center' width='100'><textarea rows='0' name='Genes",k,"3' cols='7'>",list.genes[k,1],"</textarea></td>
						<td align='center' width='200'><select size='1' name='Genes",k,"4'>",list.names[k,1],"</select></td>
						</tr>",sep="")
				
				
				
					write(list.terms, file = paste(wd,"/",results.dir,"/html/Genes ",nom," - ",taxoname," Annotations - Specificity Level ",i,".htm",sep=""), append = TRUE, sep = " ")
				}
				}
				
				end.page <- paste("</table>
						<p>
						..........................................................................................................................</p>
						<p>
						This file was produced on ",date()," with <a target='_blank' href='http://corneliu.henegar.info/FunNet.htm' style='text-decoration: none'>
						<font size='3'>FunNet ",.funnet.version,"</font></a> based on GO and KEGG annotations updated on ",annot.date,".</p>
						</body>
				
						</html>",sep="")
				
				write(end.page, file = paste(wd,"/",results.dir,"/html/Genes ",nom," - ",taxoname," Annotations - Specificity Level ",i,".htm",sep=""), append = TRUE, sep = " ")
						
		
		}
	
	}
	
	
}else{
	
	if(!is.null(annot.matrix)){
				if(nrow(annot.matrix) > 1){
					print.data <- print.data[order(print.data[,"exp.nr"],decreasing=TRUE),]
					annot.matrix <- annot.matrix[rownames(print.data),]
				}
				
				# list (as a string) of genes for each term
				
				list.genes <- matrix("",nrow(annot.matrix),1)
				list.names <- matrix("",nrow(annot.matrix),1)
				kegg.genes <- NULL
				
	
				for(k in 1:nrow(annot.matrix)){
					ligne <- annot.matrix[k,]
					a <- names(ligne[ligne == 1])
					
					"" -> b
					"" -> d
					for(j in 1:length(a)){
						paste(b,a[j],'\n',sep="") -> b
						paste(d,"<option value='",locus.name[locus.name[,1] == a[j],2],"'>",locus.name[locus.name[,1] == a[j],2],"</option>",sep="") -> d
					
					}
					b -> list.genes[k,1]
					d -> list.names[k,1]
					rm(a,b,d)
				}
				
				if(!go){	# for KEGG annotations prepare genes to be placed on the maps
					
					if(org=="HS"){org<-"hsa"}
					if(org=="MM"){org<-"mmu"}
					if(org=="SC"){org<-"sce"}
					kegg.genes <- matrix("",nrow(annot.matrix),1)
					
					for(k in 1:nrow(annot.matrix)){
						ligne <- annot.matrix[k,]
						a <- names(ligne[ligne == 1])
										
						"" -> b
						for(j in 1:length(a)){	
							paste(b,a[j],'%09',bgcolor,',black/',sep="") -> b
										
						}
						b -> kegg.genes[k,1]
						rm(a,b)
					}				
				}
				
				# routine for saving the results of treatment for GO annotations as HTML files
				
					wd <- getwd()
					exp.term.lnk <- ""
					exp.id.lnk <- ""
				
				
					head <- paste("<html>
				
					<head>
					<meta http-equiv='Content-Type' content='text/html; charset=windows-1252'>
					<title>Genes ",nom," - ",taxoname," Annotations</title>
					</head>
				
					<body link='#0000FF' vlink='#0000FF' alink='#0000FF'>
				
					<table id='table0' style='width: 700px; border-collapse: collapse' cellPadding='0' border='0'>
						<tr>
							<td style='font-weight: 700; font-size: 10pt; vertical-align: middle; color: navy; font-style: normal; font-family: Arial, sans-serif; white-space: nowrap; height: 38pt; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
							<span lang='en-us' style='LETTER-SPACING: 2pt'><font size='3'>Genes ",nom," -
							",taxoname," Annotations</font></span></td>
						</tr>
						<tr>
							<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
							&nbsp;</td>
						</tr>
						<tr>
							<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
							<font face='Times New Roman' size='3'>
							<b><span lang='en-us'>Used notations:</span></b></font></td>
						</tr>
						<tr>
							<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
							&nbsp;</td>
						</tr>
						<tr>
							<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
							<ul>
								<li><font face='Times New Roman' size='3'><span lang='en-us'><b>List Hits</b> - the number of genes
								annotated by the considered ",taxoname," category or annotation cluster within
								the analyzed list of target genes</span> </font> </li>
							</ul>
							</td>
						</tr>
						<tr>
							<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
							<ul>
								<li><font face='Times New Roman' size='3'><span lang='en-us'><b>List </b></span><b>Total</b><span lang='en-us'>
								- the number of genes within the analyzed list of target genes
								having at least one ",taxoname," annotation</span> </font> </li>
							</ul>
							</td>
						</tr>
						<tr>
							<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
							<ul>
								<li><font face='Times New Roman' size='3'><b>Population<span lang='en-us'> Hits</span></b><span lang='en-us'>
								- the number of genes, available on the entire microarray, annotated
								by the considered ",taxoname," category or annotation cluster</span>
								</font>
								</li>
							</ul>
							</td>
						</tr>
						<tr>
							<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
							<ul>
								<li><font face='Times New Roman' size='3'><b>Population Total</b> - the number of genes available on the
								entire microarray and having at least one ",taxoname," annotation </font>
								</li>
							</ul>
							</td>
						</tr>
						<tr>
							<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
							<ul>
								<li><font face='Times New Roman' size='3'><span lang='en-us'><b>P-value</b> - the significance p-value of
								the gene enrichment of the considered ",taxoname," category or annotation
								cluster, calculated with a unilateral Fisher exact test</span>
								</font> </li>
							</ul>
							</td>
						</tr>
				
					</table>
					<p>",sep="")
	
					write(head, file = paste(wd,"/",results.dir,"/html/Genes ",nom," - ",taxoname," Annotations.htm",sep=""), append = TRUE, sep = " ")
					
						table.head <- paste("
						<p>
					..........................................................................................................................</p>
						<p><b>Terms</b></p>
						<table border='1' id='table2' style='border-collapse: collapse' bordercolor='#808080'>
							<tr>
								<td align='center' width='45'><b>Rank</b></td>
								<td align='center' width='70'><b>ID</b></td>
								<td align='center' width='200'><b>Term</b></td>
								<td align='center' width='45'>
								<p style='margin-top: 0; margin-bottom: 0'><b>List </b></p>
								<p style='margin-top: 0; margin-bottom: 0'><b>Hits</b></td>
								<td align='center' width='45'>
								<p style='margin-top: 0; margin-bottom: 0'><b>List </b></p>
								<p style='margin-top: 0; margin-bottom: 0'><b>Total</b></td>
								<td align='center' width='70'>
								<p style='margin-top: 0; margin-bottom: 0'><b>Population </b></p>
								<p style='margin-top: 0; margin-bottom: 0'><b>Hits</b></td>
								<td align='center' width='70'>
								<p style='margin-top: 0; margin-bottom: 0'><b>Population </b></p>
								<p style='margin-top: 0; margin-bottom: 0'><b>Total</b></td>
								<td align='center' width='100'><b>P-value</b></td>
								<td align='center' width='100'><b>Genes ID's</b></td>
								<td align='center' width='200'><b>Genes Names</b></td>
						</tr>",sep="")
					
					write(table.head, file = paste(wd,"/",results.dir,"/html/Genes ",nom," - ",taxoname," Annotations.htm",sep=""), append = TRUE, sep = " ")
					
					for(k in 1:nrow(print.data)){
						if(go){	# for GO annotations
					
							exp.term.lnk <- paste("<a target='_blank' href='http://www.ebi.ac.uk/ego/DisplayGoTerm?id=",rownames(print.data)[k],"' style='text-decoration: none'>
										<font size='3'>",terms.name[rownames(print.data)[k],2],"</font></a>",sep="")
							exp.id.lnk <- paste("<a target='_blank' href='http://www.ebi.ac.uk/ego/DisplayGoTerm?id=",rownames(print.data)[k],"' style='text-decoration: none'>
										<font size='3'>",rownames(print.data)[k],"</font></a>",sep="")
						}else{	# for KEGG annotations

							exp.term.lnk <- paste("<a target='_blank' href='http://www.genome.jp/kegg-bin/mark_pathway_www?@",org,rownames(print.data)[k],"/default%3dwhite/reference%3d%23FFFF33/",kegg.genes[k,1],"' style='text-decoration: none'>
										<font size='3'>",terms.name[rownames(print.data)[k],2],"</font></a>",sep="")
							exp.id.lnk <- paste("<a target='_blank' href='http://www.genome.jp/kegg-bin/mark_pathway_www?@",org,rownames(print.data)[k],"/default%3dwhite/reference%3d%23FFFF33/",kegg.genes[k,1],"' style='text-decoration: none'>
										<font size='3'>",rownames(print.data)[k],"</font></a>",sep="")
						
						}
					
					
					
						list.terms <- paste("
							<tr>
							<td align='center' width='45'><font size='3'>",k,"</font></td>
							<td align='center' width='70'>
							",exp.id.lnk,"</td>
							<td align='center' width='200'>
							",exp.term.lnk,"</td>
							<td align='center' width='45'><font size='3'>",print.data[k,"exp.nr"],"</font></td>
							<td align='center' width='45'><font size='3'>",print.data[k,"exp.total"],"</font></td>
							<td align='center' width='70'><font size='3'>",print.data[k,"pop.hits"],"</font></td>
							<td align='center' width='70'><font size='3'>",print.data[k,"pop.total"],"</font></td>
							<td align='center' width='70'><font size='3'>",format(print.data[k,"pval"], digits = 3, scientific = TRUE),"</font></td>
							<td align='center' width='100'><textarea rows='0' name='Genes",k,"3' cols='7'>",list.genes[k,1],"</textarea></td>
							<td align='center' width='200'><select size='1' name='Genes",k,"4'>",list.names[k,1],"</select></td>
							</tr>",sep="")
					
					
					
						write(list.terms, file = paste(wd,"/",results.dir,"/html/Genes ",nom," - ",taxoname," Annotations.htm",sep=""), append = TRUE, sep = " ")
					
					}
					
					end.page <- paste("</table>
							<p>
							..........................................................................................................................</p>
							<p>
							This file was produced on ",date()," with <a target='_blank' href='http://corneliu.henegar.info/FunNet.htm' style='text-decoration: none'>
							<font size='3'>FunNet ",.funnet.version,"</font></a> based on GO and KEGG annotations updated on ",annot.date,".</p>
							</body>
					
							</html>",sep="")
					
					write(end.page, file = paste(wd,"/",results.dir,"/html/Genes ",nom," - ",taxoname," Annotations.htm",sep=""), append = TRUE, sep = " ")
							
			
			}
		
	
}



	rm()
}




#############################################################################################
#
# 17. Function .profil.plot.two() -> Plot functional profiles as back-to-back bar plots (up and down-regulated themes)
#
#############################################################################################


.profil.plot.two <- function(up.annot=NULL, down.annot=NULL, terms.name=NULL, taxoname=NULL, height=NA, width=NA, 
				compute.dim=TRUE, dev ="png", extra=""){

	

	if(!is.null(down.annot) & !is.null(terms.name) & !is.null(up.annot)){
		
		if(!is.null(taxoname)){
		
			rownames(terms.name) <- terms.name[,1]
			
			graph.name <- NULL
			if(taxoname == "bp"){graph.name <- "GO Biological Process"}
			if(taxoname == "cc"){graph.name <- "GO Cellular Component"}
			if(taxoname == "mf"){graph.name <- "GO Molecular Function"}
			if(taxoname == "kegg"){graph.name <- "KEGG"}
			
			if(taxoname == "GO Biological Process"){graph.name <- "GO Biological Process"}
			if(taxoname == "GO Cellular Component"){graph.name <- "GO Cellular Component"}
			if(taxoname == "GO Molecular Function"){graph.name <- "GO Molecular Function"}
			if(taxoname == "KEGG"){graph.name <- "KEGG"}
			
			if(is.list(up.annot$annot.matrix) & is.list(down.annot$annot.matrix)){

				for(i in 1:min(length(up.annot$annot.matrix),length(down.annot$annot.matrix))){

					if(!is.null(up.annot$annot.matrix[[i]]) & !is.null(down.annot$annot.matrix[[i]])){
						
						up.matrix <- as.data.frame(up.annot$annot.matrix[[i]])
						down.matrix <- as.data.frame(down.annot$annot.matrix[[i]])
						up.matrix <- up.matrix[!(rownames(up.matrix) %in% c("GO:0008150","GO:0000004","GO:0007582","GO:0005575","GO:0008372","GO:0003674","GO:0005554")),]
						down.matrix <- down.matrix[!(rownames(down.matrix) %in% c("GO:0008150","GO:0000004","GO:0007582","GO:0005575","GO:0008372","GO:0003674","GO:0005554")),]
						
						x.up <- apply(up.matrix,1,sum)
						y.up <- apply(up.matrix,2,sum)
						y.up <- length(y.up[y.up != 0])
						x.down <- apply(down.matrix,1,sum)
						y.down <- apply(down.matrix,2,sum)
						y.down <- length(y.down[y.down != 0])

						wdth <- max(length(x.up), length(x.down))

						left <- sort(x.up*100/y.up, decreasing=TRUE)
						right <- sort(x.down*100/y.down, decreasing=TRUE)

						if(length(left)>12){left <- left[1:12] 
									wdth <- 12}
						if(length(right)>12){right <- right[1:12]
									wdth <- 12}
						if(compute.dim){
							max.left_right <- max(c(nchar(terms.name[names(left),2]),nchar(terms.name[names(right),2]),nchar("Transcriptional domain coverage (%)")),na.rm=TRUE)
							width <- 20*max.left_right
							if(width < 600){width <- 600}
						
							height <- 300 + 43*wdth
							#if(height < 450){height <- 450}
						}
						
						if(dev == "ps"){
							CairoPS(file=paste(extra,taxoname," Level ",i,".ps",sep=""), onefile = FALSE, width=width/100, height=height/100, horizontal=FALSE,
								paper="special", print.it=TRUE)
						}else if(dev == "pdf"){
							CairoPDF(file=paste(extra,taxoname," Level ",i,".pdf",sep=""), width=width/100, height=height/100, paper="special")
						}else{
							CairoPNG(file=paste(extra,taxoname," Level ",i,".png",sep=""), height=0.85*height, width=0.85*width, res=600)
						}
	

						par(mfrow=c(1,2), bg="white", las=1, adj=1, mar =c(6,1,5,0.5), ask=F)

						barplot(-left, col=rainbow(n=length(left), start=0, end=1/6), axisnames=FALSE, 
								beside=TRUE, axes=FALSE, horiz = TRUE, xlim=c(-100,0), 
								space=1, width=0.8, ylim=c(0,(0.4+1.7*wdth)))
						title(main = "Up-regulated Transcripts", font.main = 1, cex = 1, line=1)

						for(j in 1:length(left)){	
							text(y=1.8*j + 0.2 -(0.2*(j-1)), x = 0, labels = terms.name[names(left)[j],2], 
								pos=2, cex=1, offset = 0.7 )
						}
						

						for(j in 1:length(left)){
							if(left[j] <= 85){
								text(y=1*j + 0.2 +(0.6*(j-1)), x = -left[j], labels = paste(round(left[j],digits=1),"%",sep=""), 
									pos=2, cex=0.8, offset = 0.7, col="blue")
							}else{
								text(y=1*j + 0.2 +(0.6*(j-1)), x = -left[j], labels = paste(round(left[j],digits=1),"%",sep=""), 
									pos=4, cex=0.8, offset = 0.7, col="blue")
							}
						}	

						axis(side=1, at=c(-100,-80,-60,-40, -20, 0),  labels = c(100,80,60,40, 20, 0), 
							tick = TRUE, line = NA,
						     	pos = 0, outer = FALSE, font = NA,
						     	lty = "solid", lwd = 1, col = NULL, hadj = NA, padj = NA)


						par(mar=c(6,0.5,5,1), adj=0)

						barplot(right, col=sort(rainbow(n=length(right), start=2/6, end=3/6)), axisnames=FALSE, 
								beside=TRUE, horiz = TRUE, axes=FALSE, xlim=c(0,100), 
								space=1, width=0.8, ylim=c(0,(0.4+1.7*wdth)))
						title(main = "Down-regulated Transcripts", font.main = 1, cex = 1, line=1)

						for(j in 1:length(right)){	
							text(y=1.8*j + 0.2 -(0.2*(j-1)), x = 0, labels = terms.name[names(right)[j],2], 
								pos=4, cex=1, offset = 0.7 )
						}
						
						for(j in 1:length(right)){	
							if(right[j] <= 85){
								text(y=1*j + 0.2 +(0.6*(j-1)), x = right[j], labels = paste(round(right[j],digits=1),"%",sep=""), 
									pos=4, cex=0.8, offset = 0.7, col="blue")
							}else{
								text(y=1*j + 0.2 +(0.6*(j-1)), x = right[j], labels = paste(round(right[j],digits=1),"%",sep=""), 
									pos=2, cex=0.8, offset = 0.7, col="blue")
							}
						}


						axis(side=1, at=c(100,80,60,40, 20, 0),  labels = c(100,80,60,40, 20, 0), tick = TRUE, 
							line = NA,
						     	pos = 0, outer = FALSE, font = NA,
						     	lty = "solid", lwd = 1, col = NULL, hadj = NA, padj = NA)

						par(mfrow=c(1,1), new=TRUE, mar=c(6,1,2,1), las=1, adj=0.5, xaxp=c(0,201,0.8*wdth+1))
						title(main = graph.name, xlab="Transcriptional domain coverage (%)", font.main = 2, 
							cex.main = 1.5, sub="", font.sub=2, cex.lab=1)
						if(dev == "pdf"){
							abline(v=(50 + 1.2*900/width), lty="dotted")
						}else{
							abline(v=(50 + 0.7*1100/width), lty="dotted")
						}
						
						dev.off()
						
					}
				}
			}else{
				
				if(!is.null(up.annot$annot.matrix) & !is.null(down.annot$annot.matrix)){

					up.matrix <- as.data.frame(up.annot$annot.matrix)
					down.matrix <- as.data.frame(down.annot$annot.matrix)
					up.matrix <- up.matrix[!(rownames(up.matrix) %in% c("GO:0008150","GO:0000004","GO:0007582","GO:0005575","GO:0008372","GO:0003674","GO:0005554")),]
					down.matrix <- down.matrix[!(rownames(down.matrix) %in% c("GO:0008150","GO:0000004","GO:0007582","GO:0005575","GO:0008372","GO:0003674","GO:0005554")),]
					
					if(is.vector(up.matrix)){
						up.matrix <- t(as.matrix(up.matrix))
					}

					if(is.vector(down.matrix)){
						down.matrix <- t(as.matrix(down.matrix))
					}

					x.up <- apply(up.matrix,1,sum)
					y.up <- apply(up.matrix,2,sum)
					y.up <- length(y.up[y.up != 0])
					x.down <- apply(down.matrix,1,sum)
					y.down <- apply(down.matrix,2,sum)
					y.down <- length(y.down[y.down != 0])

					wdth <- max(length(x.up), length(x.down))

					left <- sort(x.up*100/y.up, decreasing=TRUE)
					right <- sort(x.down*100/y.down, decreasing=TRUE)

					if(length(left)>12){left <- left[1:12] 
								wdth <- 12}
					if(length(right)>12){right <- right[1:12]
								wdth <- 12}
					if(compute.dim){
						max.left_right <- max(c(nchar(terms.name[names(left),2]),nchar(terms.name[names(right),2]),nchar("Transcriptional domain coverage (%)")),na.rm=TRUE)
						width <- 20*max.left_right
						if(width < 600){width <- 600}
					
						height <- 300 + 43*wdth
						#if(height < 450){height <- 450}
					}
						
								
					if(dev == "ps"){
						CairoPS(file=paste(extra,taxoname,".ps",sep=""), onefile = FALSE, width=width/100, height=height/100, horizontal=FALSE,
							paper="special", print.it=TRUE)
					}else if(dev == "pdf"){
						CairoPDF(file=paste(extra,taxoname,".pdf",sep=""), width=width/100, height=height/100,	paper="special")
					}else{
						CairoPNG(file=paste(extra,taxoname,".png",sep=""), height=0.85*height, width=0.85*width, res=600)
					}

					par(mfrow=c(1,2), bg="white", las=1, adj=1, mar =c(6,1,5,0.5), ask=F)

					barplot(-left, col=rainbow(n=length(left), start=0, end=1/6), axisnames=FALSE, 
							beside=TRUE, axes=FALSE, horiz = TRUE, xlim=c(-100,0), 
							space=1, width=0.8, ylim=c(0,(0.4+1.7*wdth)))
					title(main = "Up-regulated Transcripts", font.main = 1, cex = 1, line=1)

					for(j in 1:length(left)){	
						text(y=1.8*j + 0.2 -(0.2*(j-1)), x = 0, labels = terms.name[names(left)[j],2], 
							pos=2, cex=1, offset = 0.7 )
					}


					for(j in 1:length(left)){
						if(left[j] <= 85){
							text(y=1*j + 0.2 +(0.6*(j-1)), x = -left[j], labels = paste(round(left[j],digits=1),"%",sep=""), 
								pos=2, cex=0.8, offset = 0.7, col="blue")
						}else{
							text(y=1*j + 0.2 +(0.6*(j-1)), x = -left[j], labels = paste(round(left[j],digits=1),"%",sep=""), 
								pos=4, cex=0.8, offset = 0.7, col="blue")
						}
					}

					axis(side=1, at=c(-100,-80,-60,-40, -20, 0),  labels = c(100,80,60,40, 20, 0), 
						tick = TRUE, line = NA,
						pos = 0, outer = FALSE, font = NA,
						lty = "solid", lwd = 1, col = NULL, hadj = NA, padj = NA)


					par(mar=c(6,0.5,5,1), adj=0)

					barplot(right, col=sort(rainbow(n=length(right), start=2/6, end=3/6)), axisnames=FALSE, 
							beside=TRUE, horiz = TRUE, axes=FALSE, xlim=c(0,100), 
							space=1, width=0.8, ylim=c(0,(0.4+1.7*wdth)))
					title(main = "Down-regulated Transcripts", font.main = 1, cex = 1, line=1)

					for(j in 1:length(right)){	
						text(y=1.8*j + 0.2 -(0.2*(j-1)), x = 0, labels = terms.name[names(right)[j],2], 
							pos=4, cex=1, offset = 0.7 )
					}

					for(j in 1:length(right)){	
						if(right[j] <= 85){
							text(y=1*j + 0.2 +(0.6*(j-1)), x = right[j], labels = paste(round(right[j],digits=1),"%",sep=""), 
								pos=4, cex=0.8, offset = 0.7, col="blue")
						}else{
							text(y=1*j + 0.2 +(0.6*(j-1)), x = right[j], labels = paste(round(right[j],digits=1),"%",sep=""), 
								pos=2, cex=0.8, offset = 0.7, col="blue")
						}
					}
					
					axis(side=1, at=c(100,80,60,40, 20, 0),  labels = c(100,80,60,40, 20, 0), tick = TRUE, 
						line = NA,
						pos = 0, outer = FALSE, font = NA,
						lty = "solid", lwd = 1, col = NULL, hadj = NA, padj = NA)

					par(mfrow=c(1,1), new=TRUE, mar=c(6,1,2,1), las=1, adj=0.5, xaxp=c(0,201,0.8*wdth+1))
					title(main = graph.name, xlab="Transcriptional domain coverage (%)", font.main = 2, 
						cex.main = 1.5, sub="", font.sub=2, cex.lab=1)
					if(dev == "pdf"){
						abline(v=(50 + 1.2*900/width), lty="dotted")
					}else{
						abline(v=(50 + 0.7*1100/width), lty="dotted")
					}

					dev.off()

				}
			}



		}
	
	}else{
		print("No data available for ploting functional profiles\\!")
	}
	
	
}



#############################################################################################
#
# 17. Function .profil.plot.one() -> Plot functional profiles as bar plots (one list of themes)
#
#############################################################################################


.profil.plot.one <- function(genes.annot=NULL, terms.name=NULL, taxoname=NULL, height=NA, width=NA, 
				compute.dim=TRUE, dev ="png", extra=""){

	

	if(!is.null(genes.annot) & !is.null(terms.name)){
		
		if(!is.null(taxoname)){
		
			rownames(terms.name) <- terms.name[,1]
			
			graph.name <- NULL
			if(taxoname == "bp"){graph.name <- "GO Biological Process"}
			if(taxoname == "cc"){graph.name <- "GO Cellular Component"}
			if(taxoname == "mf"){graph.name <- "GO Molecular Function"}
			if(taxoname == "kegg"){graph.name <- "KEGG"}
			
			if(taxoname == "GO Biological Process"){graph.name <- "GO Biological Process"}
			if(taxoname == "GO Cellular Component"){graph.name <- "GO Cellular Component"}
			if(taxoname == "GO Molecular Function"){graph.name <- "GO Molecular Function"}
			if(taxoname == "KEGG"){graph.name <- "KEGG"}
			
			if(is.list(genes.annot$annot.matrix)){

				for(i in 1:length(genes.annot$annot.matrix)){

					if(!is.null(genes.annot$annot.matrix[[i]])){

						genes.matrix <- as.data.frame(genes.annot$annot.matrix[[i]])
						genes.matrix <- genes.matrix[!(rownames(genes.matrix) %in% c("GO:0008150","GO:0000004","GO:0007582","GO:0005575","GO:0008372","GO:0003674","GO:0005554")),]
							
						if(is.vector(genes.matrix)){
							genes.matrix <- t(as.matrix(genes.matrix))
						}
						
						x.genes <- apply(genes.matrix,1,sum)
						y.genes <- apply(genes.matrix,2,sum)
						y.genes <- length(y.genes[y.genes != 0])
						

						wdth <- length(x.genes)

						right <- sort(x.genes*100/y.genes, decreasing=TRUE)

						
						if(length(right)>12){right <- right[1:12]
									wdth <- 12}
						
						if(compute.dim){
							max.right <- max(c(nchar(terms.name[names(right),2]),nchar("Transcriptional domain coverage (%)")))
							width <- 10.5*max.right
							if(width < 300){width <- 300}
						
							height <- 300 + 43*wdth
							#if(height < 450){height <- 450}
						}
						
						if(dev == "ps"){
							CairoPS(file=paste(extra,taxoname," Level ",i,".ps",sep=""), onefile = FALSE, width=width/100, height=height/100, horizontal=FALSE,
								paper="special", print.it=TRUE)
						}else if(dev == "pdf"){
							CairoPDF(file=paste(extra,taxoname," Level ",i,".pdf",sep=""), width=width/100, height=height/100, paper="special")
						}else{
							CairoPNG(file=paste(extra,taxoname," Level ",i,".png",sep=""), height=0.85*height, width=0.85*width, res=600)
						}
	

						par(mfrow=c(1,1), bg="white", las=1, adj=0, mar =c(6,1,5,1), ask=F)

						barplot(right, col=sort(rainbow(n=length(right), start=3/6, end=4/6)), axisnames=FALSE, 
							beside=TRUE, horiz = TRUE, axes=FALSE, xlim=c(0,100), 
							space=1, width=0.8, ylim=c(0,(0.4+1.7*wdth)))
						
						par(adj=0.5)
						title(main = graph.name, xlab="Transcriptional domain coverage (%)", font.main = 2, 
							cex.main = 1.5, sub="", font.sub=2, cex.lab=1)

						for(j in 1:length(right)){	
							text(y=1.8*j + 0.2 -(0.2*(j-1)), x = 0, labels = terms.name[names(right)[j],2], 
								pos=4, cex=1, offset = 0.7 )
						}
						
						for(j in 1:length(right)){	
							if(right[j] <= 85){
								text(y=1*j + 0.2 +(0.6*(j-1)), x = right[j], labels = paste(round(right[j],digits=1),"%",sep=""), 
									pos=4, cex=0.8, offset = 0.7, col="red")
							}else{
								text(y=1*j + 0.2 +(0.6*(j-1)), x = right[j], labels = paste(round(right[j],digits=1),"%",sep=""), 
									pos=2, cex=0.8, offset = 0.7, col="red")
							}
						}

						axis(side=1, at=c(100,80,60,40, 20, 0),  labels = c(100,80,60,40, 20, 0), tick = TRUE, 
							line = NA,	pos = NA, outer = FALSE, font = NA,
						     	lty = "solid", lwd = 1, col = NULL, hadj = NA, padj = NA)
						
						dev.off()
						
					}
				}
			}else{
				
				if(!is.null(genes.annot$annot.matrix)){

					genes.matrix <- as.data.frame(genes.annot$annot.matrix)
					genes.matrix <- genes.matrix[!(rownames(genes.matrix) %in% c("GO:0008150","GO:0000004","GO:0007582","GO:0005575","GO:0008372","GO:0003674","GO:0005554")),]
					
					if(is.vector(genes.matrix)){
						genes.matrix <- t(as.matrix(genes.matrix))
					}
					
					x.genes <- apply(genes.matrix,1,sum)
					y.genes <- apply(genes.matrix,2,sum)
					y.genes <- length(y.genes[y.genes != 0])
					
					wdth <- length(x.genes)

					right <- sort(x.genes*100/y.genes, decreasing=TRUE)

					if(length(right)>12){right <- right[1:12]
								wdth <- 12}
					if(compute.dim){
						max.right <- max(c(nchar(terms.name[names(right),2]),nchar("Transcriptional domain coverage (%)")))
						width <- 10.5*max.right
						if(width < 300){height <- 300}
					
						height <- 300 + 43*wdth
						#if(height < 450){height <- 450}
					}
								
					if(dev == "ps"){
						CairoPS(file=paste(extra,taxoname,".ps",sep=""), onefile = FALSE, width=width/100, height=height/100, horizontal=FALSE,
							paper="special", print.it=TRUE)
					}else if(dev == "pdf"){
						CairoPDF(file=paste(extra,taxoname,".pdf",sep=""), width=width/100, height=height/100,	paper="special")
					}else{
						CairoPNG(file=paste(extra,taxoname,".png",sep=""), height=0.85*height, width=0.85*width, res=600)
					}

					par(mfrow=c(1,1), bg="white", las=1, adj=0, mar =c(6,1,5,1), ask=F)

					barplot(right, col=sort(rainbow(n=length(right), start=3/6, end=4/6)), axisnames=FALSE, 
							beside=TRUE, horiz = TRUE, axes=FALSE, xlim=c(0,100), 
							space=1, width=0.8, ylim=c(0,(0.4+1.7*wdth)))
					
					par(adj=0.5)
					title(main = graph.name, xlab="Transcriptional domain coverage (%)", font.main = 2, 
						cex.main = 1.5, sub="", font.sub=2, cex.lab=1)

					for(j in 1:length(right)){	
						text(y=1.8*j + 0.2 -(0.2*(j-1)), x = 0, labels = terms.name[names(right)[j],2], 
							pos=4, cex=1, offset = 0.7 )
					}

					for(j in 1:length(right)){	
						if(right[j] <= 85){
							text(y=1*j + 0.2 +(0.6*(j-1)), x = right[j], labels = paste(round(right[j],digits=1),"%",sep=""), 
								pos=4, cex=0.8, offset = 0.7, col="red")
						}else{
							text(y=1*j + 0.2 +(0.6*(j-1)), x = right[j], labels = paste(round(right[j],digits=1),"%",sep=""), 
								pos=2, cex=0.8, offset = 0.7, col="red")
						}
					}

					axis(side=1, at=c(100,80,60,40, 20, 0),  labels = c(100,80,60,40, 20, 0), tick = TRUE, 
						line = NA,
						pos = NA, outer = FALSE, font = NA,
						lty = "solid", lwd = 1, col = NULL, hadj = NA, padj = NA)
					
					dev.off()

				}
			}



		}
	
	}else{
		print("No data available for ploting functional profiles\\!")
	}
	
	
}


#############################################################################################
#
# 18. Function .central.plot.two() -> Plot centralities of functional themes as back-to-back bar plots (up and down-regulated themes)
#
#############################################################################################


.central.plot.two <- function(up.annot=NULL, down.annot=NULL, clusters=NULL, taxoname=NULL, height=NA, width=NA, 
				compute.dim=TRUE, dev ="png", extra="", annot.clust.method="umilds"){
				
	graph.name <- NULL
	if(!is.null(taxoname)){
		if(taxoname == "bp"){graph.name <- "GO Biological Process"}
		if(taxoname == "cc"){graph.name <- "GO Cellular Component"}
		if(taxoname == "mf"){graph.name <- "GO Molecular Function"}
		if(taxoname == "kegg"){graph.name <- "KEGG"}

		if(taxoname == "GO Biological Process"){graph.name <- "GO Biological Process"}
		if(taxoname == "GO Cellular Component"){graph.name <- "GO Cellular Component"}
		if(taxoname == "GO Molecular Function"){graph.name <- "GO Molecular Function"}
		if(taxoname == "KEGG"){graph.name <- "KEGG"}
	}

	if(!is.null(down.annot) & !is.null(clusters) & !is.null(up.annot)){
		
		up.matrix <- as.data.frame(up.annot)
		down.matrix <- as.data.frame(down.annot)
		up.matrix <- up.matrix[!(rownames(up.matrix) %in% c("GO:0008150","GO:0000004","GO:0007582","GO:0005575","GO:0008372","GO:0003674","GO:0005554")),]
		down.matrix <- down.matrix[!(rownames(down.matrix) %in% c("GO:0008150","GO:0000004","GO:0007582","GO:0005575","GO:0008372","GO:0003674","GO:0005554")),]
		
		
		terms.name <- clusters$terms.name
		rownames(terms.name) <- terms.name[,1]

		adj.matrix <- clusters$annot.proximity

		rownames(adj.matrix) <- rownames(adj.matrix)
		colnames(adj.matrix) <- colnames(adj.matrix)

		centrality <- matrix(NA,nrow(adj.matrix),4)
		rownames(centrality) <- rownames(adj.matrix)
		colnames(centrality) <- c("degree","scaled_degree","betweenness","scaled_betweenness")

		if(annot.clust.method %in% c("umilds","spectral")){
			centrality[,1] <- degree(adj.matrix, gmode="graph", diag=FALSE, rescale=FALSE)
			centrality[,2] <- degree(adj.matrix, gmode="graph", diag=FALSE, rescale=TRUE)
			centrality[,3] <- betweenness(adj.matrix, gmode="graph", diag=FALSE, cmode="undirected", rescale=FALSE)
			centrality[,4] <- betweenness(adj.matrix, gmode="graph", diag=FALSE, cmode="undirected", rescale=TRUE)
		}else if(annot.clust.method=="ucknn"){
			centrality[,1] <- degree(adj.matrix, gmode="digraph", diag=FALSE, rescale=FALSE)
			centrality[,2] <- degree(adj.matrix, gmode="digraph", diag=FALSE, rescale=TRUE)
			centrality[,3] <- betweenness(adj.matrix, gmode="digraph", diag=FALSE, cmode="directed", rescale=FALSE)
			centrality[,4] <- betweenness(adj.matrix, gmode="digraph", diag=FALSE, cmode="directed", rescale=TRUE)
		}
		
		for(i in 1:nrow(centrality)){
			if(is.nan(centrality[i,4])){
				centrality[i,4] <- 0
			}
		}

		if(sum(centrality[,4])==0){
			centrality[,4] <- 100/nrow(centrality)
		}else{
			centrality[,4] <- 100*centrality[,4]
		}
		
		
		x.up <- centrality[rownames(up.matrix),4]
		names(x.up) <- rownames(up.matrix)
		x.down <- centrality[rownames(down.matrix),4]
		names(x.down) <- rownames(down.matrix)
		
		x <- c(x.up, x.down)

		if(length(x) > 25){
			x <- sort(x,decreasing=T)[1:25]
		}

		up.color <- NULL
		if(length(x.up[names(x.up) %in% names(x)]) > 0){
			x.up <- sort(x.up[names(x.up) %in% names(x)],decreasing=T)

			up.color <- rainbow(n=length(x.up), start=0, end=1/6)
			names(up.color) <- names(x.up)
		}

		down.color <- NULL
		if(length(x.down[names(x.down) %in% names(x)]) > 0){
			x.down <- sort(x.down[names(x.down) %in% names(x)],decreasing=T)

			down.color <- rainbow(n=length(x.down), start=2/6, end=3/6)
			names(down.color) <- names(x.down)
		}

		couleurs <- c(up.color,down.color)

		length.clusters <- NULL

		for(i in 1:length(clusters$gene.cluster)){
			length.clusters <- c(length.clusters, length(clusters$gene.cluster[[i]]))
		}

		names(length.clusters) <- 1:length(length.clusters)
		length.clusters <- sort(length.clusters,decreasing=T)

		# prepare panels

		bar.nr <- 0
		module.label <- NULL
		bar.clusters.left <- NULL
		bar.label.left <- NULL
		bar.clusters.right <- NULL
		bar.label.right <- NULL


		for(i in 1:length(length.clusters)){

			a <- clusters$id.cluster[[as.integer(names(length.clusters)[i])]]


			if(length(a[a %in% names(x)]) >0){

				a <- a[a %in% names(x)]
				a <- sort(x[a],decreasing=T)
				left <- sort(a[names(a) %in% names(x.up)], decreasing=T)
				right <- sort(a[names(a) %in% names(x.down)], decreasing=T)

				if(is.null(bar.clusters.left) && is.null(bar.clusters.right)){
					bar.clusters.left <- list(NULL)
					bar.clusters.right <- list(NULL)

					valeur.left <- left
					couleur.left <- couleurs[names(left)]
					bar.label.left <- terms.name[names(left),2]

					valeur.right <- right
					couleur.right <- couleurs[names(right)]
					bar.label.right <- terms.name[names(right),2]

					if(length(left) > length(right)){				
						valeur.right <- c(valeur.right,as.vector(matrix(NA,length(left)-length(right),1)))
						couleur.right <- c(couleur.right,as.character(matrix("#00000000",length(left)-length(right),1)))
						bar.label.right <- c(bar.label.right,as.character(matrix("",length(left)-length(right),1)))
					}else if(length(right) > length(left)){
						valeur.left <- c(valeur.left,as.vector(matrix(NA,length(right)-length(left),1)))
						couleur.left <- c(couleur.left,as.character(matrix("#00000000",length(right)-length(left),1)))
						bar.label.left <- c(bar.label.left,as.character(matrix("",length(right)-length(left),1)))
					}

					bar.clusters.left[[1]] <- list(valeur=valeur.left,couleur=couleur.left)
					bar.clusters.right[[1]] <- list(valeur=valeur.right,couleur=couleur.right)

					module.label <- c(as.character(matrix("",max(length(left),length(right)),1)),"Module 1")
				}else{
					valeur.left <- c(as.vector(matrix(NA,bar.nr,1)), left)
					couleur.left <- c(as.character(matrix("#00000000",bar.nr,1)), couleurs[names(left)])
					bar.label.left <- c(bar.label.left,"", terms.name[names(left),2])

					valeur.right <- c(as.vector(matrix(NA,bar.nr,1)), right)
					couleur.right <- c(as.character(matrix("#00000000",bar.nr,1)), couleurs[names(right)])
					bar.label.right <- c(bar.label.right, "", terms.name[names(right),2])


					if(length(left) > length(right)){				
						valeur.right <- c(valeur.right,as.vector(matrix(NA,length(left)-length(right),1)))
						couleur.right <- c(couleur.right,as.character(matrix("#00000000",length(left)-length(right),1)))
						bar.label.right <- c(bar.label.right,as.character(matrix("",length(left)-length(right),1)))
					}else if(length(right) > length(left)){
						valeur.left <- c(valeur.left,as.vector(matrix(NA,length(right)-length(left),1)))
						couleur.left <- c(couleur.left,as.character(matrix("#00000000",length(right)-length(left),1)))
						bar.label.left <- c(bar.label.left,as.character(matrix("",length(right)-length(left),1)))
					}

					bar.clusters.left[[i]] <- list(valeur=valeur.left,couleur=couleur.left)
					bar.clusters.right[[i]] <- list(valeur=valeur.right,couleur=couleur.right)			
					module.label <- c(module.label,as.character(matrix("",max(length(left),length(right)),1)),
								paste("Module ",i,sep=""))
				}

				bar.nr <- bar.nr + max(length(left),length(right)) + 1
			}
		}
		
		if(compute.dim){
			max.left_right <- max(c(nchar(bar.label.left),nchar(bar.label.right),nchar("Betweenness centrality (%)")))
			width <- 18.5*max.left_right
			if(width < 600){width <- 600}
		
			height <- 300 + 50*length(x)
			#if(height < 450){height <- 450}
			
		}

		if(dev == "ps"){
			CairoPS(file=paste(extra,taxoname," Centrality.ps",sep=""), onefile = FALSE, width=width/100, height=height/100, horizontal=FALSE,
				paper="special", print.it=TRUE)
		}else if(dev == "pdf"){
			CairoPDF(file=paste(extra,taxoname," Centrality.pdf",sep=""), width=width/100, height=height/100, paper="special")
		}else{
			CairoPNG(file=paste(extra,taxoname," Centrality.png",sep=""), height=height, width=width, res=600)
		}

		par(mfrow=c(1,2), bg="white", las=1, adj=1, mar =c(6,1,5,0.5), ask=F)

		# plot left panel

		if(!is.null(bar.clusters.left)){

			barplot(-bar.clusters.left[[1]]$valeur, col=bar.clusters.left[[1]]$couleur, axisnames=FALSE, beside=TRUE, horiz = TRUE, axes=FALSE, 
				xlim=c(-100,0), space=1, width=0.8, ylim=c(0,1.55*bar.nr))

			if(length(bar.clusters.left)>=2){

				for(i in 2:length(bar.clusters.left)){

					barplot(-bar.clusters.left[[i]]$valeur, col=bar.clusters.left[[i]]$couleur, axisnames=FALSE, beside=TRUE, 
						horiz = TRUE, axes=FALSE, xlim=c(-100,0), space=1, width=0.8, ylim=c(0,1.55*bar.nr), add=T)
				}
			}
		

			title(main = "Up-regulated Transcripts", font.main = 1, cex = 1, line=1)

			for(i in 1:length(bar.label.left)){	
				text(y=1.8*i + 0.2 -(0.2*(i-1)), x = 0, labels = bar.label.left[i], pos=2, cex=0.8, offset = 0.7 )
			}

			for(i in 1:length(bar.clusters.left)){
				for(j in 1:length(bar.clusters.left[[i]]$valeur)){
					if(!is.na(bar.clusters.left[[i]]$valeur[j]) && bar.clusters.left[[i]]$valeur[j] <= 85){
						text(y=1*j + 0.2 +(0.6*(j-1)), x = -(bar.clusters.left[[i]]$valeur[j]), labels = paste(round(bar.clusters.left[[i]]$valeur[j],
							digits=1),"%",sep=""), pos=2, cex=0.8, offset = 0.7, col="blue")
					}else{
						text(y=1*j + 0.2 +(0.6*(j-1)), x = -(bar.clusters.left[[i]]$valeur[j]), labels = paste(round(bar.clusters.left[[i]]$valeur[j],
							digits=1),"%",sep=""), pos=4, cex=0.8, offset = 0.7, col="blue")
					}
				}
			}

			axis(side=1, at=c(-100,-80,-60,-40, -20, 0), labels = c(100,80,60,40, 20, 0), tick = TRUE, line = NA,
			     pos = NA, outer = FALSE, font = NA,
			     lty = "solid", lwd = 1, col = NULL, hadj = NA, padj = NA)

		}

		# plot right panel

		par(mar=c(6,0.5,5,1), adj=0)

		if(!is.null(bar.clusters.right)){

			barplot(bar.clusters.right[[1]]$valeur, col=bar.clusters.right[[1]]$couleur, axisnames=FALSE, beside=TRUE, horiz = TRUE, 
				axes=FALSE, xlim=c(0,100), space=1, width=0.8, ylim=c(0,1.55*bar.nr))

			if(length(bar.clusters.right)>=2){

				for(i in 2:length(bar.clusters.right)){

					barplot(bar.clusters.right[[i]]$valeur, col=bar.clusters.right[[i]]$couleur, axisnames=FALSE, beside=TRUE, 
						horiz = TRUE, axes=FALSE, xlim=c(0,100), space=1, width=0.8, ylim=c(0,1.55*bar.nr), add=T)
				}
			}
		

			title(main = "Down-regulated Transcripts", font.main = 1, cex = 1, line=1)

			for(i in 1:length(bar.label.right)){	
				text(y=1.8*i + 0.2 -(0.2*(i-1)), x = 0, labels = bar.label.right[i], pos=4, cex=0.8, offset = 0.7 )
			}

			for(i in 1:length(bar.clusters.right)){
				for(j in 1:length(bar.clusters.right[[i]]$valeur)){
					if(!is.na(bar.clusters.right[[i]]$valeur[j]) && bar.clusters.right[[i]]$valeur[j] <= 85){
						text(y=1*j + 0.2 +(0.6*(j-1)), x = bar.clusters.right[[i]]$valeur[j], labels = paste(round(bar.clusters.right[[i]]$valeur[j],
							digits=1),"%",sep=""), pos=4, cex=0.8, offset = 0.7, col="blue")
					}else{
						text(y=1*j + 0.2 +(0.6*(j-1)), x = bar.clusters.right[[i]]$valeur[j], labels = paste(round(bar.clusters.right[[i]]$valeur[j],
							digits=1),"%",sep=""), pos=2, cex=0.8, offset = 0.7, col="blue")
					}
				}
			}
			

			axis(side=1, at=c(100,80,60,40, 20, 0), labels = c(100,80,60,40, 20, 0), tick = TRUE, line = NA,
			     pos = NA, outer = FALSE, font = NA,
			     lty = "solid", lwd = 1, col = NULL, hadj = NA, padj = NA)

		}

		# plot title and module labels


		par(mfrow=c(1,1), new=TRUE, mar=c(6,1,2,1), las=1, adj=0.5, xaxp=c(0,201,11))
		
		if(dev == "pdf"){
			abline(v=(50 + 1.2*900/width), lty="dotted",col="gray")
		}else{
			abline(v=(50 + 0.7*1100/width), lty="dotted",col="gray")
		}
		

		par(font=2)
		for(i in 1:length(module.label)){	
			text(y=1.8*(i-0.3) + 0.2 -(0.2*(i-1)), x = 50.7, labels = module.label[i], pos=NULL, cex=0.8, offset = 0)
		}
		par(font=1)

		title(main = graph.name, xlab="Betweenness centrality (%)", font.main = 2, 
			cex.main = 1.5, sub="", font.sub=2, cex.lab=1)
		dev.off()
	}else{
		print("No data available for ploting functional clusters\\!")
	}

}


#############################################################################################
#
# 18. Function .central.plot.one() -> Plot centralities of functional themes as bar plots
#
#############################################################################################


.central.plot.one <- function(genes.annot=NULL, clusters=NULL, taxoname=NULL, height=NA, width=NA, 
				compute.dim=TRUE, dev ="png", extra="", annot.clust.method="umilds"){
				
	graph.name <- NULL
	if(!is.null(taxoname)){
		if(taxoname == "bp"){graph.name <- "GO Biological Process"}
		if(taxoname == "cc"){graph.name <- "GO Cellular Component"}
		if(taxoname == "mf"){graph.name <- "GO Molecular Function"}
		if(taxoname == "kegg"){graph.name <- "KEGG"}

		if(taxoname == "GO Biological Process"){graph.name <- "GO Biological Process"}
		if(taxoname == "GO Cellular Component"){graph.name <- "GO Cellular Component"}
		if(taxoname == "GO Molecular Function"){graph.name <- "GO Molecular Function"}
		if(taxoname == "KEGG"){graph.name <- "KEGG"}
	}

	if(!is.null(genes.annot) & !is.null(clusters)){
		genes.matrix <- as.data.frame(genes.annot)
		genes.matrix <- genes.matrix[!(rownames(genes.matrix) %in% c("GO:0008150","GO:0000004","GO:0007582","GO:0005575","GO:0008372","GO:0003674","GO:0005554")),]
		
		
		
		terms.name <- clusters$terms.name
		rownames(terms.name) <- terms.name[,1]

		adj.matrix <- clusters$annot.proximity

		rownames(adj.matrix) <- rownames(adj.matrix)
		colnames(adj.matrix) <- colnames(adj.matrix)

		centrality <- matrix(NA,nrow(adj.matrix),4)
		rownames(centrality) <- rownames(adj.matrix)
		colnames(centrality) <- c("degree","scaled_degree","betweenness","scaled_betweenness")

		if(annot.clust.method %in% c("umilds","spectral")){
			centrality[,1] <- degree(adj.matrix, gmode="graph", diag=FALSE, rescale=FALSE)
			centrality[,2] <- degree(adj.matrix, gmode="graph", diag=FALSE, rescale=TRUE)
			centrality[,3] <- betweenness(adj.matrix, gmode="graph", diag=FALSE, cmode="undirected", rescale=FALSE)
			centrality[,4] <- betweenness(adj.matrix, gmode="graph", diag=FALSE, cmode="undirected", rescale=TRUE)
		}else if(annot.clust.method=="ucknn"){
			centrality[,1] <- degree(adj.matrix, gmode="digraph", diag=FALSE, rescale=FALSE)
			centrality[,2] <- degree(adj.matrix, gmode="digraph", diag=FALSE, rescale=TRUE)
			centrality[,3] <- betweenness(adj.matrix, gmode="digraph", diag=FALSE, cmode="directed", rescale=FALSE)
			centrality[,4] <- betweenness(adj.matrix, gmode="digraph", diag=FALSE, cmode="directed", rescale=TRUE)
		}
		
		for(i in 1:nrow(centrality)){
			if(is.nan(centrality[i,4])){
				centrality[i,4] <- 0
			}
		}

		if(sum(centrality[,4])==0){
			centrality[,4] <- 100/nrow(centrality)
		}else{
			centrality[,4] <- 100*centrality[,4]
		}
		
		
		x <- centrality[rownames(genes.matrix),4]

		if(length(x) > 25){
			x <- sort(x,decreasing=T)[1:25]
		}

		genes.color <- NULL
		x <- sort(x[names(x) %in% names(x)],decreasing=T)
		genes.color <- sort(rainbow(n=length(x), start=3/6, end=4/6))
		names(genes.color) <- names(x)
		

		couleurs <- genes.color

		length.clusters <- NULL

		for(i in 1:length(clusters$gene.cluster)){
			length.clusters <- c(length.clusters, length(clusters$gene.cluster[[i]]))
		}

		names(length.clusters) <- 1:length(length.clusters)
		length.clusters <- sort(length.clusters,decreasing=T)

		# prepare panels

		bar.nr <- 0
		module.label <- NULL
		bar.clusters.right <- NULL
		bar.label.right <- NULL


		for(i in 1:length(length.clusters)){

			a <- clusters$id.cluster[[as.integer(names(length.clusters)[i])]]


			if(length(a[a %in% names(x)]) >0){

				a <- a[a %in% names(x)]
				a <- sort(x[a],decreasing=T)
				right <- sort(a[names(a) %in% names(x)], decreasing=T)

				if(is.null(bar.clusters.right)){
					bar.clusters.right <- list(NULL)

					valeur.right <- right
					couleur.right <- couleurs[names(right)]
					bar.label.right <- terms.name[names(right),2]

					bar.clusters.right[[1]] <- list(valeur=valeur.right,couleur=couleur.right)
					module.label <- c(as.character(matrix("",length(right),1)),"Module 1")
				}else{
					valeur.right <- c(as.vector(matrix(NA,bar.nr,1)), right)
					couleur.right <- c(as.character(matrix("#00000000",bar.nr,1)), couleurs[names(right)])
					bar.label.right <- c(bar.label.right, "", terms.name[names(right),2])

					bar.clusters.right[[i]] <- list(valeur=valeur.right,couleur=couleur.right)			
					module.label <- c(module.label,as.character(matrix("",length(right),1)),
								paste("Module ",i,sep=""))
				}

				bar.nr <- bar.nr + length(right) + 1
			}
		}
		
		if(compute.dim){
			max.right <- max(c(nchar(bar.label.right),nchar("Betweenness centrality (%)")))
			width <- 9.5*max.right
			if(width < 300){height <- 300}
		
			height <- 300 + 50*length(x)
			if(height < 450){height <- 450}
			
		}

		if(dev == "ps"){
			CairoPS(file=paste(extra,taxoname," Centrality.ps",sep=""), onefile = FALSE, width=width/100, height=height/100, horizontal=FALSE,
				paper="special", print.it=TRUE)
		}else if(dev == "pdf"){
			CairoPDF(file=paste(extra,taxoname," Centrality.pdf",sep=""), width=width/100, height=height/100, paper="special")
		}else{
			CairoPNG(file=paste(extra,taxoname," Centrality.png",sep=""), height=height, width=width, res=600)
		}

		par(mfrow=c(1,1), bg="white", las=1, adj=0, mar =c(6,1,5,1), ask=F)

		# plot modules


		if(!is.null(bar.clusters.right)){

			barplot(bar.clusters.right[[1]]$valeur, col=bar.clusters.right[[1]]$couleur, axisnames=FALSE, beside=TRUE, horiz = TRUE, 
				axes=FALSE, xlim=c(0,100), space=1, width=0.8, ylim=c(0,1.55*bar.nr))

			if(length(bar.clusters.right)>=2){

				for(i in 2:length(bar.clusters.right)){

					barplot(bar.clusters.right[[i]]$valeur, col=bar.clusters.right[[i]]$couleur, axisnames=FALSE, beside=TRUE, 
						horiz = TRUE, axes=FALSE, xlim=c(0,100), space=1, width=0.8, ylim=c(0,1.55*bar.nr), add=T)
				}
			}
		
			for(i in 1:length(bar.label.right)){	
				text(y=1.8*i + 0.2 -(0.2*(i-1)), x = 0, labels = bar.label.right[i], pos=4, cex=0.8, offset = 0.7 )
			}

			for(i in 1:length(bar.clusters.right)){
				for(j in 1:length(bar.clusters.right[[i]]$valeur)){
					if(!is.na(bar.clusters.right[[i]]$valeur[j]) && bar.clusters.right[[i]]$valeur[j] <= 85){
						text(y=1*j + 0.2 +(0.6*(j-1)), x = bar.clusters.right[[i]]$valeur[j], labels = paste(round(bar.clusters.right[[i]]$valeur[j],
							digits=1),"%",sep=""), pos=4, cex=0.8, offset = 0.7, col="blue")
					}else{
						text(y=1*j + 0.2 +(0.6*(j-1)), x = bar.clusters.right[[i]]$valeur[j], labels = paste(round(bar.clusters.right[[i]]$valeur[j],
							digits=1),"%",sep=""), pos=2, cex=0.8, offset = 0.7, col="blue")
					}
				}
			}

			axis(side=1, at=c(100,80,60,40, 20, 0), labels = c(100,80,60,40, 20, 0), tick = TRUE, line = NA,
			     pos = NA, outer = FALSE, font = NA,
			     lty = "solid", lwd = 1, col = NULL, hadj = NA, padj = NA)

		}

		# plot title and module labels


		par(mfrow=c(1,1), new=TRUE, mar=c(6,1,2,1), las=1, adj=0.5, xaxp=c(0,101,11))
		
		par(font=2)
		for(i in 1:length(module.label)){	
			text(y=1.8*(i-0.3) + 0.2 -(0.2*(i-1)), x = 50.7, labels = module.label[i], pos=NULL, cex=0.8, offset = 0)
		}
		par(font=1)

		title(main = graph.name, xlab="Betweenness centrality (%)", font.main = 2, 
			cex.main = 1.5, sub="", font.sub=2, cex.lab=1)
		dev.off()
	}else{
		print("No data available for ploting functional clusters\\!")
	}

}


#############################################################################################
#
# 18. Function .clusters.plot() -> Experimental: plots functional clusters as back-to-back bar plots (up and down-regulated themes)
#
#############################################################################################

# experimental!!!

.clusters.plot <- function(up.annot=NULL, down.annot=NULL, clusters=NULL, taxoname=NULL, height=700, width=900, 
				dev ="jpeg", extra=""){

	if(!is.null(down.annot) & !is.null(clusters) & !is.null(up.annot)){
		up.matrix <- up.annot$annot.matrix
		down.matrix <- down.annot$annot.matrix

		x.up <- apply(up.matrix,1,sum)
		y.up <- apply(up.matrix,2,sum)
		y.up <- length(y.up[y.up != 0])
		x.down <- apply(down.matrix,1,sum)
		y.down <- apply(down.matrix,2,sum)
		y.down <- length(y.down[y.down != 0])
		x <- c(x.up * 100 / y.up, x.down * 100 / y.down)

		if(length(x) > 25){
			x <- sort(x,decreasing=T)[1:25]
		}

		up.color <- NULL
		if(length(x.up[names(x.up) %in% names(x)]) > 0){
			x.up <- sort(x.up[names(x.up) %in% names(x)],decreasing=T)

			up.color <- rainbow(n=length(x.up), start=0, end=1/6)
			names(up.color) <- names(x.up)
		}

		down.color <- NULL
		if(length(x.down[names(x.down) %in% names(x)]) > 0){
			x.down <- sort(x.down[names(x.down) %in% names(x)],decreasing=T)

			down.color <- rainbow(n=length(x.down), start=3/6, end=4/6)
			names(down.color) <- names(x.down)
		}

		couleurs <- c(up.color,down.color)

		length.clusters <- NULL

		for(i in 1:length(clusters$gene.cluster)){
			length.clusters <- c(length.clusters, length(clusters$gene.cluster[[i]]))
		}

		names(length.clusters) <- 1:length(length.clusters)
		length.clusters <- sort(length.clusters,decreasing=T)

		# prepare panels

		bar.nr <- 0
		module.label <- NULL
		bar.clusters.left <- NULL
		bar.label.left <- NULL
		bar.clusters.right <- NULL
		bar.label.right <- NULL


		for(i in 1:length(length.clusters)){

			a <- clusters$id.cluster[[as.integer(names(length.clusters)[i])]]


			if(length(a[a %in% names(x)]) >0){

				a <- a[a %in% names(x)]
				a <- sort(x[a],decreasing=T)
				left <- sort(a[names(a) %in% names(x.up)], decreasing=T)
				right <- sort(a[names(a) %in% names(x.down)], decreasing=T)

				if(is.null(bar.clusters.left) && is.null(bar.clusters.right)){
					bar.clusters.left <- list(NULL)
					bar.clusters.right <- list(NULL)

					valeur.left <- left
					couleur.left <- couleurs[names(left)]
					bar.label.left <- clusters$terms.name[names(left),2]

					valeur.right <- right
					couleur.right <- couleurs[names(right)]
					bar.label.right <- clusters$terms.name[names(right),2]

					if(length(left) > length(right)){				
						valeur.right <- c(valeur.right,as.vector(matrix(NA,length(left)-length(right),1)))
						couleur.right <- c(couleur.right,as.character(matrix("#00000000",length(left)-length(right),1)))
						bar.label.right <- c(bar.label.right,as.character(matrix("",length(left)-length(right),1)))
					}else if(length(right) > length(left)){
						valeur.left <- c(valeur.left,as.vector(matrix(NA,length(right)-length(left),1)))
						couleur.left <- c(couleur.left,as.character(matrix("#00000000",length(right)-length(left),1)))
						bar.label.left <- c(bar.label.left,as.character(matrix("",length(right)-length(left),1)))
					}

					bar.clusters.left[[1]] <- list(valeur=valeur.left,couleur=couleur.left)
					bar.clusters.right[[1]] <- list(valeur=valeur.right,couleur=couleur.right)

					module.label <- c(as.character(matrix("",max(length(left),length(right)),1)),"Module 1")
				}else{
					valeur.left <- c(as.vector(matrix(NA,bar.nr,1)), left)
					couleur.left <- c(as.character(matrix("#00000000",bar.nr,1)), couleurs[names(left)])
					bar.label.left <- c(bar.label.left,"", clusters$terms.name[names(left),2])

					valeur.right <- c(as.vector(matrix(NA,bar.nr,1)), right)
					couleur.right <- c(as.character(matrix("#00000000",bar.nr,1)), couleurs[names(right)])
					bar.label.right <- c(bar.label.right, "", clusters$terms.name[names(right),2])


					if(length(left) > length(right)){				
						valeur.right <- c(valeur.right,as.vector(matrix(NA,length(left)-length(right),1)))
						couleur.right <- c(couleur.right,as.character(matrix("#00000000",length(left)-length(right),1)))
						bar.label.right <- c(bar.label.right,as.character(matrix("",length(left)-length(right),1)))
					}else if(length(right) > length(left)){
						valeur.left <- c(valeur.left,as.vector(matrix(NA,length(right)-length(left),1)))
						couleur.left <- c(couleur.left,as.character(matrix("#00000000",length(right)-length(left),1)))
						bar.label.left <- c(bar.label.left,as.character(matrix("",length(right)-length(left),1)))
					}

					bar.clusters.left[[i]] <- list(valeur=valeur.left,couleur=couleur.left)
					bar.clusters.right[[i]] <- list(valeur=valeur.right,couleur=couleur.right)			
					module.label <- c(module.label,as.character(matrix("",max(length(left),length(right)),1)),
								paste("Module ",i,sep=""))
				}

				bar.nr <- bar.nr + max(length(left),length(right)) + 1
			}
		}



		par(mfrow=c(1,2), bg="white", las=1, adj=1, mar =c(6,1,5,0.5), ask=F)

		# plot left panel

		if(!is.null(bar.clusters.left)){

			barplot(-bar.clusters.left[[1]]$valeur, col=bar.clusters.left[[1]]$couleur, axisnames=FALSE, beside=TRUE, horiz = TRUE, axes=FALSE, 
				xlim=c(-100,0), space=1, width=0.8, ylim=c(0,1.55*bar.nr))

			if(length(bar.clusters.left)>=2){

				for(i in 2:length(bar.clusters.left)){

					barplot(-bar.clusters.left[[i]]$valeur, col=bar.clusters.left[[i]]$couleur, axisnames=FALSE, beside=TRUE, 
						horiz = TRUE, axes=FALSE, xlim=c(-100,0), space=1, width=0.8, ylim=c(0,1.55*bar.nr), add=T)
				}
			}
		}

		title(main = "Upregulated Transcripts", font.main = 1, cex = 1, line=1)

		for(i in 1:length(bar.label.left)){	
			text(y=1.8*i + 0.2 -(0.2*(i-1)), x = 0, labels = bar.label.left[i], pos=2, cex=0.8, offset = 0.7 )
		}

		axis(side=1, at=c(-100,-80,-60,-40, -20, 0), labels = c(100,80,60,40, 20, 0), tick = TRUE, line = NA,
		     pos = NA, outer = FALSE, font = NA,
		     lty = "solid", lwd = 1, col = NULL, hadj = NA, padj = NA)



		# plot right panel

		par(mar=c(6,0.5,5,1), adj=0)

		if(!is.null(bar.clusters.right)){

			barplot(bar.clusters.right[[1]]$valeur, col=bar.clusters.right[[1]]$couleur, axisnames=FALSE, beside=TRUE, horiz = TRUE, 
				axes=FALSE, xlim=c(0,100), space=1, width=0.8, ylim=c(0,1.55*bar.nr))

			if(length(bar.clusters.right)>=2){

				for(i in 2:length(bar.clusters.right)){

					barplot(bar.clusters.right[[i]]$valeur, col=bar.clusters.right[[i]]$couleur, axisnames=FALSE, beside=TRUE, 
						horiz = TRUE, axes=FALSE, xlim=c(0,100), space=1, width=0.8, ylim=c(0,1.55*bar.nr), add=T)
				}
			}
		}

		title(main = "Downregulated Transcripts", font.main = 1, cex = 1, line=1)

		for(i in 1:length(bar.label.right)){	
			text(y=1.8*i + 0.2 -(0.2*(i-1)), x = 0, labels = bar.label.right[i], pos=4, cex=0.8, offset = 0.7 )
		}

		axis(side=1, at=c(100,80,60,40, 20, 0), labels = c(100,80,60,40, 20, 0), tick = TRUE, line = NA,
		     pos = NA, outer = FALSE, font = NA,
		     lty = "solid", lwd = 1, col = NULL, hadj = NA, padj = NA)



		# plot title and module labels


		par(mfrow=c(1,1), new=TRUE, mar=c(6,1,2,1), las=1, adj=0.5, xaxp=c(0,201,11))

		abline(v=51, lty="dotted",col="gray")

		par(font=2)
		for(i in 1:length(module.label)){	
			text(y=1.8*(i-0.3) + 0.2 -(0.2*(i-1)), x = 50.7, labels = module.label[i], pos=NULL, cex=0.8, offset = 0)
		}
		par(font=1)

		title(main = "KEGG interaction modules", xlab="Transcript space coverage (%)", font.main = 2, 
			cex = 2, sub="", font.sub=2)
	}else{
		print("No data available for ploting functional clusters\\!")
	}

}



#############################################################################################
#
# 18. Function .pvalues() -> Calculate the p-value for each annotating term using an exact Fisher test
#
#############################################################################################

#Array DATAS contains all the data necessary to build the contingency table for calculating the Fisher test

.pvalues <- function(datas){
	a <- as.matrix(datas[,1])
	b <- as.matrix(datas[,2]) - as.matrix(datas[,1])
	cc <- as.matrix(datas[,3]) - as.matrix(datas[,1])
	d <- as.matrix(datas[,4]) + as.matrix(datas[,1]) - as.matrix(datas[,2]) - as.matrix(datas[,3])
	p <- matrix(0,length(t(a)),1)
	for(i in 1:length(t(a))){
		p[i] <- as.double(matrix(fisher.test(matrix(c(a[i,1],b[i,1],cc[i,1],d[i,1]),2,2),alternative = "g"),1,1))
	}	

# test for p-values > 1 (to avoid errors during Storey q-values calculation)	

	for(i in 1:length(p)){
		if(p[i] > 1){
			p[i] <- 1
		}
	}

	return(p)
	rm()
}





#############################################################################################
#
# 19. Function .fdr.adjust() -> P-values adjustment using FDR of Benjamini & Hochberg (1995) modified by Paciorek (2004)
#
#############################################################################################



.fdr.adjust <- function(pvalues, qlevel=0.05, method="original", adjust=NULL){	


	x.fdr <- fdr(pvals=pvalues,qlevel=qlevel,method=method,adjustment.method=adjust)	 
	y.fdr <- as.vector(matrix(1,length(pvalues),1))

	if(!is.null(x.fdr)){
		for(i in 1:length(x.fdr)){ y.fdr[x.fdr[i]] <- pvalues[x.fdr[i]]}
	}

	return(y.fdr)

}






# ECML-2006 paper routines used to calculate mutual information coefficient (MI) between 2 genes and performs permutations of gene expression data in order 
#   to estimate MI coefficient distribution, and also to compute distance matrix between annotations



#############################################################################################
#
# 20. Function .sil.h() -> Calculate the silhouette for a partition of ucknn clusters (asymetrical distances between objects)
#
#############################################################################################

.sil.h <- function(memb, annot.distance){
	
	clust.list <- list(NULL)

	for(i in 1:length(unique(memb))){
		clust.vect <- memb[memb == i]
		clust.list[[i]] <- names(clust.vect)
	}
	
	sil.list <- NULL
	for(i in 1:length(clust.list)){
		sil.clust <- NULL
		for(j in 1:length(clust.list[[i]])){
			if(length(clust.list[[i]]) == 1){
				sil.clust <- c(sil.clust,0)
			}else{
				x <- mean(annot.distance[clust.list[[i]], clust.list[[i]][j]])
			
				y <- NULL
				for(l in 1:length(clust.list)){
					if(l != i){
						y <- c(y,mean(annot.distance[clust.list[[l]], clust.list[[i]][j]]))
					}
				}
				y <- min(y,na.rm=TRUE)
				sil.clust <- c(sil.clust,(y - x)/max(x,y,na.rm=TRUE))
			
				rm(x,y)
			}
		}
		sil.list <- c(sil.list,mean(sil.clust))
		
	}
	
	return(sil.list)
	


}



#############################################################################################
#
# 21. Function .annot.citers.rank() -> Calculate the citers rank vector starting from an ensemble of functional annotations (vectors of gene identifiers)
#
#############################################################################################

.annot.citers.rank <- function(annot.distance, k){

	# starting from the matrix of unilateral distances between annotations compute ranks from distances by columns
	
	annot.rank.matrix <- matrix(NA,nrow(annot.distance),ncol(annot.distance))
	colnames(annot.rank.matrix) <- colnames(annot.distance)
	rownames(annot.rank.matrix) <- rownames(annot.distance)
	
	for(i in 1:ncol(annot.distance)){
		
		annot.rank.matrix[,i] <- rank(annot.distance[,i],ties.method="min")
	
	}
	
	
	citers.rank <- matrix(NA,nrow(annot.rank.matrix),1)
	
	for(i in 1:nrow(annot.rank.matrix)){
		citers.rank[i] <- sum(annot.rank.matrix[i,]  <= k) # results in the number of elements not in the sum of their values
	}
	
	citers.rank <- as.vector(citers.rank)
	names(citers.rank) <- rownames(annot.rank.matrix)
	if(length(citers.rank[citers.rank != 1]) > 0){
		citers.rank <- citers.rank[citers.rank != 1] # get rid of those who are refering none
	
		citers.rank <- rank(citers.rank,ties.method="min")
		citers.rank <- rank(abs(citers.rank - max(citers.rank,na.rm=TRUE)),ties.method="min")

		citers.rank <- sort(citers.rank)

	}else{
		citers.rank <- as.vector(NULL)
	}
	return(citers.rank)
}



#############################################################################################
#
# 22. Function .presumptive.proximity() -> Transforms a symetrical distance matrix into a proximity matrix using a gaussian kernel
#
#############################################################################################

.presumptive.proximity <- function(annot.distance, sigma=NA){

	if(is.na(sigma)){
		sigma <- median(annot.distance[lower.tri(annot.distance,diag=FALSE)])
	}
				
	annot.proximity <- exp(-(annot.distance^2)/(2*(sigma^2)))
	diag(annot.proximity) <- 0
	return(annot.proximity)
}



#############################################################################################
#
# 22. Function .annotation.distance() -> Calculate the distance matrix between an ensemble of functional annotations (vectors of gene identifiers)
#
#############################################################################################

.annotation.distance <- function(annot.matrix, coexp.matrix, measure = "unilat.pond.norm.mean"){
	
	# annotations on rows and genes on columns
	
	annot.distance <- matrix(NA, nrow(annot.matrix), nrow(annot.matrix))
	annot.names <- rownames(annot.matrix)
	colnames(annot.distance) <- annot.names
  	rownames(annot.distance) <- annot.names

	# for each annotation in the list 
	  
  	for(i in 1:nrow(annot.matrix)){
  	
  	# for each annotation in the list except the actual one
  	
  		for(j in 1:nrow(annot.matrix)){
  			if(j != i){  				
  				annot.distance[j,i] <- .presumptive.distance(annot1 = names(annot.matrix[i,][annot.matrix[i,] == 1]), 
  						annot2 = names(annot.matrix[j,][annot.matrix[j,] == 1]), coexp.matrix = coexp.matrix, 
  						measure = measure)
  			}  	
  		}  	
  	}
  	
  	diag(annot.distance) <- 0
  	
  	return(annot.distance)

}



#############################################################################################
#
# 24. Function .presumptive.distance() -> Calculate distance between two annotations (vectors of gene identifiers)
#
#############################################################################################

.presumptive.distance <- function(annot1, annot2, coexp.matrix, measure = "norm.sum"){

	distance <- 1	# corresponding to the maximum distance 1 the minimum being 0

	# total number of strongly co-regulated genes (MI threshold) divised by the number of genes annoted by the two annotations

	if(measure == "norm.sum"){	# does not take into account the "real" distances between genes
 		if(length(annot1) >0 & length(annot2) >0){
			strong1 <- NULL  # vectors of strongly co-regulated genes in the two annotations
			strong2 <- NULL
			
			# for each gene in annot1 search the closest gene in annot2
			try(afin.matrix <- as.matrix(coexp.matrix[annot1, annot2]))
			#print(paste(length(annot1),":",length(annot2)))
			#print(paste(nrow(afin.matrix),":",ncol(afin.matrix)))


			try(strong1 <- apply(afin.matrix,1,sum))
			try(strong1 <- names(strong1[strong1 > 0]))

			try(strong2 <- apply(afin.matrix,2,sum))
			try(strong2 <- names(strong2[strong2 > 0]))

			strong <- unique(c(strong1,strong2))

			# caculate a normalzed distance as the sum of "strong" genes rapported as the sum of total genes (annot1+2) then normalize to unit
			if(length(strong) > 0){
				try(distance <- (1 - (length(strong)/length(unique(c(annot1,annot2))))))
				#print(paste("distance:",distance))
			}
		}
		
	# mean of min distances between strongly related genes within two annotations normalized to the maximum of this mean and ponderated by the previous measure accounting only gene numbers
	
	}else if(measure == "pond.norm.mean"){	# take into account the "real" distances between genes ponderated by the previous measure (norm.sum)
		
		strong <- NULL  # vector of strongly co-regulated genes in the two annotations
		gene.dist <- NULL	# vector of min distances among genes which will serve to calculate mean(min) distances among genes strongly correlated 		
		  
		# for each gene in annot1 search the closest gene in annot2
		  
		for(i in 1:length(annot1)){
		   	d.vector <- c(coexp.matrix[annot2, annot1[i]], coexp.matrix[annot1[i], annot2])
		    
		    	# if there is at least one gene in annot2 which correlates stronger than the threshold then count the gene
		    
		    	if(max(as.numeric(d.vector), na.rm=TRUE) > 0){
		      		strong <- c(strong, annot1[i])
		      		gene.dist <- c(gene.dist, max(as.numeric(d.vector), na.rm=TRUE))
		    	}  
		}
		  
		# for each gene in annot2 search the closest gene in annot1
		  
		for(i in 1:length(annot2)){
		    	d.vector <- c(coexp.matrix[annot1, annot2[i]], coexp.matrix[annot2[i], annot1])
		    
		    	if(max(as.numeric(d.vector), na.rm=TRUE) > 0){
		      		strong <- c(strong, annot2[i])
		      		gene.dist <- c(gene.dist, max(as.numeric(d.vector), na.rm=TRUE))
		   	}  
		}
				
		# caculate a normalzed distance as the sum of "strong" genes rapported as the sum of total genes (annot1+2) then normalize to unit
		if(length(strong) > 0){
		  	distance <- (1 - ((length(strong)/(length(annot1) + length(annot2))) * (mean(gene.dist) / 2 )))
		}	
	
	# equivalent of "norm.sum" for the unilateral case
	
	}else if(measure == "unilat.norm.sum"){
		
		strong <- NULL  # vector of strongly co-regulated genes in the two annotations
		  		
		  
		# for each gene in annot1 search the closest gene in annot2
		  
		for(i in 1:length(annot1)){
		   	d.vector <- c(coexp.matrix[annot2, annot1[i]], coexp.matrix[annot1[i], annot2])
		    
		    	# if there is at least one gene in annot2 which correlates stronger than the threshold then count the gene
		    
		    	if(max(as.numeric(d.vector), na.rm=TRUE) > 0){
		      		strong <- c(strong, annot1[i])
		    	}  
  		}
		
		# caculate a normalized distance as the sum of "strong" genes related to annot2 rapported to the sum of total genes (annot1) in order to normalize to unit
		if(length(strong) > 0){
		  	distance <- (1 - (length(strong)/length(annot1)))
		}
	
	# equivalent of "pond.norm.mean" for the unilateral case
	
	}else if(measure == "unilat.pond.norm.mean"){
		
		strong <- NULL  # vector of strongly co-regulated genes in the two annotations
				gene.dist <- NULL	# vector of min distances among genes which will serve to calculate mean(min) distances among genes strongly correlated 		
				  
				# for each gene in annot1 search the closest gene in annot2
				  
				for(i in 1:length(annot1)){
				   	d.vector <- c(coexp.matrix[annot2, annot1[i]], coexp.matrix[annot1[i], annot2])
				    
				    	# if there is at least one gene in annot2 which correlates stronger than the threshold then count the gene
				    
				    	if(max(as.numeric(d.vector), na.rm=TRUE) > 0){
				      		strong <- c(strong, annot1[i])
				      		gene.dist <- c(gene.dist, max(as.numeric(d.vector), na.rm=TRUE))
				    	}  
		}
		
		# caculate a normalzed distance as the sum of "strong" genes rapported as the sum of total genes (annot1+2) then normalize to unit
				if(length(strong) > 0){
				  	distance <- (1 - (((length(strong)/length(annot1)) * (mean(gene.dist) / 2))))
				}	
	
		
	}

	return(distance)

}




#############################################################################################
#
# 25. Function .uc.knn() -> Aglomerative kNN clustering of genomic annotations (ECML 2006)
#
#############################################################################################

.uc.knn <- function(annot.matrix, coexp.matrix, alpha = 0.05, taxoname, annot.prox.measure = "unilat.pond.norm.mean",splits = 10, normal = TRUE, terms.name){
  	
  	cat(paste("\n\t",taxoname," annotations UC-KNN clustering started...",sep=""))
 # getting significantly enriched annotations data
  	
	rownames(terms.name) <- terms.name[,1]
	   	 
 # calculate the annotation matrix of citers
    
  	annot.distance <- .annotation.distance(annot.matrix = annot.matrix, coexp.matrix = coexp.matrix, measure = annot.prox.measure)

    	class.same <- 0 # count the cases in which kNN give same result as simple
	class.differ <- 0 # count the cases in kNN result differs from simple
	part.list <- list(NULL) # the list of cluster partitions
	k.list <- list(NULL) # the list with k of KNN

  for(k in 1:(ncol(annot.distance)-1)){
  # find the best citers for kNN in ranked order
  	best.citers <- .annot.citers.rank(annot.distance, k = k) # k is the kNN
  	
  # cluster annotations    
		
	
	if(length(best.citers) == 1){
		x <- as.vector(matrix(1,nrow(annot.distance),1))
		names(x) <- rownames(annot.distance)
		part.list[[length(part.list)+1]] <- x
		k.list[[length(k.list)+1]] <- k
		rm(x)
	}else if(length(best.citers) >= 2){		
		if(length(best.citers) > 2){
			lim <- (length(best.citers)-1)
		}else{
			lim <- 2
		}
		for(i in 2:lim){
			clust.part <- matrix(0, i, ncol(annot.distance)-i) # the future cluster matrix
			
			rownames(clust.part) <- names(best.citers[1:i])
			
			colnames(clust.part) <- colnames(annot.distance)[!colnames(annot.distance) %in% rownames(clust.part)]
			
	
			for(j in 1:ncol(clust.part)){
				test <- as.vector(annot.distance[,colnames(clust.part)[j]])
				names(test) <- rownames(annot.distance)
				test <- rank(test,ties.method="min")
				test <- test[test <= k]
				
				test.matrix <- annot.distance[rownames(clust.part),c(names(test),colnames(clust.part)[j])]
				
				test.best <- matrix(NA,nrow(test.matrix),1)
				test.best <- as.vector(test.best)
				names(test.best) <- rownames(test.matrix)
	
				for(l in 1:nrow(test.matrix)){
					test.best[l] <- sum(test.matrix[l,]) # find the best seed to cluster with
				}
				
				best.seed <- names(test.best[test.best == min(test.best, na.rm=TRUE)])
				clust.part[best.seed,j] <- 1 # cluster to best seed
				
				# count if simple differs from kNN
				simple.matrix <- as.matrix(annot.distance[rownames(clust.part),colnames(clust.part)[j]])
				
				if(simple.matrix[best.seed,1] != min(simple.matrix[,1],na.rm=TRUE)){
					class.differ <- class.differ + 1
				}else{
					class.same <- class.same + 1
				}
			}
			# add seeds to their clusters
			for(j in 1:nrow(clust.part)){
				x <- matrix(0,nrow(clust.part),1)
				rownames(x) <- rownames(clust.part)
				colnames(x) <- rownames(clust.part)[j]
				x[j,1] <- 1
				clust.part <- cbind(x,clust.part)
				
			}
			# build a "cutree" like vector of clusters for silhouette calculation
			
			part.vector <- as.vector(NULL)
			
			for(j in 1:nrow(clust.part)){
				x <- clust.part[j,]
				x <- x[x == 1]
				x <- x +j -1
				part.vector <- c(part.vector, x)
			}
			
			part.list[[length(part.list)+1]] <- part.vector
			k.list[[length(k.list)+1]] <- k
		}
	}
 }
  	
  	
  	.sil.hc <- as.vector(NULL)	# vector of silhouettes

  # calculate cluster silhouettes
	
  	for(i in 2:length(part.list)){ 
    	# cut the tree 
    		memb <- part.list[[i]]
    		sil <- mean(.sil.h(memb, annot.distance))
    		.sil.hc <- c(.sil.hc, sil)

  	}
	

	names(.sil.hc) <- 2:length(part.list)	# choose the max silhouette
	#cat(.sil.hc)
	# find the best partition 
	.sil.hc <- .sil.hc[!is.nan(.sil.hc)]
	best.index <- as.numeric(names(.sil.hc[.sil.hc == max(.sil.hc, na.rm=TRUE)]))
	cat(paste("\n\t\t\tBest index values = ",length(best.index),sep=""))
	
	optim.clust <- length(unique(part.list[[best.index[1]]]))
	optim.k <- k.list[[best.index[1]]]
	
	for(i in 1:length(best.index)){
		optim.clust <- length(unique(part.list[[best.index[i]]]))
		optim.k <- k.list[[best.index[i]]]
		cat(paste("\n\t\tOptimum clusters number = ",optim.clust,sep=""))
		cat(paste("\n\t\tOptimum kNN = ",optim.k,sep=""))
	}
	
	best.partition <- part.list[[best.index[1]]]
	cluster.length <- as.vector(NULL)
	best.index <- best.index[1]
	sil.part <- mean(.sil.h(best.partition, annot.distance))
	sil.cluster <- .sil.h(best.partition, annot.distance)
	
	id.cluster <- list(NULL)	
	term.cluster <- list(NULL)
	gene.cluster <- list(NULL)
	
	for(i in 1:(optim.clust)){
		
		id.cluster[[i]] <- as.character(names(best.partition[best.partition == i]))
		cluster.length <- c(cluster.length,length(id.cluster[[i]]))
	
	}
	
	if(length(id.cluster) > 0){
		
		for(i in 1:length(id.cluster)){
			term.cluster[[i]] <- terms.name[id.cluster[[i]],2]
			
			gene.vect.matrix <- annot.matrix[id.cluster[[i]],]
			gene.vect <- NULL
			
			if(is.matrix(gene.vect.matrix)){
			
				for(j in 1:nrow(gene.vect.matrix)){
					gene.vect <- c(gene.vect, names(gene.vect.matrix[j,][gene.vect.matrix[j,] == 1]))			
				}
			}else{
				gene.vect <- c(gene.vect, names(gene.vect.matrix[gene.vect.matrix == 1]))	
			}
			
			gene.cluster[[i]] <- unique(gene.vect)
		}	
		
		
		gene.connect <- matrix(0,nrow(coexp.matrix),length(gene.cluster)+4)
		rownames(gene.connect) <- rownames(coexp.matrix)
		colnames(gene.connect) <- c("annotation_module","assigned_module",as.character(1:length(gene.cluster)),"total_intra","total_net")
			
		for(i in 1:length(gene.cluster)){
			gene.connect[as.character(gene.cluster[[i]]),"annotation_module"] <- i
			
		}
		
		gene.connect[,"assigned_module"] <- gene.connect[,"annotation_module"]
		
		for(i in 1:nrow(gene.connect)){
				
			for(j in 1:length(gene.cluster)){
				x <- coexp.matrix[rownames(gene.connect)[i],as.character(gene.cluster[[j]])]
				gene.connect[i,as.character(j)] <- length(x[x > 0])
					
			}
			
			gene.connect[i,"total_intra"] <- sum(gene.connect[i,as.character(1:length(gene.cluster))])
			
			y <- coexp.matrix[rownames(gene.connect)[i],]
			gene.connect[i,"total_net"] <- length(y[y > 0])
		}
		
		for(i in 1:nrow(gene.connect)){
			if(gene.connect[i,"assigned_module"] == 0){
				x <- gene.connect[i,as.character(1:length(gene.cluster))]
				names(x) <- 1:length(x)
				gene.connect[i,"assigned_module"] <- as.numeric(names(x)[x == max(x)][1])
			}
		}
		
		for(i in 1:length(gene.cluster)){
			lngth <- sqrt(sum(gene.connect[,as.character(i)] * gene.connect[,as.character(i)]))
			gene.connect[,as.character(i)] <- gene.connect[,as.character(i)]/lngth
			rm(lngth)
		}
		
		# total_intra
		lngth <- sqrt(sum(gene.connect[,"total_intra"] * gene.connect[,"total_intra"]))
		gene.connect[,"total_intra"] <- gene.connect[,"total_intra"]/lngth
		rm(lngth)
		
		# total_net
		lngth <- sqrt(sum(gene.connect[,"total_net"] * gene.connect[,"total_net"]))
		gene.connect[,"total_net"] <- gene.connect[,"total_net"]/lngth
		rm(lngth)

		names(best.partition) <- terms.name[names(best.partition),2]
		annot.proximity <- 1 - annot.distance
		
		clusters <- list(annot.distance=annot.distance, annot.proximity=annot.proximity, sil=.sil.hc, sil.part=sil.part,sil.cluster=sil.cluster,
					cluster.length=cluster.length,best.index=best.index,id.cluster=id.cluster,term.cluster=term.cluster,
					terms.name=terms.name, gene.cluster=gene.cluster,best.partition=best.partition,gene.connect=gene.connect)
		  		
	}else{
		clusters <- NULL
		cat("\n\t\tNo clusters were created...")
	}
	#cat(paste("\n\t\t\tSame = ",class.same,sep=""))
	#cat(paste("\n\t\t\tDiffer = ",class.differ,sep=""))
	#cat(paste("\n\t\t\tDiffer = ",class.differ*100/class.same,sep=""))
	return(clusters)

}




#############################################################################################
#
# 26. Function fdr() -> FDR correction routine available from Chris Paciorek website
#
#############################################################################################


# File:      fdr.r 
# Date:      5/10/04
# Version:   0.1.3  
# Author:    Chris Paciorek - please contact the author with bug
#            reports: paciorek AT alumni.cmu.edu
# License:   GPL version 2 or later
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# Purpose:   implement False Discovery Rate (FDR) functions for multiple testing, following the Ventura et al. reference below
# Usage:     source('fdr.r');  fdr(my.pvals)
# References:
#             Ventura, V., C.J. Paciorek, and J.S. Risbey.  2004.  Controlling the proportion of falsely-rejected hypotheses when conducting multiple tests with climatological data.  Journal of Climate, in press.  Also Carnegie Mellon University, Department of Statistics technical report 775 (www.stat.cmu.edu/tr/tr775/tr775.html).
#             Benjamini, Y, and Y. Hochberg. 1995. Controlling the false discovery rate: a practical and powerful approach to multiple testing.  JRSSB 57:289-300.
#             Benjamini, Y. and D. Yekutieli.  2001.  The control of the false discovery rate in multiple testing under dependency. Annals of Statistics 29:1165-1188.
#             Benjamini, Y., A. Krieger, and D. Yekutieli.  2001.  Two staged linear step up FDR controlling procedure.  Technical Report, Department of Statistics and Operations Research, Tel Aviv University.  URL: http://www.math.tau.ac.il/~ybenja/Papers.html
#             Storey, J. 2002.  A direct approach to false discovery rates.  JRSSB 64: 479--498.
#

fdr <- function(pvals,qlevel=0.05,method="original",adjustment.method=NULL,adjustment.args=NULL){
#
# Description:
#
#    This is the main function designed for general usage for determining significance based on the FDR approach.
#
# Arguments:
#
#   pvals (required):  a vector of pvals on which to conduct the multiple testing
#
#   qlevel: the proportion of false positives desired
#
#   method: method for performing the testing.  'original' follows Benjamini & Hochberg (1995); 'general' is much more conservative, requiring no assumptions on the p-values (see Benjamini & Yekutieli (2001)).  We recommend using 'original', and if desired, using 'adjustment.method="mean" ' to increase power.
#
#   adjustment.method: method for increasing the power of the procedure by estimating the proportion of alternative p-values, one of "mean", the modified Storey estimator that we suggest in Ventura et al. (2004), "storey", the method of Storey (2002), or "two-stage", the iterative approach of Benjamini et al. (2001)
#
#   adjustment.args: arguments to adjustment.method; see .prop.alt() for description, but note that for "two-stage", qlevel and fdr.method are taken from the qlevel and method arguments to fdr()
#
# Value:
#
#   NULL if no significant tests, or a vector of the indices of the significant tests
#
# Examples:
#
#   signif <- fdr(pvals,method="original",adjustment.method="mean")
#
  n <- length(pvals)

  a <- 0   # initialize proportion of alternative hypotheses
  if(!is.null(adjustment.method)){
    if(adjustment.method=="two-stage"){  # set things up for the "two-stage" estimator
      qlevel <- qlevel/(1+qlevel)  # see Benjamini et al. (2001) for proof that this controls the FDR at level qlevel
      adjustment.args$qlevel <- qlevel
      adjustment.args$fdr.method <- method
      cat(paste('Adjusting cutoff using two-stage method, with method ',method,' and qlevel ',round(qlevel,4),'\n',sep=""))
    }
    if(adjustment.method=="mean" & is.null(adjustment.args)){
      adjustment.args <- list(edf.lower=0.8,num.steps=20)  # default arguments for "mean" method of Ventura et al. (2004)
      cat(paste('Adjusting cutoff using mean method, with edf.lower=0.8 and num.steps=20\n',sep=""))
    }
    a <- .prop.alt(pvals,adjustment.method,adjustment.args)
  }
  if(a==1){    # all hypotheses are estimated to be alternatives
    return(1:n)
  } else{      # adjust for estimate of a; default is 0
    qlevel <- qlevel/(1-a)
  }

  return(.fdr.master(pvals,qlevel,method))
}

.fdr.master <- function(pvals,qlevel=0.05,method="original"){
#
# Description:
#
#    This is an internal function that performs various versions of the FDR procedure, but without the modification described in section 4 of our J of Climate paper.
#
# Arguments:
#
#   pvals (required):  a vector of pvals on which to conduct the multiple testing
#
#   qlevel: the proportion of false positives desired
#
#   method: one of 'original', the original method of Benjamini & Hochberg (1995), or 'general', the method of Benjamini & Yekutieli (2001), which requires no assumptions about the p-values, but which is much more conservative.  We recommend 'original' for climatological data, and suspect it works well generally for spatial data.
#
# Value:
#
#   NULL if no significant tests, or a vector of the indices of the significant tests
#
  n <- length(pvals)
  if(method=="general"){
    qlevel <- qlevel/sum(1/(1:n))  # This requires fewer assumptions but is much more conservative
  } else{
    if(method!="original"){
      stop(paste("No method of type: ",method,sep=""))
    }
  }
  return(.fdr.basic(pvals,qlevel))
}


.fdr.basic <- function(pvals,qlevel=0.05){
#
# Description:
#
#    This is an internal function that performs the basic FDR of Benjamini & Hochberg (1995).
#
# Arguments:
#
#   pvals (required):  a vector of pvals on which to conduct the multiple testing
#
#   qlevel: the proportion of false positives desired
#
# Value:
#
#   NULL if no significant tests, or a vector of the indices of the significant tests
#
  n <- length(pvals)
  sorted.pvals <- sort(pvals)
  sort.index <- order(pvals)
  indices <- (1:n)*(sorted.pvals<=qlevel*(1:n)/n)
  num.reject <- max(indices)
  if(num.reject){
    indices <- 1:num.reject
    return(sort(sort.index[indices]))  
  } else{
    return(NULL)
  }
}


.storey <- function(edf.quantile,pvals){
#
# Description:
#
#    This is an internal function that calculates the basic Storey (2002) estimator of a, the proportion of alternative hypotheses.
#
# Arguments:
#
#   edf.quantile (required):  the quantile of the empirical distribution function at which to estimate a
# 
#   pvals (required):  a vector of pvals on which to estimate a
#
# Value:
#
#   estimate of a, the number of alternative hypotheses
#
#
  if(edf.quantile >=1 | edf.quantile <=0){
    stop('edf.quantile should be between 0 and 1')
  }
  a <- (mean(pvals<=edf.quantile)-edf.quantile)/(1-edf.quantile)
  if(a>0){
    return(a)
  } else{
    return(0)
  }
}


.prop.alt <- function(pvals,adjustment.method="mean",adjustment.args=list(edf.lower=0.8,num.steps=20)){
#
# Description:
#
#    This is an internal function that calculates an estimate of a, the proportion of alternative hypotheses, using one of several methods.
#
# Arguments:
#
#   pvals (required):  a vector of pvals from which to estimate a
#
#   adjustment.method: method for  estimating the proportion of alternative p-values, one of "mean", the modified Storey estimator suggested in Ventura et al. (2004); ".storey", the method of Storey (2002); or "two-stage", the iterative approach of Benjamini et al. (2001)
#
#   adjustment.args: arguments to adjustment.method;
#      for "mean", specify edf.lower, the smallest quantile at which to estimate a, and num.steps, the number of quantiles to use - the approach uses the average of the Storey (2002) estimator for the num.steps quantiles starting at edf.lower and finishing just less than 1
#      for ".storey", specify edf.quantile, the quantile at which to calculate the estimator
#      for "two-stage", the method uses a standard FDR approach to estimate which p-values are significant; this number is the estimate of a; therefore the method requires specification of qlevel, the proportion of false positives and fdr.method ('original' or 'general'), the FDR method to be used.  We do not recommend 'general' as this is very conservative and will underestimate a.
#  
# Value:
#
#   estimate of a, the number of alternative hypotheses
#
# Examples:
#
#   a <- .prop.alt(pvals,adjustment.method="mean")
#
  n <- length(pvals)
  if(adjustment.method=="two-stage"){
    if(is.null(adjustment.args$qlevel) | is.null(adjustment.args$fdr.method)){
      stop("adjustment.args$qlevel or adjustment.args$fdr.method not specified.  Two-stage estimation of the number of alternative hypotheses requires specification of the FDR threshold and FDR method ('original' or 'general')")
    }
    return(length(.fdr.master(pvals,adjustment.args$qlevel,method=adjustment.args$fdr.method))/n)
  }

  if(adjustment.method==".storey"){
    if(is.null(adjustment.args$edf.quantile)){
      stop("adjustment.args$edf.quantile not specified. Using Storey's method for estimating  the number of alternative hypotheses requires specification of the argument of the p-value EDF at which to do the estimation (a number close to one is recommended)")
    }
    return(.storey(adjustment.args$edf.quantile,pvals))
  }

  if(adjustment.method=="mean"){
    if(is.null(adjustment.args$edf.lower) | is.null(adjustment.args$num.steps)){
      stop("adjustment.args$edf.lower or adjustment.args$num.steps is not specified. Using the method of Ventura et al. (2004) for estimating  the number of alternative hypotheses requires specification of the lowest quantile of the p-value EDF at which to do the estimation (a number close to one is recommended) and the number of steps between edf.lower and 1, starting at edf.lower, at which to do the estimation")
    }
    if(adjustment.args$edf.lower >=1 | adjustment.args$edf.lower<=0){
      stop("adjustment.args$edf.lower must be between 0 and 1");
    }
    if(adjustment.args$num.steps<1 | adjustment.args$num.steps%%1!=0){
      stop("adjustment.args$num.steps must be an integer greater than 0")
    }
    stepsize <- (1-adjustment.args$edf.lower)/adjustment.args$num.steps
    edf.quantiles <- matrix(seq(from=adjustment.args$edf.lower,by=stepsize,len=adjustment.args$num.steps),nr=adjustment.args$num.steps,nc=1)
    a.vec <- apply(edf.quantiles,1,.storey,pvals)
    return(mean(a.vec))
  }
}




# The folowing functions are useful for gene co-expression network analysis and are the courtesy of 
# Steve Horvath, Bin Zhang, Jun Dong, Andy Yip
# Zhang B, Horvath S (2005) A General Framework for Weighted Gene Co-Expression Network Analysis. 
# Statistical Applications in Genetics and Molecular Biology.


#####################################################################################################
################################################################################################################################
# A) Assessing scale free topology and choosing the parameters of the adjacency function
#    using the scale free topology criterion (Zhang and Horvath 05)
#################################################################################################################


# ===================================================
#For hard thresholding, we use the .signum (step) function
if(exists(".signum") ) rm(.signum); 
.signum=function(corhelp,tau1){
	adjmat1= as.matrix(abs(corhelp)>=tau1)
	dimnames(adjmat1) <- dimnames(corhelp)
	diag(adjmat1) <- 0
	adjmat1
}

# ===================================================
# For soft thresholding, one can use the .sigmoid function 
# But we have focused on the power adjacency function in the tutorial...
if (exists(".sigmoid") ) rm(.sigmoid); 
.sigmoid=function(ss,mu1=0.8,alpha1=20){
	1/(1+exp(-alpha1*(ss-mu1)))
}







#This function is useful for speeding up the connectivity calculation.
#The idea is to partition the adjacency matrix into consecutive baches of a given size.
#In principle, the larger the batch size the faster is the calculation. But smaller batchsizes require #less memory...
# Input: gene expression data set where *rows* correspond to microarray samples and columns correspond to genes.  
if (exists(".SoftConnectivity") ) rm(.SoftConnectivity);
.SoftConnectivity=function(datE, power=6,batchsize=1500) {
	no.genes=dim(datE)[[2]]
	no.samples=dim(datE)[[1]]
	if (no.genes<no.samples | no.genes<10 | no.samples<5 ){
		print("Error: Something seems to be wrong. Make sure that the input data frame has genes as rows and array samples as columns. ") 
	} else {
		sum1=function(x) sum(x,na.rm=T)
		k=rep(NA,no.genes)
		ad1 <- NULL
		no.batches=as.integer(no.genes/ batchsize)
		if (no.batches>0) {
			for (i in 1:no.batches) {
			print(paste("batch number = ", i))
			index1=c(1:batchsize)+(i-1)* batchsize
			ad1=abs(cor(datE[,index1], datE,use="p",method="spearman"))^power
			ad1[is.na(ad1)]=0
			k[index1]=apply(ad1,1,sum1)
			} # end of for (i in 1:no.batches
		} # end of if (no.batches>0)...
		if (no.genes-no.batches*batchsize>0 ) {
			restindex=c((no.batches*batchsize+1):no.genes)
			ad1=abs(cor(datE[,restindex], datE,use="p",method="spearman"))^power
			ad1[is.na(ad1)]=0
			k[restindex]=apply(ad1,1,sum1)
		} # end of if
	} # end of else statement
	k
} # end of function





# ===================================================
# The function .PickHardThreshold can help one to estimate the cut-off value 
# when using the .signum (step) function.
# The first column lists the threshold ("cut"),
# the second column lists the corresponding p-value based on the Fisher Transform 
# of the co-expression. 
# The third column reports the resulting scale free topology fitting index R^2.
# The fourth column reports the slope of the fitting line, it shoud be negative for 
# biologically meaningul networks.
# The fifth column reports the fitting index for the truncated exponential model. 
# Usually we ignore it.
# The remaining columns list the mean, median and maximum resulting connectivity.
# To pick a hard threshold (cut) with the scale free topology criterion:
# aim for high scale free R^2 (column 3), high connectivity (col 6) and negative slope 
# (around -1, col 4).
# The output is a list with 2 components. The first component lists a sugggested cut-off
# while the second component contains the whole table.
# The removeFirst option removes the first point (k=0, P(k=0)) from the regression fit.
# no.breaks specifies how many intervals used to estimate the frequency p(k) i.e. the no. of points in the 
# scale free topology plot.

if (exists(".PickHardThreshold")) rm(.PickHardThreshold);
.PickHardThreshold <- function(datExpr1,RsquaredCut=0.85, cutvector=seq(0.1,0.9,by=0.05),
				removeFirst=FALSE,no.breaks=10, coexp.method="spearman"){
	no.genes <- dim(datExpr1)[[2]]
	no.samples= dim(datExpr1)[[1]]
	colname1=c("Cut","p-value", "scale law R^2", "slope="  ,"truncated R^2","mean(k)","median(k)","max(k)")
	datout=data.frame(matrix(NA,nrow=length(cutvector),ncol=length(colname1) ))
	names(datout)=colname1
	datout[,1]=cutvector
	for (i in c(1:length(cutvector) ) ){
		cut1=cutvector[i]
		datout[i,2]=2*(1-pt(sqrt(no.samples-1)*cut1/sqrt(1-cut1^2),no.samples-1))
	}
	if(exists("fun1")) rm(fun1)

	fun1=function(x,coexp.method) {
		if(coexp.method %in% c("spearman","pearson","kendall")){
			corx <- abs(cor(x,datExpr1,use="p", method=coexp.method))
		}else if(coexp.method == "euclid"){
			corx <- 1 - (as.matrix(dist(rbind(x,t(datExpr1))))/max(dist(t(datExpr1)),na.rm=TRUE))
			corx <- corx[1,2:ncol(corx)]
		}
		out1=rep(NA, length(cutvector) )
		for (j in c(1:length(cutvector))) {out1[j]=sum(corx>cutvector[j])}
		out1
	} # end of fun1

	datk=t(apply(datExpr1,2,function(x) fun1(x,coexp.method=coexp.method)))
	for (i in c(1:length(cutvector) ) ){
		nolinkshelp <- datk[,i]-1
		cut2=cut(nolinkshelp,no.breaks)
		binned.k=tapply(nolinkshelp,cut2,mean)
		freq1=as.vector(tapply(nolinkshelp,cut2,length)/length(nolinkshelp))
		# The following code corrects for missing values etc
		breaks1=seq(from=min(nolinkshelp),to=max(nolinkshelp),length=no.breaks+1)
		hist1=hist(nolinkshelp,breaks=breaks1,equidist=F,plot=FALSE,right=TRUE)
		binned.k2=hist1$mids
		binned.k=ifelse(is.na(binned.k),binned.k2,binned.k)
		binned.k=ifelse(binned.k==0,binned.k2,binned.k)
		freq1=ifelse(is.na(freq1),0,freq1)
		xx= as.vector(log10(binned.k))
		if(removeFirst) {freq1=freq1[-1]; xx=xx[-1]}
		plot(xx,log10(freq1+.000000001),xlab="log10(k)",ylab="log10(p(k))" )
		lm1= lm(as.numeric(log10(freq1+.000000001))~ xx )
		lm2=lm(as.numeric(log10(freq1+.000000001))~ xx+I(10^xx) )
		datout[i,3]=summary(lm1)$adj.r.squared 
		datout[i,4]=summary(lm1)$coefficients[2,1]  
		datout[i,5]=summary(lm2)$adj.r.squared
		datout[i,6]=mean(nolinkshelp)
		datout[i,7]= median(nolinkshelp)
		datout[i,8]= max(nolinkshelp) 
	} 
	datout=signif(datout,3) 
	print(data.frame(datout));
	# the cut-off is chosen as smallest cut with R^2>RsquaredCut 
	ind1=datout[,3]>RsquaredCut
	indcut=NA
	indcut=ifelse(sum(ind1)>0,min(c(1:length(ind1))[ind1]),indcut)
	# this estimates the cut-off value that should be used. 
	# Don't trust it. You need to consider slope and mean connectivity as well!
	cut.estimate=cutvector[indcut][[1]]
	list(estimate=cut.estimate, tablou=data.frame(datout));
	#write(paste("Cut estimate: ",cut.estimate,"\n",sep=""),file = paste(getwd(),"/",results.dir,"/scale_free_hard_threshold.txt",sep=""),append=FALSE)
	#write.table(datout,file = paste(getwd(),"/",results.dir,"/scale_free_hard_threshold.txt",sep=""),append=TRUE, col.names=TRUE,sep="\t")
} # end of function











# ===========================================================
# The function .PickSoftThreshold allows one to estimate the power parameter when using
# a soft thresholding approach with the use of the power function AF(s)=s^Power
# The function .PickSoftThreshold allows one to estimate the power parameter when using
# a soft thresholding approach with the use of the power function AF(s)=s^Power
# The removeFirst option removes the first point (k=1, P(k=1)) from the regression fit.
if (exists(".PickSoftThreshold")) rm(.PickSoftThreshold);
.PickSoftThreshold=function(datExpr1,RsquaredCut=0.85, powervector=c(seq(1,10,by=1),seq(12,20,by=2)),
				removeFirst=FALSE,no.breaks=10,coexp.method="spearman") {
	no.genes <- ncol(datExpr1)
	no.samples <- nrow(datExpr1)
	colname1=c("Power", "scale law R^2" ,"slope", "truncated R^2","mean(k)","median(k)","max(k)")
	datout <- as.data.frame(matrix(0,nrow=length(powervector),ncol=length(colname1) ))
	colnames(datout)=colname1
	datout[,1]=powervector
	
	if(exists("fun1")) rm(fun1)
	
	fun1=function(x,coexp.method) {
		if(coexp.method %in% c("spearman","pearson","kendall")){
			corx <- abs(cor(x,datExpr1,use="p", method=coexp.method))
		}else if(coexp.method == "euclid"){
			corx <- 1 - (as.matrix(dist(rbind(x,t(datExpr1))))/max(dist(t(datExpr1)),na.rm=TRUE))
			corx <- corx[1,2:ncol(corx)]
		}
		out1=rep(NA, length(powervector) )
		for (j in c(1:length(powervector))) {out1[j]=sum(corx^powervector[j])}
		out1
	} # end of fun1

	
	datk=t(apply(datExpr1,2,function(x) fun1(x,coexp.method=coexp.method)))
	for (i in c(1:length(powervector) ) ){
		nolinkshelp <- datk[,i]-1
		cut2=cut(nolinkshelp,no.breaks)
		binned.k=tapply(nolinkshelp,cut2,mean)
		freq1=as.vector(tapply(nolinkshelp,cut2,length)/length(nolinkshelp))
		# The following code corrects for missing values etc
		breaks1=seq(from=min(nolinkshelp),to=max(nolinkshelp),length=no.breaks+1)
		hist1=hist(nolinkshelp,breaks=breaks1,equidist=F,plot=FALSE,right=TRUE)
		binned.k2=hist1$mids
		binned.k=ifelse(is.na(binned.k),binned.k2,binned.k)
		binned.k=ifelse(binned.k==0,binned.k2,binned.k)
		freq1=ifelse(is.na(freq1),0,freq1)

		xx= as.vector(log10(binned.k))
		if(removeFirst) {freq1=freq1[-1]; xx=xx[-1]}
		plot(xx,log10(freq1+.000000001),xlab="log10(k)",ylab="log10(p(k))" )
		lm1= lm(as.numeric(log10(freq1+.000000001))~ xx )
		lm2=lm(as.numeric(log10(freq1+.000000001))~ xx+I(10^xx) )
		datout[i,2]=summary(lm1)$adj.r.squared 
		datout[i,3]=summary(lm1)$coefficients[2,1]  
		datout[i,4]=summary(lm2)$adj.r.squared
		datout[i,5]=mean(nolinkshelp)
		datout[i,6]= median(nolinkshelp)
		datout[i,7]= max(nolinkshelp) 
	} 
	datout=signif(datout,3) 
	print(data.frame(datout));
	# the cut-off is chosen as smallest cut with R^2>RsquaredCut 
	ind1=datout[,2]>RsquaredCut
	indcut=NA
	indcut=ifelse(sum(ind1)>0,min(c(1:length(ind1))[ind1]),indcut)
	# this estimates the power value that should be used. 
	# Don't trust it. You need to consider slope and mean connectivity as well!
	power.estimate=powervector[indcut][[1]]
	list(estimate=power.estimate, tablou=data.frame(datout));
	#write(paste("Power estimate: ",power.estimate,"\n",sep=""),file = paste(getwd(),"/",results.dir,"/scale_free_soft_threshold.txt",sep=""),append=FALSE)
	#write.table(datout,file = paste(getwd(),"/",results.dir,"/scale_free_soft_threshold.txt",sep=""),append=TRUE, col.names=TRUE,sep="\t")
}






# ===================================================
# The function .ScaleFreePlot creates a plot for checking scale free topology
# when truncated1=T is specificed, it provides the R^2 measures for the following
# degree distributions: a) scale free topology, b) log-log R^2 and c) truncated exponential R^2

# The function .ScaleFreePlot creates a plot for checking scale free topology
if(exists(".ScaleFreePlot")) rm(.ScaleFreePlot) ; 
.ScaleFreePlot=function(kk,no.breaks=10,AF1="" ,truncated1=FALSE, removeFirst=FALSE){
	cut1=cut(kk,no.breaks)
	binned.k=tapply(kk,cut1,mean)
	freq1=tapply(kk,cut1,length)/length(kk)
	# The following code corrects for missing values etc
	breaks1=seq(from=min(kk),to=max(kk),length=no.breaks+1)
	hist1=hist(kk,breaks=breaks1,equidist=F,plot=FALSE,right=TRUE)
	binned.k2=hist1$mids
	binned.k=ifelse(is.na(binned.k),binned.k2,binned.k)
	binned.k=ifelse(binned.k==0,binned.k2,binned.k)
	freq1=ifelse(is.na(freq1),0,freq1)
	plot(log10(binned.k),log10(freq1+.000000001),xlab="log10(k)",ylab="log10(p(k))" )
	xx= as.vector(log10(binned.k))
	if(removeFirst) {freq1=freq1[-1]; xx=xx[-1]}
	lm1=lm(as.numeric(log10(freq1+.000000001))~ xx )
	lines(xx,predict(lm1),col=1)
	OUTPUT=data.frame(ScaleFreeRsquared=round(summary(lm1)$adj.r.squared,2),Slope=round(lm1$coefficients[[2]],2))
	if (truncated1==TRUE) { 
		lm2=lm(as.numeric(log10(freq1+.000000001))~ xx+I(10^xx) );
		OUTPUT=data.frame(ScaleFreeRsquared=round(summary(lm1)$adj.r.squared,2),Slope=round(lm1$coefficients[[2]],2),
		TruncatedRsquared=round(summary(lm2)$adj.r.squared,2))
		print("the red line corresponds to the truncated exponential fit")
		lines(xx,predict(lm2),col=2);
		title(paste(AF1, 
		", scale free R^2=",as.character(round(summary(lm1)$adj.r.squared,2)),
		", slope=", round(lm1$coefficients[[2]],2),
	", trunc.R^2=",as.character(round(summary(lm2)$adj.r.squared,2))))} else { 
		title(paste(AF1, ", scale R^2=",as.character(round(summary(lm1)$adj.r.squared,2)) , 
		", slope=", round(lm1$coefficients[[2]],2)))
	}
	OUTPUT
} # end of function 









#################################################################################################################
################################################################################################################################
# B) Computing the topological overlap matrix 
#################################################################################################################
#################################################################################################################



# ===================================================
#The function .TOMdist computes a dissimilarity 
# based on the topological overlap matrix (Ravasz et al)
# Input: an Adjacency matrix with entries in [0,1]
if(exists(".TOMdist")) rm(.TOMdist);
.TOMdist=function(adjmat1, maxADJ=FALSE) {
	diag(adjmat1)=0;
	adjmat1[is.na(adjmat1)]=0;
	maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) ); 
	
	if (maxh1>1 | minh1 < -1 ) {
		print(paste("ERROR: the adjacency matrix contains entries that are larger than 1 or smaller than 0!!!, max=",maxh1,", min=",minh1)) 
	} else { 
		if (  max(c(as.dist(abs(adjmat1-t(adjmat1)))))>0   ) {
			print("ERROR: non-symmetric adjacency matrix!!!") 
		} else { 
			kk=apply(adjmat1,2,sum)
			maxADJconst=1
			if (maxADJ==TRUE) maxADJconst=max(c(as.dist(adjmat1 ))) 
			Dhelp1=matrix(kk,ncol=length(kk),nrow=length(kk))
			denomTOM= pmin(as.dist(Dhelp1),as.dist(t(Dhelp1)))   +as.dist(maxADJconst-adjmat1); 
			gc();gc();
			numTOM=as.dist(adjmat1 %*% adjmat1 +adjmat1);
			#TOMmatrix=numTOM/denomTOM
			# this turns the TOM matrix into a dissimilarity 
			out1=1-as.matrix(numTOM/denomTOM) 
			diag(out1)=1
			out1
		}
	}
}



# ===================================================
# This function computes a TOMk dissimilarity
# which generalizes the topological overlap matrix (Ravasz et al)
# Input: an Adjacency matrix with entries in [0,1]
# WARNING:  ONLY FOR UNWEIGHTED NETWORKS, i.e. the adjacency matrix contains binary entries...
# This function is explained in Yip and Horvath (2005)
# http://www.genetics.ucla.edu/labs/horvath/GTOM/
if(exists(".TOMkdist")) rm(.TOMkdist);
.TOMkdist = function(adjmat1,k=1){
    maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) );
    if (k!=round(abs(k))) {
        stop("k must be a positive integer!!!", call.=TRUE);}
    if (maxh1>1 | minh1 < 0 ){
        print(paste("ERROR: entries of the adjacency matrix must be between inclusively 0 and 1!!!, max=",maxh1,", min=",minh1))}
    else {
if (  max(c(as.dist(abs(adjmat1-t(adjmat1)))))>0   ) {print("ERROR: non-symmetric adjacency matrix!!!") } else { 

        B <- adjmat1;
        if (k>=2) {
            for (i in 2:k) {
                diag(B) <- diag(B) + 1;
                B = B %*% adjmat1;}}   # this gives the number of paths with length at most k connecting a pair
        B <- (B>0);   # this gives the k-step reachability from a node to another
        diag(B) <- 0;   # exclude each node being its own neighbor
        B <- B %*% B   # this gives the number of common k-step-neighbor that a pair of nodes share

        Nk <- diag(B);
        B <- B +adjmat1;   # numerator
        diag(B) <- 1;
        denomTOM=outer(Nk,Nk,FUN="pmin")+1-adjmat1;
        diag(denomTOM) <- 1;
        1 - B/denomTOM   # this turns the TOM matrix into a dissimilarity
}}
}


# IGNORE THIS function...
# The function .TOMdistROW computes the TOM distance of a gene (node)
# with that of all other genes in the network.
# WhichRow is an integer that specifies which row of the adjacency matrix
# corresponds to the gene of interest.
# Output=vector of TOM distances.
if (exists(".TOMdistROW") ) rm(.TOMdistROW) 
.TOMdistROW=function(WhichRow=1, adjmat1, maxADJ=FALSE) {
	diag(adjmat1)=0;
	maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) ); 
	if (maxh1>1 | minh1 < 0 ) {print(paste("ERROR: the adjacency matrix contains entries that are larger than 1 or smaller than 0!!!, 
						max=",maxh1,", min=",minh1)) } else { 
	kk=apply(adjmat1,2,sum)
	numTOM=adjmat1[WhichRow,] %*% adjmat1 +adjmat1[WhichRow,]; 
	numTOM[WhichRow]=1
	maxADJconst=1
	if (maxADJ==TRUE) maxADJconst=max(c(as.dist(adjmat1 ))) 
		denomTOM=pmin(kk[WhichRow],kk)+maxADJconst-adjmat1[WhichRow,]; denomTOM[WhichRow]=1
		#TOMmatrix=numTOM/denomTOM
		# this turns the TOM matrix into a dissimilarity 
		1-numTOM/denomTOM 
	}
}

