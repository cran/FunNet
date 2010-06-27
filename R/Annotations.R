.packageName <- "FunNet"


.MiseEnFormeKegg <- function (ligne, Fichier){ 
	vec1 <- ligne
	LL <-as.character(vec1[1]) 
	if (!identical(all.equal( as.numeric(LL), as.numeric(NA)),TRUE) ){   
		Kmap <- as.character(vec1[-1])
		Kmap <- Kmap[Kmap != ""]  
		NbKmap <- length(Kmap)

		LocusLink  <- gl(1, NbKmap, label = LL) 

		LLcorrespondance <- cbind.data.frame(LocusLink, Kmap )  
		write.table(LLcorrespondance, Fichier,col.names =FALSE,row.names=FALSE,append = TRUE, sep= "\t",quote = FALSE)  


	}

}


.MiseEnFormeKeggSC <- function(ligne, Fichier) { 
	vec1 <- ligne
	LL <-as.character(vec1[1]) 

	Kmap <- as.character(vec1[-1])
	Kmap <- Kmap[Kmap != ""]  
	NbKmap <- length(Kmap)

	LocusLink  <- gl(1, NbKmap, label = LL) 

	LLcorrespondance <- cbind.data.frame(LocusLink, Kmap)  
	write.table(LLcorrespondance, Fichier, col.names = FALSE, row.names=FALSE,append = TRUE, sep= "\t",quote = FALSE)  
}




.LoadKEGG <- function(espece){

	DIResp <- espece

	download.file(paste("ftp://ftp.genome.ad.jp/pub/kegg/pathway/organisms/",espece,"/",espece,"_gene_map.tab",sep=""), paste(getwd(),"/Annotations/KEGG/",DIResp,"/tmp.txt",sep=""), mode = "w" )

	FichierFinalNonQuote <- file(paste(getwd(),"/Annotations/KEGG/",DIResp,"/ll_keggNonQuote.txt",sep=""), "w+")

	FichierFinal <- file(paste(getwd(),"/Annotations/KEGG/",DIResp,"/ll_kegg.txt",sep=""), "w")

	tempo <- file(paste(getwd(),"/Annotations/KEGG/",DIResp,"/tmp.txt",sep=""), "r") 
	tempo2 <- file(paste(getwd(),"/Annotations/KEGG/",DIResp,"/tmp2.txt",sep=""), "w+")



	NbCol <- count.fields(file=tempo)
	MaxNbCol <- max(NbCol)


	seek(tempo, where = 0)


	DTtemp <- read.table ( file = tempo ,na.strings = "",fill=TRUE,colClasses="character",sep= "\t",header = FALSE,quote = "",comment.char = "")

	write.table (DTtemp, tempo2 ,col.names =FALSE,row.names=FALSE,quote =FALSE, sep = "\t")
	DTtemp <- read.table ( file = tempo2 ,fill=TRUE,colClasses="character",col.names =c(1: (MaxNbCol)),sep= "", quote = "",comment.char = "")



	close(tempo)
	close(tempo2)
	unlink(paste(getwd(),"/Annotations/KEGG/",DIResp,"/tmp.txt",sep=""))
	unlink (paste(getwd(),"/Annotations/KEGG/",DIResp,"/tmp2.txt",sep=""))





	if(espece != "sce"){
	
		suppressWarnings (apply (DTtemp, 1, FUN = ".MiseEnFormeKegg" ,FichierFinalNonQuote))
		annot.base[[espece]]$KEGG.file.annot <<- read.table (file = FichierFinalNonQuote ,na.strings = "",fill=TRUE,colClasses="character",sep= "\t",header = FALSE,quote = "",comment.char = "")
		
		write.table (annot.base[[espece]]$KEGG.file.annot, FichierFinal,col.names =FALSE,row.names=FALSE,append = TRUE, sep= "\t")
		colnames(annot.base[[espece]]$KEGG.file.annot) <<- c("GeneID","KEGG")
	
	}else if(espece == "sce"){
	
		suppressWarnings (apply(DTtemp, 1, FUN = ".MiseEnFormeKeggSC" ,FichierFinalNonQuote))
		annot.base$sce$KEGG.file.annot <<- read.table(file = FichierFinalNonQuote, na.strings = "",fill=TRUE,colClasses="character",sep= "\t",header = FALSE,quote = "",comment.char = "")
		
		DTLL <- .GeneInfo()
		DTLL <- as.matrix(DTLL[DTLL[,1] == species[species[,"name"]=="sce","taxid"],])
		rownames(DTLL) <- as.character(DTLL[,4])
		annot.base$sce$KEGG.file.annot <<- as.matrix(annot.base$sce$KEGG.file.annot)
		rownames(annot.base$sce$KEGG.file.annot) <<- as.character(annot.base$sce$KEGG.file.annot[,1])
		x <- DTLL[DTLL[,4] %in% as.character(annot.base$sce$KEGG.file.annot[,1]),2]
		annot.base$sce$KEGG.file.annot <<- annot.base$sce$KEGG.file.annot[annot.base$sce$KEGG.file.annot[,1] %in% names(x),]

		y <- NULL

		for(i in 1:nrow(annot.base$sce$KEGG.file.annot)){				

			y <- rbind(y,c(x[annot.base$sce$KEGG.file.annot[i,1]],annot.base$sce$KEGG.file.annot[i,2]))
		}
		colnames(y) <- as.character(matrix("",1,ncol(y)))
		annot.base$sce$KEGG.file.annot <<- y
		colnames(annot.base$sce$KEGG.file.annot) <<- c("GeneID","KEGG")
		
		write.table (annot.base$sce$KEGG.file.annot, FichierFinal, col.names =FALSE, row.names=FALSE, append = TRUE, sep= "\t") 
	}


	close (FichierFinalNonQuote)
	unlink(paste(getwd(),"/Annotations/KEGG/",DIResp,"/ll_keggNonQuote.txt",sep=""))
	close (FichierFinal)


}



.Kegg <- function (species) { 


	dir.create(paste(getwd(),"/Annotations/KEGG",sep=""))


	download.file ( "ftp://ftp.genome.ad.jp/pub/kegg/pathway/map_title.tab" ,paste(getwd(), "/Annotations/KEGG/Temp_kegg_terms.txt",sep=""), mode = "w" )
	temp <- file (paste(getwd(),"/Annotations/KEGG/Temp_kegg_terms.txt",sep=""),"r")
	Fichkeggterms <- file (paste(getwd(),"/Annotations/KEGG/kegg_terms.txt",sep=""),"w")

	KEGG.terms.name <<- read.table (file = temp ,na.strings = "",fill=TRUE,colClasses="character",sep= "\t",header = FALSE,quote = "",comment.char = "")
	write.table (KEGG.terms.name, Fichkeggterms ,col.names =FALSE,row.names=FALSE, sep = "\t")
	colnames(KEGG.terms.name) <<- c("KEGG","Name")

	close(temp)
	close(Fichkeggterms)
	unlink(paste(getwd(),"/Annotations/KEGG/Temp_kegg_terms.txt",sep=""))

	for(i in 1:nrow(species)){	
		dir.create(paste(getwd(),"/Annotations/KEGG/",species[i,"name"],sep=""))
		.LoadKEGG(species[i,"name"])	
	}

}



.MakeGoTerms <- function (DTGoTermsAndIds) {


	DTGoTermsAndIds <- DTGoTermsAndIds[substr(DTGoTermsAndIds[,1],1,1) != "!" ,]
	dimens <- nrow(DTGoTermsAndIds)
	DTGoTermsAndIds <- cbind.data.frame(DTGoTermsAndIds, row.names = c(1: dimens)) 
	GO.terms.name <<- DTGoTermsAndIds[,1:2] 


	FichierFinal <- file(paste(getwd(),"/Annotations/GO/go_terms.txt",sep=""), "w") 
	write.table (GO.terms.name, FichierFinal, col.names =FALSE, row.names=FALSE, sep = "\t")

	colnames(GO.terms.name) <<- c("GO","Name")

	close(FichierFinal)

}



.MakeGoHierarchy <- function () {


	NomFichier  <- "go_daily-termdb-tables.tar.gz"

	download.file(paste ("http://archive.geneontology.org/latest-termdb/", NomFichier,sep = ""),paste (getwd(),"/go_daily-termdb-tables.tar.gz",sep=""),mode="wb")

	system("gunzip go_daily-termdb-tables.tar.gz")
	system("tar -xf go_daily-termdb-tables.tar")


	FichTerm2term <- file(paste(getwd(),"/go_","daily","-termdb-tables/term2term.txt",sep=""), "r") 
	FichTerm <- file(paste(getwd(),"/go_","daily","-termdb-tables/term.txt",sep=""), "r") 
	FichFinal <- file(paste(getwd(),"/Annotations/GO/go_hierarchy.txt",sep=""), "w") 


	DTterm2term <- read.table(file = FichTerm2term ,na.strings = "",fill=TRUE,colClasses="character",col.names =c(1:5),sep= "\t",header = FALSE, quote = "",comment.char = "")
	DTterm <- read.table(file = FichTerm ,na.strings = "",fill=TRUE,colClasses="character",col.names =c(1:7),sep= "\t",header = FALSE, quote = "",comment.char = "")
	GO.terms.hierarchy <<- cbind(DTterm[DTterm2term[,4],4],DTterm[ DTterm2term[,3],4])


	write.table(GO.terms.hierarchy, FichFinal, col.names =FALSE, row.names=FALSE, sep= "\t")
	colnames(GO.terms.hierarchy) <<- c("GO child","GO parent")


	close(FichTerm2term)
	close(FichTerm)
	close (FichFinal)
	unlink(paste (getwd(),"/go_daily-termdb-tables.tar",sep=""))
	unlink(paste(getwd(),"/go_","daily","-termdb-tables",sep=""), recursive = TRUE)

	return(GO.terms.hierarchy)

}




.TroisiemeFiltrage <- function (DTGoTermsAndIds,DTLoc2GoFiltre1, FichOntology ,initialeOntology) {

	OntDTGoTermsAndIds <- DTGoTermsAndIds[ DTGoTermsAndIds[,3] == initialeOntology,]
	DTOntology <- DTLoc2GoFiltre1 [ DTLoc2GoFiltre1[,2] %in% OntDTGoTermsAndIds[,1],]

	DTOntology <- unique.data.frame(DTOntology) 

	OrdreLignes <- do.call("order", DTOntology[,c( 1,2)])  
	DTOntology <- DTOntology [ OrdreLignes,]

	dimens <- nrow(DTOntology)
	DTOntology <- cbind.data.frame(DTOntology, row.names = c(1: dimens))

	write.table(DTOntology, FichOntology ,col.names =FALSE,row.names=FALSE,sep="\t")
	return(DTOntology)
}



.GeneInfo <- function(){

	LLTemp <- file(paste(getwd(),"/Annotations/gene_info.txt",sep=""), "r")

	DTLLbrut <- read.table(file = LLTemp,na.strings = "-",fill=TRUE,colClasses="character",sep= "\t",header = FALSE, quote = "",comment.char = "#")

	DTLL <- cbind(DTLLbrut[,1],DTLLbrut[,2],DTLLbrut[,3],DTLLbrut[,4],DTLLbrut[,7],DTLLbrut[,8],DTLLbrut[,9])

	close(LLTemp)
	return(DTLL) 


}




.MakeCorrespondanceLLGO <- function (DTLoc2Go, DTGoTermsAndIds,DTGOHierarchy, espece) {  

	DIResp <- espece
	DTLoc2Go <- DTLoc2Go[DTLoc2Go[,1] == species[species[,"name"]==espece,"taxid"],2:3]



	FichBioProcess <- file(paste(getwd(),"/Annotations/GO/",DIResp,"/biological_process.txt",sep= ""), "w") 
	FichCellComp <- file(paste(getwd(),"/Annotations/GO/",DIResp,"/cellular_component.txt",sep= ""), "w") 
	FichMolFunc <- file(paste(getwd(),"/Annotations/GO/",DIResp,"/molecular_function.txt",sep= ""), "w")



	DTLoc2GoFiltre1 <- DTLoc2Go[ DTLoc2Go[,2] %in% DTGoTermsAndIds[,1],1:2]



	annot.base[[espece]]$GO.DIR.BP.file.annot <<- .TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichBioProcess ,"P")
	colnames(annot.base[[espece]]$GO.DIR.BP.file.annot) <<- c("GeneID","GO")
	annot.base[[espece]]$GO.DIR.CC.file.annot <<- .TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichCellComp ,"C")
	colnames(annot.base[[espece]]$GO.DIR.CC.file.annot) <<- c("GeneID","GO")
	annot.base[[espece]]$GO.DIR.MF.file.annot <<- .TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichMolFunc ,"F")
	colnames(annot.base[[espece]]$GO.DIR.MF.file.annot) <<- c("GeneID","GO")


	close(FichBioProcess)
	close(FichCellComp)
	close (FichMolFunc)


}



.MakeLocusNames <- function( espece, DTLL) {  

	DTLLesp <- DTLL[DTLL[,1] == species[species[,"name"]==espece,"taxid"],2:ncol(DTLL)]
	DIResp <- espece



	DTLLesp <- DTLLesp[, c(1,6,2)]
	DTLLesp <- unique.data.frame (DTLLesp) 
	DTLLesp <- cbind.data.frame(DTLLesp, row.names = DTLLesp[,1]) 


	FichierFinal <- file(paste(getwd(),"/Annotations/LL/",DIResp,"/locus_name.txt",sep= ""), "w") 

	write.table(DTLLesp, FichierFinal ,col.names =FALSE,row.names=FALSE, sep="\t")


	annot.base[[espece]]$locus.name <<- DTLLesp
	colnames(annot.base[[espece]]$locus.name) <<- c("GeneID","Name","Symbol")

	close(FichierFinal)

}


.GeneToUnigene <- function(espece, DTLL=NULL){

	UGTemp <- file(paste(getwd(),"/Annotations/gene2unigene.txt",sep=""), "r")
	UGLLbrut <- read.table(file = UGTemp,na.strings = "-",fill=TRUE,colClasses="character",sep= "\t",header = FALSE, quote = "",comment.char = "#")
	close(UGTemp)


	if(!is.null(DTLL)){
		DTLLesp <- DTLL[DTLL[,1] == species[species[,"name"]==espece,"taxid"],2:ncol(DTLL)]
		DIResp <- espece
	}



	DTLLesp <- as.character(DTLLesp[, 1])
	DTLLesp <- unique(DTLLesp)

	UGLLesp <- NULL

	if(espece != "sce"){
		UGLLesp <- UGLLbrut[UGLLbrut[,1] %in% DTLLesp,]
		UGLLesp <- unique.data.frame(UGLLesp) 
		UGLLesp <- cbind.data.frame(UGLLesp, row.names = c(1:nrow(UGLLesp))) 
	}

	FichierFinal <- file(paste(getwd(),"/Annotations/LL/",DIResp,"/unigene.txt",sep= ""), "w") 

	write.table(UGLLesp, FichierFinal ,col.names =FALSE,row.names=FALSE, sep="\t")


	if(espece != "sce"){
		annot.base[[espece]]$unigene <<- UGLLesp
		colnames(annot.base[[espece]]$unigene) <<- c("GeneID","UniGene")
	}else if(espece == "sce"){
		annot.base$sce$orf <<- DTLL[DTLL[,1] == species[species[,"name"]=="sce","taxid"],c(2,3)]
		annot.base$sce$orf <<- cbind.data.frame(annot.base$sce$orf, row.names = as.character(annot.base$sce$orf[,1]))
		colnames(annot.base$sce$orf) <<- c("GeneID","ORF")
	}

	close(FichierFinal)

}






.goAndLL <- function(species){ 


	dir.create(paste(getwd(),"/Annotations/GO",sep=""))



	download.file("http://www.geneontology.org/doc/GO.terms_and_ids" ,paste(getwd(),"/Annotations/GO/GoTermsTmp.txt",sep=""), mode = "w")

	FichGoTermsAndIds <- file(paste(getwd(),"/Annotations/GO/GoTermsTmp.txt",sep=""), "r") 

	DTGoTermsAndIds <- read.table(file = FichGoTermsAndIds ,na.strings = "",fill=TRUE,colClasses="character",col.names =c(1:3),sep= "\t",header = FALSE, quote = "",comment.char = "")
	close(FichGoTermsAndIds)
	unlink(paste(getwd(),"/Annotations/GO/GoTermsTmp.txt",sep=""))



	download.file("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz" ,paste(getwd(),"/Annotations/GO/gene2go.gz",sep=""), mode = "wb")

	Gene2GOTempgz <- gzfile(paste(getwd(),"/Annotations/GO/gene2go.gz",sep=""),"rb") 

	Gene2GOTemp <- file(paste(getwd(),"/Annotations/GO/gene2go.txt",sep=""), "w+")
	vecTemp <- readLines(Gene2GOTempgz)
	cat(vecTemp , file = Gene2GOTemp, sep = "\n")


	DTLoc2Go <- read.table(file = Gene2GOTemp,na.strings = "",fill=TRUE ,colClasses="character",col.names =c(1:8),sep= "\t",header = FALSE, quote = "",comment.char = "#")
	close(Gene2GOTemp)
	unlink(paste(getwd(),"/Annotations/GO/gene2go.txt",sep=""))
	close(Gene2GOTempgz)
	unlink(paste(getwd(),"/Annotations/GO/gene2go.gz",sep=""))



	.MakeGoTerms(DTGoTermsAndIds) 

	DTGoHierarchy <- .MakeGoHierarchy()


	DTLL <- .GeneInfo()

	for(i in 1:nrow(species)){
		dir.create(paste(getwd(),"/Annotations/GO/",species[i,"name"],sep=""))
		.MakeCorrespondanceLLGO(DTLoc2Go,DTGoTermsAndIds, DTGoHierarchy , species[i,"name"])	
	}


	dir.create(paste(getwd(),"/Annotations/LL",sep=""))

	for(i in 1:nrow(species)){
		dir.create(paste(getwd(),"/Annotations/LL/",species[i,"name"],sep=""))
		.MakeLocusNames(species[i,"name"], DTLL)
		.GeneToUnigene(species[i,"name"], DTLL)
	}

}




########################################################################

annotations <- function(cust.specs=NULL){   

	species <<- cust.specs
	annot.base <<- list(NULL)

	if(is.null(species)){
		species <<- cbind(c("hsa","mmu","rno","sce","gga"),c("9606","10090","10116","559292","9031"))
		colnames(species) <<- c("name","taxid")
		rownames(species) <<- species[,"taxid"]
	}
	
	GO.DIR.BP.file.annot <- NULL
	GO.DIR.CC.file.annot <- NULL
	GO.DIR.MF.file.annot <- NULL
	KEGG.file.annot <- NULL
	locus.name <- NULL
	unigene <- NULL
	orf <- NULL
	
	for(i in 1:nrow(species)){
		if(species[i,"name"] != "sce"){
			annot.base[[i]] <<- list(GO.DIR.BP.file.annot,GO.DIR.CC.file.annot,GO.DIR.MF.file.annot,KEGG.file.annot,locus.name,unigene)
			names(annot.base[[i]]) <<- c("GO.DIR.BP.file.annot","GO.DIR.CC.file.annot","GO.DIR.MF.file.annot","KEGG.file.annot","locus.name","unigene")
		}else if(species[i,"name"] == "sce"){
			annot.base[[i]] <<- list(GO.DIR.BP.file.annot,GO.DIR.CC.file.annot,GO.DIR.MF.file.annot,KEGG.file.annot,locus.name,orf)
			names(annot.base[[i]]) <<- c("GO.DIR.BP.file.annot","GO.DIR.CC.file.annot","GO.DIR.MF.file.annot","KEGG.file.annot","locus.name","orf")
		}	
	}
	
	names(annot.base) <<- species[,"name"]
	rm(GO.DIR.BP.file.annot,GO.DIR.CC.file.annot,GO.DIR.MF.file.annot,KEGG.file.annot,locus.name,unigene,orf)

	dir.create(paste(getwd(),"/Annotations", sep= ""))
	
	download.file(paste("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz",sep="") ,paste( getwd(),"/Annotations/Temp.gz",sep=""), mode = "wb" )
	
	LLTempgz <- gzfile(paste(getwd(),"/Annotations/Temp.gz",sep=""),"rb")
	LLTemp <- file(paste(getwd(),"/Annotations/gene_info.txt",sep=""), "w+")
	vecTemp <- readLines(LLTempgz)
	cat(vecTemp, file = LLTemp, sep = "\n")	
	rm(vecTemp)
	close(LLTempgz)
	close(LLTemp)
	unlink(paste(getwd(),"/Annotations/Temp.gz",sep=""))
	
	download.file(paste("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2unigene",sep="") ,paste( getwd(),"/Annotations/gene2unigene.txt",sep=""), mode = "wb" )
		
	.Kegg(species=species)
	.goAndLL(species=species)
	
	unlink(paste(getwd(),"/Annotations/gene_info.txt",sep=""))
	unlink(paste(getwd(),"/Annotations/gene2unigene.txt",sep=""))
	annot.date <- format(Sys.time(), "%Y %b %d")
		
	save(GO.terms.hierarchy,GO.terms.name,KEGG.terms.name,annot.base,annot.date,species,file="sysdata.rda",compress=TRUE) 
	
	unlink(paste(getwd(),"/Annotations", sep= ""), recursive=TRUE)
	
	q(save="no")

}


