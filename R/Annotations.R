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




.LoadKEGG <- function(espece = espece){

	if(espece == "hsa"){ 
		DIResp <- "HS"
	}else if(espece == "mmu"){   
		DIResp <- "MM"
	}else if(espece == "rno"){   
		DIResp <- "RN"
	}else if(espece == "sce"){
		DIResp <- "SC"
	}


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





	if(espece == "hsa" || espece == "mmu" || espece == "rno"){		
		suppressWarnings (apply (DTtemp, 1, FUN = ".MiseEnFormeKegg" ,FichierFinalNonQuote))
	}else if(espece == "sce"){
		suppressWarnings (apply (DTtemp, 1, FUN = ".MiseEnFormeKeggSC" ,FichierFinalNonQuote))
	}

	if(espece == "hsa"){ 

		HS.KEGG.file.annot <<- read.table ( file = FichierFinalNonQuote ,na.strings = "",fill=TRUE,colClasses="character",sep= "\t",header = FALSE,quote = "",comment.char = "")
		write.table (HS.KEGG.file.annot, FichierFinal,col.names =FALSE,row.names=FALSE,append = TRUE, sep= "\t")
		colnames(HS.KEGG.file.annot) <<- c("GeneID","KEGG")

	}else if(espece == "mmu"){

		MM.KEGG.file.annot <<- read.table ( file = FichierFinalNonQuote ,na.strings = "",fill=TRUE,colClasses="character",sep= "\t",header = FALSE,quote = "",comment.char = "")
		write.table (MM.KEGG.file.annot, FichierFinal,col.names =FALSE,row.names=FALSE,append = TRUE, sep= "\t")
		colnames(MM.KEGG.file.annot) <<- c("GeneID","KEGG")

	}else if(espece == "rno"){

		RN.KEGG.file.annot <<- read.table ( file = FichierFinalNonQuote ,na.strings = "",fill=TRUE,colClasses="character",sep= "\t",header = FALSE,quote = "",comment.char = "")
		write.table (RN.KEGG.file.annot, FichierFinal,col.names =FALSE,row.names=FALSE,append = TRUE, sep= "\t")
		colnames(RN.KEGG.file.annot) <<- c("GeneID","KEGG")

	}else if(espece == "sce"){

		SC.KEGG.file.annot <<- read.table ( file = FichierFinalNonQuote, na.strings = "",fill=TRUE,colClasses="character",sep= "\t",header = FALSE,quote = "",comment.char = "")

		DTLL <- .GeneInfo()
		DTLL <- as.matrix(DTLL[DTLL[,1] == "4932",])
		rownames(DTLL) <- as.character(DTLL[,4])
		SC.KEGG.file.annot <<- as.matrix(SC.KEGG.file.annot)
		rownames(SC.KEGG.file.annot) <<- as.character(SC.KEGG.file.annot[,1])
		x <- DTLL[DTLL[,4] %in% as.character(SC.KEGG.file.annot[,1]),2]
		SC.KEGG.file.annot <<- SC.KEGG.file.annot[SC.KEGG.file.annot[,1] %in% names(x),]

		y <- NULL

		for(i in 1:nrow(SC.KEGG.file.annot)){				

			y <- rbind(y,c(x[SC.KEGG.file.annot[i,1]],SC.KEGG.file.annot[i,2]))
		}
		colnames(y) <- as.character(matrix("",1,ncol(y)))
		SC.KEGG.file.annot <<- y
		colnames(SC.KEGG.file.annot) <<- c("GeneID","KEGG")

		write.table (SC.KEGG.file.annot, FichierFinal, col.names =FALSE, row.names=FALSE, append = TRUE, sep= "\t") 

	}


	close (FichierFinalNonQuote)
	unlink(paste(getwd(),"/Annotations/KEGG/",DIResp,"/ll_keggNonQuote.txt",sep=""))
	close (FichierFinal)


}



.Kegg <- function () { 


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


	dir.create(paste(getwd(),"/Annotations/KEGG/HS",sep=""))
	.LoadKEGG("hsa")


	dir.create(paste(getwd(),"/Annotations/KEGG/MM",sep=""))
	.LoadKEGG("mmu")

	dir.create(paste(getwd(),"/Annotations/KEGG/RN",sep=""))
	.LoadKEGG("rno")


	dir.create(paste(getwd(),"/Annotations/KEGG/SC",sep=""))
	.LoadKEGG("sce")
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



.MakeGoHierarchy <- function (DATE="") {


	if(DATE == ""){
		DATE <- Sys.time()
		DATE <- substr(as.character(DATE),1,7)  
		DATE <- gsub("-","",DATE)   
	}

	NomFichier  <- paste("go_", DATE , "-termdb-tables.tar.gz" ,sep = "")

	download.file(paste ("http://godatabase-archive.stanford.edu/latest/", NomFichier,sep = ""),paste (getwd(),"/go_DATE-termdb-tables.tar.gz",sep=""),mode="wb")

	system("gunzip go_DATE-termdb-tables.tar.gz")
	system("tar -xf go_DATE-termdb-tables.tar")


	FichTerm2term <- file(paste(getwd(),"/go_",DATE,"-termdb-tables/term2term.txt",sep=""), "r") 
	FichTerm <- file(paste(getwd(),"/go_",DATE,"-termdb-tables/term.txt",sep=""), "r") 
	FichFinal <- file(paste(getwd(),"/Annotations/GO/go_hierarchy.txt",sep=""), "w") 


	DTterm2term <- read.table(file = FichTerm2term ,na.strings = "",fill=TRUE,colClasses="character",col.names =c(1:5),sep= "\t",header = FALSE, quote = "",comment.char = "")
	DTterm <- read.table(file = FichTerm ,na.strings = "",fill=TRUE,colClasses="character",col.names =c(1:6),sep= "\t",header = FALSE, quote = "",comment.char = "")
	GO.terms.hierarchy <<- cbind(DTterm[DTterm2term[,4],4],DTterm[ DTterm2term[,3],4])


	write.table(GO.terms.hierarchy, FichFinal, col.names =FALSE, row.names=FALSE, sep= "\t")
	colnames(GO.terms.hierarchy) <<- c("GO child","GO parent")


	close(FichTerm2term)
	close(FichTerm)
	close (FichFinal)
	unlink(paste (getwd(),"/go_DATE-termdb-tables.tar",sep=""))
	unlink(paste(getwd(),"/go_",DATE,"-termdb-tables",sep=""), recursive = TRUE)

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




.MakeCorrespondanceLLGO <- function (DTLoc2Go, DTGoTermsAndIds,DTGOHierarchy, espece = espece) {  

	if(espece == "hs"){
		DIResp <- "HS"
		DTLoc2Go <- DTLoc2Go[DTLoc2Go[,1] == "9606",2:3]
	}else if(espece == "mm"){
		DIResp <- "MM"
		DTLoc2Go <- DTLoc2Go[DTLoc2Go[,1] == "10090",2:3]
	}else if(espece == "rn"){
		DIResp <- "RN"
		DTLoc2Go <- DTLoc2Go[DTLoc2Go[,1] == "10116",2:3]
	}else if(espece == "sc"){
		DIResp <- "SC"
		DTLoc2Go <- DTLoc2Go[DTLoc2Go[,1] == "4932",2:3]
	}



	FichBioProcess <- file(paste(getwd(),"/Annotations/GO/",DIResp,"/biological_process.txt",sep= ""), "w") 
	FichCellComp <- file(paste(getwd(),"/Annotations/GO/",DIResp,"/cellular_component.txt",sep= ""), "w") 
	FichMolFunc <- file(paste(getwd(),"/Annotations/GO/",DIResp,"/molecular_function.txt",sep= ""), "w")



	DTLoc2GoFiltre1 <- DTLoc2Go[ DTLoc2Go[,2] %in% DTGoTermsAndIds[,1],1:2]



	if(espece == "hs"){
		HS.GO.DIR.BP.file.annot <<- .TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichBioProcess ,"P")
		colnames(HS.GO.DIR.BP.file.annot) <<- c("GeneID","GO")
		HS.GO.DIR.CC.file.annot <<- .TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichCellComp ,"C")
		colnames(HS.GO.DIR.CC.file.annot) <<- c("GeneID","GO")
		HS.GO.DIR.MF.file.annot <<- .TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichMolFunc ,"F")
		colnames(HS.GO.DIR.MF.file.annot) <<- c("GeneID","GO")
	}else if(espece == "mm"){
		MM.GO.DIR.BP.file.annot <<- .TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichBioProcess ,"P")
		colnames(MM.GO.DIR.BP.file.annot) <<- c("GeneID","GO")
		MM.GO.DIR.CC.file.annot <<- .TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichCellComp ,"C")
		colnames(MM.GO.DIR.CC.file.annot) <<- c("GeneID","GO")
		MM.GO.DIR.MF.file.annot <<- .TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichMolFunc ,"F")
		colnames(MM.GO.DIR.MF.file.annot) <<- c("GeneID","GO")
	}else if(espece == "rn"){
		RN.GO.DIR.BP.file.annot <<- .TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichBioProcess ,"P")
		colnames(RN.GO.DIR.BP.file.annot) <<- c("GeneID","GO")
		RN.GO.DIR.CC.file.annot <<- .TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichCellComp ,"C")
		colnames(RN.GO.DIR.CC.file.annot) <<- c("GeneID","GO")
		RN.GO.DIR.MF.file.annot <<- .TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichMolFunc ,"F")
		colnames(RN.GO.DIR.MF.file.annot) <<- c("GeneID","GO")
	}else if(espece == "sc"){
		SC.GO.DIR.BP.file.annot <<- .TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichBioProcess ,"P")
		colnames(SC.GO.DIR.BP.file.annot) <<- c("GeneID","GO")
		SC.GO.DIR.CC.file.annot <<- .TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichCellComp ,"C")
		colnames(SC.GO.DIR.CC.file.annot) <<- c("GeneID","GO")
		SC.GO.DIR.MF.file.annot <<- .TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichMolFunc ,"F")
		colnames(SC.GO.DIR.MF.file.annot) <<- c("GeneID","GO")

	}


	close(FichBioProcess)
	close(FichCellComp)
	close (FichMolFunc)


}



.MakeLocusNames <- function( espece, DTLL) {  

	if(espece == "hs"){
		DTLLesp <- DTLL[DTLL[,1] == "9606",2:ncol(DTLL)]
		DIResp <- "HS"
	}else if(espece == "mm"){
		DTLLesp <- DTLL[DTLL[,1] == "10090",2:ncol(DTLL)]
		DIResp <- "MM"
	}else if(espece == "rn"){
		DTLLesp <- DTLL[DTLL[,1] == "10116",2:ncol(DTLL)]
		DIResp <- "RN"
	}else if(espece == "sc"){
		DTLLesp <- DTLL[DTLL[,1] == "4932",2:ncol(DTLL)]
		DIResp <- "SC"
	}



	DTLLesp <- DTLLesp[, c(1,6,2)]
	DTLLesp <- unique.data.frame (DTLLesp) 
	DTLLesp <- cbind.data.frame(DTLLesp, row.names = DTLLesp[,1]) 


	FichierFinal <- file(paste(getwd(),"/Annotations/LL/",DIResp,"/locus_name.txt",sep= ""), "w") 

	write.table(DTLLesp, FichierFinal ,col.names =FALSE,row.names=FALSE, sep="\t")


	if(espece == "hs"){
		HS.locus.name <<- DTLLesp
		colnames(HS.locus.name) <<- c("GeneID","Name","Symbol")
	}else if(espece == "mm"){
		MM.locus.name <<- DTLLesp
		colnames(MM.locus.name) <<- c("GeneID","Name","Symbol")
	}else if(espece == "rn"){
		RN.locus.name <<- DTLLesp
		colnames(RN.locus.name) <<- c("GeneID","Name","Symbol")
	}else if(espece == "sc"){
		SC.locus.name <<- DTLLesp
		colnames(SC.locus.name) <<- c("GeneID","Name","Symbol")
	}

	close(FichierFinal)

}


.GeneToUnigene <- function(espece, DTLL=NULL){

	UGTemp <- file(paste(getwd(),"/Annotations/gene2unigene.txt",sep=""), "r")
	UGLLbrut <- read.table(file = UGTemp,na.strings = "-",fill=TRUE,colClasses="character",sep= "\t",header = FALSE, quote = "",comment.char = "#")
	close(UGTemp)


	if(espece == "hs" & !is.null(DTLL)){
		DTLLesp <- DTLL[DTLL[,1] == "9606",2:ncol(DTLL)]
		DIResp <- "HS"
	}else if(espece == "mm" & !is.null(DTLL)){
		DTLLesp <- DTLL[DTLL[,1] == "10090",2:ncol(DTLL)]
		DIResp <- "MM"
	}else if(espece == "rn" & !is.null(DTLL)){
		DTLLesp <- DTLL[DTLL[,1] == "10116",2:ncol(DTLL)]
		DIResp <- "RN"
	}else if(espece == "sc" & !is.null(DTLL)){
		DTLLesp <- DTLL[DTLL[,1] == "4932",2:ncol(DTLL)]
		DIResp <- "SC"
	}



	DTLLesp <- as.character(DTLLesp[, 1])
	DTLLesp <- unique(DTLLesp)

	UGLLesp <- NULL

	if(espece %in% c("hs","mm","rn")){
		UGLLesp <- UGLLbrut[UGLLbrut[,1] %in% DTLLesp,]
		UGLLesp <- unique.data.frame(UGLLesp) 
		UGLLesp <- cbind.data.frame(UGLLesp, row.names = c(1:nrow(UGLLesp))) 
	}

	FichierFinal <- file(paste(getwd(),"/Annotations/LL/",DIResp,"/unigene.txt",sep= ""), "w") 

	write.table(UGLLesp, FichierFinal ,col.names =FALSE,row.names=FALSE, sep="\t")


	if(espece == "hs"){
		HS.unigene <<- UGLLesp
		colnames(HS.unigene) <<- c("GeneID","UniGene")
	}else if(espece == "mm"){
		MM.unigene <<- UGLLesp
		colnames(MM.unigene) <<- c("GeneID","UniGene")
	}else if(espece == "rn"){
		RN.unigene <<- UGLLesp
		colnames(RN.unigene) <<- c("GeneID","UniGene")
	}else if(espece == "sc"){
		SC.orf <<- DTLL[DTLL[,1] == "4932",c(2,3)]
		SC.orf <<- cbind.data.frame(SC.orf, row.names = as.character(SC.orf[,1]))
		colnames(SC.orf) <<- c("GeneID","ORF")
	}

	close(FichierFinal)

}






.goAndLL <- function(DATE=""){ 


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

	DTGoHierarchy <- .MakeGoHierarchy(DATE=DATE)


	DTLL <- .GeneInfo()


	dir.create(paste(getwd(),"/Annotations/GO/HS",sep=""))
	.MakeCorrespondanceLLGO(DTLoc2Go,DTGoTermsAndIds, DTGoHierarchy , "hs")


	dir.create(paste(getwd(),"/Annotations/GO/MM",sep=""))
	.MakeCorrespondanceLLGO (DTLoc2Go, DTGoTermsAndIds, DTGoHierarchy , "mm")

	dir.create(paste(getwd(),"/Annotations/GO/RN",sep=""))
	.MakeCorrespondanceLLGO (DTLoc2Go, DTGoTermsAndIds, DTGoHierarchy , "rn")


	dir.create(paste(getwd(),"/Annotations/GO/SC",sep=""))
	.MakeCorrespondanceLLGO (DTLoc2Go, DTGoTermsAndIds, DTGoHierarchy , "sc")




	dir.create(paste(getwd(),"/Annotations/LL",sep=""))


	dir.create(paste(getwd(),"/Annotations/LL/HS",sep=""))
	.MakeLocusNames("hs", DTLL)
	.GeneToUnigene("hs", DTLL)

	dir.create(paste(getwd(),"/Annotations/LL/MM",sep=""))
	.MakeLocusNames("mm", DTLL)
	.GeneToUnigene("mm", DTLL)

	dir.create(paste(getwd(),"/Annotations/LL/RN",sep=""))
	.MakeLocusNames("rn", DTLL)
	.GeneToUnigene("rn", DTLL)

	dir.create(paste(getwd(),"/Annotations/LL/SC",sep=""))
	.MakeLocusNames("sc", DTLL)
	.GeneToUnigene("sc", DTLL)

}




########################################################################

annotations <- function(date.annot=""){   

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
		
	.Kegg()
	.goAndLL(DATE=date.annot)
	
	unlink(paste(getwd(),"/Annotations/gene_info.txt",sep=""))
	unlink(paste(getwd(),"/Annotations/gene2unigene.txt",sep=""))
	annot.date <- format(Sys.time(), "%Y %b %d")
		
	save(GO.terms.hierarchy,GO.terms.name,HS.GO.DIR.BP.file.annot,HS.GO.DIR.CC.file.annot,HS.GO.DIR.MF.file.annot,
		HS.KEGG.file.annot,HS.locus.name,HS.unigene,KEGG.terms.name,MM.GO.DIR.BP.file.annot,MM.GO.DIR.CC.file.annot,MM.GO.DIR.MF.file.annot,
		MM.KEGG.file.annot,MM.locus.name,MM.unigene,RN.GO.DIR.BP.file.annot,RN.GO.DIR.CC.file.annot,RN.GO.DIR.MF.file.annot,
		RN.KEGG.file.annot,RN.locus.name,RN.unigene,SC.GO.DIR.BP.file.annot,SC.GO.DIR.CC.file.annot,SC.GO.DIR.MF.file.annot,
		SC.KEGG.file.annot,SC.locus.name,SC.orf,annot.date,file="sysdata.rda",compress=TRUE) 
	
	unlink(paste(getwd(),"/Annotations", sep= ""), recursive=TRUE)
	
	q(save="no")

}


