  \name{FunNet}
          \alias{FunNet}
          \title{Integrative Functional Analysis of Transcriptional Networks}
          \description{
            FunNet is an integrative tool for analyzing gene co-expression networks built from 
            microarray expression data. The analytic model implemented in this library
            involves two abstraction layers: transcriptional and functional (biological roles). 
            A functional profiling technique using Gene Ontology & KEGG annotations is 
            applied to extract a list of relevant biological themes from 
            microarray expression profiling data. Afterwards multiple-instance 
            representations are built to relate significant themes to their 
            transcriptional instances (i.e. the two layers of the model). An adapted
            non-linear dynamical system model is used to quantify the proximity of relevant 
            genomic themes based on the similarity of the expression profiles of their gene instances.
            Eventually an unsupervised multiple-instance clustering procedure, relying on 
            the two abstraction layers, is used to identify the structure of the co-expression
            network composed from modules of functionally related transcripts. Functional 
            and transcriptional maps of the co-expression network are provided separately
            together with detailed information on the network centrality of related transcripts
            and genomic themes.
          }
          \usage{
          FunNet(wd="", org="hsa", two.lists=TRUE, up.frame=NULL, down.frame=NULL,
	  	genes.frame=NULL, restrict=FALSE, ref.list=NULL, logged=FALSE,
	  	discriminant=FALSE, go.bp=TRUE, go.cc=TRUE, go.mf=TRUE, kegg=TRUE,
	  	annot.method="specificity", annot.details=TRUE, 
	  	direct=FALSE, enriched=TRUE, fdr=NA, build.annot.net=TRUE,
	  	coexp.matrix=NULL, coexp.method="spearman", estimate.th=FALSE, 
	  	hard.th=NA, soft.th=NA, topological = FALSE, keep.sign=FALSE, level=NA, 
	  	annot.clust.method="umilds", annot.prox.measure="unilat.pond.norm.mean",
	  	test.recovery=FALSE, test.robust=FALSE, replace.annot=NA, random.annot=FALSE,
	  	build.gene.net=FALSE, gene.clust.method="hclust", gene.net.details=FALSE,
	  	gene.clusters=NA, alpha=0.05, RV=0.90, sigma=NA, keep.rdata=FALSE, zip=TRUE)
          }
          \arguments{
            \item{wd}{sets the working directory where the expression data files are to be found and 
            		where results are to be stored.}
            \item{org}{indicates the biological species to which analyzable transcript expression data is related; 
            		currently only four possibilities are available with FunNet: "hsa" for human expression data,
            		"mmu" for mouse (Mus Musculus) expression data, "rna" for rat (Rattus Norvegicus) expression data,
            		"sce" for yeast (Saccharomyces Cerevisiae) expression data, and "gga" for chicken (Gallus gallus) data. Default value is "hsa".}
            \item{up.frame}{a dataframe containing expression data for up-regulated genes for the case where two lists 
            		of genes are analyzed comparatively (i.e. up- vs. down-regulated genes). The first column should 
            		contain the GeneID's of the analyzed genes (text data), while the other columns should provide expression 
            		measurements (numeric data) in the analyzed microarray experiments (one column for each array). 
            		The default value is NULL.}
	    \item{down.frame}{a dataframe containing expression data for down-regulated genes for the case where two lists 
            		of genes are analyzed comparatively (i.e. up- vs. down-regulated genes). The first column should 
            		contain the GeneID's of the analyzed genes (text data), while the other columns should provide expression 
            		measurements (numeric data) in the analyzed microarray experiments (one column for each array). 
            		The default value is NULL.}
	    \item{genes.frame}{a dataframe containing expression data for analyzed genes for the case where one list
            		of genes is analyzed. The first column should contain the GeneID's of the analyzed genes (text data), 
            		while the other columns should provide expression measurements (numeric data) in the analyzed 
            		microarray experiments (one column for each array). The default value is NULL.}
	    \item{ref.list}{a dataframe containing a single column providing GeneID's (text data) to be used as the reference list 
	    		of analyzed genes (i.e. all genes spotted on the microarrays used for expression profiling). The default 
	    		value is NULL.}
            \item{two.lists}{possible values are TRUE if a discriminatory functional analysis of two lists of transcripts
            		is required (e.g. significantly up-regulated transcripts versus down-regulated transcripts) or 
            		FALSE if only one list of transcripts is to be analyzed. Please see the provided datasets for the required 
            		format of the data files. The default value of this parameter is TRUE.}
            \item{restrict}{possible values are TRUE if a reference list of transcripts is provided for the statistical
            		significance calculation of the transcript enrichment of the biological annotations or FALSE if 
            		such a restriction is not imposed and the transcript enrichment significance is therefore estimated 
            		with regards of the whole genome. The purpose of the reference list is to correct the significance of 
            		the enrichment calculations for those situations in which expression data is not available for the 
            		whole genome but only for a fraction of it, either because of microarray processing errors which 
            		limits the number of transcripts available for analysis, or for the case of dedicated microarrays, 
            		designed to scan only a fraction of the genome. The transcripts should be identified only by their 
            		EntrezGene ID number. The default value for this parameter is FALSE.}
            \item{logged}{possible values are TRUE or FALSE. This parameter indicates whether the expression measurements 
            		are provided as log fold changes or as simple fold changes. This parameter is used exclusively in 
            		the filtering and pretreatment of the provided expression data.}
            \item{discriminant}{possible values are TRUE or FALSE. This parameter indicate whether a discriminant annotation 
            		of two lists of genes should be performed (if TRUE) or alternatively only an independent annotation of 
            		each list of genes (if FALSE). In the case of a discriminant analysis the enrichment of each list is 
            		computed with regards to their union, while otherwise the provided reference list is used for enrichment 
            		computation. The default value for this parameter is FALSE.}
            \item{go.bp}{possible values are TRUE or FALSE. This parameter indicates whether GO Biological Process annotations
            		should be used. Default value is TRUE.}
            \item{go.cc}{possible values are TRUE or FALSE. This parameter indicates whether GO Cellular Component annotations
            		should be used. Default value is TRUE.}
            \item{go.mf}{possible values are TRUE or FALSE. This parameter indicates whether GO Molecular Function annotations
            		should be used. Default value is TRUE.}
            \item{kegg}{possible values are TRUE or FALSE. This parameter indicates whether KEGG annotations should be used.
            		Default value is TRUE.}
            \item{annot.method}{indicates the type of ontological annotation to be used (for GO annotations only). Three values 
            		are possible: "specificity", "terminological", "decorrelated". "Specificity" will direct the enrichment 
            		computation with regards to the precision of the annotation. In this situation "specificity" levels 
            		are identified in reference to the most precise annotations available (i.e. direct annotations) which 
            		constitute the first specificity level. Further levels of annotation precision are identified by 
            		following annotation inheritance in the GO lattice. "Terminological" will indicate that the enrichment 
            		computations should be performed in reference to the terminological levels of GO (i.e. levels of 
            		conceptual precision within the ontology). "Decorrelated" will use a specific decorrelating technique 
            		to reduce informational redundancy within GO. A single list of enriched categories is provided in this 
            		last case, while in the previous two a hierarchy of lists, one for each level of annotation specificity 
            		or conceptual precision, is provided. The default value for this parameter is "specificity".}
            \item{annot.details}{indicates whether detailed results reflecting enrichment computation should be stored as HTML 
            		files. The default value for this parameter is TRUE.}
            \item{direct}{indicates whether the terminological annotation (see \code{annot.method} parameter) should consider the 
            		precision of annotating categories or not. If FALSE (default) a conventional terminological annotation 
            		is performed. Otherwise the enrichment computation considers for each ontological level only  
            		categories with a similar annotation precision.}
            \item{enriched}{if FALSE it allows to use all available annotations regardless of their enrichment significance. 
            		Provided only for experimental purposes. The default value for this parameter is TRUE.}
            \item{fdr}{if numeric it uses the Storey FDR approach for enrichment p-values correction to a False Discovery Rate equal 
            		to the FDR in percentages. The default value for this parameter is NA (no FDR correction is performed in this case.}
            \item{build.annot.net}{indicates whether a functional interaction network should be computed. The default value for 
            		this parameter is TRUE.}
            \item{coexp.matrix}{allows to specify another co-expression matrix than the one computed by default. The default 
	    		value is NULL.}
            \item{coexp.method}{indicates the co-expression measure to be used to compare gene expression profiles. Four values 
            		are possible: "spearman", "pearson", "kendall", "euclid". The first three are conventional correlation 
            		coefficients while the last one indicates an Euclidean distance between expression profiles.}
            \item{estimate.th}{indicates whether an estimation of the co-expression threshold based on a theoretical scale-free 
            		distribution (power law) should be performed. Both a "hard" (i.e. the adjustment quality for various discrete 
            		values of the used co-expression measure) and a "soft" (i.e. the adjustment quality for various power law 
            		exponents) threshold are computed (for details please see Zhang & Horvath's paper indicated as reference). 
            		Based on these calculations the values of the hard or soft threshold must be manually chosen and indicated 
            		in the respective parameters. The default value for this parameter is FALSE.}
            \item{hard.th}{indicates the value of the discrete hard threshold for co-expression significance. Either the hard or the soft
            		threshold must be specified except when their estimation is requested. The default value is NA.}
            \item{soft.th}{indicates the value of the discrete soft threshold for co-expression significance. Either the hard or the soft
            		threshold must be specified except when their estimation is requested. The default value is NA.}
            \item{topological}{indicates whether a topological measure of similarity between gene expression profiles should be computed 
            		based on the co-expression network build with the co-expression measure indicated by the \code{coexp.method} 
            		parameter. The default value for this parameter is FALSE.}
            \item{keep.sign}{indicates whether the sign of the correlation coeficient should be considered when computing gene co-expression. 
            		The default value for this parameter is FALSE.}
            \item{level}{indicates the level of terminological or specificity annotation to be used for building functional interaction 
            		networks. The default value for this parameter is NA which implies that the most specific (i.e. the first) level
            		will be used. It has no impact on KEGG and "decorrelated" GO annotation.}
            \item{annot.clust.method}{indicates the multiple-instance clustering algorithm to be used to explore the structure of 
            		functional interaction modules. The two possibilities "umilds" and "ucknn" are detailed in the respective 
            		references provided below. The default value is "umilds" which performs a dynamical estimation of the proximity 
            		between annotating categories based on the co-expression of their annotated transcripts and then uses a spectral 
            		clustering technique to explore the modular structure of the network.}
            \item{annot.prox.measure}{indicates the type of measures to be used for computing the proximity between annotating categories. 
            		It has five values: "dynamical", "unilat.pond.norm.mean", "unilat.norm.sum", "norm.sum" and "pond.norm.mean". The 
            		first one is to be used with the "umilds" algorithm (see \code{annot.clust.method} parameter), while the others are 
            		symetrical or asymetrical proximity measures to be used with the "ucknn" algorithm. By default "umilds" algorithm 
            		uses the dynamical proximity estimation so there is no need to indicate it except for experimental purposes.}
            \item{test.recovery}{this parameter is purely experimental and therefore it should not be modified. The default value is 
            		FALSE.}
            \item{test.robust}{this parameter is purely experimental and therefore it should not be modified. The default value is 
            		FALSE.}
            \item{replace.annot}{this parameter is purely experimental and therefore it should not be modified. The default value is 
            		NA.}
            \item{random.annot}{this parameter is purely experimental and therefore it should not be modified. The default value is 
            		FALSE.}
            \item{build.gene.net}{indicates whether a simple gene co-expression network should be computed. If TRUE the results
            		of this computation does not take into account functional assignment of analyzed genes. The default value 
            		is FALSE.}
            \item{gene.clust.method}{this parameter is experimental at this stage and should not be modified.}
            \item{gene.net.details}{indicates whether detailed results containing network centrality information in relation to 
            		the analyzed transcripts should be provided. The default value is FALSE.}
            \item{gene.clusters}{allows to specify a predefined number of gene clusters to be considered in the analysis of 
            		the gene co-expression network. By default the number of clusters (i.e. modules) is computed by optimizing 
            		the Silhouette of the resulting gene clustering partitions.}
            \item{alpha}{indicates the threshold of p-values significance (alpha) resulting from statistical calculations concerning
            		transcript enrichment of biological annotations. Default value is 0.05.}
            \item{RV}{allows to control the strength of the co-inertia level used to evaluate the convergence of the dynamical model 
            		on which relies the dynamical computation of proximity between annotating categories. It is used in 
            		combination with a co-inertia analysis and a Mantel test (see \code{ade4} package). The default value is 0.9.}
            \item{sigma}{this parameter is experimental at this stage and therefore it should not be modified.}
            \item{keep.rdata}{this parameter allows to remove the temporary RData files in the end of the computations. The default value is FALSE}
            \item{zip}{this parameter allows to create a ZIP archive for storing the results in the end of the computations. The default value is TRUE}
          }
          \details{FunNet can be used with the currently available R distributions (tested with distributions posterior to 2.5.0),
          		either with Microsoft Windows operating environments (tested with Windows XP), Mac OS (tested with OS X Leopard) 
          		or, better, with a Linux operating environment (tested with Debian 4.0 and OpenSuse 10.3). Please be aware that 
          		FunNet analysis implies a lot of computations and therefore high processing power and good stability of 
          		the operating system are absolute requirements.

          Together with the FunNet algorithm this package provide also:\cr
	    1. GO and KEGG annotations (as of February 2008) automatically extracted from their respective web resources\cr\cr
	    2. The routine for the automated extraction and update of the functional annotations from their respective
			web resources. The use of this routine is simple: \code{annotations(date.annot = "")}. Under common
			circumstances these routine will provide up-to-date annotations, stored into environmental variables,
			directly formatted for FunNet's use. Some errors may be seen when using this routine related to a
			lack of availability of the GO annotations for the current month. In case of extraction errors, 
			explained most usually by a delay in updating GO web servers, the release date can be expressly 
			indicated (see \code{\link{annotations}}).\cr\cr
	    3. Four test data sets (see examples below and the dedicated man pages). Two of these datasets are related to 
			adipose tissue expression profiling in obese subjects at baseline and after a bariatric surgery. 
			The other two are yeast datasets related to the cell cycle and DNA repairing processes induced by
			irradiation.\cr

	    The format of the data should be respected in order to perform a successful analysis. The only transcript 
			identification system acceptable for FunNet analysis is EntrezGene GeneID's. The transcript 
			expression data should be organized in dataframes within one row for each transcript. The 
			first column contains the transcript identifiers for each transcript and the rest of them 
			the expression level of that transcript in each of the available microarray samples. 
			See the provided test data for more details.\cr

	    The results of the FunNet analysis of transcript expression data are stored as HTML, tab separated text or R 
			data files in a "Results" subfolder of the working folder. For each type of available biological
			annotations and for each list of transcript expression data to be analyzed (one or two), FunNet 
			provides a ranked list with the significantly enriched annotating categories, as well as network 
			structures as text files designed to be imported in Cytoscape for graphical analysis. 
			Detailed findings on the terminological composition and transcript enrichment significance of the 
			resulting functional clusters, as well as various network centrality measures are equally provided.}
          \references{
	          1. Prifti E, Zucker JD, Clement K, Henegar C. FunNet: an integrative tool for exploring transcriptional interactions. 
	  			Bioinformatics. 2008 Nov 15;24(22):2636-8.
	            
	          2. Henegar C, Tordjman J, Achard V, Lacasa D, Cremer I, Guerre-Millo M, Poitou C, Basdevant A,		
	  			Stich V, Viguerie N, Langin D, Bedossa P, Zucker J-D, Clement K. Adipose tissue transcriptomic 
	  			signature highlights the pathologic relevance of extracellular matrix in human obesity. 
	  			Genome Biology 2008, 9(1):R14.
	  
	  	  3. Henegar C, Clement K, and Zucker JD (2006). Unsupervised multiple-instance learning for functional profiling 
	  			of genomic data. Lecture Notes in Computer Science: ECML 2006. 
	  			Springer Berlin / Heidelberg, 4212/2006 : 186-197.
	  
	  	  4. Henegar C, Cancello R, Rome S, Vidal H, Clement K, Zucker JD. Clustering biological annotations and gene 
	  			expression data to identify putatively co-regulated biological processes. J Bioinform Comput Biol.
	  			2006 Aug;4(4):833-52.
	  
	  	  5. Cancello R, Henegar C, Viguerie N, Taleb S, Poitou C, Rouault C, Coupaye M, Pelloux V, Hugol D, Bouillot 
	  			JL, Bouloumie A, Barbatelli G, Cinti S, Svensson PA, Barsh GS, Zucker JD, Basdevant A, Langin D,
	  			Clement K. Reduction of macrophage infiltration and chemoattractant gene expression changes in 
	  			white adipose tissue of morbidly obese subjects after surgery-induced weight loss. 
	  			Diabetes 2005; 54(8):2277-86.
	  
	  	  6. Zhang B, Horvath S. A general framework for weighted gene co-expression network analysis. Stat Appl 
	  			Genet Mol Biol 4 (2005) Article17.
	  
	  	  7. FunNet websites: \url{http://corneliu.henegar.info/FunNet.htm}, \url{http://www.funnet.ws},
	  			\url{http://www.funnet.info}
          }
          \seealso{
            \code{\link{cluster}, \link{annotations}, link{FunCluster}}.
          }
          \concept{gene co-expression interaction networks}
          \concept{functional interactions from gene expression data}
          \concept{gene expression functional profiling}
          \concept{functional profiling of cDNA microarray data}
          \concept{clustering functional annotations and gene expression data}
          \concept{identifying co-regulated biological processes from gene expression data and functional annotations}
          \concept{Gene Ontology functional profiling of microarray data}
          \concept{KEGG functional profiling of microarray data}
          \examples{
          \dontrun{
          ## load bypass data (see Diabetes and Genome Biology papers for details)
          data(bypass)

          ## or load adipose tissue expression profiling data (see Genome Biology paper for details)
          data(obese)

          ## most common use
          data(obese)
          FunNet(org="HS", two.lists=TRUE, up.frame=up.frame, down.frame=down.frame,
	  	genes.frame=NULL, restrict=TRUE, ref.list=ref.list, logged=TRUE,
	  	discriminant=TRUE, go.bp=TRUE, go.cc=TRUE, go.mf=TRUE, kegg=TRUE,
	  	annot.method="specificity", annot.details=TRUE, 
	  	direct=FALSE, enriched=TRUE, fdr=NA, build.annot.net=TRUE,
	  	coexp.matrix=NULL, coexp.method="spearman", estimate.th=FALSE, 
	  	hard.th=0.8, soft.th=NA, topological = FALSE, keep.sign=FALSE, level=1, 
	  	annot.clust.method="umilds", annot.prox.measure="dynamical",
	  	test.recovery=FALSE, test.robust=FALSE, replace.annot=NA, random.annot=FALSE,
	  	build.gene.net=TRUE, gene.clust.method="hclust", gene.net.details=TRUE,
	  	gene.clusters=NA, alpha=0.05, RV=0.90, sigma=NA, keep.rdata=FALSE, zip=TRUE)                   
          }}
          \keyword{cluster}
