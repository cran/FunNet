\name{FunNet.R-package}
\alias{FunNet.R-package}
\docType{package}
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
\details{
          Together with the FunNet algorithm this package provides also:\cr
	    1. GO and KEGG annotations automatically extracted from their respective web resources and updated on a regular
			basis\cr\cr
	    2. The routine for the automated extraction and update of the functional annotations from their respective
			web resources. The use of this routine is simple: \code{annotations()}. Under common
			circumstances these routine will provide up-to-date annotations, stored into environmental variables,
			directly formatted for FunNet's use.\cr\cr
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
			resulting functional clusters, as well as various network centrality measures are equally provided.
}
\author{
          Corneliu Henegar \email{corneliu@henegar.info}
}
\references{
      1. Prifti E, Zucker JD, Clement K, Henegar C. Interactional and functional centrality in transcriptional co-expression 			networks. Bioinformatics. 2010 Oct 19. [Epub ahead of print].
       
      2. Prifti E, Zucker JD, Clement K, Henegar C. FunNet: an integrative tool for exploring transcriptional interactions. 
			Bioinformatics. 2008 Nov 15;24(22):2636-8.
          
      3. Henegar C, Tordjman J, Achard V, Lacasa D, Cremer I, Guerre-Millo M, Poitou C, Basdevant A,
      			Stich V, Viguerie N, Langin D, Bedossa P, Zucker J-D, Clement K. Adipose tissue transcriptomic 
			signature highlights the pathologic relevance of extracellular matrix in human obesity. 
			Genome Biology 2008, 9(1):R14.

	4. Henegar C, Clement K, and Zucker JD (2006). Unsupervised multiple-instance learning for functional profiling 
			of genomic data. Lecture Notes in Computer Science: ECML 2006. 
			Springer Berlin / Heidelberg, 4212/2006 : 186-197.

	5. Henegar C, Cancello R, Rome S, Vidal H, Clement K, Zucker JD. Clustering biological annotations and gene 
			expression data to identify putatively co-regulated biological processes. J Bioinform Comput Biol.
			2006 Aug;4(4):833-52.

	6. Cancello R, Henegar C, Viguerie N, Taleb S, Poitou C, Rouault C, Coupaye M, Pelloux V, Hugol D, Bouillot 
			JL, Bouloumie A, Barbatelli G, Cinti S, Svensson PA, Barsh GS, Zucker JD, Basdevant A, Langin D,
			Clement K. Reduction of macrophage infiltration and chemoattractant gene expression changes in 
			white adipose tissue of morbidly obese subjects after surgery-induced weight loss. 
			Diabetes 2005; 54(8):2277-86.

	7. Zhang B, Horvath S. A general framework for weighted gene co-expression network analysis. Stat Appl 
			Genet Mol Biol 4 (2005) Article17.

	8. FunNet websites: \url{http://corneliu.henegar.info/FunNet.htm}, \url{http://www.funnet.ws},
			\url{http://www.funnet.info}
}
\keyword{ cluster }

\seealso{
            \code{\link{cluster}}.
}
\examples{
          \dontrun{
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
	  	test.recovery=FALSE, test.robust=FALSE, replace.annot=NA, 
	  	build.gene.net=TRUE, gene.clust.method="hclust", gene.net.details=TRUE,
	  	gene.clusters=NA, alpha=0.05, RV=0.90, sigma=NA, keep.rdata=FALSE, zip=TRUE)                   

          ## updating annotations
          annotations()
          }
}
