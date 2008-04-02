\name{FunNet.R-package}
\alias{FunNet.R-package}
\docType{package}
\title{Functional Analysis of Gene Co-expression Networks}
\description{
          FunNet is an integrative tool for analyzing gene co-expression networks from 
	  microarray data involving two abstraction levels: genes and biological themes. 
	  A functional profiling technique using Gene Ontology & KEGG annotations is 
	  applied to extract a list of significantly enriched biological themes from 
	  microarray expression profiling data. Afterwards multiple-instance 
	  representations are built to relate significant themes to their 
	  transcriptional instances. A specific non-linear dynamical system model is 
	  used to extract information about the proximity of relevant genomic themes 
	  from the similarity of the expression profiles of their gene instances.
	  Eventually an unsupervised multiple-instance clustering procedure is used 
	  to identify the structure of the co-expression network composed from modules 
	  of functionally related transcripts. A functional and a transcriptional image
	  of the co-expression network are provided separately together with detailed 
	  information on the network centrality of related transcripts and genomic themes.
}
\details{
          Together with the FunNet algorithm this package provides also:\cr
	            1. GO and KEGG annotations automatically extracted from their respective web resources and updated on a regular
	            		basis\cr\cr
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
	            		resulting functional clusters, as well as various network centrality measures are equally provided.
}
\author{
          Corneliu Henegar \email{corneliu@henegar.info}
}
\references{
          1. Henegar C, Tordjman J, Achard V, Lacasa D, Cremer I, Guerre-Millo M, Poitou C, Basdevant A,		
	            	Stich V, Viguerie N, Langin D, Bedossa P, Zucker J-D, Clement K. Adipose tissue transcriptomic 
	            	signature highlights the pathologic relevance of extracellular matrix in human obesity. In positive
	            	revision in Genome Biology journal (October 2007).
	            	
          2. Henegar C, Cancello R, Rome S, Vidal H, Clement K, Zucker JD. Clustering biological annotations and gene 
	            	expression data to identify putatively co-regulated biological processes. J Bioinform Comput Biol.
	            	2006 Aug;4(4) :833-52.
                    
          3. Cancello R, Henegar C, Viguerie N, Taleb S, Poitou C, Rouault C, Coupaye M, Pelloux V, Hugol D, Bouillot JL, 
          		Bouloumie A, Barbatelli G, Cinti S, Svensson PA, Barsh GS, Zucker JD, Basdevant A, Langin D, 
          		Clement K. Reduction of macrophage infiltration and chemoattractant gene expression changes in 
          		white adipose tissue of morbidly obese subjects after surgery-induced weight loss. 
          		Diabetes 2005; 54(8):2277-86.

          4. FunNet website: \url{http://corneliu.henegar.info/FunNet.htm}
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
		discriminant=TRUE, annot.method="specificity", annot.details=TRUE, 
		direct=FALSE, enriched=TRUE, fdr=FALSE, build.annot.net=TRUE,
		coexp.matrix=NULL, coexp.method="spearman", estimate.th=FALSE, 
		hard.th=0.8, soft.th=NA, topological = FALSE, keep.sign=FALSE,
		level=1, annot.clust.method="umilds", annot.prox.measure="dynamical",
		test.recovery=FALSE, test.robust=FALSE, replace.annot=NA, 
		build.gene.net=TRUE, gene.clust.method="hclust", gene.net.details=TRUE,
		gene.clusters=NA, alpha=0.05, RV=0.90, sigma=NA)

          ## the following example forces the use of a previous GO release 
          (e.g. October 2007) for updating annotations
          annotations(date.annot = "200710")
          }
}
