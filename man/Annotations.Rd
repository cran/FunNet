  \name{Annotations}
          \alias{annotations}
          \title{Integrative Functional Analysis of Transcriptional Networks}
          \description{
            This routine is close to a similar one belonging to the package FunCluster. The actual one, provided with the package FunNet, 
            performs a slightly different automated extraction and update of the Gene Ontology & KEGG annotations which are needed for 
            FunNet analysis. The difference relies in the fact that additionally to the gene annotation data, extracted from the NCBI ressources,
            this routine provides also the ontological lattice of GO required by FunNet enrichment computation routine for specificity, 
            terminological or decorrelated annotation.
          }
          \usage{
          annotations(cust.specs=NULL)
          }
          \arguments{
            \item{cust.specs}{this parameter allows to specify a custom list of organisms for extracting GO and KEGG annotations 
			from respective databases. It is experimental at this moment and therefore no specific documentation is 
			provided at this time about its use (although an example of usage is available within the "Annotations.R" 
			script included in the package).}
          }
          \details{For details concerning the utilization of the annotations data by the FunNet package please see FunNet help 
          		or man page \code{help(FunNet)}.
          		The "Annotations" routine is allowing the automated extraction and update of the functional
          		annotations from their respective web resources. Under common circumstances this routine will
          		provide up-to-date annotations, stored into environmental variables and directly formatted for
          		FunNet use. Some errors may be seen when using this routine to update annotations provided within the package
          		in relation to frequent changes in the format of annotation files provided by the NCBI. Please 
          		update your package with the last version available as this should correct such issues in most cases.\cr
          		
          		The annotation data is provided as a unique R data archive which should be manualy loaded into R 
          		after loading the FunNet package in order to be able to use the updated annotations for your analysis.\cr
          		
            		Important note for Microsoft Windows users: in order to use this routine you will need additional
            		software for handling TAR and GZIP archives. This software is freely available for Windows under the 
            		GNU license.\cr
            		
            		For TAR packages please go to: \cr
            		\url{http://gnuwin32.sourceforge.net/packages/tar.htm}. 
            		For GZIP you can go to: \cr
            		\url{http://gnuwin32.sourceforge.net/packages/gzip.htm}.
            		The TAR and GZIP executable and their dependencies (DLL's) should be placed somewhere into the PATH 
            		(like "C:/Windows" for example) in order to be available for R calls.
            }
          \references{
      1. Prifti E, Zucker JD, Clement K, Henegar C. Interactional and functional centrality in transcriptional co-expression 
      		networks. Bioinformatics. 2010 Oct 19. [Epub ahead of print]
       
      2. Prifti E, Zucker JD, Clement K, Henegar C. FunNet: an integrative tool for exploring transcriptional interactions. 
			Bioinformatics. 2008 Nov 15;24(22):2636-8.
          
      3. Henegar C, Tordjman J, Achard V, Lacasa D, Cremer I, Guerre-Millo M, Poitou C, Basdevant A,			Stich V, Viguerie N, Langin D, Bedossa P, Zucker J-D, Clement K. Adipose tissue transcriptomic 
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

          \seealso{
            \code{\link{FunNet}}.
          }
          \examples{
          \dontrun{
          annotations()
          }}
          \keyword{cluster}
