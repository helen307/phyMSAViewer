#' Produces the phylogenetic tree graph from sequence file
#'
#' A function that produces a phylogenetic tree based on the fasta file uploaded
#' by the users.
#'
#' @param mySeqs A AAStringSet of sequences, length of each sequence varies.
#'
#' @return Returns a graph of a phylogenetic tree calculated on the sequence dataset
#'         along with its full MSA plot.
#' @examples
#' \dontrun{
#' # Use the mySeqs dataset for this function:
#' seqToPhy(mySeqs)
#' # Access the phylogenetic tree with MSA plot directly.
#' }
#'
#'
#' @references
#'   1. Code for retrieving Uniprot data were borrowed from: http://rforbiochemists.blogspot.com/2016/12/drawing-simple-phylogenetic-tree-of.html
#'   2. Pages, H., Aboyoun, P., Gentleman, R., & DebRoy, S. (2016). Biostrings: String objects representing biological sequences, and matching algorithms. R package version, 2(0), 10-18129.
#'   3. Paradis, E., Claude, J., & Strimmer, K. (2004). APE: analyses of phylogenetics and evolution in R language. Bioinformatics, 20(2), 289-290.
#'   4. Yu, G., Smith, D. K., Zhu, H., Guan, Y., & Lam, T. T. Y. (2017). ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods in Ecology and Evolution, 8(1), 28-36.
#'
#' @export
#'
#' @importFrom Biostrings readAAStringSet
#' @importFrom ape nj
#' @import msa
#' @import ggtree
#'

seqToPhy <- function(mySeqs){

  # perform msa with MUSCLE
  to_align <- msa::msa(mySeqs)

  # Build alignment
  my_align <- msa::msaConvert(to_align, type="seqinr::alignment")
  write.fasta(as.list(my_align$seq),
              my_align$nam,
              file.out="msa.fasta")

  # pair-wise distance
  dis <- dist.alignment(my_align, "identity")

  # neighbor-joining algo
  myTree <- ape::nj(dis)

  # final plot: phy + msa
  ggtree::msaplot(p=ggtree(myTree) + geom_tiplab(hjust=1,vjust=-1),
                  fasta="msa.fasta")
}

