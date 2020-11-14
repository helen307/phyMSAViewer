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
#' # Use the mySeqs dataset for this function:
#' uniprotToPhy <- function(mySeqs)
#' # Access the phylogenetic tree with MSA plot directly.
#' uniprotToPhy
#'
#' @export
#'
#' @importFrom Biostrings readAAStringSet
#' @importFrom ape nj
#' @import msa
#' @import ggtree
#'
seqToPhy <- function(mySeqs){

  # perform multiple sequence alignment
  myAln <- msa::msa(mySeqs)

  # Build tree
  myAln2 <- msa::msaConvert(myAln, type="seqinr::alignment")

  # write into fasta
  write.fasta(as.list(myAln2$seq),myAln2$nam,file.out="msa.fasta")

  d <- dist.alignment(myAln2, "identity")
  myTree <- ape::nj(d) # neighbor-joining

  ggtree::msaplot(p=ggtree(myTree), fasta="msa.fasta", window=c(50, 200))
}
