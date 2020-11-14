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

  # perform msa with MUSCLE
  to_align <- msa::msa(mySeqs,method="Muscle")

  # Build alignment
  my_align <- msa::msaConvert(to_align, type="seqinr::alignment")

  # write into fasta
  write.fasta(as.list(my_align$seq),
              my_align$nam,
              file.out="msa.fasta")

  # pair-wise distance
  dis <- dist.alignment(my_align, "identity")

  # neighbor-joining algo
  myTree <- ape::nj(dis)

  # final plot: phy + msa
  ggtree::msaplot(p=ggtree(myTree),
                  fasta="msa.fasta",
                  window=c(50, 200))
}
