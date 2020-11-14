#' Produces the phylogenetic tree graph
#'
#' A function that produces a phylogenetic tree based on the Uniprot ID entered
#' by the users.
#'
#' @param ID A character string that includes the Uniprot IDs separated by "OR"
#'
#' @return Returns a graph of a phylogenetic tree calculated on the Uniprot IDs
#'         entered along with its full MSA plot.
#' @examples
#' # Use the following UniprotID for this function: P19838, Q00653, Q01201
#' uniprotToPhy <- function(ID="AC=P19838 OR AC=Q00653 OR AC=Q01201")
#' # Access the phylogenetic tree with MSA plot directly.
#' uniprotToPhy
#'
#' @export
#'
#' @importFrom Biostrings readAAStringSet
#' @importFrom ape nj
#' @import msa
#' @import ggtree
#' @import seqinr
#'

uniprotToPhy <- function(ID){

  uniID <- ID
  if(!is.character(uniID) | length(uniID) < 1){
    return("Please provide a valid string of Uniprot ID")
  }

  mybank <- seqinr::choosebank(bank = "swissprot")
  seq1 <- seqinr::query("relSeq", uniID)
  # N.B. protein info NOT returned in the same order as requested


  # get sequence
  seq2 <- seqinr::getSequence(seq1)

  # write into fasta file
  write.fasta(sequences = seq2,
              names = getName(seq1),
              nbchar = 80, file.out = "seqs.fasta")

  # read sequence from the fasta file
  mySeqs <- Biostrings::readAAStringSet("seqs")   # from package Biostrings

  # perform multiple sequence alignment
  to_align <- msa::msa(mySeqs, method="Muscle")

  # Build tree
  my_align <- msa::msaConvert(to_align,
                              type="seqinr::alignment")

  # write into fasta
  write.fasta(as.list(my_align$seq),
              my_align$nam,
              file.out="msa.fasta")

  d <- dist.alignment(my_align, "identity")
  myTree <- ape::nj(d) # neighbor-joining

  ggtree::msaplot(p=ggtree(myTree),
                  fasta="msa.fasta",
                  window=c(50, 200))
}
