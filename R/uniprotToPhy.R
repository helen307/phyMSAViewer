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
  ### step 4.1 Choose a bank
  seqinr::choosebank()
  # UniProt Knowledgebase Release 2016_08 of 07-Sep-2016 Last Updated: Oct  4, 2016

  ### step 4.2 Make the query
  mybank <- seqinr::choosebank(bank = "swissprot")
  rel_seq <- seqinr::query("relSeq", uniID)
  # N.B. protein info NOT returned in the same order as requested


  # get sequence
  rel_seqs <- seqinr::getSequence(rel_seq)

  # write into fasta file
  write.fasta(sequences = rel_seqs,
              names = getName(rel_seq),
              nbchar = 80, file.out = "relseqs.fasta")

  # read sequence from the fasta file
  mySeqs <- Biostrings::readAAStringSet("relseqs.fasta")   # from package Biostrings

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
