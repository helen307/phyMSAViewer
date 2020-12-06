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
#' \dontrun{
#' # Use the following UniprotID for this
#' # function: P19838, Q00653, Q01201
#' uniprotToPhy("AC=P19838 OR AC=Q00653 OR AC=Q01201")
#' # Access the phylogenetic tree with MSA plot directly.
#' }
#'
#'
#'@references
#'   1. Code for retrieving Uniprot data were borrowed from: http://rforbiochemists.blogspot.com/2016/12/drawing-simple-phylogenetic-tree-of.html
#'   2. Pages, H., Aboyoun, P., Gentleman, R., & DebRoy, S. (2016). Biostrings: String objects representing biological sequences, and matching algorithms. R package version, 2(0), 10-18129.
#'   3. Paradis, E., Claude, J., & Strimmer, K. (2004). APE: analyses of phylogenetics and evolution in R language. Bioinformatics, 20(2), 289-290.
#'   4. Yu, G., Smith, D. K., Zhu, H., Guan, Y., & Lam, T. T. Y. (2017). ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods in Ecology and Evolution, 8(1), 28-36.
#'   5. Charif, D., Lobry, J. R., Necsulea, A., Palmeira, L., Penel, S., Perriere, G., & Penel, M. S. (2020). Package ‘seqinr’.
#'
#' @export
#'
#' @importFrom Biostrings readAAStringSet
#' @importFrom ape nj
#' @import msa
#' @import ggtree
#' @import seqinr
#' @import stringr
#'
#'
uniprotToPhy <- function(ID){

  uniID <- ID

  if (!grepl("[A-Z][0-9]{5}", uniID)) {
    print("Please provide a valid string of Uniprot ID")
  } else if ((nchar(uniID) == 6 || nchar(uniID) == 9) &&
           stringr::str_count(uniID, "OR") == 0) {
    # ensure the IDs are separated by 'OR's
    print("Please use ORs to separate multiple Uniprot IDs")
  } else if ((stringr::str_count(uniID, " OR ") + 1) != stringr::str_count(uniID, "AC=")){
    # ensure each ID is led by a 'AC='
    print("Please use AC= to lead all the Uniprot IDs")
  } else {
    print("Valid input!")
  }

  mybank <- seqinr::choosebank(bank = "swissprot")
  seq1 <- seqinr::query("relSeq", uniID)

  # get sequence
  seq2 <- seqinr::getSequence(seq1)
  seqinr::write.fasta(sequences = seq2,
              names = getName(seq1),
              nbchar = 80, file.out = "seqs.fasta")

  # read sequence from the fasta file
  mySeqs <- Biostrings::readAAStringSet("seqs.fasta")   # from package Biostrings

  # perform multiple sequence alignment
  to_align <- msa::msa(mySeqs, method="Muscle")

  # Build tree

  my_align <- msa::msaConvert(to_align,
                              type="seqinr::alignment")
  seqinr::write.fasta(as.list(my_align$seq),
              my_align$nam,
              file.out="msa.fasta")
  # pair-wise
  d <- dist.alignment(my_align, "identity")
  # neighbor-joining
  myTree <- ape::nj(d)

  # final plot
  ggtree::msaplot(p=ggtree(myTree) + geom_tiplab(hjust=1,vjust=-1),
                  fasta="msa.fasta")
}
