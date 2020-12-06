#' Produces a fasta file that contains the desired Uniprot IDs' sequences.
#'
#' A function that produces a fasta file based on the Uniprot ID entered
#' by the users.
#'
#' @param ID A character string that includes the Uniprot IDs separated by "OR"
#'
#' @return Creates a fasta file called seqs.fasta in the main directory.
#'
#' @examples
#' # Use the following UniprotID for this function: P19838, Q00653, Q01201
#' \dontrun{
#' uniprotToFasta("AC=P19838 OR AC=Q00653 OR AC=Q01201")
#' }
#'
#' # A file called seqs.fasta will appear in the main directory.
#'
#' @references
#'   1. Brennan, P. (1970, January 01). Drawing a simple phylogenetic tree of the human rel homology domain family. Retrieved December 04, 2020, from http://rforbiochemists.blogspot.com/2016/12/drawing-simple-phylogenetic-tree-of.html
#'   2. Code for retrieving Uniprot data were borrowed from: http://rforbiochemists.blogspot.com/2016/12/drawing-simple-phylogenetic-tree-of.html
#'   3. Pages, H., Aboyoun, P., Gentleman, R., & DebRoy, S. (2016). Biostrings: String objects representing biological sequences, and matching algorithms. R package version, 2(0), 10-18129.
#'   4. Paradis, E., Claude, J., & Strimmer, K. (2004). APE: analyses of phylogenetics and evolution in R language. Bioinformatics, 20(2), 289-290.
#'   5. Yu, G., Smith, D. K., Zhu, H., Guan, Y., & Lam, T. T. Y. (2017). ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods in Ecology and Evolution, 8(1), 28-36.
#'   6. Charif, D., Lobry, J. R., Necsulea, A., Palmeira, L., Penel, S., Perriere, G., & Penel, M. S. (2020). Package ‘seqinr’.
#'
#' @export
#'
#' @import seqinr
#' @import stringr
#'
uniprotToFasta <- function(ID){

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
}
