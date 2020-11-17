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
#'   Code for retrieving Uniprot data were borrowed from: http://rforbiochemists.blogspot.com/2016/12/drawing-simple-phylogenetic-tree-of.html
#'
#' @export
#'
#' @import seqinr
#'
uniprotToFasta <- function(ID){

  uniID <- ID
  if(!is.character(uniID) | length(uniID) < 1){
    return("Please provide a valid string of Uniprot ID")
  }

  mybank <- seqinr::choosebank(bank = "swissprot")
  seq1 <- seqinr::query("relSeq", uniID)

  # get sequence
  seq2 <- seqinr::getSequence(seq1)
  seqinr::write.fasta(sequences = seq2,
                      names = getName(seq1),
                      nbchar = 80, file.out = "seqs.fasta")
}
