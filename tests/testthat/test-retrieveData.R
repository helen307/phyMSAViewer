context("uniprotToPhy function")
library(phyMSAViewer)
library(seqinr)
library(testthat)
library(stringr)


test_that("seq1 works", {
  user_input <- "AC=P19838 OR AC=Q00653 OR AC=Q01201"
  cnt <- stringr::str_count(uniID, "OR") + 1

  mybank <- seqinr::choosebank(bank = "swissprot")
  seq1 <- seqinr::query("relSeq", user_input)

  expect_equal(as.integer(as.character(seq1)[3]), cnt)
})

test_that("seq2 works", {
  user_input <- "AC=P19838 OR AC=Q00653 OR AC=Q01201"
  cnt <- 3

  mybank <- seqinr::choosebank(bank = "swissprot")
  seq1 <- seqinr::query("relSeq", user_input)
  seq2 <- seqinr::getSequence(seq1)

  expect_equal(length(seq2), cnt)
})


test_that("read works", {
  user_input <- "AC=P19838 OR AC=Q00653 OR AC=Q01201"
  cnt <- 3

  mybank <- seqinr::choosebank(bank = "swissprot")
  seq1 <- seqinr::query("relSeq", user_input)
  seq2 <- seqinr::getSequence(seq1)
  seqinr::write.fasta(sequences = seq2,
                      names = getName(seq1),
                      nbchar = 80, file.out = "seqs.fasta")
  mySeqs <- Biostrings::readAAStringSet("seqs.fasta")   # from package Biostrings



  tot_seq2 <- length(seq2[[1]]) + length(seq2[[2]]) + length(seq2[[3]])
  tot_mySeqs <- length(mySeqs[[1]]) + length(mySeqs[[2]]) + length(mySeqs[[3]])

  expect_equal(tot_seq2, tot_mySeqs)
})
