#' @useDynLib BinaryDosageCompress
#' @importFrom Rcpp sourceCpp
#' @exportPattern "^[[:alpha:]]+"
NULL



#' Get all paths to the .vcf and .vcf.gz files in the given directories
#'
#' This function reads an array of directories and returns an array of all paths to the .vcf and .vcf.gz files in the directories
#'
#' @param dirs an array of directories to search the .vcf and .vcf.gz files
#' @return an array of all paths to the .vcf and .vcf.gz files in the directories
#' @export
findVcfFiles <- function(dirs) {
  all_files <- c()

  # Loop through each directory
  for (dir in dirs) {
    # Check if the directory exists
    if (!dir.exists(dir)) {
      warning(paste("Directory does not exist:", dir))
      next
    }

    # List .vcf and .gz files in the directory
    files_in_dir <- list.files(path = dir, pattern = "\\.vcf$|\\.vcf.gz$", full.names = TRUE)
    all_files <- c(all_files, files_in_dir)
  }

  return(all_files)
}

#' Read SNP with position array and index
#'
#' This function reads a dataframe of sample ID, DS and GP, given an array of positions for each SNP in the bd file and the index for the SNP in the file
#'
#' @param file a bd file to read the SNP
#' @param positions an array of positions for each SNP in the bd file
#' @param index the index for the SNP in the file
#' @return a dataframe of sample ID, DS and GP
#' @export
readSNPwithPositionVector <- function(filename, positions, index) {
  return(readSNPwithPosition(filename, positions[index]))
}

#' Read SNP's Dosage with position array and index
#'
#' This function reads an array of Dosage data, given an array of positions for each SNP in the bd file and the index for the SNP in the file
#'
#' @param file a bd file to read the SNP
#' @param positions an array of positions for each SNP in the bd file
#' @param index the index for the SNP in the file
#' @return an array of Dosage data
#' @export
readSNPonlyDSwithPositionVector <- function(filename, positions, index) {
  return(readSNPonlyDSwithPosition(filename, positions[index]))
}
