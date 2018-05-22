
#' Backsolve with band 1 upper Cholesky.
#'
#' Backsolve with band 1 upper Cholesky.
#' @param U An upper triangular square matrix with non-zero entries only on the
#'   main diagonal and the first superdiagonal.
#' @param z A vector with as many elements as the number of rows of U.
#' @return A vector.
#' @importFrom Matrix as
#' @export
band1_backsolve <- function(U, z) {
  if (nrow(U) != ncol(U)) stop("U must be square.")
  if (nrow(U) != length(z)) stop("The dimensions of U and z must match.")
  if (!(class(U) %in% c("dgCMatrix", "dgRMatrix", "dgTMatrix", "dsCMatrix",
                        "dtCMatrix", "nsparseMatrix"))) {
    U <- as(U, "CsparseMatrix")
  }
  return(band1_backsolve_cpp(U, z))
}
