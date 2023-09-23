runifm <- function(n, m=NULL, a=0, b=1) {
  # Generate a random uniform matrix
  # If m = NULL, an n x n matrix is generated, else n x m
  if(is.null(m)) return(matrix(runif(n^2, a=a, b=b), nrow = n))
  return(matrix(runif(n*m, a=a, b=b), nrow = n))
}

rpsdm <- function(n) {
  # for any matrix A, A'A is positive semi-definite
  A <- runifm(n)
  return(A%*%t(A))
}

rpdm <- function(n) {
  # Given a positive semo-definite matrix, make all
  # eigen values positive to have a psd.
  A <- rpsdm(n)
  return(A + diag(runif(1), n))
}

is.sqm <- function(mat) {
  return(nrow(mat)==ncol(mat))
}

is.sym <- function(mat) {
  if(!is.sqm(mat)) stop("matrix must be square")
  return(all(mat == t(mat)))
}

is.psd <- function(mat) {
  # Check if all eigen values are greater than 
  # or equal to zero
  if(!is.sym(mat)) stop("matrix must be symmetric")
  return(all(eigen(A)$values >= 0))
}

is.pd <- function(mat) {
  # Check if all eigen values are greater than zero
  if(!is.sym(mat)) stop("matrix must be symmetric")
  return(min(eigen(mat)$values > 0))
}

make_sym <-function(mat) {
  if(!is.sqm(mat)) stop("matrix must be square")
  if(is.sym(mat)) return(mat)
  return((mat + t(mat))/2)
}

make_psd <- function(mat) {
  if(is.psd(mat)) return(mat)
  svdfact <- svd(mat)
  u <- svdfact$u
  v <- svdfact$v
  sigma <- svdfact$d
  pos_sigma <- abs(d) 
  return(u %*% diag(pos_d) %*% t(v))
}

make_pd <- function(mat) {
  if(is.pd(mat)) return(mat)
  if(!is.sym(mat)) {
    mat <- make_sym(mat)
  lambda_min <- min_eig(mat)
  if(lambda_min <= 0) {
    epsilon <- abs(lambda_min) + 1e-6
    mat <- mat + epsilon + diag(1, nrow(mat))
  }
}

min_eig <- function(mat) {
  return(min(eigen(mat)$values))
}
