toeplitz <- function(n, use_sparse = F) {
  if (n < 0) stop("n < 1");
  
  if (n == 1) return(matrix(data = c(2)));
  
  d <- diag(2, nrow = n, ncol = n);
  l <- diag(-1, n - 1, n);
  u <- diag(-1, n, n - 1);
  m <- matrix(data = 0, ncol = n, nrow = n);
  m[1:n, 2:n] <- u;
  m[2:n, 1:n] <- m[2:n, 1:n] + l;
  m <- m + d;
  if (use_sparse) m <- Matrix(m);
  return(m);
}

#'@export
matrix_K <- function(n, use_sparse) toeplitz(n, use_sparse);

#'@export
matrix_C <- function(n, use_sparse) {
  m_C <- toeplitz(n, use_sparse);
  m_C[1, n] <- -1;
  m_C[n, 1] <- -1;
  return(m_C);
}

#'@export
matrix_T <- function(n, use_sparse) {
  m_T <- toeplitz(n);
  m_T[1, 1] <- 1;
  
  return(m_T);
}

#'@export
matrix_B <- function(n, use_sparse) {
  m_B <- toeplitz(n, use_sparse);
  m_B[1, 1] <- 1;
  m_B[n, n] <- 1;
  return(m_B);
}


#'@export
KTCB <- function(n, use_sparse) {
  if (n < 1) stop("n < 1");
  
  m_K <- matrix_K(n, use_sparse);
  m_T <- matrix_T(n, use_sparse);
  m_C <- matrix_C(n, use_sparse);
  m_B <- matrix_B(n, use_sparse);
  
  return(
    list(
      "K" = m_K,
      "T" = m_T,
      "C" = m_C,
      "B" = m_B
    )
  );
}
