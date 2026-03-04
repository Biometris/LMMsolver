#
# calculations of integrals of the form \int B_i(x) B_j(x) dx
#
bernstein <- function(p, xi, k) {
  choose(p, k)*xi^k*(1-xi)^(p-k)
}

bernstein_matrix <- function(p, x) {
  x <- seq(0, 1, length = p+1)
  V <- matrix(0, nrow=p+1, ncol=p+1)
  for (i in seq_len(p+1)) {
    for (j in seq_len(p+1)) {
      V[i,j] <- bernstein(p, x[i],j-1)
    }
  }
  V
}

bernstein_gram_coef <- function(p) {
  G <- matrix(0, nrow = p + 1, ncol = p + 1)
  for (i in 0:p) {
    for (j in 0:p) {
      G[i + 1, j + 1] <-
        choose(p, i) * choose(p, j) /
        ((2*p+1)*choose(2 * p , i + j ))
    }
  }

  G
}

# p: degree of B-splines basis
# q: Number of B-splines
GramMatrix <- function(p, q) {

  ## --- Reference points on one element [0,1]
  x_ref <- seq(0, 1, length = p + 1)

  ## --- Basis evaluations on reference points
  knots <- PsplinesKnots(xmin = 0, xmax = 1, degree = p, nseg = 1)
  B_spline_vals <- as.matrix(Bsplines(knots, x = x_ref))
  B_bernstein_vals <- bernstein_matrix(p, x_ref)

  ## --- Bernstein Gram matrix (analytic)
  G_bernstein <- bernstein_gram_coef(p)

  ## --- CHANGE OF BASIS: B-spline -> Bernstein
  ##     B_spline(x) = Bernstein(x) %*% C
  C_bspline_to_bernstein <- solve(B_bernstein_vals, B_spline_vals)

  ## --- INTEGRATE in Bernstein basis and TRANSFORM BACK
  ##     M_elem = C^T G C
  M_elem <- t(C_bspline_to_bernstein) %*%
    G_bernstein %*%
    C_bspline_to_bernstein

  ## --- Global assembly
  nseg <- q - p
  M <- matrix(0, q, q)
  for (k in 1:nseg) {
    idx <- k:(k + p)
    M[idx, idx] <- M[idx, idx] + M_elem
  }

  M
}

integrate_bspline <- function(a, p, xmin = 0, xmax = 1) {

  q <- length(a)
  nseg <- q - p
  h <- (xmax - xmin) / nseg

  ## reference points for degree-p polynomials
  x_ref <- seq(0, 1, length = p + 1)

  ## B-spline and Bernstein evaluations on [0,1]
  knots <- PsplinesKnots(
    xmin = 0, xmax = 1, degree = p, nseg = 1
  )
  B_spline_vals <- as.matrix(Bsplines(knots, x_ref))
  B_bernstein_vals <- bernstein_matrix(p, x_ref)

  ## change of basis: B-spline -> Bernstein
  C_bspline_to_bernstein <- solve(B_bernstein_vals, B_spline_vals)

  ## integral on one segment
  int <- 0
  for (k in 1:nseg) {
    a_loc <- a[k:(k + p)]
    c_loc <- C_bspline_to_bernstein %*% a_loc
    int <- int + sum(c_loc)
  }

  ## scale
  int * h / (p + 1)
}

# p: degree of B-splines basis
# q: number of B-splines
ortho_int_condition <- function(p, q)
{
  x <- rep(NA,q)
  nseg <- q - p
  dx <- 1/nseg
  for (i in 1:q) {
    b <- rep(0, q)
    b[i] <- 1
    x[i] <- integrate_bspline(b, p=p, xmin = 0, xmax=1)/dx
  }
  x
}

integrate_x_bspline <- function(a, p, xmin = 0, xmax = 1) {

  q <- length(a)
  nseg <- q - p
  h <- (xmax - xmin) / nseg

  ## reference points for degree-p polynomials
  x_ref <- seq(0, 1, length = p + 1)

  ## B-spline and Bernstein evaluations on [0,1]
  knots <- PsplinesKnots(
    xmin = 0, xmax = 1, degree = p, nseg = 1
  )
  B_spline_vals <- as.matrix(Bsplines(knots, x_ref))
  B_bernstein_vals <- bernstein_matrix(p, x_ref)

  ## change of basis: B-spline -> Bernstein
  C_bspline_to_bernstein <- solve(B_bernstein_vals, B_spline_vals)

  ## first moment on one segment
  moment <- 0
  for (k in 1:nseg) {
    a_loc <- a[k:(k + p)]
    c_loc <- C_bspline_to_bernstein %*% a_loc

    ## ∫_0^1 x b_k^{(p)}(x) dx = (k+1)/((p+1)(p+2))
    bern_moment <- sum(c_loc * (1:(p + 1))) / ((p + 1) * (p + 2))

    ## shift from local to global x
    x_left <- xmin + (k - 1) * h
    moment <- moment + h * (x_left * sum(c_loc) / (p + 1) + h * bern_moment)
  }

  moment
}

ortho_x_int_condition <- function(p, q, xmin = 0, xmax = 1) {
  x <- numeric(q)
  nseg <- q - p
  dx <- (xmax - xmin) / nseg
  for (i in 1:q) {
    b <- rep(0, q)
    b[i] <- 1
    x[i] <- integrate_x_bspline(b, p = p, xmin = xmin, xmax = xmax) / dx
  }
  x
}

# helper
cross3 <- function(u, v) {
  c(
    u[2]*v[3] - u[3]*v[2],
    u[3]*v[1] - u[1]*v[3],
    u[1]*v[2] - u[2]*v[1]
  )
}

# Returns a q x (q-2) matrix M, orthogonal to
# int_f(x) and int_x*f(x)
ortho_x_diff_matrix <- function(p, q) {
  w0 <- ortho_int_condition(p = p, q = q)
  w1 <- ortho_x_int_condition(p = p, q = q)
  M <- matrix(0, nrow = q, ncol = q - 2)
  nseg <- q - p
  for (i in 1:(q-2)) {
    # solve locally a 3x3 model:
    wloc0 <- c(w0[i],w0[i+1],w0[i+2])
    wloc1 <- c(w1[i],w1[i+1],w1[i+2])
    a <- cross3(wloc0, wloc1)
    M[i,i] <- a[1]
    M[i+1,i] <- a[2]
    M[i+2,i] <- a[3]
  }
  spam::as.spam(M*nseg)
}




Gram_analytic <- function(p) {

  if (p == 0) {
    # piecewise constant
    M <- matrix(1, nrow = 1, ncol = 1)

  } else if (p == 1) {
    # linear B-splines
    M <- (1/6) * matrix(
      c(
        2, 1,
        1, 2
      ),
      nrow = 2, byrow = TRUE
    )

  } else if (p == 2) {
    # quadratic B-splines
    M <- (1/120) * matrix(
      c(
        6, 13,  1,
        13, 54, 13,
        1, 13,  6
      ),
      nrow = 3, byrow = TRUE
    )

  } else if (p == 3) {
    # cubic P-splines (extended knots)
    M <- matrix(
      c(
        1/252,     43/1680,  1/84,     1/5040,
        43/1680,   33/140,   311/1680, 1/84,
        1/84,      311/1680, 33/140,   43/1680,
        1/5040,    1/84,     43/1680,  1/252
      ),
      nrow = 4, byrow = TRUE
    )

  } else {
    stop("Analytic Gram matrix implemented only for p = 0,1,2,3")
  }

  M
}

# Returns a q x (q-1) matrix M, see IWSM-abstract
ortho_diff_matrix <- function(p, q) {
  # calculate the weights w_j
  w <- ortho_int_condition(p = p, q = q)
  M <- matrix(0, nrow = q, ncol = q - 1)
  for (i in 1:(q - 1)) {
    M[i, i]     <-  w[i + 1]
    M[i+1, i]   <- -w[i]
  }
  spam::as.spam(M)
}


