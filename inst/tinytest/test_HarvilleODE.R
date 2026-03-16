
##
## Mask used in tests
##
Mask <- rbind(
  c(1,1),
  c(2,2),
  c(2,1),
  c(3,3)
)

##
## Test matrix <-> vector conversion
##
C <- matrix(
  c(4,1,0,
    1,5,0,
    0,0,6),
  3,3, byrow = TRUE
)

theta <- mat2vec(C, Mask)
C2 <- vec2mat(theta, Mask)

expect_equal(C2, C)


##
## Test eta -> xi -> eta roundtrip
##
eta <- c(0.3, -0.4, 0.2, 0.5)

xi <- eta_to_xi(eta, Mask)
eta2 <- xi_to_eta(xi, Mask)

expect_equal(eta2, eta, tolerance = 1e-10)


##
## Test xi -> theta consistency
##
xi <- matrix(
  c(2,0,0,
    1,3,0,
    0,0,1),
  3,3, byrow = TRUE
)

theta <- xi_to_theta(xi, Mask)
thetaMat <- vec2mat(theta, Mask)

expect_equal(thetaMat, tcrossprod(xi))


##
## Test theta -> xi inversion
##
xi <- eta_to_xi(c(0.2,0.1,0.3,-0.2), Mask)

theta <- xi_to_theta(xi, Mask)
thetaMat <- vec2mat(theta, Mask)

xi2 <- theta_to_xi(thetaMat)

expect_equal(xi2, xi, tolerance = 1e-10)


##
## Test rotation helper
##
xi <- matrix(
  c(2,0,0,
    1,3,0,
    0,0,1),
  3,3, byrow = TRUE
)

M <- diag(c(1,2,3))

rot <- rotate_theta(M, xi)

expect_equal(rot, crossprod(xi, M %*% xi))


##
## Full transformation pipeline test
##
eta <- rnorm(4)

xi <- eta_to_xi(eta, Mask)
theta <- xi_to_theta(xi, Mask)

thetaMat <- vec2mat(theta, Mask)

xi2 <- theta_to_xi(thetaMat)
eta2 <- xi_to_eta(xi2, Mask)

expect_equal(eta2, eta, tolerance = 1e-10)


##
## Additional edge case: zero correlations
##
eta <- c(0, 0, 0, 0)

xi <- eta_to_xi(eta, Mask)
theta <- xi_to_theta(xi, Mask)

thetaMat <- vec2mat(theta, Mask)

expect_true(isSymmetric(thetaMat))
expect_true(all(diag(thetaMat) > 0))
