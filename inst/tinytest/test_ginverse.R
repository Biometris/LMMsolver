library(Matrix)
library(spam)

## ------------------------------------------------------------------
## Test data
## ------------------------------------------------------------------

lev <- as.character(1:5)
K <- diag(1, 5)
dimnames(K) <- list(lev, lev)

K_matrix <- K
K_Matrix <- Matrix::Matrix(K, sparse = TRUE)
K_spam <- spam::as.spam(K)

data_ok <- data.frame(
  id = factor(c("1", "2", "3", "4", "5", "1"))
)


## ==================================================================
## as.ginverse()
## ==================================================================

## Valid matrix
g <- as.ginverse(
  list(id = K_matrix),
  levels = list(id = lev)
)

expect_true(inherits(g, "ginverse"))
expect_identical(names(g), "id")
expect_identical(attr(g, "levels"), list(id = lev))
expect_equal(attr(g, "tol"), 1e-10)


## Valid Matrix object
g <- as.ginverse(
  list(id = K_Matrix),
  levels = list(id = lev)
)

expect_true(inherits(g, "ginverse"))
expect_true(inherits(g$id, "Matrix"))
expect_identical(attr(g, "levels"), list(id = lev))


## Valid spam object
g <- as.ginverse(
  list(id = K_spam),
  levels = list(id = lev)
)

expect_true(inherits(g, "ginverse"))
expect_true(inherits(g$id, "spam"))
expect_identical(attr(g, "levels"), list(id = lev))


## precisionMatrices is not a list
expect_error(
  as.ginverse(K_matrix, levels = list(id = lev)),
  "precisionMatrices must be a list"
)


## precisionMatrices is an unnamed list
expect_error(
  as.ginverse(list(K_matrix), levels = list(id = lev)),
  "precisionMatrices must be a named list"
)


## levels is not a list
expect_error(
  as.ginverse(
    list(id = K_matrix),
    levels = lev
  ),
  "'levels' must be a named list"
)


## levels is an unnamed list
expect_error(
  as.ginverse(
    list(id = K_matrix),
    levels = list(lev)
  ),
  "'levels' must be a named list"
)


## levels has different names
expect_error(
  as.ginverse(
    list(id = K_matrix),
    levels = list(animal = lev)
  ),
  "'levels' must have the same names as 'precisionMatrices'"
)


## Invalid matrix class
expect_error(
  as.ginverse(
    list(id = 1:5),
    levels = list(id = lev)
  ),
  "'id' must be a matrix, Matrix, or spam object"
)


## Non-square matrix
K_rect <- matrix(1, nrow = 5, ncol = 4)

expect_error(
  as.ginverse(
    list(id = K_rect),
    levels = list(id = lev)
  ),
  "'id' must be square"
)


## levels is NULL
expect_error(
  as.ginverse(
    list(id = K_matrix),
    levels = list(id = NULL)
  ),
  "'levels\\[\\[id\\]\\]' must not be NULL"
)


## Wrong number of levels
expect_error(
  as.ginverse(
    list(id = K_matrix),
    levels = list(id = as.character(1:4))
  ),
  "'levels\\[\\[id\\]\\]' must have length 5"
)


## Duplicate levels
expect_error(
  as.ginverse(
    list(id = K_matrix),
    levels = list(id = c("1", "2", "3", "4", "4"))
  ),
  "'levels\\[\\[id\\]\\]' contains duplicated levels"
)


## levels is not character
expect_error(
  as.ginverse(
    list(id = K_matrix),
    levels = list(id = 1:5)
  ),
  "'levels\\[\\[id\\]\\]' must be a character vector"
)


## Missing levels argument
expect_error(
  as.ginverse(list(id = K_matrix)),
  "argument \"levels\" is missing"
)


## Custom tolerance
g <- as.ginverse(
  list(id = K_matrix),
  levels = list(id = lev),
  tol = 1e-8
)

expect_equal(attr(g, "tol"), 1e-8)


## Multiple matrices
K2 <- diag(1, 3)
lev2 <- c("a", "b", "c")

g <- as.ginverse(
  list(id = K_matrix, animal = K2),
  levels = list(id = lev, animal = lev2)
)

expect_identical(names(g), c("id", "animal"))
expect_identical(
  attr(g, "levels"),
  list(id = lev, animal = lev2)
)


## ==================================================================
## checkGinverseAgainstData()
## ==================================================================

g <- as.ginverse(
  list(id = K_matrix),
  levels = list(id = lev)
)


## Valid case
out <- LMMsolver:::checkGinverseAgainstData(
  g,
  data_ok,
  random_terms = "id"
)

expect_true(is.list(out))
expect_identical(names(out), "id")
expect_equal(dim(out$id), c(5, 5))


## Valid case with reordering
data_reordered <- data.frame(
  id = factor(
    c("3", "1", "5", "2"),
    levels = c("3", "1", "5", "2")
  )
)

K2 <- matrix(
  1:9,
  nrow = 3,
  dimnames = list(c("1", "2", "3"), c("1", "2", "3"))
)

#g2 <- as.ginverse(
#  list(id = K2),
#  levels = list(id = c("1", "2", "3"))
#)

#out <- LMMsolver:::checkGinverseAgainstData(
#  g2,
#  data_reordered,
#  random_terms = "id"
#)

#expect_identical(
#  attr(out$id, "dimnames"),
#  list(c("3", "1", "2"), c("3", "1", "2"))
#)


## ginverse has wrong class
expect_error(
  LMMsolver:::checkGinverseAgainstData(
    list(id = K_matrix),
    data_ok,
    random_terms = "id"
  ),
  "ginverse must be created with as.ginverse\\(\\)"
)


## ginverse has no levels attribute
g_no_levels <- list(id = K_matrix)
class(g_no_levels) <- c("ginverse", "list")

expect_error(
  LMMsolver:::checkGinverseAgainstData(
    g_no_levels,
    data_ok,
    random_terms = "id"
  ),
  "ginverse does not contain level information"
)


## ginverse not present in random effects
expect_error(
  LMMsolver:::checkGinverseAgainstData(
    g,
    data_ok,
    random_terms = "animal"
  ),
  "ginverse 'id' not present in random effects"
)


## Data column missing
expect_error(
  LMMsolver:::checkGinverseAgainstData(
    g,
    data.frame(animal = factor(lev)),
    random_terms = "id"
  ),
  "Column 'id' not found in data"
)

## ginverse levels are NULL
g_bad <- g
attr(g_bad, "levels") <- list(id = NULL)

expect_error(
  LMMsolver:::checkGinverseAgainstData(
    g_bad,
    data_ok,
    random_terms = "id"
  ),
  "No level information available for ginverse 'id'"
)

## Number of levels does not match matrix dimensions
g_bad <- g
attr(g_bad, "levels") <- list(id = as.character(1:4))

expect_error(
  LMMsolver:::checkGinverseAgainstData(
    g_bad,
    data_ok,
    random_terms = "id"
  ),
  "Number of levels for ginverse 'id' does not match matrix dimensions"
)


## Duplicated ginverse levels
g_bad <- g
attr(g_bad, "levels") <- list(
  id = c("1", "2", "3", "4", "4")
)

expect_error(
  LMMsolver:::checkGinverseAgainstData(
    g_bad,
    data_ok,
    random_terms = "id"
  ),
  "ginverse 'id' contains duplicated levels"
)


## Data contains a level absent from ginverse
data_missing <- data.frame(
  id = factor(c("1", "2", "3", "6"))
)

#expect_error(
#  LMMsolver:::checkGinverseAgainstData(
#    g,
#    data_missing,
#   random_terms = "id"
#  ),
#  "ginverse 'id' missing levels: 5"
#)


## spam survives the check as spam
g_spam <- as.ginverse(
  list(id = K_spam),
  levels = list(id = lev)
)

out <- LMMsolver:::checkGinverseAgainstData(
  g_spam,
  data_ok,
  random_terms = "id"
)

expect_true(inherits(out$id, "spam"))


