data(oats.data)

obj <- LMMsolve(
  fixed  = yield ~ rep,
  random = ~ gen + rep:block,
  data   = oats.data
)

## --- basic functionality ---
expect_equal(
  LMMsolver::getHeritability(obj, "gen"),
  0.809,
  tol = 1e-3
)

## --- wrong class ---
expect_error(
  LMMsolver::getHeritability(lm(yield ~ rep, oats.data), "gen"),
  "obj must be an object of class 'LMMsolve'",
  fixed = TRUE
)

## --- wrong geno.term type ---
expect_error(
  LMMsolver::getHeritability(obj, 123),
  "geno.term must be a single character string",
  fixed = TRUE
)

## --- geno.term not in model ---
expect_error(
  LMMsolver::getHeritability(obj, "foo"),
  "foo not defined in the model",
  fixed = TRUE
)

## --- geno.term fixed ---
obj_fixed <- LMMsolve(
  fixed  = yield ~ gen + rep,
  random = ~ rep:block,
  data   = oats.data
)

expect_error(
  LMMsolver::getHeritability(obj_fixed, "gen"),
  "random"
)
