data(oats.data)

obj <- LMMsolve(
  fixed  = yield ~ rep,
  random = ~ gen + rep:block,
  data   = oats.data
)

## --- basic functionality ---
expect_equal(
  getHeritability(obj, "gen"),
  0.809,
  tol = 1e-3
)

## --- wrong class ---
expect_error(
  getHeritability(lm(yield ~ rep, oats.data), "gen"),
  "gen not defined in the model",
  fixed = TRUE
)

## --- wrong geno.term type ---
expect_error(
  getHeritability(obj, 123),
  "123 not defined in the model",
  fixed = TRUE
)

## --- geno.term not in model ---
expect_error(
  getHeritability(obj, "foo"),
  "not defined"
)

## --- geno.term fixed ---
obj_fixed <- LMMsolve(
  fixed  = yield ~ gen + rep,
  random = ~ rep:block,
  data   = oats.data
)

expect_error(
  getHeritability(obj_fixed, "gen"),
  "random"
)
