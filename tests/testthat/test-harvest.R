# BADARG validation RED test deferred to Task 7 where harvest() first exists.
# Calling .Call("C_rk_calibrate") before R_init_leafblower registers the symbol
# crashes R with an opaque symbol-not-found error — not a named assertion failure.
# Task 7 adds the real RED test: expect_error(harvest(df, tgt, min_weight=5, max_weight=5), "min_weight")
