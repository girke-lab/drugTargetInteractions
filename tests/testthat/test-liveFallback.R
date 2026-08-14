## Tests for the live-API fallback mechanism (.dtiLiveOrCached) and the
## classed failure condition it keys on (.dtiSignalApiFailure).
##
## Entirely network-free: every case is driven by a synthetic expression,
## which is the point - the mechanism's whole job is to behave correctly
## when an upstream API misbehaves, and that is not something to leave to
## whether a real service happens to be down during a test run.
##
## The regression these guard against: the transport helpers degrade to
## NULL and emit a *warning* on an HTTP 500/timeout rather than erroring,
## so an error-only tryCatch let the empty live result through and never
## reached the fixture. Every REST source was affected; only callers that
## stop outright (biomaRt) were ever covered.

fx <- "ensembl_paralogs_nlrp3.rds"   # any shipped fixture with >0 rows

quietly <- function(expr) suppressWarnings(suppressMessages(expr))

test_that("a hard error falls back to the shipped fixture", {
    out <- quietly(.dtiLiveOrCached(stop("boom"), fixture = fx, label = "X"))
    expect_true(isTRUE(attr(out, "dtiCached")))
    expect_gt(nrow(out), 0L)
    expect_false(is.na(attr(out, "dtiCachedDate")))
})

test_that("a dtiApiFailure warning falls back too", {
    ## The case that was silently broken: warn + degrade to empty, exactly
    ## as .dtiApiGET() does on an HTTP 500.
    out <- quietly(.dtiLiveOrCached({
        .dtiSignalApiFailure("drugTargetInteractions API GET failed: HTTP 500")
        data.frame()
    }, fixture = fx, label = "X"))
    expect_true(isTRUE(attr(out, "dtiCached")))
    expect_gt(nrow(out), 0L)
})

test_that("an unrelated warning does NOT trigger a fallback", {
    ## Catching all warnings would over-trigger; the live result stands.
    out <- quietly(.dtiLiveOrCached({
        warning("some unrelated warning")
        data.frame(a = 1)
    }, fixture = fx, label = "X"))
    expect_false(isTRUE(attr(out, "dtiCached")))
    expect_identical(nrow(out), 1L)
})

test_that("a genuinely empty live result is not mistaken for an outage", {
    ## 0 rows is a correct answer for e.g. a gene with no paralogs, so it
    ## must not be used as the fallback trigger.
    out <- quietly(.dtiLiveOrCached(data.frame(), fixture = fx, label = "X"))
    expect_false(isTRUE(attr(out, "dtiCached")))
    expect_identical(nrow(out), 0L)
})

test_that("a successful live call is returned untouched", {
    out <- quietly(.dtiLiveOrCached(data.frame(a = 1:3), fixture = fx,
                                    label = "X"))
    expect_false(isTRUE(attr(out, "dtiCached")))
    expect_identical(out$a, 1:3)
})

test_that("the fallback reports which path was taken, and why", {
    expect_message(
        suppressWarnings(.dtiLiveOrCached({
            .dtiSignalApiFailure("HTTP 503 Service Unavailable")
            NULL
        }, fixture = fx, label = "Ensembl")),
        "live call to 'Ensembl' failed .*503.*using shipped cached result")
})

test_that(".dtiSignalApiFailure signals a warning of the expected class", {
    cnd <- tryCatch(.dtiSignalApiFailure("nope"), warning = function(w) w)
    expect_s3_class(cnd, "dtiApiFailure")
    expect_s3_class(cnd, "warning")
    expect_identical(conditionMessage(cnd), "nope")
    ## still an ordinary warning to any caller not looking for the class
    expect_warning(.dtiSignalApiFailure("nope"), "nope")
})

test_that("the transport helpers signal the class on a real transport failure", {
    ## Unroutable host: fails fast without depending on any live service.
    ## Catch `warning` rather than `condition` - httr2 signals its own
    ## httr2_perform instrumentation condition first, which is not the
    ## one under test (and which tryCatch's class-specific handler in
    ## .dtiLiveOrCached() likewise ignores).
    cnd <- tryCatch(
        .dtiApiGET("http://127.0.0.1:9/definitely-not-listening",
                   timeout = 2L, maxTries = 1L),
        warning = function(w) w)
    expect_s3_class(cnd, "dtiApiFailure")
})

test_that("hardStop = TRUE still errors rather than warning", {
    expect_error(
        .dtiApiGET("http://127.0.0.1:9/definitely-not-listening",
                   timeout = 2L, maxTries = 1L, hardStop = TRUE))
})
