context("mergeReplicates methods")

test_that("mergeReplicates on a matrix returns a matrix", {
    valueData <- matrix(1:100, ncol = 20,
                        dimnames = list(NULL, letters[1:20]))

    test1 <- rep(FALSE, 20)
    test1[c(2, 3, 5)] <- TRUE
    test2 <- rep(FALSE, 20)
    test2[c(11, 16, 20)] <- TRUE
    test3 <- rep(FALSE, 20)
    test3[c(9, 14)] <- TRUE

    LL <- IRanges::LogicalList(pt1 = test1, pt2 = test2, pt3 = test3)
    nCOLS <- sum(apply(as.matrix(LL), 2, function(x) !any(x)), length(LL))
    mergedObj <- mergeReplicates(valueData, replicates = LL,
                                  simplify = BiocGenerics::mean)

    expect_true(is(mergedObj, "matrix"))
    expect_identical(ncol(mergedObj), nCOLS)

})
