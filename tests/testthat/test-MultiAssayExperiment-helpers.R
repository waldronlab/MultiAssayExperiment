test_that("MultiAssayExperiment anyReplicated returns same order", {
    m <- matrix(0, 3, 3, dimnames=list(letters[1:3], letters[1:3]))
    m2 <- matrix(0, 0, 0)
    m3 <- matrix(0, 1, 1, dimnames=list("d", "d"))
    obs <- MultiAssayExperiment(list(m=m, m2=m2, m3=m3))
     # anyReplicated(obs)
})
