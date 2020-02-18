

counts <- c()
r <- 10
for(i in rbinom(100, r, .75)){
	counts <- c(counts, i, r - i)
}
rCounts <- counts[seq(1, length(counts), 2)]
aCounts <- counts[seq(2, length(counts), 2)]

genoLLHsum(counts[1:2], 4, .01)
genoLLHsum(counts[1:2], 6, .01)

genoEM(rCounts, aCounts, 2, .01, 1000, .001, FALSE)
genoEM(rCounts, aCounts, 3, .01, 1000, .001, FALSE)
genoEM(rCounts, aCounts, 4, .01, 1000, .001, TRUE)
genoEM(rCounts, aCounts, 5, .01, 1000, .001, FALSE)
genoEM(rCounts, aCounts, 6, .01, 1000, .001, TRUE)
genoEM(rCounts, aCounts, 7, .01, 1000, .001, FALSE)
genoEM(rCounts, aCounts, 8, .01, 1000, .001, FALSE)


