

counts <- c()
r <- 10
for(i in rbinom(100, r, .75)){
	counts <- c(counts, i, r - i)
}


genoLLHsum(counts[1:2], 4, .01)
genoLLHsum(counts[1:2], 6, .01)

genoEM(counts, 2, .01, 10000, .00)
genoEM(counts, 3, .01, 1000, .001)
genoEM(counts, 4, .01, 1000, .001)
genoEM(counts, 5, .01, 1000, .001)
genoEM(counts, 6, .01, 1000, .001)
genoEM(counts, 7, .01, 1000, .001)
genoEM(counts, 8, .01, 1000, .001)

