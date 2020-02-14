

counts <- c()
r <- 10
for(i in rbinom(100, r, .75)){
	counts <- c(counts, i, r - i)
}
# counts <- rep(c(20,20), 100)

testPloidy(counts, 2, .01)
testPloidy(counts, 4, .01)
testPloidy(counts, 5, .01)
testPloidy(counts, 6, .01)

testPloidy(counts, 4, .01) -
testPloidy(counts, 6, .01)

testPloidy(counts[1:2], 4, .01) -
testPloidy(counts[1:2], 6, .01)

for(i in seq(1, length(counts), 2)){
	print(counts[i] / r)
}
summary(counts[seq(1,length(counts),2)]/r)

4/6
5/6

exp(testPloidy(counts, 4, .01))
exp(testPloidy(counts, 6, .01))

testPloidy(c(3,1), 4, .01)
testPloidy(c(3,1), 6, .01)

27/40
tetra <- c()
sexa <- c()
for(j in 1:1000){
	counts <- c()
	r <- 100
	for(i in rbinom(20, r, .75)){
		counts <- c(counts, i, r - i)
	}

	tetra <- c(tetra, testPloidy(counts, 4, .01))
	sexa <- c(sexa, testPloidy(counts, 6, .01))
}
summary(tetra-sexa)

summary(tetra)
summary(sexa)

genoLLHsum(counts[1:2], 4, .01)
genoLLHsum(counts[1:2], 6, .01)

genoEM(counts, 2, .01, 10)
genoEM(counts, 3, .01, 10)
genoEM(counts, 4, .01, 10)
genoEM(counts, 5, .01, 10)
genoEM(counts, 6, .01, 10)
genoEM(counts, 7, .01, 10)
genoEM(counts, 8, .01, 10)


