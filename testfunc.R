# testing functions during development

counts <- c()
r <- 10
for(i in rbinom(100, r, .75)){
	counts <- c(counts, i, r - i)
}
rCounts <- counts[seq(1, length(counts), 2)]
aCounts <- counts[seq(2, length(counts), 2)]

genoEM(rCounts, aCounts, 2, rep(1, length(rCounts)), rep(.01, length(rCounts)), 1000, .001, FALSE)
genoEM(rCounts, aCounts, 3, rep(1, length(rCounts)), rep(.01, length(rCounts)), 1000, .001, FALSE)
genoEM(rCounts, aCounts, 4, rep(1, length(rCounts)), rep(.01, length(rCounts)), 1000, .001, TRUE)
genoEM(rCounts, aCounts, 5, rep(1, length(rCounts)), rep(.01, length(rCounts)), 1000, .001, FALSE)
genoEM(rCounts, aCounts, 6, rep(1, length(rCounts)), rep(.01, length(rCounts)), 1000, .001, TRUE)
genoEM(rCounts, aCounts, 7, rep(1, length(rCounts)), rep(.01, length(rCounts)), 1000, .001, FALSE)
genoEM(rCounts, aCounts, 8, rep(1, length(rCounts)), rep(.01, length(rCounts)), 1000, .001, FALSE)



# pull out read counts from sample sturgeon data

genosFiles <- dir(path = "S:\\Eagle Fish Genetics Lab\\Tom\\sturgeon ploidy\\Sacramento_parentage_genos",
			   pattern = "\\.genos$")

genosFiles <- paste0("S:\\Eagle Fish Genetics Lab\\Tom\\sturgeon ploidy\\Sacramento_parentage_genos\\",
				 genosFiles)

refCounts <- matrix(nrow = 0, ncol = 325)
altCounts <- matrix(nrow = 0, ncol = 325)
for(f in genosFiles){
	rReads <- c() # ref
	aReads <- c() # alt
	mNames <-c() # locus names
	gFile <- file(f, "r")
	# get sample name
	line <- readLines(gFile, n = 1)
	sName <- gsub("\\.fastq$", "", strsplit(line, ",")[[1]][1])

	line <- readLines(gFile, n = 1) # read first marker line
	while(length(line) > 0){
		sep <- strsplit(line, ",")[[1]]
		mNames <- c(mNames, sep[1])
		rReads <- c(rReads, as.numeric(gsub("^[ACTG-]=", "", sep[2])))
		aReads <- c(aReads, as.numeric(gsub("^[ACTG-]=", "", sep[3])))
		line <- readLines(gFile, n = 1)
	}
	close(gFile)

	# save data, use names to make sure all in same order
	names(rReads) <- mNames
	names(aReads) <- mNames
	refCounts <- rbind(refCounts, rReads)
	altCounts <- rbind(altCounts, aReads)
	rownames(refCounts)[nrow(refCounts)] <- sName
	rownames(altCounts)[nrow(altCounts)] <- sName


}

# save(refCounts, altCounts, file = "sturgeonData.rda")
load("sturgeonData.rda")

fpTest <- funkyPloid(refCounts[1:2,], altCounts[1:2,], ploidy = c(4,6), maxRep = 10000, maxDiff = .0001)

mB <- refCounts[1,] + altCounts[1,] > 0
genoEM(refCounts[1,mB], altCounts[1,mB],
              4, rep(1, ncol(refCounts)), rep(.01, ncol(refCounts)), 100,
              .00, TRUE)
genoEM(refCounts[1,mB], altCounts[1,mB],
              6, rep(1, ncol(refCounts)), rep(.01, ncol(refCounts)), 100,
              .00, TRUE)

genoEM(refCounts[1,1], altCounts[1,1],
              4, rep(1, ncol(refCounts)), rep(.01, ncol(refCounts)), 100,
              .001, TRUE)
genoEM(refCounts[1,1], 0,
              4, rep(1, ncol(refCounts)), rep(.01, ncol(refCounts)), 100,
              .001, TRUE)
tllr <- c()
for(i in which(mB)){
	tllr <- c(tllr, genoEM(refCounts[1,i], altCounts[1,i],
	              4, rep(1, ncol(refCounts)), rep(.01, ncol(refCounts)), 100,
	              .001, FALSE))
}
refCounts[1, which(mB)[which(is.nan(tllr))]]
altCounts[1, which(mB)[which(is.nan(tllr))]]

genoEM(1000, 1200,
              4, rep(1, ncol(refCounts)), rep(.01, ncol(refCounts)), 100,
              .001, TRUE)
system.time(
fpTest <- funkyPloid(refCounts, altCounts, ploidy = c(4,5,6), maxRep = 10000, maxDiff = .0001)
)

fpTest$Ind <- gsub("^i09[0-9]_[0-9]+_P745[0-9]_WSTG2017-", "", fpTest$Ind)
head(fpTest)
sum(fpTest[,3] == 0)
# 48
sum(fpTest[,4] == 0)
# 19
sum(fpTest[,5] == 0)
# 29

write.table(fpTest, "fpTest.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

totC <- refCounts + altCounts

hist(rowSums(totC))
sort(rowSums(totC))
totC <- totC[rowSums(totC) > 40000,]

hist(totC[,8])

apply(totC, 1, mean)
var(totC)


totCProp <- t(apply(totC, 1, function(x) x / sum(x)))
totCProp <- apply(totCProp, 2, function(p) log(p / (1 - p)))
hist(totCProp[,1])
hist(totCProp[,2])
hist(totCProp[,3])
hist(totCProp[,4])
hist(totCProp[,5])
hist(totCProp[,6])
hist(totCProp[,7])
hist(totCProp[,8])
hist(totCProp[,9])
hist(totCProp[,10])

# fitting some models
library(MGLM)

str(totC)

mnModel <- MGLMfit(totC, dist = "MN")
show(mnModel)

dmModel <- MGLMfit(totC, dist = "DM")
show(dmModel)

# Distribution: Dirichlet Multinomial
# Log-likelihood: -179494.8
# BIC: 360421.8
# AIC: 359639.6
# LRT test p value: <0.0001
# Iterations: 7

class(dmModel)
coef(dmModel)
str(dmModel)

dmModel@estimate
save.image("after_DM_fit.RData")


# training / test data sets from Stuart

# pull out read counts from sample sturgeon data

genosFiles <- dir(path = "S:\\Eagle Fish Genetics Lab\\Tom\\sturgeon ploidy\\ploidy_genos",
			   pattern = "\\.genos$")

genosFiles <- paste0("S:\\Eagle Fish Genetics Lab\\Tom\\sturgeon ploidy\\ploidy_genos\\",
				 genosFiles)

refCounts <- matrix(nrow = 0, ncol = 325)
altCounts <- matrix(nrow = 0, ncol = 325)
sName_true <- c()
for(f in genosFiles){
	rReads <- c() # ref
	aReads <- c() # alt
	mNames <-c() # locus names
	gFile <- file(f, "r")
	# get sample name
	line <- readLines(gFile, n = 1)
	# not pulling name from .genos file
	sName <- gsub("S:\\\\Eagle Fish Genetics Lab\\\\Tom\\\\sturgeon ploidy\\\\ploidy_genos\\\\", "", f)
	sName_true <- c(sName_true, gsub("\\.fastq$", "", strsplit(line, ",")[[1]][1]))

	line <- readLines(gFile, n = 1) # read first marker line
	while(length(line) > 0){
		sep <- strsplit(line, ",")[[1]]
		mNames <- c(mNames, sep[1])
		rReads <- c(rReads, as.numeric(gsub("^[ACTG-]=", "", sep[2])))
		aReads <- c(aReads, as.numeric(gsub("^[ACTG-]=", "", sep[3])))
		line <- readLines(gFile, n = 1)
	}
	close(gFile)

	# save data, use names to make sure all in same order
	names(rReads) <- mNames
	names(aReads) <- mNames
	refCounts <- rbind(refCounts, rReads)
	altCounts <- rbind(altCounts, aReads)
	rownames(refCounts)[nrow(refCounts)] <- sName
	rownames(altCounts)[nrow(altCounts)] <- sName


}

# save.image("sturgeonData2.rda")
load("sturgeonData2.rda")


system.time(
fp_456_new <- funkyPloid(refCounts, altCounts, ploidy = c(4,5,6), maxIter = 10000, maxDiff = .0001)
)
head(fp_456_new)
fp_456_new[,1]
ploid <- gsub("(.+-)|(_._.+)", "", fp_456_new[,1])
ploid2 <- gsub("(^.+((WSTGUCD20)|(20-BLCJ))-)|(-[0-9]+$)", "", sName_true)
all <- cbind(ploid2, ploid, sName_true, fp_456_new, stringsAsFactors = FALSE)
plot(log(fp_456_new$LLR_4 + .0001), log(fp_456_new$LLR_6 + .0001), col = c("red", "blue", "green")[as.factor(ploid)])
plot(fp_456_new$LLR_4, fp_456_new$LLR_6)

write.table(all, "summary.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(fp_456_new, "summary_for_Stuart.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


sName_true
toPlot <- c("WSTGUCD20-12N_2_10X.genos",
"WSTGUCD20-8N_1_10X.genos",
"WSTGUCD20-8N_5_10X.genos",
"WSTGUCD20-UNK_C_10X.genos")
pdf("oddInds.pdf")
for(i in toPlot){

		hist(refCounts[i,] / (refCounts[i,] + altCounts[i,]), breaks = 40, xaxt = "n", main = i)
		axis(side=1, at=c(seq(0,1,.25), seq(0,1,1/6)), labels=round(c(seq(0,1,.25), seq(0,1,1/6)),2), cex.axis = .9
			)

}
dev.off()
i <- "WSTG20-BLCJ-UNK_A_10X.genos"

funkyPloid(refCounts, altCounts, ploidy = c(4,5,6), maxIter = 10000, maxDiff = .0001)
gprops <- genoProps(refCounts, altCounts, ploidy = 4, maxIter = 10000, maxDiff = .0000001)

head(gprops)
max(gprops$numIter)
min(gprops$Loci)

gprops[1:20,]

