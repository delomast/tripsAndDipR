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
load("after_DM_fit.RData")



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
gprops4 <- genoProps(refCounts, altCounts, ploidy = 4, maxIter = 10000, maxDiff = .0000001)
gprops6 <- genoProps(refCounts, altCounts, ploidy = 6, maxIter = 10000, maxDiff = .0000001)

indsLowMinLLR <- fp_456_new[grepl("10X", fp_456_new$Ind) & apply(fp_456_new[,3:5], 1, median) < 2000,1]

apply(gprops4[grepl("8N", gprops4$Ind),5:9],2,mean)

gprops4[gprops4$Ind %in% indsLowMinLLR,5:9]


head(gprops)
max(gprops$numIter)
min(gprops$Loci)

gprops[1:20,]

########################
## looking at Monte Carlo
########################

load("after_DM_fit.RData")


ploidy <- 6
alpha <- dmModel@estimate

# just making up some genotype proportions
# based on extended HWE (no double reductions)
totalGenoProbs <- data.frame()
obs <- 100
# for(i in 1:length(alpha)){
for(i in 1:100){
	reF <- runif(1)
	genoprobs <- rep(NA, ploidy + 1)
	for(p in 0:ploidy){
		genoprobs[p+1] <- choose(ploidy, p) * (reF^p) * ((1 - reF)^(ploidy - p))
	}
	# posterior with a uniform prior
	genoprobs <- as.vector(rmultinom(1, obs, genoprobs))
	totalGenoProbs <- rbind(totalGenoProbs, genoprobs)
}

randSamples <- rPloidySamples(1000, 20000, 4, alpha, eps = NULL, h = NULL, totalGenoProbs,
				genotypePropsAreKnown = FALSE)
str(randSamples)

mcError <- funkyPloid(randSamples$counts, randSamples$counts_alt, ploidy = c(4,5,6), maxIter = 10000, maxDiff = .0001)

table(mcError[,3])




##########################

# testing different input formats

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

fp_sep <- funkyPloid(refCounts, altCounts, ploidy = c(4,5,6), maxIter = 10000, maxDiff = .000001)

twocol <- data.frame()
for(i in 1:ncol(refCounts)){
	twocol <- rbind(twocol, refCounts[,i], altCounts[,i])
}
twocol <- as.data.frame(t(twocol))

fp_together <- funkyPloid(twocol, ploidy = c(4,5,6), maxIter = 10000, maxDiff = .000001)

identical(fp_sep[,2:ncol(fp_sep)], fp_together[,2:ncol(fp_together)])
all.equal(fp_sep, fp_together)


i <- rownames(refCounts)[grepl("1LF12_OE5Ex2E4C", rownames(refCounts))]
hist(refCounts[i,] / (refCounts[i,] + altCounts[i,]), breaks = 40, xaxt = "n", main = i)
axis(side=1, at=c(seq(0,1,.25), seq(0,1,1/6)), labels=round(c(seq(0,1,.25), seq(0,1,1/6)),2), cex.axis = .9
	)

plot(refCounts[i,], altCounts[i,])
sum(refCounts[i,], altCounts[i,])

genoProps(refCounts[i,,drop = FALSE], altCounts[i,,drop = FALSE], ploidy = 4, maxIter = 10000, maxDiff = .000001)
genoProps(refCounts[i,,drop = FALSE], altCounts[i,,drop = FALSE], ploidy = 6, maxIter = 10000, maxDiff = .000001)
sum(log(1 / (1+r+a)))

fourSixUnif <- cbind(
	genoProps(refCounts, altCounts, ploidy = 4, maxIter = 10000, maxDiff = .000001)$LLH,
genoProps(refCounts, altCounts, ploidy = 6, maxIter = 10000, maxDiff = .000001)$LLH,
	apply(refCounts + altCounts, 1, function(x){
		sum(-log(1+x))
	})

)
head(fourSixUnif)
fourSixUnif <- t(apply(fourSixUnif, 1, function(x) max(x) - x))


toPlot <- c("3HF6_OE5ExO610",
"3HF2_OE5ExO610",
"1HF17_6F6BxO610",
"3HF5_OE5ExO610",
"1LF12_OE5Ex2E4C",
"3HF8_OE5ExO610",
"1LF15_OE5Ex2E4C",
"1LF11_6F6BxO610",
"1HF23_OE5Ex2E4C",
"3LF10_OE5ExO610",
"1HF14_6F6BxO610",
"1LF47_6F6BxO610",
"3LF9_OE5ExO610",
"1HF45_OE5Ex2E4C",
"1HF41_6F6BxO610",
"1HF48_6158x2E4C",
"1HF19_OE5Ex2E4C"
)

pdf("oddInds.pdf")
for(p in toPlot){
	i <- rownames(refCounts)[grepl(p, rownames(refCounts))]
	r <- refCounts[i,]
	a <- altCounts[i,]
	# b <- r > 5 & a > 5
	# r <- r[b]
	# a <- a[b]
	hist(r / (r + a), breaks = 40, xaxt = "n", main = i)
	axis(side=1, at=c(seq(0,1,.25), seq(0,1,1/6)), labels=round(c(seq(0,1,.25), seq(0,1,1/6)),2), cex.axis = .9)
}
dev.off()

pdf("goodPlots.pdf")
for(p in c("1HF36_6158x2E4C", "DAM-2808")){
	i <- rownames(refCounts)[grepl(p, rownames(refCounts))]
	r <- refCounts[i,]
	a <- altCounts[i,]
	hist(r / (r + a), breaks = 40, xaxt = "n", main = i)
	axis(side=1, at=c(seq(0,1,.25), seq(0,1,1/6)), labels=round(c(seq(0,1,.25), seq(0,1,1/6)),2), cex.axis = .9)
}
dev.off()



# looking at known ploidy
ref8N <- refCounts[grepl("8N", rownames(refCounts)),]
alt8N <- altCounts[grepl("8N", rownames(altCounts)),]

fpInd <- funkyPloid(ref8N, alt8N, ploidy = c(4,6,8), maxIter = 10000, maxDiff = .000001)
inds8n <- fpInd[fpInd$LLR_6 > 0,1]
ref8N <- refCounts[inds8n,]
alt8N <- altCounts[inds8n,]
fpInd <- funkyPloid(ref8N, alt8N, ploidy = c(4,6,8), maxIter = 10000, maxDiff = .000001)
head(fpInd)
gp4Ind <- genoProps(ref8N, alt8N, ploidy = 4, maxIter = 10000, maxDiff = .000001)
gp8Ind <- genoProps(ref8N, alt8N, ploidy = 8, maxIter = 10000, maxDiff = .000001)
head(gp4Ind)
head(gp8Ind)

i <- "WSTG20-BLCJ-8N_1_10X.genos"
hist(refCounts[i,] / (refCounts[i,] + altCounts[i,]), breaks = 40, xaxt = "n", main = i)
axis(side=1, at=c(seq(0,1,.25), seq(0,1,1/8)), labels=round(c(seq(0,1,.25), seq(0,1,1/8)),2), cex.axis = .9
	)

fpLoc <- funkyPloid(t(ref8N), t(alt8N), ploidy = c(4,8), maxIter = 10000, maxDiff = .000001)
round(fpLoc[,3:5], 2)
head(fpLoc)
gp4Loc <- genoProps(t(ref8N), t(alt8N), ploidy = 4, maxIter = 10000, maxDiff = .000001)
gp8Loc <- genoProps(t(ref8N), t(alt8N), ploidy = 8, maxIter = 10000, maxDiff = .000001)
head(gp4Loc)
head(gp8Loc)

i <- "Atr_10322-43"
hist(t(ref8N)[i,] / (t(ref8N)[i,] +
t(alt8N)[i,])
)

fpInd <- funkyPloid(ref8N, alt8N, ploidy = c(4,6,8), maxIter = 10000, maxDiff = .000001, noise = TRUE)
gp4Ind <- genoProps(ref8N, alt8N, ploidy = 4, maxIter = 10000, maxDiff = .000001, noise = TRUE)
gp8Ind <- genoProps(ref8N, alt8N, ploidy = 8, maxIter = 10000, maxDiff = .000001, noise = TRUE)

fpInd2 <- funkyPloid(ref8N, alt8N, ploidy = c(4,5,6), maxIter = 10000, maxDiff = .000001, noise = TRUE)



multTau <- function(pars, r, a, ploidy){
	if(any(pars <= 0 | pars >= 1)) return (Inf)
	b <- r + a > 0
	# ncat <- ploidy + 1
	ncat <- ploidy + 2
	tau <- pars[1:ncat]
	mixW <- pars[(ncat + 1):(2 * ncat)]
	mixW <- mixW / sum(mixW)
# 	return(
# 		-llh_calc_BB(r[b], a[b], tau, mixW,
#               ploidy, rep(1, sum(b)), rep(.01, sum(b)))
# 	)
	return(
		-llh_calc_BB_noise(r[b], a[b], tau, mixW,
              ploidy, rep(1, sum(b)), rep(.01, sum(b)))
	)
}

oneTau <- function(pars, r, a, ploidy){
	if(any(pars <= 0 | pars >= 1)) return (Inf)
	b <- r + a > 0
	tau <- pars[1]
	mixW <- pars[2:length(pars)]
	mixW <- mixW / sum(mixW)
	return(
		-llh_calc_BB(r[b], a[b], rep(tau, length(mixW)), mixW,
              ploidy, rep(1, sum(b)), rep(.01, sum(b)))
	)
# 	return(
# 		-llh_calc_BB_noise(r[b], a[b], rep(tau, length(mixW)), mixW,
#               ploidy, rep(1, sum(b)), rep(.01, sum(b)))
# 	)
}

pToTest <- 6
c("1HF36_6158x2E4C", "DAM-2808")
p <- "DAM-2808"
p <- "1HF36_6158x2E4C"
i <- rownames(refCounts)[grepl(p, rownames(refCounts))]

# r2 <- optim(rep(.01, pToTest + 2), oneTau, method = "Nelder-Mead", r = refCounts[i,], a = altCounts[i,], ploidy = pToTest,
# 		  control = list(maxit=10000))

r2 <- optim(rep(.01, (pToTest+2)*2), multTau, method = "Nelder-Mead", r = refCounts[i,], a = altCounts[i,], ploidy = pToTest,
		  control = list(maxit=10000))

mixW <- r2$par[(pToTest+1):length(r2$par)]
mixW / sum(mixW)


#  No noise
# 6 : 1758.296
# 4: 1824

# 4: 1746
# 6: 1809

# noise
# 6 : 1767
# 4: 1828

# 4: 1763
# 6: 1792

bpoints <- c()
binpoints <- c()
for(k in 0:pToTest){
	# s <- (1 - r2$par[1]) / r2$par[1]
	s <- (1 - r2$par[k+1]) / r2$par[k+1]
	p <- k / pToTest
	p <- p * .99 + (1 - p) * .01
	bpoints <- c(bpoints, rbeta(1000, p * s, (1-p) * s))
	binpoints <- c(binpoints, rbeta(1000, p * 1e9, (1-p) * 1e9))

}

hist(bpoints, breaks = 40, freq = FALSE, add = FALSE, col = "blue")
hist(binpoints, breaks = 40, freq = FALSE, add = TRUE, col = "red")
hist(refCounts[i,] / (refCounts[i,] + altCounts[i,]), breaks = 40, freq = FALSE, add = TRUE, col = "gray")

hist(refCounts[i,] / (refCounts[i,] + altCounts[i,]), breaks = 40, freq = FALSE, add = FALSE, col = "gray")
hist(bpoints, breaks = 40, freq = FALSE, add = TRUE, col = "blue")


pToTest <- 4
totalGenoProbs2 <- totalGenoProbs
totalGenoProbs2[,2] <- totalGenoProbs2[,2] + 1000
mcCounts <- rPloidySamples(2, 10000, pToTest, rep(10,100), genotypeCounts = totalGenoProbs)
mcCounts <- rPloidySamples(2, 10000, pToTest, rep(10,100), genotypeCounts = totalGenoProbs2)

r2 <- optim(rep(.02, pToTest + 2), oneTau, method = "Nelder-Mead", r = mcCounts$counts[1,],
		  a = mcCounts$counts_alt[1,], ploidy = pToTest, control = list(maxit=10000))

bpoints <- c()
binpoints <- c()
for(k in 0:pToTest){
	s <- (1 - r2$par[1]) / r2$par[1]
	p <- k / pToTest
	p <- p * .99 + (1 - p) * .01
	bpoints <- c(bpoints, rbeta(1000, p * s, (1-p) * s))
	binpoints <- c(binpoints, rbeta(1000, p * 1e9, (1-p) * 1e9))

}

hist(bpoints, breaks = 40, freq = FALSE, add = FALSE, col = "blue")
hist(binpoints, breaks = 40, freq = FALSE, add = TRUE, col = "red")
hist(mcCounts[[1]][1,] / (mcCounts[[1]][1,] + mcCounts[[2]][1,]), breaks = 40, freq = FALSE, add = TRUE, col = "gray")


# Beta-binomial mixture models

r <- refCounts[i,]
a <- altCounts[i,]
b <- r+a > 0

BBpolyEM(r[b], a[b], 4, h = rep(1, sum(b)), eps = rep(.01, sum(b)), noise = TRUE,
	    mdiff = .0001, maxrep = 30)

BBpolyEM(r[b], a[b], 4, h = rep(1, sum(b)), eps = rep(.01, sum(b)), noise = FALSE,
	    mdiff = .0001, maxrep = 30)

gp4 <- genoProps(refCounts, altCounts, ploidy = 4, h = NULL, eps = NULL,
				   maxIter = 100, maxDiff = .001, model = "BB_noise")
gp5 <- genoProps(refCounts, altCounts, ploidy = 5, h = NULL, eps = NULL,
				   maxIter = 100, maxDiff = .001, model = "BB_noise")
gp6 <- genoProps(refCounts, altCounts, ploidy = 6, h = NULL, eps = NULL,
				   maxIter = 100, maxDiff = .001, model = "BB_noise")

llr <- cbind(gp4$LLH, gp5$LLH, gp6$LLH)
noise <- cbind(gp4$noise, gp5$noise, gp6$noise)
llr <- t(apply(llr, 1, function(x) round(max(x) - x, 3)))

tolook <- cbind(rownames(refCounts), llr)
write.table(tolook, "tolook.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(noise, "noise.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

debugonce(genoProps)

# BB with known

gp4 <- genoProps(refCounts, altCounts, ploidy = 4, h = NULL, eps = NULL,
				   maxIter = 100, maxDiff = .001, model = "BB_noise")
gp5 <- genoProps(refCounts, altCounts, ploidy = 5, h = NULL, eps = NULL,
				   maxIter = 100, maxDiff = .001, model = "BB_noise")
gp6 <- genoProps(refCounts, altCounts, ploidy = 6, h = NULL, eps = NULL,
				   maxIter = 100, maxDiff = .001, model = "BB_noise")

llr <- cbind(gp4$LLH, gp5$LLH, gp6$LLH)
noise <- cbind(gp4$noise, gp5$noise, gp6$noise)
llr <- t(apply(llr, 1, function(x) round(max(x) - x, 3)))
tolook <- cbind(rownames(refCounts), llr)
write.table(tolook, "tolook2.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(noise, "noise2.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


gp4 <- genoProps(refCounts[1:2,], altCounts[1:2,], ploidy = 4, h = NULL, eps = NULL,
				   maxIter = 100, maxDiff = .001, model = "BB_noise")
