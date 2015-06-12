library(devtools)
reload()
source("R/functions.R")
source("R/S3_helpers.R")

load(system.file("data", "b3_wang_sdist.Rdata", package = "SimGTF"))
fwang_sdist = wang_sdist
fwang_sdist$b_1 = sample(c(-1,1),nrow(fwang_sdist),replace=TRUE)*(abs(fwang_sdist$b_1) + log(2)) 
summary(fwang_sdist)
fwang_sdist = subset(fwang_sdist, b_0 > -15 & b_0 < -7)
summary(fwang_sdist)

load(system.file("data", "Ensembl_80_Homo_sapien_GTFSimData.Rdata", package="SimGTF"))



N1 = N2 = 3 #number of samples per treatment
phex = PhenoExpression(
	trans_dist = fwang_sdist,
	sample_names = paste("count",1:(N1+N2),sep="_"),
	treatments = c(rep(1,N1),rep(2,N2)),
	E_lib_size = rep(6*10^6, N1+N2),
	de_prob = .2,
	min_transcripts = 1,
	max_transcripts = 2
)

sim = DisjointGenomeSimulation(sim_data, phex)
gene_counts = geneLocusCounts(sim)
disjoint_counts = disjointExonCounts(sim)
gene_counts2 = geneLocusExactTest(gene_counts, phex$treatments)
disjoint_counts2 = disjointExactTest(disjoint_counts, phex$treatments)


gres = GeneLocusResults(sim)
gcounts = counts(gres)
gde = actuallyDE(gres)

