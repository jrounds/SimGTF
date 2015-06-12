library(data.table)
library(devtools)
wd ="/media/sf_work/15spring/working_SimGTF/SimGTF"
setwd(wd)
install()
library(SimGTF)

load(system.file("data", "b3_wang_sdist.Rdata", package = "SimGTF"))
fwang_sdist = wang_sdist
fwang_sdist$b_1 = sample(c(-1,1),nrow(fwang_sdist),replace=TRUE)*(abs(fwang_sdist$b_1) + log(2)) 
summary(fwang_sdist)
fwang_sdist = subset(fwang_sdist, b_0 > -15 & b_0 < -7)
summary(fwang_sdist)

load(system.file("data", "Ensembl_80_Homo_sapien_GTFSimData.Rdata", package="SimGTF")) 
ens80_sim_data = sim_data



max_transcripts = c(1,6,11)
i = 1
results = vector("list", length(max_transcripts))


for(i in seq_along(min_transcripts)){
	cat(i, "\n")
	ntrans = max_transcripts[i]
	de_prob = 1 - (1-1/2)^(1/ntrans)
	N1 = N2 = 3 #number of samples per treatment
	phex = PhenoExpression(
		trans_dist = fwang_sdist,
		sample_names = paste("count",1:(N1+N2),sep="_"),
		treatments = c(rep(1,N1),rep(2,N2)),
		E_lib_size = rep(6*10^6, N1+N2),
		de_prob = de_prob,
		min_transcripts = 1,
		max_transcripts = ntrans
	)
	
	sim = DisjointGenomeSimulation(ens80_sim_data, phex, seed=i)
	results[[i]] = GeneLocusResults(sim)
	sim = NULL
	gc()
}

save(results, file ="a3_demo_results.Rdata")

tpr = lapply(results, TPR)
fpr = lapply(results, FPR)
pdf("a3_Example_ROC.pdf")
plot.ROC(fpr,tpr, legend.labels=paste("# Transcripts", min_transcripts))
dev.off()
i = 1
plot(fpr[[i]], tpr[[i]], xlab="False Positive Rate", ylab="True Positive Rate", col=i, type="l", lwd=2)
for(i in 2:length(results))
	lines(fpr[[i]], tpr[[i]], col=i, lwd=2)

