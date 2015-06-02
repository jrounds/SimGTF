#############################################################################################################################################################################
# Jeremiah Rounds
# roundsjeremiah@gmail.com
#######################################################################################

#######################################################################################
# This script demonstrates processing a GTF into a data structure suitable to simulation
#######################################################################################
library(data.table)
wd ="/media/sf_work/15spring/working_SimGTF/SimGTF/examples"
setwd(wd)
source("../R/functions.R")
data_dir = "../data"

if(FALSE){
	#LONG OP NEEDS MEMORY RUN ONCE
	setwd("/noback/jrounds/SimGTF/examples")
	file = paste0(data_dir,"/Homo_sapiens.GRCh38.80.gtf")

	gtfdt = readGTF(file)

	gtf = gtfdt$full
	save(gtf, file = "gtf.Rdata")

	gtf_exons = getPCExons(gtf)
	sim_data = gtfToSimData(gtf_exons)
	 
	save(sim_data, file= "sim_data.Rdata")


}
load("sim_data.Rdata")

	library(data.table)
	library(GenomicRanges)
	

	#######################################################################################
	# sample distribution for genes note the +log(2)
	#######################################################################################

	load("../data/b3_wang_sdist.Rdata")
	fwang_sdist = wang_sdist
	fwang_sdist$b_1 = sample(c(-1,1),nrow(fwang_sdist),replace=TRUE)*(abs(fwang_sdist$b_1) + log(2)) 
	summary(fwang_sdist)
	fwang_sdist = subset(fwang_sdist, b_0 > -15 & b_0 < -7)
	summary(fwang_sdist)
	
	sim_data$gene_dist = fwang_sdist
	
	#######################################################################################
	# Sim genes
	#######################################################################################
	source("../R/functions.R")
			set.seed(1)

			n1 = n2 = 3
			treatments = c(rep(1,n1),rep(2,n2))
			E_lib_size = rep(6*10^6, length(treatments))
			#E_lib_coverages = rep(7, length(treatments))
			count_names = paste("count",1:length(treatments),sep="_")
			names(treatments) = names(E_lib_size) = count_names
			de_prob = .2
			min_transcripts = 1
			max_transcripts = 2
			group_transcripts_id = "gene_locus_id"
			zero_prob = 0
			
			
		base_sim = simDisjointGenomeNonParametric(
								sim_data,
								treatment = treatments, 
								E_lib_size = E_lib_size, 
								de_prob  = de_prob,
								min_transcripts = min_transcripts,
								max_transcripts = max_transcripts) 
		colSums(base_sim$trans[,count_names,with=FALSE])

		ls(base_sim)
		#trans is the first sim
		#disjoint is the exon level counts on the genome
		#trans_ex is trans with exons but the same exon may appear in multiple trans
		#for now though since max_transcripts = 1 it will be the case that in gene_locus_id 
		#disjoint_id is unique
		source("../R/functions.R")
		gene_counts = geneLocusCounts(base_sim$trans, count_names)
		disjoint_counts = geneLocusCounts(base_sim$trans_ex, count_names, c("gene_locus_id", "disjoint_id"))
		source("../R/functions.R")
		colSums(gene_counts[,count_names,with=FALSE])
		colSums(disjoint_counts[,count_names,with=FALSE])




		#gene$gene_total = rowSums(gene[,count_name,with=FALSE])
		gene_counts2 = geneLocusExactTest(gene_counts, treatments)
		disjoint_counts2 = disjointExactTest(disjoint_counts, treatments)
		
		
		
		
