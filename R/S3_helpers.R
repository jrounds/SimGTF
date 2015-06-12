PhenoExpression = function(sample_names, E_lib_size, treatments, trans_dist, de_prob=.20, error = 0, min_transcripts=1, max_transcripts=1){
	names(E_lib_size) = sample_names
	names(treatments) = sample_names
	this = environment()
	class(this) = c("PhenoExpression", class(this))
	return(this)
	


}

GTFSimData = function(gtf_exons, ignore.strand=TRUE){
	e = gtfToSimData(gtf_exons, ignore.strand)
	class(e) = c("GTFSimData", class(e))
	return(e)
}

DisjointGenomeSimulation = function(sim_data, phex, seed = NULL){
	if(!is.null(seed))
		set.seed(seed)
	e = simDisjointGenomeNonParametric(sim_data, phex)
	e$seed = seed
	class(e) = c("DisjointGenomeSimulation", class(e))
	return(e)
}
geneLocusCounts = function(this){
	UseMethod("geneLocusCounts", this)
}
geneLocusCounts.DisjointGenomeSimulation = function(base_sim){
	gene_counts = geneLocusCounts2(base_sim$trans, base_sim$pheno$sample_names)
	class(gene_counts) = c("geneLocusCounts", class(gene_counts))
	return(gene_counts)

}
disjointExonCounts = function(this){
	UseMethod("disjointExonCounts", this)

}
disjointExonCounts.DisjointGenomeSimulation = function(base_sim){
	disjoint_counts = geneLocusCounts2(base_sim$trans_ex, base_sim$pheno$sample_names, 
				c("gene_locus_id", "disjoint_id"))
	class(disjoint_counts) = c("disjointExonCounts", class(disjoint_counts))
	return(disjoint_counts)
	
}

GeneLocusResults = function(sim){
	gene_counts = geneLocusCounts(sim)
	gene_counts2 = geneLocusExactTest(gene_counts, phex$treatments)
	res = allTestMetrics(gene_counts2)
	
	e = new.env()
	e$counts = gene_counts2
	e$fdr_pvalue = "gene_locus_exact"
	e$raw_pvalue = "gene_locus_exact_raw"
	e$group = "gene_locus_id"
	e$results = res
	e$pheno = sim$pheno
	e$pheno$trans_dist = NULL
	class(e) = c("GTFSimResults", class(e))
	return(e)
	
}
DisjointExonResults = function(sim){
	disjoint_counts = disjointExonCounts(sim)
	disjoint_counts2 = disjointExactTest(disjoint_counts, phex$treatments)
	res = allTestMetrics(disjoint_counts2)
	
	e = new.env()
	e$counts = gene_counts2
	e$fdr_pvalue = "disjoint_exact"
	e$raw_pvalue = "disjoint_exact_raw"
	e$group = "disjoint_id"
	e$results = res
	e$pheno = sim$pheno
	e$pheno$trans_dist =  NULL
	class(e) = c("GTFSimResults", class(e))
	return(e)
	
}

counts = function(this){
	UseMethod("counts", this)
}
counts.GTFSimResults = function(e){
	 count_names = e$sim$pheno$sample_names
	 m = as.matrix(e$counts[,count_names, with=FALSE])
	 rownames(m) = e$counts[[e$group]]
	 return(m)

}
actuallyDE = function(this){
	UseMethod("actuallyDE", this)
}
actuallyDE.GTFSimResults = function(e){

	 m = e$counts[["actually_de"]]
	 names(m) = e$counts[[e$group]]
	 return(m)

}
extract = function(this, select){
	UseMethod("extract", this)
}
extract.GTFSimResults = function(e, select){
		#print(select)
		if(select == "pvalue")
			select = e$raw_pvalue
		if(select == "adjusted.pvalue")
			select = e$fdr_pvalue
		v = e$counts[[select]]
		names(v) = e$counts[[e$group]]
		return(v)
}

TPR = function(this){
	UseMethod("TPR", this)
}
TPR.GTFSimResults = function(e){
	a = e$results[[e$fdr_pvalue]]$alpha
	tpr = e$results[[e$fdr_pvalue]]$tpr
	names(tpr) = a
	return(tpr)
	
}

FPR = function(this){
	UseMethod("FPR", this)
}
FPR.GTFSimResults = function(e){
	a = e$results[[e$fdr_pvalue]]$alpha
	fpr = e$results[[e$fdr_pvalue]]$fpr
	names(fpr) = a
	return(fpr)
	
}
plot.GTFSimResults = function(e, type=""){
	if(type == "FDR"){
		plotFDR(e$results, which=e$fdr_pvalue)
		return()
	}
	if(type == "TPR"){
		plotTPR(e$results, which=e$fdr_pvalue)
		return()
	}
	if(type == "ROC"){
		plotROC(e$results, which=e$fdr_pvalue)
		return()
	}
	plotFDR(e$results, which=e$fdr_pvalue)
	plotTPR(e$results, which=e$fdr_pvalue)
	plotROC(e$results, which=e$fdr_pvalue)
	
}

FPRList = function(results){
	fpr = lapply(results, FPR)
	attr(fpr, "var_name") = "False Positive Rate"
	class(fpr) = c("MetricList", class(fpr))
	return(fpr)
}
TPRList = function(results){
	tpr = lapply(results, TPR)
	attr(tpr, "var_name") = "True Positive Rate"
	class(tpr) = c("MetricList", class(tpr))
	return(tpr)
}

plot.ROC = function(fpr, tpr, legend.labels=NULL){
	i = 1
	plot(fpr[[i]], tpr[[i]], xlab="False Positive Rate", ylab="True Positive Rate", col=i, type="l", lwd=2)
	for(i in 2:length(fpr))
		lines(fpr[[i]], tpr[[i]], col=i, lwd=2)
	if(!is.null(legend.labels)){
		legend("bottomright", legend.labels, col=1:length(fpr), lty=1, lwd=2)
	}

}





if(FALSE){  #I found this S3 example instructional...
setHasBreakfast <- function(elObjeto, newValue)
        {
                print("Calling the base setHasBreakfast function")
                UseMethod("setHasBreakfast",elObjeto)
                print("Note this is not executed!")
        }

setHasBreakfast.default <- function(elObjeto, newValue)
        {
                print("You screwed up. I do not know how to handle this object.")
                return(elObjeto)
        }


setHasBreakfast.NorthAmerican <- function(elObjeto, newValue)
        {
                print("In setHasBreakfast.NorthAmerican and setting the value")
                elObjeto$hasBreakfast <- newValue
                return(elObjeto)
        }
	}