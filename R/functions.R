#' setTerm
#'
#' @param color what color to change a linux terminal to
setTerm = function(color){
	system(paste("setterm -term linux -background", color))
}

upSimGTF = function(){
	system("scp -rp /media/sf_work/15spring/working_SimGTF/SimGTF falconer.stat.purdue.edu:/noback/jrounds/SimGTF")
	
}
#' progressFunction
#'
#' @param nrow total length of progress 
#' @param foo function that will be procssed during progress
progressFunction = function(nrow, foo){
	prog = txtProgressBar(0,nrow, style=3)
	prog_count = 0
	foo2 = function(sd){
		prog_count <<- prog_count + nrow(sd)
		setTxtProgressBar(prog, prog_count)
		return(foo(sd))
	}
	return(foo2)
}
#' readGTF
#' 
#' @param file source file to read (e.g. unzipped Ensembl gtf
#' @param ... other arguments to read.table
readGTF = function(file,...){
	x = read.table(file, stringsAsFactors=FALSE, sep="\t", ...)
	cat("Read complete... \nProcessing...\n")
	c9 = x[,9]
	c9 = gsub("; ",";", c9,fixed=TRUE)
	c9 = strsplit(c9,";",fixed=TRUE)
	hmm = lapply(c9,function(c) strsplit(c, " ",fixed=TRUE))
	f = unlist(lapply(hmm, function(h) 
		lapply(h, function(h0){
			if(length(h0) > 0)
				h0[[1]]
			})
	))
	fields = unique(f)
	blank = as.list(rep(NA, length(fields)))
	names(blank) = fields
	reprocess = lapply(hmm, function(h){
		ret = blank
		for(h0 in h){
			if(length(h0) > 1)
				ret[[h0[[1]]]] = h0[[2]]
		}
		return(ret)

	})
	data = rbindlist(reprocess)
	colnames(x) = c("seqname","source","feature","start","end","score","strand","frame","attribute")
	full = cbind(x[,1:8], data)
	full = as.data.table(full)
	#full
	#return(as.list(environment()))
	return(full)
}



#######################################################################################
# Functions for building requisite data structures from a gtf
#######################################################################################


#######################################################################################
# processing into subexons
#######################################################################################
#' fixTypes
#'
#' @param d data.table on which to fix up types to a standard representatation for genomic data
fixTypes = function(d){
	if(class(d$strand) != "integer" && class(d$strand) != "numeric"){
		s = integer(length=nrow(d))
		m = d$strand == "+"
		s[m] = 1
		m = d$strand == "-"
		s[m] = -1
		d[, strand := as.integer(s)]
	}
	n = d$seqnames
	d$seqnames = NULL
	d$seqnames = as.character(n)
	setkeyv(d, c("seqnames","strand","start"))
	return(d)
}
#' as.GRanges
#'
#' @param dt data.table to turn into a GRanges object
as.GRanges = function(dt){
	GRanges(seqname=dt$seqname, strand=dt$strand, 
		ranges=IRanges(start=dt$start, end=dt$end))
}
#' disjoinExon
#' 
#' @param dt disjoint_transcripts from a gtf_data object
#' @param ignore.strand should disjoint counting bins be caculated by ignoring strand (as ooposed to over-lapping across strands?)
#' @param id_lead leading string for a disjoint exon key field.
#' @return disjoint genome exon data.table
disjoinExons = function(dt, ignore.strand=FALSE, id_lead="dis"){
	gr_exons = as.GRanges(dt)
	gr_exons = disjoin(gr_exons, ignore.strand)
	dt_exons = as.data.table(as.data.frame(gr_exons))
	dt_exons$seqnames = as.character(dt_exons$seqnames)
	setkeyv(dt_exons,c("seqnames","strand","start","end"))
	dt_exons = fixTypes(dt_exons)
	dt_exons$disjoint_id = paste(id_lead,1:nrow(dt_exons),sep="_")
	dt_exons$source = "disjoinExons"
	dt_exons$feature = "disjoint_exon"
	if(ignore.strand)
		dt_exons$strand = NA
	setkeyv(dt_exons,"disjoint_id")
	class(dt_exons) = c("Disjoint_Exons", class(dt_exons)) 
	dt_exons
}


#' createDisjointMap
#' 
#' @param dt Object of class "Transcripts"
#' @param disjoint Object from disjoinExons()
#' @param ignore.strand was disjoint created ignoring strand?
#' @seealso disjoinExons
#' @return a map between the transcript_id of dt to  to the disjoint exons of disjoint
createDisjointMap = function(dt, disjoint, ignore.strand=FALSE){
	setkey(dt, exon_id)
	udt = unique(dt)
	gr_exons = as.GRanges(disjoint)
	gr_gene = as.GRanges(udt)
	o = findOverlaps(gr_gene, gr_exons, ignore.strand=ignore.strand)
	gene_sub_map = as.data.table(as.data.frame(o))
	setnames(gene_sub_map,c("exon_row", "disjoint_row"))
	gene_sub_map
	gene_sub_map$exon_id = udt$exon_id[gene_sub_map$exon_row]
	gene_sub_map$disjoint_id = disjoint$disjoint_id[gene_sub_map$disjoint_row]
	gene_sub_map$exon_row = gene_sub_map$disjoint_row = NULL
	setkey(gene_sub_map, exon_id)
	#class(gene_sub_map) = c("Map_To_Disjoint", class(gene_sub_map))
	return(gene_sub_map)

}
#' createDisjointTranscripts
#'
#' @param dt Object of class "Transcripts"
#' @param disjoint Object from disjoinExons()
#' @param map Object from createDisjointMap()
#' @seealso disjointExons, createDisjointMap
#' @return exons in dt turned into sub-exon disjoint counting bins organized by transcript_id 
createDisjointTranscripts = function(dt, disjoint, map){
	setkey(map, exon_id)
	setkey(disjoint, disjoint_id)
	setkey(dt, exon_id)
	dtt = dt[map, allow.cartesian=TRUE]
	dtt$start = NULL
	dtt$end = NULL
	dtt$strand = NULL
	dtt$seqname = NULL
	dtt$source = NULL
	dtt$feature = NULL
	setkey(dtt, disjoint_id)
	dtt = dtt[disjoint]
	#class(dtt) = c("Disjoint_Transcripts", class(dtt))
	dtt
	
	
}


#dt is from createDisjointTranscripts

#' createDisjointGenes
#'
#' @param dt Disjoint_Transcripts object from createDisjointTranscripts
#' @return reduces a gene_id to its disjoint exon bins
createDisjointGenes = function(dt){
	setkeyv(dt, c("gene_id","disjoint_id"))
	dtg = unique(dt)
	drop = grep("transcript",names(dtg), value=TRUE)
	dtg[,drop:=NULL,with=FALSE]
	return(dtg)
	
	
}

#dg is from createDisjointGenes

#'  createDisjointGeneLoci
#'
#' @param dg Object from createDisjointGenes
#' @return finds overlapping gene regions and merges them into a gene loci.
createDisjointGeneLoci = function(dg){
	#genes have to be merged into non-overlapping loci
	
	map = dg[,c("gene_id","disjoint_id"),with=FALSE]
	setkey(map, "disjoint_id")
	
	gid = unique(dg$gene_id)
	gid = as.list(gid)
	names(gid) = unlist(gid)
	op = progressFunction(nrow(dg),function(.SD){
		did = .SD$disjoint_id
		thisgid = unique(.SD$gene_id)
		gid[[thisgid]] <<- unique(map[did]$gene_id)
		return(NULL)
	})
	dg[,op(.SD), by=gene_id, .SDcols=c("disjoint_id", "gene_id")]
	cat("\n")
	lengths = unlist(lapply(gid, length))
	table(lengths)
	lengths = sort(lengths,dec=TRUE)
	
	#not quite done yet
	for(i in seq_along(lengths)){
		thisgid = names(lengths)[i]
		thisgid = gid[[thisgid]]
		#if(length(thisgid) == 7)
		#	browser()
		for(g in thisgid)
			gid[[g]] = unique(c(gid[[g]], thisgid))
	
	}
	lengths = unlist(lapply(gid, length))
	table(lengths)
	lengths = sort(lengths,dec=TRUE)
	head(lengths)
	gene_locus = unlist(lapply(gid, function(g) paste(sort(g),collapse="+")))
	names(gene_locus) = names(gid)
	gene_locus_names = sort(unique(unlist(gene_locus)))
	gene_locus_id = paste("gene_locus",1:length(gene_locus_names),sep="_")
	names(gene_locus_id) = gene_locus_names
	
	
	mapped = gene_locus_id[gene_locus[dg$gene_id]]
	names(mapped) = NULL
	dg$gene_locus_id = mapped
	

	return(dg)
	
}
#' appendGeneLocus
#'
#' @param from data.table with a gene_locus_id field
#' @param to data.table to have a gene_locus_id_field
#' @return The _to_ input with a gene_locus_id field
appendGeneLocus = function(from,to){
	dg = from
	dt = to
	setkey(dt, gene_id)
	setkey(dg, gene_id)
	udg = unique(dg)
	dt$gene_locus_id = udg[dt,gene_locus_id]
	return(dt)
}
#' getPCExons
#'
#' @param gtf Transcripts data object
#' @param keep_seqnames which chromosomes to keep in the gtf
#' @return gtf with only protein coding exon data for keep_seqname chromosomes.
getPCExons = function(gtf, keep_seqnames = c(1:22,"X","Y")){
	exons = subset(gtf, feature=="exon" & 
		gene_biotype == "protein_coding" & 
			(seqname %in% keep_seqnames)
			)
	return(exons)
}
#' gtfToSimData
#' @param gtf_exons Transcript data object from gtf
#' @param ignore.strand Should disjoint simulation data be constructed by ignoring strand?
#' @return all objects necessary to simulate disjoint exons
gtfToSimData = function(gtf_exons, ignore.strand=TRUE){
	require(data.table)
	require(GenomicRanges)
	gtf_exons = copy(gtf_exons)
	



	disjoint = disjoinExons(gtf_exons, ignore.strand,"pc_d")
	disjoint$source = "SimGTFScript"

	disjoint_map= createDisjointMap(gtf_exons, disjoint, ignore.strand)
	#######################################################################################
	#and finally a convenience data structure for transcripts but this is modified at the end
	#######################################################################################
	disjoint_transcripts = createDisjointTranscripts(gtf_exons,disjoint, disjoint_map)



	#######################################################################################
	# convenience data structure for genes
	#######################################################################################


	disjoint_genes = createDisjointGenes(disjoint_transcripts)


	length(unique(disjoint_genes$exon_id))


	#shared disjoint_id between genes?
	ngenes = disjoint_genes[,list(ngenes = length(unique(gene_id))), by=disjoint_id]
	ngenes = ngenes$ngenes
	summary(ngenes)
	table(ngenes)

	disjoint_genes = createDisjointGeneLoci(disjoint_genes)
	#save(disjoint_genes, file= "a2_Drosophilia_disjoint_genes.Rdata")

	length(unique(disjoint_genes$gene_locus_id))

	length(unique(disjoint_genes$gene_id))


	#######################################################################################
	# For simulation truth we need a disjoint_genes and transcripts without overlap
	#######################################################################################

	disjoint_transcripts = appendGeneLocus(disjoint_genes, disjoint_transcripts)
	return(environment())
	
	


}


#######################################################################################
# Parse GTF Functions ^^^^^^^^^^^^^^
# Now starting
# Sim data functions
#######################################################################################

simTranscriptSwitchingCounts = function(trans, treatments, lib_sizes){

	n = names(treatments)
	treatmentsp = paste("prop",treatments,sep="_")
	names(treatmentsp) = n

	require(MASS)
	ntrans = nrow(trans)
	for(c in names(treatmentsp)){
		E = trans[[treatmentsp[c]]]*lib_sizes[c]
		trans[[c]] = rnegbin(ntrans,E,1/trans$dispersion)
	}
	
	return(trans)
}
simDisjointTranscriptSwitchingCounts = function(dt,trans, counts, treatments, verbose=TRUE){
	n = colnames(dt)
	tr1_trans = trans
	
	op = function(tr= 1, var="tr1_transcript_id"){
			setkey(dt, "transcript_id")
			#ignore var names
			tr2_trans = copy(trans)
			setkeyv(tr2_trans,var)
			tr2_dt = dt[tr2_trans]
			n = colnames(tr2_dt)
			drop = grep("^i.",n,value=TRUE)
			tr2_dt[,drop:= NULL,with=FALSE]  #dtt now has a repeated obs for each disjoint exon
			zero = counts[treatments == tr]
			tr2_dt[, zero:= 0, with=FALSE]
			return(tr2_dt)
	}
	
	tr1_dt = op(1,"tr1_transcript_id")
	tr2_dt = op(2,"tr2_transcript_id")
	
	#######################################################################################
	# recombine and manip
	#######################################################################################
	dtt = rbind(tr1_dt, tr2_dt, fill=TRUE)
	i = which(! colnames(tr1_dt) %in% colnames(tr2_dt))
		
	if(verbose){
		cat("Before allocation...\n")
		print(colSums(trans[,counts,with=FALSE]))
	}
	colSums(tr1_dt[,counts, with=FALSE])
	colSums(tr2_dt[,counts, with=FALSE])
	
	dtt[, transcript_width := sum(width), by = transcript_id]
	bquote({keep = .(names(dtt))})
	 keep = c("gene_id",  "gene_name",  "transcript_id", "exon_number", "exon_id", 
    "disjoint_id", "width", "gene_locus_id", "dispersion","prop_1", 
	"prop_2", "actually_de", "transcript_width")

	#######################################################################################
	# need to compress treatment 1 and treatment 2 in the case of not differentially expressed
	#######################################################################################
	op = progressFunction(nrow(dtt), function(.SD){
		ret = colSums(.SD[,counts,with=FALSE])
		names(ret) = counts
		ret = as.list(ret)
		for(c0 in keep){
			uc0 = unique(.SD[[c0]])
			if(length(uc0) == 0)
				next
			if(length(uc0) == 1){
				ret[[c0]] = uc0
			}else{
				browser()
			}
		}
		
		return(as.list(ret))
	})
	cat("Compressing disjoint transcript regions...\n")
	dtt0 = dtt[, op(.SD), by=c("disjoint_id","transcript_id")]
	cat("\n")

	#######################################################################################
	# There is something really f'd going on with data.table so I am hacking in
	# a loop over the transcripts
	setkey(dtt0, transcript_id)
	allt = unique(dtt0$transcript_id)
	t0 = "ENST00000000233"
	work = lapply(allt, function(t0){
		sd = copy(dtt0[t0])
		p = sd$width/unique(sd$transcript_width)
		
		for(c in counts){
			
			sd[[c]]= as.numeric(rmultinom(1,sd[[c]][1], p))
		}
		return(sd)
	
	})
	dtt1 = rbindlist(work)
	
	
	#######################################################################################
	# Have to fix up transcripts for treatments
	#######################################################################################

	
	if(verbose){
		cat("After allocation...\n")
		print(colSums(dtt1[,counts,with=FALSE]))
	}
	return(dtt1)
}

simWithForceDESwitching = function(sim_data,  treatments, E_lib_sizes,	
				de_prob=.20, min_transcripts=1, max_transcripts=Inf,  force_de = NULL,only_sim_trans=FALSE)
				{
	if(any(is.null(names(treatments))) || any(names(treatments) != names(E_lib_sizes))){
		stop("Count names should be in treatments and lib_sizes.")
	}
	dt = copy(sim_data$disjoint_transcripts)
	sdist = as.data.frame(copy(sim_data$gene_dist))
	
	fragment_size = 150  #not really used but for back of envelope calculations
	count_names = names(treatments)
	zero_prob = 0

	#######################################################################################
	# I need to filter dt to only have genes with 2 or more transcripts
	#######################################################################################
	op = function(.SD){
		if(length(unique(.SD$transcript_id)) == 1)
			return(NULL)
		return(.SD)
	}
	fdt = dt[,op(.SD),by=gene_locus_id]
	#if(!all(force_de %in% unique(fdt[[group_transcripts_id]])))
	#	stop("force_de has a gene with 1 transcript.")
		
	#######################################################################################
	# We need to know how big the transcripts are in order to make sure they have a proper
	# size differentially
	#######################################################################################
	w = dt[,list(width = sum(width) ), by = c("gene_locus_id","transcript_id")]
	
	#######################################################################################
	# I need treatment 1 transcripts
	#######################################################################################
	opmin = function(.SD){
		avail = which.min(.SD$width)
		avail = sample(avail, 1)
		#keep = sample(avail, 1)
		return(.SD[avail])
	}
	t1 = w[, opmin(.SD), by=gene_locus_id]
	ut1 = unique(t1$transcript_id)
	#w2 = subset(w, !(w$transcript_id %in% ut1))
	w2 = w
	opmax = function(.SD){
		avail = which.max(.SD$width)
		avail = sample(avail, 1)
		#keep = sample(avail, 1)
		return(.SD[avail])
	}
	t2 = w2[,opmax(.SD), by=gene_locus_id]
	all(sort(unique(t2[[group_transcripts_id]])) == sort( unique(t1[[group_transcripts_id]])))

	#######################################################################################
	# unique 
	#######################################################################################
	setkeyv(t1, group_transcripts_id)
	setkeyv(t2, group_transcripts_id)
	ut1 = unique(t1)
	ut2 = unique(t2)
	
	#######################################################################################
	# I need to assign these each parameters
	#######################################################################################
	ut1row = sample(1:nrow(sdist), nrow(ut1), replace=TRUE)
	ut2row = sample(1:nrow(sdist), nrow(ut2), replace=TRUE)
	
	ut1_2 = cbind(ut1, sdist[ut1row,])
	ut2_2 = cbind(ut2, sdist[ut1row,])
	setkeyv(ut1_2, group_transcripts_id)
	setkeyv(ut2_2, group_transcripts_id)
	setkey(t1, transcript_id)
	setkey(t2, transcript_id)
	g1 = ut1_2$transcript_id[1]
	g2 = ut2_2$transcript_id[1]
	t1[g1]
	t2[g2]
	#######################################################################################
	# I need to assign differential expression
	#######################################################################################
	de = sample(c(FALSE,TRUE), nrow(ut1_2), replace=TRUE,prob=c(1-de_prob, de_prob))
	names(de) = ut1_2[[group_transcripts_id]]
	de[force_de] = TRUE
	p =  exp(ut1_2$b_0)
	p = p*ut1_2$width/sum(p*ut1_2$width)
	ut1_2$prop_1 = p
	ut1_2$prop_2 = p
	p =  exp(ut1_2$b_0)
	p = p*ut2_2$width/sum(p*ut2_2$width)
	ut2_2$prop_1 = p
	ut2_2$prop_2 = p
	ut1_2$actually_de = de
	
	#######################################################################################
	# ut1 divides into 2 cases.  One case the transcript gets used twice
	# another case the transcript gets used once and changes in another treatment
	#######################################################################################
	ut1_twice = subset(ut1_2, !actually_de) #use twice
	ut1_twice$tr2_transcript_id = ut1_twice$transcript_id
	ut1_twice$tr1_transcript_id = ut1_twice$transcript_id
	ut1_once = subset(ut1_2, actually_de)
	ut1_once$prop_2 =  NULL
	
	ut2_once = subset(ut2_2, de)
	#the two onces need to be joined
	ut2_once$prop_1 = NULL
	setkeyv(ut1_once, group_transcripts_id)
	setkeyv(ut2_once, group_transcripts_id)
	
	utde = ut1_once[ut2_once]
	utde$tr2_transcript_id = utde$i.transcript_id
	utde$tr1_transcript_id = utde$transcript_id
	n = grep("i.",colnames(utde), fixed=TRUE,invert=TRUE)
	utde = as.data.table(as.data.frame(utde)[,n])
	
	expr_t = rbind(ut1_twice, utde)
	sum(expr_t$prop_1)
	sum(expr_t$prop_2)
	
	ret = list()
	ret$force_de = force_de
	ret$de_prob = de_prob
	trans = simTranscriptSwitchingCounts(expr_t, treatments, E_lib_sizes)	
	ret$trans = trans
	if(only_sim_trans)
		return(as.environment(ret))
	trans_ex = simDisjointTranscriptSwitchingCounts(dt,trans,count_names,treatments)
	disjoint = extractDisjointGenomeSwitchingCounts(trans_ex, count_names)
	
	


	ret$trans_ex = trans_ex
	ret$disjoint = disjoint
	return(as.environment(ret))
		
				
}
 

 
 
appendSimInfoMixture = function(e){
	return(appendSimInfo(e))
}
reSimCountsFromTrans = function(original){
	#original inputs
	treatment = original$treatment
	E_lib_size = original$E_lib_size	
	de_prob= original$de_prob
	min_transcripts  = original$min_transcripts
	max_transcripts = original$max_transcripts
	group_transcripts_id = original$group_transcripts_id
	zero_prob = original$zero_prob
	fragment_size = original$fragment_size
	only_sim_trans = original$only_sim_trans
	mixture = original$mixture
	trans = as.list(as.data.frame(original$trans))
	trans[[1]] = original$trans[[1]] # tricking R into copying this dt
	trans = as.data.table(as.data.frame(trans))
	ntrans = length(trans$transcript_id) #shouldnt change
	total_transcript_width = original$total_transcript_width  # shouldnt change
	
	#caclualted values
	mixture$prior = mixture$prior/sum(mixture$prior)
	count_name = names(treatment)
	E_coverage = E_lib_size*fragment_size/total_transcript_width
	names(E_coverage) = names(E_lib_size)

	#new simulations
	trans = simTranscriptCounts(trans, treatment, E_lib_size)
	if(!only_sim_trans){
		trans_ex = simDisjointTranscriptCounts(dt,trans,count_name)
		disjoint = extractDisjointGenomeCounts(trans_ex, count_name)
	}
	#confirm trans is the same except for counts
	e = environment()
	class(e) = c(class(e), "transcript_simulation")
	return(appendSimInfo(e))
}


#######################################################################################
# The 3rd simulation strat!  
# A non-parametric simulation
#######################################################################################
extractTranscriptData = function(dt,min_transcripts, max_transcripts, group_transcripts_id,zero_prob, force_de = NULL){
	#transcript data extraction
	by =  unique(c(group_transcripts_id,"gene_id","gene_name",
							"transcript_id","transcript_name"))
	op = function(.SD){
	list(seqnames=unique(.SD$seqnames), strand=unique(.SD$strand),
				min.start=min(.SD$start),max.end=max(.SD$end),width=sum(.SD$width))
	
	}
	trans = dt[,op(.SD), by =by] 
	
	#dropping group_transcript with zero_prob
	keep_names = unique(trans[[group_transcripts_id]])
	keep = runif(length(keep_names)) > zero_prob
	keep = keep_names[keep]
	setkeyv(trans, group_transcripts_id)
	trans = trans[keep]
	
	#making transcripts fit within in min and max numbers
	op = function(.SD){

		tid = .SD$transcript_id
		nt = length(tid)
		ret =list()
		ret$ntranscripts = nt
		ret$transcript_id = tid
		ret$nkept_transcripts = length(tid)
		if(nt < min_transcripts)
			return(NULL)
		if(nt == min_transcripts)
			return(ret)
		#now in the nt > min_transcripts case
		if(nt <= max_transcripts)
			return(ret)
		
		#now in the nt > min_transcripts and nt > max_transcripts
		ret$transcript_id = sample(tid, max_transcripts)
		ret$nkept_transcripts = length(ret$transcript_id)
		return( ret)
		
	
	}
	keep = trans[,op(.SD), by=group_transcripts_id, .SDcols=c("transcript_id")]
	setkeyv(trans, c(group_transcripts_id,"transcript_id"))
	setkeyv(keep, c(group_transcripts_id, "transcript_id"))
	trans2 = trans[keep]

	return(trans2)
}
simTranscriptParamsNonParametric = function(trans, sdist, de_prob, force_de = NULL){
	i1 = sample(1:nrow(sdist), nrow(trans), replace=TRUE)
	sparam = sdist[i1,]
	de = sample(c(FALSE,TRUE), nrow(trans), replace=TRUE, prob =c(1-de_prob,de_prob))

	if(!is.null(force_de)){	
		tid = trans$transcript_id
		gid = trans$gene_locus_id
		names(tid) = gid
		names(de) = tid
		de[tid[force_de]] = TRUE
	
	}
	sparam$actually_de = de
	
	#the non-de genes need to be assigned one of these values at random in b_0
	# I don't know why I picked "t" this is really "pi"
	t1 = exp(sparam$b_0)
	t2 = exp(sparam$b_0 + sparam$b_1)
	t = t1
	swap = runif(length(t)) >= .5
	t[swap] = t2[swap]
	de = sparam$actually_de
	t1[!de] = t2[!de] = t[!de]
	#these now need to be normalized to sum to 1 in each treatment
	t1 = t1/sum(t1)
	prop_de = sum(t1[de])
	prop_nde = sum(t1[!de])
	#t2 needs to be made to sum to 1 but the not de need to continue to have the same total proportion
	t2[de] = t2[de]/sum(t2[de])*prop_de
	t2[!de] = t2[!de]/sum(t2[!de])*prop_nde
	sum(t1)
	sum(t2)

	new_b_0 = log(t1)
	new_b_1 = log(t2) - new_b_0

	sparam$b_0 = new_b_0
	sparam$b_1 = new_b_1
	sparam$prop_1 = t1
	sparam$prop_2 = t2
	trans1 = cbind(trans,sparam)
	#pairs(sparam)
	head(trans1)
	return(trans1)


}

#based on simTranscriptNonParametric
simTranscriptParamsWithForceDE = function(trans, sdist, de_prob, force_de = NULL){
	i1 = sample(1:nrow(sdist), nrow(trans), replace=TRUE)
	sparam = sdist[i1,]
	sparam$actually_de = sample(c(FALSE,TRUE), nrow(trans), replace=TRUE, prob =c(1-de_prob,de_prob))
	if(!is.null(force_de)){
		m = trans$gene_id %in% force_de
		sparam$actually_de[m] = TRUE
	}
	#the non-de genes need to be assigned one of these values at random in b_0
	# I don't know why I picked "t" this is really "pi"
	t1 = exp(sparam$b_0)
	t2 = exp(sparam$b_0 + sparam$b_1)
	t = t1
	swap = runif(length(t)) >= .5
	t[swap] = t2[swap]
	de = sparam$actually_de
	t1[!de] = t2[!de] = t[!de]
	#these now need to be normalized to sum to 1 in each treatment
	t1 = t1/sum(t1)
	prop_de = sum(t1[de])
	prop_nde = sum(t1[!de])
	#t2 needs to be made to sum to 1 but the not de need to continue to have the same total proportion
	t2[de] = t2[de]/sum(t2[de])*prop_de
	t2[!de] = t2[!de]/sum(t2[!de])*prop_nde
	sum(t1)
	sum(t2)

	new_b_0 = log(t1)
	new_b_1 = log(t2) - new_b_0

	sparam$b_0 = new_b_0
	sparam$b_1 = new_b_1
	sparam$prop_1 = t1
	sparam$prop_2 = t2
	trans1 = cbind(trans,sparam)
	#pairs(sparam)
	head(trans1)
	return(trans1)


}
simTranscriptCounts = function(trans, treatments, lib_sizes){
	if(all(unique(treatments) %in% c("1","2")))	{
		n = names(treatments)
		treatments = paste("prop",treatments,sep="_")
		names(treatments) = n
	} else {
		stop("treatments must be 1 or 2")
	}
	require(MASS)
	ntrans = nrow(trans)
	for(c in names(treatments)){
		E = trans[[treatments[c]]]*lib_sizes[c]
		trans[[c]] = rnegbin(ntrans,E,1/trans$dispersion)
	}
	return(trans)
}
simDisjointTranscriptCounts = function(dt,trans, counts,  verbose=TRUE){
	n = colnames(dt)
	setkey(trans, transcript_id)
	setkey(dt, transcript_id)
	if(verbose){
		cat("Before allocation...\n")
		print(colSums(trans[,counts,with=FALSE]))
	}
	#browser()
	dtt = dt[trans]
	n = names(dtt)
	drop = grep("^i.",n,value=TRUE)
	dtt[,drop:= NULL,with=FALSE]  #dtt now has a repeated obs for each disjoint exon
	foo = progressFunction(nrow(dtt),function(.SD){
		ret = list()
		p = .SD$width/sum(.SD$width)
		
		for(c in counts){
			
			ret[[c]] = as.integer(rmultinom(1, .SD[[c]][1], p))
		}
		return(ret)
	})
	dtt[, counts := foo(.SD), by="transcript_id", with = FALSE, .SDcols=c(counts, "width")]
	cat("\n")
	
	#######################################################################################
	#
	if(verbose){
		cat("After allocation...\n")
		print(colSums(dtt[,counts,with=FALSE]))
	}
	return(dtt)
}

extractDisjointGenomeCounts = function(trans_ex,counts,verbose=TRUE){
	if(verbose){
		cat("Before genome feature counting...\n")
		print(colSums(trans_ex[,counts,with=FALSE]))
	}
	keep = c("width")
	op = progressFunction(nrow(trans_ex),function(.SD){
		ret = list()		
		for(k in keep)
			ret[[k]] = unique(.SD[[k]])
			#ret[[k]] = .SD[[k]][1]
			#
		for(c in counts)
			ret[[c]] = sum(.SD[[c]])
		ret$actually_de = any(.SD$actually_de)

		

		return(ret)
	})
	genome = trans_ex[,op(.SD) , by="disjoint_id", .SDcols = c(keep,counts,"actually_de")]
	cat("\n")
	if(verbose){
		cat("After genome feature counting...\n")
		print(colSums(genome[,counts,with=FALSE]))
	}

	return(genome)
}
 appendSimInfo = function(sim){

	evalq({
		#reports
		phenotype = list()
		phenotype$name = count_name
		phenotype$treatment = treatment	
		phenotype$E_lib_size = E_lib_size
		phenotype$lib_size = colSums(trans[,count_name,with=FALSE])	
		phenotype$E_transcript_coverage = E_coverage
		phenotype$transcript_coverage = phenotype$lib_size*fragment_size/total_transcript_width
		phenotype = as.data.table(as.data.frame(phenotype,stringsAsFactors=FALSE))
		phenotype
		phenotype$name = as.character(count_name)
		f = c( "min_transcripts","max_transcripts", "group_transcripts_id", 
	 "fragment_size", "total_transcript_width")
		param = lapply(f, get, envir=environment())
		names(param) = f
		param = as.data.table(as.data.frame(param))
		m = names(phenotype)[c(-1,-2)]
		for(m0 in m) param[[paste("mean",m0,sep="_")]] = mean(phenotype[[m0]])
		param
		
	}, env = sim)
	return(sim)
}
geneLocusCounts = function(trans, count_names, by=c("gene_locus_id")){
	op = progressFunction(nrow(trans),function(.SD){
		ret = list()
		ret$actually_de = any(.SD$actually_de)
		for(c0 in count_names){
			ret[[c0]] = sum(.SD[[c0]])
		}
		ret$nexpressed_transcripts = length(unique(.SD$transcript_id))
		ret$nkept_transcripts = unique(.SD$nkept_transcripts)
		return(ret)
	})
	ret = trans[,op(.SD), by=by]
	c("\n")
	return(ret)
}
#' featureExactTest
#'
#' @param gene_counts data.table of feature counts
#' @param treatments vector indicating sample groups.  names(treatments) should be feature counts in gene_counts. 
#' @param field_disp Use this field as tag.wise dispersion instead of asking edgeR for it.  NULL indicates use edgeR.
featureExactTest = function(counts,treatments, field_disp=NULL, has_sim_truth=TRUE, feature_type ="gene"){
	require(edgeR)
	field_id = paste(feature_type, "_", "id",sep="")
	ret_raw = paste(feature_type, "_exact_raw",sep="")
	ret_bh = paste(feature_type, "_exact",sep="")
	ret_tested = paste(feature_type, "_exact_tested", sep="")
	ret_total = paste(feature_type, "_total", sep="")
	#######################################################################################
	# How does edgeR do?
	#######################################################################################
	count_names = names(treatments)
	y_counts = as.matrix(copy(counts[, count_names, with=FALSE]))

	row = rowSums(y_counts)
	m = row  > 0
	y = DGEList(y_counts[m,], group=treatments)
	y = calcNormFactors(y)
	
	if(!is.null(field_disp)){
		disp = gene_counts[[field_disp]][m]
		y$tagwise.dispersion = disp
	}else{
		y = estimateDisp(y, prior.df=10)
		disp = y$tagwise.dispersion
	}
	e = exactTest(y)
	p = e@.Data[[1]]$PValue
	q = p.adjust(p,"BH")
	discovery = q < .05
	#results = list()
	#results[[field_id]]= gene_counts[[field_id]]
	ngene = nrow(counts)
	if(has_sim_truth){
		#results$diff_exp = gene_counts$diff_exp
		#results$evidently_diff_exp = gene_counts$evidently_diff_exp
		#results$actually_de = results$diff_exp
	}
	counts[[ret_total]] = row
	counts[[ret_tested]] = m
	oounts[[ret_raw]] = gene_counts[[ret_bh]] = rep(1,ngene)
	x = counts[[ret_bh]]
	x[m] = q
	counts[[ret_bh]] = x
	x = counts[[ret_raw]]
	x[m] = p
	counts[[ret_raw]] = x
	if(is.null(field_disp)){
		x = rep(as.numeric(NA), ngene)
		x[m] = disp
		counts[["tagwise.dispersion"]]= x
		
		
	}
	
	return(counts)

}
geneExactTest = function(..., feature_type="gene") featureExactTest(...,feature_type=feature_type)
disjointExactTest = function(...,feature_type="disjoint") featureExactTest(...,feature_type=feature_type)
geneLocusExactTest = function(..., feature_type="gene_locus") featureExactTest(..., feature_type=feature_type)

#' simDisjointGenomeNonParametric
#'
#' @param gtf_data Object from gtfToSimData
#' @param gene_dist data.frame from which to sample with replacement b_0, b_1, and dispersion triplets.   Column names should be those.
#' @param treatment vector of numeric 1 and 2 indicating sample assignments.
#' @param E_lib_size vector indicating expected library size
#' @param de_prob probability that a simulated transcript will be differentially expressed between treatments
#' @param min_transcripts minimum number of transcripts a gene_locus_id should have to be included in the simulation
#' @param max_transcripts maximum number of transcripts allowed per gene_locus_id (excess transcripts are removed and they are sampled at random for simulation).
#' @param group_transcripts_id what field should transcripts be grouped by for min_transcripts/max_transcripts
#' @param only_sim_transcripts Do not simulate the full disjoint exon/genome?
#' @param force_de set of gene names to force into the differentially expressed category
simDisjointGenomeNonParametric = function(gtf_data, gene_dist, treatment, E_lib_size,	
				de_prob=.50, min_transcripts=1, max_transcripts=Inf, 
				group_transcripts_id = "gene_locus_id", only_sim_trans=FALSE, force_de = NULL){
		
	if(any(is.null(names(treatment))) || any(names(treatment) != names(E_lib_size))){
		stop("Count names should be in treatments and lib_sizes.")
	}
	dt = gtf_data$disjoint_transcripts
	sdist= gene_dist
	
	fragment_size = 150  #not really used but for back of envelope calculations
	count_name = names(treatment)
	zero_prob = 0
	trans = extractTranscriptData(dt, min_transcripts, max_transcripts,
			group_transcripts_id, zero_prob)
	total_transcript_width = sum(trans$width)


	E_coverage = E_lib_size*fragment_size/total_transcript_width
	names(E_coverage) = names(E_lib_size)
	
	ntrans = length(trans$transcript_id)
	
	
	#######################################################################################
	# simulation
	#######################################################################################
	trans = simTranscriptParamsNonParametric(trans,sdist, de_prob, force_de)
	trans = simTranscriptCounts(trans, treatment, E_lib_size)

	if(!only_sim_trans){
		trans_ex = simDisjointTranscriptCounts(dt,trans,count_name)
		disjoint = extractDisjointGenomeCounts(trans_ex, count_name)
	}
	
	

	
	e = environment()
	class(e) = c(class(e), "transcript_simulation")
	return(appendSimInfo(e))				
				
}

#based on
simWithForceDE = function(force_de,dt, sdist, treatment, E_lib_size,	
				de_prob=.20, min_transcripts=1, max_transcripts=Inf, 
				group_transcripts_id = "gene_locus_id", only_sim_trans=FALSE)
				{
	if(any(is.null(names(treatment))) || any(names(treatment) != names(E_lib_size))){
		stop("Count names should be in treatments and lib_sizes.")
	}
	fragment_size = 150  #not really used but for back of envelope calculations
	count_name = names(treatment)
	zero_prob = 0
	trans = extractTranscriptData(dt, min_transcripts, max_transcripts,
			group_transcripts_id, zero_prob)
	total_transcript_width = sum(trans$width)
	

	E_coverage = E_lib_size*fragment_size/total_transcript_width
	names(E_coverage) = names(E_lib_size)
	
	ntrans = length(trans$transcript_id)
	
	
	#######################################################################################
	# simulation
	#######################################################################################
	trans = simTranscriptParamsWithForceDE(trans,sdist, de_prob, force_de)
	trans = simTranscriptCounts(trans, treatment, E_lib_size)

	if(!only_sim_trans){
		trans_ex = simDisjointTranscriptCounts(dt,trans,count_name)
		disjoint = extractDisjointGenomeCounts(trans_ex, count_name)
	}
	
	

	
	e = environment()
	class(e) = c(class(e), "transcript_simulation")
	return(appendSimInfo(e))				
				
}

