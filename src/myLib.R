fuseProject <- function(meltCount, caseID = "case") {
    project <- rep(NA, nrow(meltCount))
    project[grep("AS.*T", as.character(meltCount[,caseID]))] <- "Uterine Adenosarcoma"
    project[grep("SCC.*T", as.character(meltCount[,caseID]))] <- "Uterine Adenosarcoma"    
    project[grep("SG.*", as.character(meltCount[,caseID]))] <- "Salivary Gland ACC"
    project[grep("R.*SG.*", as.character(meltCount[,caseID]))] <- "Salivary Gland ACC"
    project[grep("AdCC.*", as.character(meltCount[,caseID]))] <- "Adenoid Cystic Carcinoma"
    project[grep("GCT.*", as.character(meltCount[,caseID]))] <- "Granular Cell Tumor"
    project[grep("WOL.*", as.character(meltCount[,caseID]))] <- "Wolffian"
    project[grep("SST.*", as.character(meltCount[,caseID]))] <- "Wolffian"
    project[grep("MCM.*", as.character(meltCount[,caseID]))] <- "Mucinus"
    project[grep("MUC.*", as.character(meltCount[,caseID]))] <- "Mucinus"
    project[grep("UNCID.*", as.character(meltCount[,caseID]))] <- "TCGA-Normals"
    project[grep("THT.*", as.character(meltCount[,caseID]))] <- "Trabecular Adenoma Thyroid"
    project[grep("LPS.*", as.character(meltCount[,caseID]))] <- "Liposarcoma"
    project[grep("BPA.*", as.character(meltCount[,caseID]))] <- "Pleomorphic Adenoma"
    project[grep("BMEC2.*", as.character(meltCount[,caseID]))] <- "Mucoepidermoid Carcinoma"
    project[grep("META.*", as.character(meltCount[,caseID]))] <- "Metaplastic Carcinoma"
    project[grep("PAP.*", as.character(meltCount[,caseID]))] <- "Papillary"
    project[grep("AMB.*", as.character(meltCount[,caseID]))] <- "AMB"
    project[grep("ONS.*", as.character(meltCount[,caseID]))] <- "Olfactory Neuroblastoma"
    project[grep("ONL.*", as.character(meltCount[,caseID]))] <- "Olfactory Neuroblastoma"
    project[grep("ONR.*", as.character(meltCount[,caseID]))] <- "Olfactory Neuroblastoma"
    project[grep("OFNB.*", as.character(meltCount[,caseID]))] <- "Olfactory Neuroblastoma"
    project[grep("STP.*", as.character(meltCount[,caseID]))] <- "Pancreatic Sem."
    project[grep("^N[0-9]{1,2}", as.character(meltCount[,caseID]))] <- "FFPE Normals"
    project[grep("^TCGA*", as.character(meltCount[,caseID]))] <- "TEAD2fs"
    project[grep("BRCA-Normal*", as.character(meltCount[,caseID]))] <- "TCGA BRCA Normals"
    project[grep("nonBRCA*", as.character(meltCount[,caseID]))] <- "TCGA non-BRCA Normals"
    project[grep("TNBC*", as.character(meltCount[,caseID]))] <- "TCGA TNBC"
    return (project)
} ## end grep fuseProject


setKnown <- function(meltCount, caseID = "case", fusionName = "fusionName") {
    known <- meltCount[,"qual"]
    known[meltCount[,caseID] == "AS4T" & meltCount[,fusionName] == "ESR1--NCOA3"] <- "Known"
    known[meltCount[,caseID] == "AS4T" & meltCount[,fusionName] == "NCOA3--ESR1"] <- "Known"
    known[meltCount[,caseID] == "AS8T" & meltCount[,fusionName] == "ESR1--NCOA2"] <- "Known"
    known[meltCount[,caseID] == "AS8T" & meltCount[,fusionName] == "NCOA2--ESR1"] <- "Known"
    known[meltCount[,caseID] == "AdCC5M1T" & meltCount[,fusionName] == "MYBL1--NFIB"] <- "Known"
    known[meltCount[,caseID] == "AdCC5M1T" & meltCount[,fusionName] == "NFIB--MYBL1"] <- "Known"
    known[meltCount[,caseID] == "AdCC5" & meltCount[,fusionName] == "MYB--NFIB"] <- "Known"
    known[meltCount[,caseID] == "AdCC5" & meltCount[,fusionName] == "NFIB--MYB"] <- "Known"
    known[meltCount[,caseID] == "AdCC11" & meltCount[,fusionName] == "ACTN1--MYBL1"] <- "Known"
    known[meltCount[,caseID] == "AdCC11" & meltCount[,fusionName] == "MYBL1--ACTN1"] <- "Known"
    known[meltCount[,caseID] == "SG14T" & meltCount[,fusionName] == "NTRK3--ETV6"] <- "Known"
    known[meltCount[,caseID] == "SG14T" & meltCount[,fusionName] == "ETV6--NTRK3"] <- "Known"
    return(known)
} ## end set known


prepStar <- function(starFusion = file, tcga = tcgaFilter, ffpe = ffpeFilter, blast = blastFilter, knownFusions = chimerDB,  header = TRUE, row.names = NULL, fill = TRUE, sep = "\t", ...) {
    ## starFusions = file containing STAR-Fusion merged star-fusion *.fusion_candidates.txt
    ## filter = list of fusions to use as a filter
    ##    hh = c("program", "case", "outputFile", "fusion_name", "JunctionReads", "SpanningFrags", "Splice_type", "LeftGene", "LeftBreakpoint", "RightGene", "RightBreakpoint", "inTCGA.Normals")
    ## "FusionName", "JunctionReadCount", "SpanningFragCount", "SpliceType", "LeftGene", "LeftBreakpoint", "RightGene", "RightBreakpoint", "LargeAnchorSupport", "LeftBreakDinuc", "LeftBreakEntropy", "RightBreakDinuc", "RightBreakEntropy"
    hh = c("program", "case", "outputFile", "FusionName", "JunctionReadCount", "SpanningFragCount", "SpliceType", "LeftGene", "LeftBreakpoint", "RightGene", "RightBreakpoint", "LargeAnchorSupport", "LeftBreakDinuc", "LeftBreakEntropy", "RightBreakDinuc", "RightBreakEntropy", "J_FFPM", "S_FFPM") ## , "inTCGA.Normals")
    ss = read.table(starFusion, header = header, row.names = row.names, fill = fill, sep = sep, ...)
    colnames(ss) = hh
    ss$fusionName = as.character(ss$FusionName)
    ss$inTCGA.Normals = ss$fusionName %in% tcga
    ss$inFFPE.Normals = ss$fusionName %in% ffpe
    ss$inBLAST = ss$fusionName %in% blast
    ss$known.fusion = ss$fusionName %in% knownFusions
    return(ss)
} ## end prepStar

prepEricScript <- function(ericscript = file, tcga = tcgaFilter, ffpe = ffpeFilter, blast = blastFilter, knownFusions = chimerDB,  header = TRUE, row.names = NULL, fill = TRUE, sep = "\t", ...) {
    ## ericScript = file containing EricScript merged ercscript *.results.filtered.tsv
    ## filter = list of fusions to use as a filter
    hh = c("program", "case", "outputFile", "GeneName1", "GeneName2", "chr1", "Breakpoint1", "strand1", "chr2", "Breakpoint2", "strand2", "EnsemblGene1", "EnsemblGene2", "crossingreads", "spanningreads", "mean.insertsize", "homology", "fusiontype", "InfoGene1", "InfoGene2", "JunctionSequence", "GeneExpr1", "GeneExpr2", "GeneExpr_Fused", "ES", "GJS", "US", "EricScore")
    ee = read.table(ericscript, header = header, row.names = row.names, fill = fill, sep = sep, ...)
    colnames(ee) = hh
    ee = ee[ee$GeneName1 != ee$GeneName2, ]
    ee$fusionName <- as.character(paste(as.character(ee$GeneName1),  as.character(ee$GeneName2), sep = "--"))
    ee$inTCGA.Normals = ee$fusionName %in% tcga
    ee$inFFPE.Normals = ee$fusionName %in% ffpe    
    ee$inBLAST = ee$fusionName %in% blast
    ee$known.fusion = ee$fusionName %in% knownFusions
    return(ee)
} ## end ericScript


prepMapsplice <- function(mapsplice = file, tcga = tcgaFilter, ffpe = ffpeFilter, blast = blastFilter, knownFusions = chimerDB,  header = TRUE,  fill = TRUE, sep = "\t", ...) {
    ## mapsplice = file containing MAPSplice merged  *.fusions_candidates.txt
    ## filter = list of fusions to use as a filter
    hh = c("program", "case", "outputFile", "doner_chrom", "acceptor_chrom", "doner_end", "acceptor_start", "id", "coverage", "strand", "rgb", "block_count", "block_size", "block_distance", "entropy", "flank_case", "flank_string", "min_mismatch", "max_mismatch", "ave_mismatch", "max_min_suffix", "max_min_prefix", "min_anchor_difference", "unique_read_count", "multi_read_count", "paired_read_count", "left_paired_read_count", "right_paired_read_count", "multiple_paired_read_count", "unique_paired_read_count", "single_read_count", "encompassing_read_pair_count", "doner_start", "acceptor_end", "doner_iosforms", "acceptor_isoforms", "obsolete1", "obsolete2", "obsolete3", "obsolete4", "minimal_doner_isoform_length", "maximal_doner_isoform_length", "minimal_acceptor_isoform_length", "maximal_acceptor_isoform_length", "paired_reads_entropy", "mismatch_per_bp", "anchor_score", "max_doner_fragment", "max_acceptor_fragment", "max_cur_fragment", "min_cur_fragment", "ave_cur_fragment", "doner_encompass_unique", "doner_encompass_multiple", "acceptor_encompass_unique", "acceptor_encompass_multiple", "doner_match_to_normal", "acceptor_match_to_normal", "doner_seq", "acceptor_seq", "match_gene_strand", "annotated_type", "fusion_type", "gene_strand", "annotated_gene_donor", "annotated_gene_acceptor")
    ma = read.table(mapsplice, header = header, fill = fill, sep = sep, ...)
    colnames(ma) = hh
    ma = ma[ma$annotated_gene_donor != ma$annotated_gene_acceptor, ]
    ma$fusionName <- as.character(paste(ma$annotated_gene_donor,  ma$annotated_gene_acceptor, sep = "--"))
    ma$inTCGA.Normals = ma$fusionName %in% tcga
    ma$inFFPE.Normals = ma$fusionName %in% ffpe    
    ma$inBLAST = ma$fusionName %in% blast
    ma$known.fusion = ma$fusionName %in% knownFusions
    return(ma)
} ## end prepMapsplice

prepDefuse <- function(defuse = file, tcga = tcgaFilter, ffpe = ffpeFilter, blast = blastFilter, knownFusions = chimerDB,  header = TRUE,  fill = TRUE, sep = "\t", ...) {
    ## defuse = file containing defuse  merged  *.fusions_candidates.txt
    ## filter = list of fusions to use as a filter
    hh = c("program", "case", "outputFile", "cluster_id", "splitr_sequence", "splitr_count",
           "splitr_span_pvalue", "splitr_pos_pvalue", "splitr_min_pvalue", "adjacent", "altsplice",
           "break_adj_entropy1", "break_adj_entropy2", "break_adj_entropy_min", "breakpoint_homology",
           "breakseqs_estislands_percident", "cdna_breakseqs_percident", "deletion", "est_breakseqs_percident",
           "eversion", "exonboundaries", "expression1", "expression2", "gene1", "gene2", "gene_align_strand1",
           "gene_align_strand2", "gene_chromosome1", "gene_chromosome2", "gene_end1", "gene_end2",
           "gene_location1", "gene_location2", "gene_name1", "gene_name2", "gene_start1", "gene_start2",
           "gene_strand1", "gene_strand2", "genome_breakseqs_percident", "genomic_break_pos1",
           "genomic_break_pos2", "genomic_strand1", "genomic_strand2", "interchromosomal",
           "interrupted_index1", "interrupted_index2", "inversion", "library_name", "max_map_count",
           "max_repeat_proportion", "mean_map_count", "min_map_count", "num_multi_map", "num_splice_variants",
           "orf", "read_through", "repeat_proportion1", "repeat_proportion2", "span_count", "span_coverage1",
           "span_coverage2", "span_coverage_max", "span_coverage_min", "splice_score", "splicing_index1",
           "splicing_index2", "upstream_gene", "downstream_gene")
    de = read.table(defuse, header = header, fill = fill, sep = sep, ...)
    colnames(de) = hh
    de = de[de$gene_name1!= de$gene_name2, ]
    de$fusionName <- as.character(paste(de$gene_name1,  de$gene_name2, sep = "--"))
    de$inTCGA.Normals = de$fusionName %in% tcga
    de$inFFPE.Normals = de$fusionName %in% ffpe    
    de$inBLAST = de$fusionName %in% blast
    de$known.fusion = de$fusionName %in% knownFusions
    return(de)
} ## end prepDefuse

prepIntegrate <- function(integrate = file, tcga = tcgaFilter, ffpe = ffpeFilter, blast = blastFilter, knownFusions = chimerDB,  skip = 1, row.names = NULL, header = TRUE, fill = TRUE, sep = "\t", ...) {
    hh = c("program", "case", "outputFile", "x3p", "x5p", "chr1", "rna_bk1", "exon_bk1",
           "chr2", "rna_bk2", "exon_bk2", "wgs_bk1", "wgs_bk2", "fusion_candidate", "reciprocal",
           "tier", "type", "en_rna", "sp_rna", "splicings", "id", "x5p_transcript", "x5p_exon",
           "x5p_exon_strand", "x5p_exon_chr", "x5p_exon_start", "x5p_exon_end", "x5p_exon_seq",
           "x5p_exon_150", "x3p_transcript", "x3p_exon", "x3p_exon_strand", "x3p_exon_chr",
           "x3p_exon_start", "x3p_exon_end", "x3p_exon_seq", "x3p_exon_150", "x5p_strand",
           "x5p_start", "x5p_end", "x3p_strand", "x3p_start", "x3p_end",
           "FUSION_ID", "TISSUE", "SPANNING_READS", "ENCOMPASSING_READS",
           "5_FPG_GENE_NAME", "5_IN_CDS?", "5_SEGMENT_TYPE", "5_SEGMENT_ID", "5_COORD_IN_SEGMENT",
           "5_FULL_AA", "5_FRAME", "3_FPG_GENE_NAME", "3_IN_CDS?", "3_SEGMENT_TYPE", "3_SEGMENT_ID",
           "3_COORD_IN_SEGMENT", "3_FULL_AA", "3_FRAME", "FPG_FRAME_DIFFERENCE", "P_VAL_CORR",
           "DRIVER_PROB", "EXPRESSION_GAIN", "5_DOMAINS_RETAINED", "3_DOMAINS_RETAINED",
           "5_DOMAINS_BROKEN", "3_DOMAINS_BROKEN", "5_PII_RETAINED", "3_PII_RETAINED", "CTF", "G",
           "H", "K", "P", "TF")
    ii = read.table(integrate, header = header, fill = fill, sep = sep, skip = skip, ...)
    colnames(ii) = hh
    ii = ii[ii$x5p != ii$x3p, ]
    ii$fusionName <- as.character(paste(ii$x5p,  ii$x3p, sep = "--"))
    ii$inTCGA.Normals = ii$fusionName %in% tcga
    ii$inFFPE.Normals = ii$fusionName %in% ffpe    
    ii$inBLAST = ii$fusionName %in% blast
    ii$known.fusion = ii$fusionName %in% knownFusions
    return(ii)
}

prepOnc <- function(oncofuse = file, tcga = tcgaFilter, ffpe = ffpeFilter, blast = blastFilter, knownFusions = chimerDB,  row.names = NULL, header = TRUE, fill = TRUE, sep = "\t", ...) {
    hh = c("program", "case", "outputFile", "SAMPLE_ID", "FUSION_ID", "TISSUE", "SPANNING_READS",
           "ENCOMPASSING_READS", "GENOMIC", "5_FPG_GENE_NAME", "5_IN_CDS?", "5_SEGMENT_TYPE", "5_SEGMENT_ID",
           "5_COORD_IN_SEGMENT", "5_FULL_AA", "5_FRAME", "3_FPG_GENE_NAME", "3_IN_CDS?", "3_SEGMENT_TYPE",
           "3_SEGMENT_ID", "3_COORD_IN_SEGMENT", "3_FULL_AA", "3_FRAME", "FPG_FRAME_DIFFERENCE", "P_VAL_CORR",
           "DRIVER_PROB", "EXPRESSION_GAIN", "5_DOMAINS_RETAINED", "3_DOMAINS_RETAINED", "5_DOMAINS_BROKEN",
           "3_DOMAINS_BROKEN", "5_PII_RETAINED", "3_PII_RETAINED", "CTF", "G", "H", "K", "P", "TF")

    onc <- read.table(oncofuse, header = header, row.names = row.names, fill = fill, sep = sep, col.names = hh, ...)
    onc$fusionName <- as.character(paste(onc$X5_FPG_GENE_NAME, onc$X3_FPG_GENE_NAME, sep = "--"))
    onc$inTCGA.Normals <- onc$fusionName %in% tcga
    onc$inFFPE.Normals = onc$fusionName %in% ffpe
    onc$inBLAST = onc$fusionName %in% blast
    onc$known.fusion = onc$fusionName %in% knownFusions
    return(onc)
} ## end prepOnc

## mapsplice header
maphh <- c("program", "case", "outputFile", "doner_chrom", "acceptor_chrom", "doner_end", "acceptor_start", "id", "coverage", "strand", "rgb", "block_count", "block_size", "block_distance", "entropy", "flank_case", "flank_string", "min_mismatch", "max_mismatch", "ave_mismatch", "max_min_suffix", "max_min_prefix", "min_anchor_difference", "unique_read_count", "multi_read_count", "paired_read_count", "left_paired_read_count", "right_paired_read_count", "multiple_paired_read_count", "unique_paired_read_count", "single_read_count", "encompassing_read_pair_count", "doner_start", "acceptor_end", "doner_iosforms", "acceptor_isoforms", "obsolete", "obsolete", "obsolete", "obsolete", "minimal_doner_isoform_length", "maximal_doner_isoform_length", "minimal_acceptor_isoform_length", "maximal_acceptor_isoform_length", "paired_reads_entropy", "mismatch_per_bp", "anchor_score", "max_doner_fragment", "max_acceptor_fragment", "max_cur_fragment", "min_cur_fragment", "ave_cur_fragment", "doner_encompass_unique", "doner_encompass_multiple", "acceptor_encompass_unique", "acceptor_encompass_multiple", "doner_match_to_normal", "acceptor_match_to_normal", "doner_seq", "acceptor_seq", "match_gene_strand", "annotated_type", "fusion_type", "gene_strand", "annotated_gene_donor", "annotated_gene_acceptor")


getNormals <- function() {
    star.normals <<- star[star$case %in% grep("Normal", star$project, value = TRUE), ]
    map.normals <<- map[map$case %in% grep("Normal", map$project, value = TRUE), ]
    ee.normals <<- ee[ee$case %in% grep("Normal", ee$project, value = TRUE), ]
    intetrate.normals <<- integrate[integrate$case %in% grep("Normal", integrate$project, value = TRUE), ]
    defuse.normals <<- defuse[defuse$case %in% grep("Normal", defuse$project, value = TRUE), ]
} ## end getNormals



getRecurrent <- function(fusion.df, caseID = "case" , fusionName = "fusionName") {
    recurrent = table(unique(fusion.df[, c(fusionName, caseID)]))
    recurrent = data.frame(recurrent > 0)    
    recurrent$n.cases = rowSums(recurrent)
    recurrent = recurrent[recurrent$n.cases > 1, ]
    return(recurrent)
} ## end getRecurrent


getCrossProject <- function(fusion.df, fusionName = "fusionName", project = "project") {
    cp = table(fusion.df[, c(fusionName, project)])
    cp = data.frame(cp > 0)
    cp$n.proj = rowSums(cp)
    cp = cp[cp$n.proj > 1, ]
#    return(cp)
}

countFusions <- function(fusion.df, fusionName = "fusionName", caseID = "case", merge.type = FALSE, type.df) {
    cn <- melt(table(fusion.df[, c(fusionName, caseID)]))
    cn$value[cn$value == 0] <- NA
    cn$qual = cn$value
    cn$qual[!is.na(cn$qual)] = "Predicted"
    cn$qual[is.na(cn$qual)] = "Not Present"
    if(merge.type) { cn = merge(cn, type.df, by = "case") }
    return(cn)
}


p <- ggplot(mpg,aes(x=class,fill=class)) + geom_bar()
ggplot_build(p)$data

## getRecurrent <- function(fusion.df, caseID = "case" , fusionName = "fusionName", project = "project") {
##    recurrent <- unique(fusion.df[ , c(caseID, fusionName, project)])
##    recurrent$mmfc <- paste(recurrent$fusionName, recurrent$case, sep = "_")
##    recurrent <- data.frame(tapply(recurrent$mmfc, recurrent$fusionName, length))
##    names(recurrent) <- "n.cases"
##    recurrent <- data.frame(recurrent[recurrent$n.cases > 1, ])
##    names(recurrent) <- "n.cases"
##    return(recurrent)
## }

string.counter<-function(strings, pattern){  
  counts<-NULL
  for(i in 1:length(strings)){
    counts[i]<-length(attr(gregexpr(pattern,strings[i])[[1]], "match.length")[attr(gregexpr(pattern,strings[i])[[1]], "match.length")>0])
  }
return(counts)
}
