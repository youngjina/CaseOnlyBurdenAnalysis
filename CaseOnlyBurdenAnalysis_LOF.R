# Usage: Rscript CaseOnlyBurdenAnalysis_LOF.R --help

# Load the packages
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(mgcv))

option_list = list(
	make_option(c("-r", "--mutation.rate.file"), type = "character", default = NULL, 
			help = "mutation rate file", metavar = "character"),
	make_option(c("-g", "--genotype.dir"), type = "character", default = NULL, 
			help = "genotype directory", metavar = "character"),
	make_option(c("-a", "--artifact.file"), type = "character", default = NULL, 
			help = "artifact file (optional)", metavar = "character"),
	make_option(c("-p", "--pli.source"), type = "character", default = "exac", 
			help = "pLI score: exac (default) or gnomad", metavar = "character"),
	make_option(c("-m", "--n.male"), type = "integer", default = NULL, 
			help = "Number of male", metavar = "number"),
	make_option(c("-f", "--n.female"), type = "integer", default = NULL, 
			help = "Number of female", metavar = "number"),
	make_option(c("-v", "--variant.type"), type = "character", default = "singleton", 
			help = "variants used as input: singleton (default) or all", metavar = "character"),
	make_option(c("-o", "--output.dir"), type = "character", default = "output", 
			help = "output directory", metavar = "character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

genotypeFile <- Sys.glob( paste( opt$genotype.dir , "*_genotypes.csv" , sep="" ) )

# X_linked.adjustment_factor
X.chr <- opt$n.male + ( opt$n.female * 2 )
auto.chr <- ( opt$n.male + opt$n.female ) * 2
X_linked.adjustment_factor <- X.chr / auto.chr

# mutation_ratio.adjustment_factor
mutation_ratio.adjustment_factor <- 0.013 / ( 1.2 * 10^(-8) ) 


##################################################################
# Observed variants
##################################################################
print( "Step 1: Observed variants" )

input.data <- read.csv( genotypeFile , header=T )

if( opt$variant.type == "singleton") {

	# Extract singletons
	singletons <- input.data %>% group_by(Variant.ID) %>% tally() %>% filter(n == 1) %>% select(Variant.ID)
	temp.observed.genotypes <- input.data %>% filter( Variant.ID %in% singletons$Variant.ID ) 
	
} else if( opt$variant.type == "all") {
	temp.observed.genotypes <- input.data
}

if( is.null(opt$artifact.file) ) {
	observed.genotypes <- temp.observed.genotypes
} else {

	# Remove artifacts
	artifacts <- data.frame( read.table( opt$artifact.file , header=F ) )
	observed.genotypes <- temp.observed.genotypes %>% filter( !( Variant.ID %in% artifacts[,1] ) )
	
}

# Calculate ratio of frameshift to nonsense
cal.ratio.mat <- observed.genotypes[ , c("Variant.ID", "Effect") ] %>% distinct()
frameshift_count <- nrow( subset( cal.ratio.mat, Effect == "frameshift_variant" ) )
nonsense_count <- nrow( subset( cal.ratio.mat, Effect == "stop_gained" | Effect == "start_lost" ) )
adj.factor <- frameshift_count / nonsense_count

# Extract X linked genes
observed.genotypes.varID <- strsplit( as.character( observed.genotypes$Variant.ID ) , "-" )
chrom <- data.frame(t(sapply(observed.genotypes.varID, `[`)))[,1]
Xlinked.genes <- as.character( gsub( "'", "", unique( observed.genotypes[ which( chrom == "X" ) , "Gene.Name" ] ) ))

# Split LOF genotypes and FrameShift genotypes
LOF.observed.genotypes <- observed.genotypes %>% filter( Effect != "frameshift_variant" )
FrameShift.observed.genotypes <- observed.genotypes %>% filter( Effect == "frameshift_variant" )

# Count the number of observed variants per gene
## LOF
LOF.observed.variant.gene.pair <- LOF.observed.genotypes[ , c("Variant.ID", "Gene.Name") ]
uniq.LOF.observed.variant.gene.pair <- LOF.observed.variant.gene.pair[ !duplicated( LOF.observed.variant.gene.pair$Variant.ID ) , ]
LOF.observed.variant.per.gene <- plyr::count( uniq.LOF.observed.variant.gene.pair, "Gene.Name" )
LOF.observed.variant.per.gene$Gene.Name <- gsub( "'", "", LOF.observed.variant.per.gene$Gene.Name )

## FrameShift
FrameShift.observed.variant.gene.pair <- FrameShift.observed.genotypes[ , c("Variant.ID", "Gene.Name") ]
uniq.FrameShift.observed.variant.gene.pair <- FrameShift.observed.variant.gene.pair[ !duplicated( FrameShift.observed.variant.gene.pair$Variant.ID ) , ]
FrameShift.observed.variant.per.gene <- plyr::count( uniq.FrameShift.observed.variant.gene.pair, "Gene.Name" )
FrameShift.observed.variant.per.gene$Gene.Name <- gsub( "'", "", FrameShift.observed.variant.per.gene$Gene.Name )


##################################################################
# [LOF] Mutation rate
##################################################################
print( "Step 2: LOF Mutation rate" )

LOF.mut.data <- fread( opt$mutation.rate.file , header=T , sep="\t" )

uniq.variant.LOF.mut.data <- LOF.mut.data[ !duplicated( LOF.mut.data$Variant.ID ) , ]
uniq.variant.LOF.mut.gene <- uniq.variant.LOF.mut.data[ , c( "Gene", "Variant.ID", "Mut.Rate" )  ]
uniq.variant.LOF.mut.gene$Gene <- gsub( "'", "", uniq.variant.LOF.mut.gene$Gene )
sum.mut.gene.Count <- uniq.variant.LOF.mut.gene %>% group_by(Gene) %>% tally()
colnames(sum.mut.gene.Count) <- c( "Gene", "Predicted.Variants" )

# Aggregate mutation rate per position into mutation rate per gene by sum
LOF.sum.mut.gene <- uniq.variant.LOF.mut.gene %>% group_by(Gene) %>% tally(Mut.Rate)
LOF.mutation_rate <- LOF.sum.mut.gene$n / mutation_ratio.adjustment_factor
names(LOF.mutation_rate) <- LOF.sum.mut.gene$Gene

expected.mutation_rate <- as.data.frame(LOF.mutation_rate)
expected.mutation_rate <- rownames_to_column(expected.mutation_rate, var = "Gene.Name")

autosome.expected.mutation_rate <- expected.mutation_rate %>% filter( !(as.character(Gene.Name) %in% Xlinked.genes) )

# Adjustment to X linked genes mutation rate
Xlinked.expected.mutation_rate <- expected.mutation_rate %>% filter( as.character(Gene.Name) %in% Xlinked.genes ) 
Xlinked.expected.mutation_rate$LOF.mutation_rate <- Xlinked.expected.mutation_rate$LOF.mutation_rate * X_linked.adjustment_factor

expected.mutation_rate.adj <- rbind( autosome.expected.mutation_rate , Xlinked.expected.mutation_rate )


##################################################################
# [FrameShift] Mutation rate
##################################################################
print( "Step 3: FrameShift Mutation rate" )

expected.mutation_rate.adj$FrameShift.mutation_rate <- expected.mutation_rate.adj$LOF.mutation_rate * adj.factor 

LOF.input.model <- merge( x = LOF.observed.variant.per.gene, y = expected.mutation_rate.adj , by = "Gene.Name" , all.x = T )

FrameShift.input.model <- merge( x = FrameShift.observed.variant.per.gene, y = expected.mutation_rate.adj , by = "Gene.Name" , all.x = T )


##################################################################
# Expected number of variants
##################################################################
print( "Step 4: Expected variants" )

LOF.Obs <- LOF.input.model[ , c( "Gene.Name", "freq" ) ]
colnames(LOF.Obs) <- c("Gene.Name", "LOF.Observed.Variant.Number")

FrameShift.Obs <- FrameShift.input.model[ , c( "Gene.Name", "freq" ) ]
colnames(FrameShift.Obs) <- c("Gene.Name", "FrameShift.Observed.Variant.Number")

# Consider all genes to calculate mutation rate
LOF.FrameShift.Obs <- merge(x=LOF.Obs, y=FrameShift.Obs, by="Gene.Name", all=T)
LOF.FrameShift.Obs[is.na(LOF.FrameShift.Obs)] <- 0

LOF.FrameShift.Obs.Mut <- merge(x=LOF.FrameShift.Obs, y=expected.mutation_rate.adj, by="Gene.Name", all.y=T)

LOF.FrameShift.Obs.Mut[ is.na( LOF.FrameShift.Obs.Mut ) ] <- 0 


##################################################################
# pLI-stratification
##################################################################
print( "Step 5: pLI-stratification" )

if( opt$pli.source == "exac" ) {

	pli.score <- read.table("/media/Dazs/youngji/Resource/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt", header=T, sep="\t")
	constraint.matrix <- pli.score[ , c("Gene.Name", "pLI", "mis_z")]
	merging.mut.exac <- merge( x=LOF.FrameShift.Obs.Mut , y=constraint.matrix, by="Gene.Name", all.x=T)

} else if( opt$pli.source == "gnomad" ) {

	pli.score <- read.table("/media/Dazs/youngji/Resource/gnomad.v2.1.1.lof_metrics.by_gene.txt", header=T, sep="\t")
	temp.constraint.matrix <- pli.score[ , c("gene", "pLI", "mis_z")]
	constraint.matrix <- temp.constraint.matrix %>% group_by(gene) %>% filter(mis_z == max(mis_z)) %>% filter(pLI == max(pLI))
	merging.mut.exac <- merge( x=LOF.FrameShift.Obs.Mut , y=constraint.matrix, by.x="Gene.Name", by.y="gene", all.x=T)

}

ssc.dnm.ihv <- read.table("/media/Dazs/youngji/Resource/SSC2_Singleton_DeNovo/Ratio_Singleton_DeNovo.txt", header=T, sep="\t")

merging.mut.pli.ssc <- merge( x=merging.mut.exac, y=ssc.dnm.ihv[ , c("gene", "n_singleton", "n.denovo", "ratio")], by.x="Gene.Name", by.y="gene", all.x=T)

merging.mut.pli.ssc$ratio[is.na(merging.mut.pli.ssc$ratio)] <- 0

pLI.decile <- mutate(merging.mut.pli.ssc, decile_rank = ntile(merging.mut.exac$pLI,10))

pLI.group_1.genotypes <- pLI.decile %>% filter( decile_rank == 1 )
pLI.group_2.genotypes <- pLI.decile %>% filter( decile_rank == 2 )
pLI.group_3.genotypes <- pLI.decile %>% filter( decile_rank == 3 )
pLI.group_4.genotypes <- pLI.decile %>% filter( decile_rank == 4 )
pLI.group_5.genotypes <- pLI.decile %>% filter( decile_rank == 5 )
pLI.group_6.genotypes <- pLI.decile %>% filter( decile_rank == 6 )
pLI.group_7.genotypes <- pLI.decile %>% filter( decile_rank == 7 )
pLI.group_8.genotypes <- pLI.decile %>% filter( decile_rank == 8 )
pLI.group_9.genotypes <- pLI.decile %>% filter( decile_rank == 9 )
pLI.group_10.genotypes <- pLI.decile %>% filter( decile_rank == 10 )
pLI.group_NA.genotypes <- pLI.decile %>% filter( is.na(decile_rank) )

adjust_group <- function( group_interval, group_prefix ) {

	pLI.stratified.Observed.Variant.Number <- base::rowSums( group_interval[ , c("LOF.Observed.Variant.Number", "FrameShift.Observed.Variant.Number")] )

	pLI.stratified.Expected.Gene.Mutation <- base::rowSums( group_interval[ , c("LOF.mutation_rate", "FrameShift.mutation_rate")] )

	pLI.stratified.Total.Variant.Number <- sum(pLI.stratified.Observed.Variant.Number)

	pLI.stratified.model <- lm( pLI.stratified.Observed.Variant.Number ~ polym( pLI.stratified.Expected.Gene.Mutation , group_interval$ratio, degree=2 ) )
	pLI.stratified.fitted.value <- fitted( pLI.stratified.model )
	
	pLI.stratified.Expected.Variant.Number <- c()
	for( i in 1 : nrow( group_interval ) ) {
		pLI.stratified.Expected.Variant.Number <- append( pLI.stratified.Expected.Variant.Number , pLI.stratified.Total.Variant.Number * ( pLI.stratified.fitted.value[ i ] / sum( pLI.stratified.fitted.value ) ) )
	}

	pLI.stratified.Enrichment <- pLI.stratified.Observed.Variant.Number / pLI.stratified.Expected.Variant.Number

	pLI.stratified.P.value <- c()
	for( i in 1 : nrow( group_interval ) ) { 
	
		suppressWarnings(
	
		pLI.stratified.P.value <- append( pLI.stratified.P.value , 
		( ppois( q = pLI.stratified.Observed.Variant.Number[i] - 1 , 
		lambda = pLI.stratified.Expected.Variant.Number[i] ,
		lower.tail = FALSE ) ) )
		
		)
	
	}
	names( pLI.stratified.P.value ) <- group_interval$Gene.Name

	pLI.stratified.Obs.Exp.matrix <- cbind( group_interval, pLI.stratified.Observed.Variant.Number, pLI.stratified.Expected.Gene.Mutation, pLI.stratified.fitted.value, pLI.stratified.Expected.Variant.Number, pLI.stratified.Enrichment, pLI.stratified.P.value )

	pLI.stratified.Obs.Exp.Pred <- merge(x=pLI.stratified.Obs.Exp.matrix, y=sum.mut.gene.Count, by.x="Gene.Name", by.y="Gene", all.x=T)

	pLI.stratified.Obs.Exp.Pred$Gene.Name <- paste( "'" , pLI.stratified.Obs.Exp.Pred$Gene , "'" , sep="" )

	pLI.stratified.sorted.anno.final.output <- pLI.stratified.Obs.Exp.Pred[ order( pLI.stratified.Obs.Exp.Pred$pLI.stratified.P.value ) , ]
	
	pLI.stratified.FINAL.output <- pLI.stratified.sorted.anno.final.output[ , c("Gene.Name", "LOF.Observed.Variant.Number", "FrameShift.Observed.Variant.Number", "pLI.stratified.Observed.Variant.Number", "LOF.mutation_rate", "FrameShift.mutation_rate", "pLI.stratified.Expected.Gene.Mutation", "pLI.stratified.Expected.Variant.Number", "pLI.stratified.fitted.value", "pLI.stratified.Enrichment", "pLI.stratified.P.value", "Predicted.Variants", "pLI", "mis_z", "n_singleton", "n.denovo", "ratio") ]
	colnames(pLI.stratified.FINAL.output) <- c("Gene.Name", "LOF.Observed.Variant.Number", "FrameShift.Observed.Variant.Number", "Observed.Variant.Number", "LOF.mutation_rate", "FrameShift.mutation_rate", "Expected.Gene.Mutation", "Expected.Variant.Number", "fitted.value", "Enrichment", "P.value", "Predicted.Variants", "pLI", "mis_z", "Singleton.Number", "DNM.Number", "Ratio.DNM.Singleton")

	write.csv( pLI.stratified.FINAL.output , paste( group_prefix, "Carlson_CaseOnly_Pvalue.csv", sep="_" ) , row.names=F , quote=F )

	pLI.stratified.burden.test.plot.input <- pLI.stratified.FINAL.output[ , c("Gene.Name", "P.value") ]
	colnames(pLI.stratified.burden.test.plot.input) <- c("Gene", "P_DOM")
	pLI.stratified.tempFile <- paste( group_prefix, ".burden.txt", sep="" )
	write.table( pLI.stratified.burden.test.plot.input, pLI.stratified.tempFile, sep="\t", quote=F, row.names=F )

	system(paste("Rscript /media/Dazs/youngji/Scripts/trapd_base_QQ.R", "--prefix", paste( group_prefix, sep="") ) )
	invisible(file.remove( pLI.stratified.tempFile ))
	
	return(pLI.stratified.burden.test.plot.input)

}

dir.create(opt$output.dir, showWarnings=F)
setwd(opt$output.dir)
sink('analysis-output.log', append=TRUE)
print(opt)
print( paste( "X_linked.adjustment_factor:", X_linked.adjustment_factor ) )
print( paste( "ratio of frameshift to nonsense:", adj.factor ) )
sink()

pLI.group_1.output <- adjust_group( pLI.group_1.genotypes , "pLI.group_1" )
pLI.group_2.output <- adjust_group( pLI.group_2.genotypes , "pLI.group_2" )
pLI.group_3.output <- adjust_group( pLI.group_3.genotypes , "pLI.group_3" )
pLI.group_4.output <- adjust_group( pLI.group_4.genotypes , "pLI.group_4" )
pLI.group_5.output <- adjust_group( pLI.group_5.genotypes , "pLI.group_5" )
pLI.group_6.output <- adjust_group( pLI.group_6.genotypes , "pLI.group_6" )
pLI.group_7.output <- adjust_group( pLI.group_7.genotypes , "pLI.group_7" )
pLI.group_8.output <- adjust_group( pLI.group_8.genotypes , "pLI.group_8" )
pLI.group_9.output <- adjust_group( pLI.group_9.genotypes , "pLI.group_9" )
pLI.group_10.output <- adjust_group( pLI.group_10.genotypes , "pLI.group_10" )
pLI.group_NA.output <- adjust_group( pLI.group_NA.genotypes , "pLI.group_NA" )

pooling.output <- rbind(
pLI.group_1.output,
pLI.group_2.output,
pLI.group_3.output,
pLI.group_4.output,
pLI.group_5.output,
pLI.group_6.output,
pLI.group_7.output,
pLI.group_8.output,
pLI.group_9.output,
pLI.group_10.output,
pLI.group_NA.output
)  


##################################################################
# Generate Final Results
##################################################################
print( "Step 6: Generate Final Results" )

pLI.stratified.tempFile <- "Final.Result.burden.txt"
write.table( pooling.output, pLI.stratified.tempFile, sep="\t", quote=F, row.names=F )

system(paste("Rscript /media/Dazs/youngji/Scripts/trapd_base_QQ.R", "--prefix", "Final.Result"))
invisible(file.remove( pLI.stratified.tempFile ))

output <- c()
file.list <- Sys.glob("pLI.group_*_Carlson_CaseOnly_Pvalue.csv")
for( inFile in file.list ) {
	input.data <- read.csv( inFile , header=T )
	output <- rbind( output, input.data )
}
pooling_stratification <- output[order(output$P.value),]

print("Top genes")
print(head(pooling_stratification))
write.csv( pooling_stratification, "Final.Result.Pvalue.csv", quote=F, row.names=F )
