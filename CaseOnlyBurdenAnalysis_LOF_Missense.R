# Usage: Rscript CaseOnlyBurdenAnalysis_LOF_Missense.R --help

# Load the packages
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(mgcv))

option_list = list(
	make_option(c("-l", "--LOF.mutation.rate.file"), type = "character", default = NULL, 
			help = "LOF mutation rate file", metavar = "character"),
	make_option(c("-r", "--Missense.mutation.rate.file"), type = "character", default = NULL, 
			help = "Missense mutation rate file", metavar = "character"),
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

LOF.observed.genotypes <- observed.genotypes %>% filter( 
Effect == "exon_loss_variant" | 
Effect == "rare_amino_acid_variant" | 
Effect == "stop_gained" | 
Effect == "stop_lost" | 
Effect == "start_lost" | 
Effect == "gene_fusion" | 
Effect == "bidirectional_gene_fusion" | 
Effect == "splice_acceptor_variant" | 
Effect == "splice_donor_variant" 
)

FrameShift.observed.genotypes <- observed.genotypes %>% filter( Effect == "frameshift_variant" )

Missense.observed.genotypes <- observed.genotypes %>% filter( 
Effect == "missense_variant" | 
Effect == "missense_variant+splice_region_variant" | 
Effect == "conservative_inframe_insertion" | 
Effect == "conservative_inframe_deletion" | 
Effect == "disruptive_inframe_insertion" | 
Effect == "disruptive_inframe_deletion" | 
Effect == "coding_sequence_variant" | 
Effect == "5_prime_UTR_truncation+exon_loss_variant" | 
Effect == "3_prime_UTR_truncation+exon_loss_variant" | 
Effect == "5_prime_UTR_premature_start_codon_gain_variant" | 
Effect == "initiator_codon_variant" | 
Effect == "initiator_codon_variant+non_canonical_start_codon"
)

# Count the number of observed variants per gene
## LOF
LOF.observed.variant.gene.pair <- LOF.observed.genotypes[ , c("Variant.ID", "Gene.Name") ]
uniq.LOF.observed.variant.gene.pair <- LOF.observed.variant.gene.pair[ !duplicated( LOF.observed.variant.gene.pair$Variant.ID ) , ]
LOF.observed.variant.per.gene <- plyr::count( uniq.LOF.observed.variant.gene.pair, "Gene.Name" )
LOF.observed.variant.per.gene$Gene.Name <- gsub( "'", "", LOF.observed.variant.per.gene$Gene.Name )

LOF.out <- strsplit( as.character( LOF.observed.genotypes$Variant.ID ) , "-" )
LOF.chrom <- data.frame(t(sapply(LOF.out, `[`)))[,1]
LOF.Xlinked.genes <- as.character( gsub( "'", "", unique( LOF.observed.genotypes[ which( LOF.chrom == "X" ) , "Gene.Name" ] ) ))

# X-linked genes
Xlinked.LOF.observed.variant.per.gene <- LOF.observed.variant.per.gene %>% filter( Gene.Name %in% LOF.Xlinked.genes )

# Autosome genes
autosome.LOF.observed.variant.per.gene <- LOF.observed.variant.per.gene %>% filter( !( Gene.Name %in% LOF.Xlinked.genes ) )


## FrameShift
FrameShift.observed.variant.gene.pair <- FrameShift.observed.genotypes[ , c("Variant.ID", "Gene.Name") ]
uniq.FrameShift.observed.variant.gene.pair <- FrameShift.observed.variant.gene.pair[ !duplicated( FrameShift.observed.variant.gene.pair$Variant.ID ) , ]
FrameShift.observed.variant.per.gene <- plyr::count( uniq.FrameShift.observed.variant.gene.pair, "Gene.Name" )
FrameShift.observed.variant.per.gene$Gene.Name <- gsub( "'", "", FrameShift.observed.variant.per.gene$Gene.Name )

FrameShift.out <- strsplit( as.character( FrameShift.observed.genotypes$Variant.ID ) , "-" )
FrameShift.chrom <- data.frame(t(sapply(FrameShift.out, `[`)))[,1]
FrameShift.Xlinked.genes <- as.character( gsub( "'", "", unique( FrameShift.observed.genotypes[ which( FrameShift.chrom == "X" ) , "Gene.Name" ] ) ))

# X-linked genes
Xlinked.FrameShift.observed.variant.per.gene <- FrameShift.observed.variant.per.gene %>% filter( Gene.Name %in% FrameShift.Xlinked.genes )

# Autosome genes
autosome.FrameShift.observed.variant.per.gene <- FrameShift.observed.variant.per.gene %>% filter( !( Gene.Name %in% FrameShift.Xlinked.genes ) )


## Missense
Missense.observed.variant.gene.pair <- Missense.observed.genotypes[ , c("Variant.ID", "Gene.Name") ]
uniq.Missense.observed.variant.gene.pair <- Missense.observed.variant.gene.pair[ !duplicated( Missense.observed.variant.gene.pair$Variant.ID ) , ]
Missense.observed.variant.per.gene <- plyr::count( uniq.Missense.observed.variant.gene.pair, "Gene.Name" )
Missense.observed.variant.per.gene$Gene.Name <- gsub( "'", "", Missense.observed.variant.per.gene$Gene.Name )

Missense.out <- strsplit( as.character( Missense.observed.genotypes$Variant.ID ) , "-" )
Missense.chrom <- data.frame(t(sapply(Missense.out, `[`)))[,1]
Missense.Xlinked.genes <- as.character( gsub( "'", "", unique( Missense.observed.genotypes[ which( Missense.chrom == "X" ) , "Gene.Name" ] ) ))

# X-linked genes
Xlinked.Missense.observed.variant.per.gene <- Missense.observed.variant.per.gene %>% filter( Gene.Name %in% Missense.Xlinked.genes )

# Autosome genes
autosome.Missense.observed.variant.per.gene <- Missense.observed.variant.per.gene %>% filter( !( Gene.Name %in% Missense.Xlinked.genes ) )


##################################################################
# [LOF] Mutation rate
##################################################################
print( "Step 2: LOF Mutation rate" )

LOF.mut.data <- fread( opt$LOF.mutation.rate.file , header=T , sep="\t" )

uniq.variant.LOF.mut.data <- LOF.mut.data[ !duplicated( LOF.mut.data$Variant.ID ) , ]
uniq.variant.LOF.mut.gene <- uniq.variant.LOF.mut.data[ , c( "Gene", "Variant.ID", "Mut.Rate" )  ]
uniq.variant.LOF.mut.gene$Gene <- gsub( "'", "", uniq.variant.LOF.mut.gene$Gene )
LOF.sum.mut.gene.Count <- uniq.variant.LOF.mut.gene %>% group_by(Gene) %>% tally()
colnames(LOF.sum.mut.gene.Count) <- c( "Gene", "LOF.Predicted.Variants" )

# Aggregate mutation rate per position into mutation rate per gene by sum
LOF.sum.mut.gene <- uniq.variant.LOF.mut.gene %>% group_by(Gene) %>% tally(Mut.Rate)
LOF.mutation_rate <- LOF.sum.mut.gene$n / mutation_ratio.adjustment_factor
names(LOF.mutation_rate) <- LOF.sum.mut.gene$Gene


##################################################################
# [FrameShift] Mutation rate
##################################################################
print( "Step 3: FrameShift Mutation rate" )

FrameShift.mutation_rate <- LOF.mutation_rate * adj.factor 


##################################################################
# [Missense] Mutation rate
##################################################################
print( "Step 4: Missense Mutation rate" )

Missense.mut.data <- fread( opt$Missense.mutation.rate.file , header=T , sep="\t" )

uniq.variant.Missense.mut.data <- Missense.mut.data[ !duplicated( Missense.mut.data$Variant.ID ) , ]
uniq.variant.Missense.mut.gene <- uniq.variant.Missense.mut.data[ , c( "Gene", "Variant.ID", "Mut.Rate" )  ]
uniq.variant.Missense.mut.gene$Gene <- gsub( "'", "", uniq.variant.Missense.mut.gene$Gene )
Missense.sum.mut.gene.Count <- uniq.variant.Missense.mut.gene %>% group_by(Gene) %>% tally()
colnames(Missense.sum.mut.gene.Count) <- c( "Gene", "Missense.Predicted.Variants" )

# Aggregate mutation rate per position into mutation rate per gene by sum
Missense.sum.mut.gene <- uniq.variant.Missense.mut.gene %>% group_by(Gene) %>% tally(Mut.Rate)
Missense.mutation_rate <- Missense.sum.mut.gene$n / mutation_ratio.adjustment_factor
names(Missense.mutation_rate) <- Missense.sum.mut.gene$Gene


##################################################################
# Merged Mutation rate
##################################################################
print( "Step 5: Merged Mutation rate" )

LOF.FrameShift.expected.mutation_rate <- cbind(LOF.mutation_rate, FrameShift.mutation_rate)
LOF.FrameShift.expected.mutation_rate <- rownames_to_column(as.data.frame(LOF.FrameShift.expected.mutation_rate))
colnames(LOF.FrameShift.expected.mutation_rate) <- c("Gene", "LOF.mutation_rate", "FrameShift.mutation_rate")

Missense.expected.mutation_rate <- rownames_to_column(as.data.frame(Missense.mutation_rate))
colnames(Missense.expected.mutation_rate) <- c("Gene", "Missense.mutation_rate")

expected.mutation_rate <- merge(x=LOF.FrameShift.expected.mutation_rate, y=Missense.expected.mutation_rate, by="Gene", all=T)
expected.mutation_rate[ is.na( expected.mutation_rate ) ] <- 0


##################################################################
# Expected number of variants
##################################################################
print( "Step 6: Expected variants" )

LOF.observed.autosome.Xlinked <- rbind( autosome.LOF.observed.variant.per.gene , Xlinked.LOF.observed.variant.per.gene )
FrameShift.observed.autosome.Xlinked <- rbind( autosome.FrameShift.observed.variant.per.gene , Xlinked.FrameShift.observed.variant.per.gene )
Missense.observed.autosome.Xlinked <- rbind( autosome.Missense.observed.variant.per.gene , Xlinked.Missense.observed.variant.per.gene )

# Merged data
LOF.input.model <- merge( x = LOF.observed.autosome.Xlinked , y = expected.mutation_rate , by.x = "Gene.Name" , by.y = "Gene" , all.y = T )
LOF.input.model[ is.na( LOF.input.model ) ] <- 0

FrameShift.input.model <- merge( x = FrameShift.observed.autosome.Xlinked , y = expected.mutation_rate , by.x = "Gene.Name" , by.y = "Gene" , all.y = T )
FrameShift.input.model[ is.na( FrameShift.input.model ) ] <- 0

Missense.input.model <- merge( x = Missense.observed.autosome.Xlinked , y = expected.mutation_rate , by.x = "Gene.Name" , by.y = "Gene" , all.y = T )
Missense.input.model[ is.na( Missense.input.model ) ] <- 0

LOF.Observed.Variant.Number <- LOF.input.model$freq
LOF.Expected.Gene.Mutation <- LOF.input.model$LOF.mutation_rate

FrameShift.Observed.Variant.Number <- FrameShift.input.model$freq
FrameShift.Expected.Gene.Mutation <- FrameShift.input.model$FrameShift.mutation_rate

Missense.Observed.Variant.Number <- Missense.input.model$freq
Missense.Expected.Gene.Mutation <- Missense.input.model$Missense.mutation_rate

LOF.Obs <- LOF.input.model[ , c( "Gene.Name", "freq" ) ]
colnames(LOF.Obs) <- c("Gene.Name", "LOF.Observed.Variant.Number")

FrameShift.Obs <- FrameShift.input.model[ , c( "Gene.Name", "freq" ) ]
colnames(FrameShift.Obs) <- c("Gene.Name", "FrameShift.Observed.Variant.Number")

Missense.Obs <- Missense.input.model[ , c( "Gene.Name", "freq" ) ]
colnames(Missense.Obs) <- c("Gene.Name", "Missense.Observed.Variant.Number")

LOF.FrameShift.Obs <- merge(x=LOF.Obs, y=FrameShift.Obs, by="Gene.Name", all=T)
LOF.FrameShift.Obs[is.na(LOF.FrameShift.Obs)] <- 0

LOF.FrameShift.Missense.Obs <- merge(x=LOF.FrameShift.Obs, y=Missense.Obs, by="Gene.Name", all=T)

LOF.FrameShift.Missense.Obs[is.na(LOF.FrameShift.Missense.Obs)] <- 0

LOF.FrameShift.Missense.Obs.Mut <- merge(x=LOF.FrameShift.Missense.Obs, y=expected.mutation_rate, by.x="Gene.Name", by.y="Gene", all.y=T)
LOF.FrameShift.Missense.Obs.Mut[ is.na( LOF.FrameShift.Missense.Obs.Mut ) ] <- 0 


Observed.Variant.Number <- base::rowSums( LOF.FrameShift.Missense.Obs.Mut[ , c("LOF.Observed.Variant.Number", "FrameShift.Observed.Variant.Number", "Missense.Observed.Variant.Number")] )

Expected.Gene.Mutation <- base::rowSums( LOF.FrameShift.Missense.Obs.Mut[ , c("LOF.mutation_rate", "FrameShift.mutation_rate", "Missense.mutation_rate")] )

Total.Variant.Number <- sum(Observed.Variant.Number)

model <- lm( Observed.Variant.Number ~ poly( Expected.Gene.Mutation , 2 ) )
fitted.value <- fitted( model )

Expected.Variant.Number <- c()
for( i in 1 : nrow( LOF.FrameShift.Missense.Obs.Mut ) ) {
	Expected.Variant.Number <- append( Expected.Variant.Number , Total.Variant.Number * ( fitted.value[ i ] / sum( fitted.value ) ) )
}

Enrichment <- Observed.Variant.Number / Expected.Variant.Number


###########################
# Significance
###########################
print( "Step 7: Significance" )

P.value <- c()
for( i in 1 : nrow( LOF.FrameShift.Missense.Obs.Mut ) ) { 
	
	suppressWarnings(
	
	P.value <- append( P.value , 
	( ppois( q = Observed.Variant.Number[i] - 1 , 
	lambda = Expected.Variant.Number[i] ,
	lower.tail = FALSE ) ) )
	
	)
	
}
names( P.value ) <- LOF.FrameShift.Missense.Obs.Mut$Gene.Name

exp.pval <- ( rank( P.value , ties.method="first" )+.5 ) / ( length( P.value )+1 )

Obs.Exp.matrix <- cbind( LOF.FrameShift.Missense.Obs.Mut, Observed.Variant.Number, Expected.Gene.Mutation, fitted.value, Expected.Variant.Number, Enrichment, P.value )

Obs.Exp.Pred <- merge( x=merge(x=Obs.Exp.matrix, y=LOF.sum.mut.gene.Count, by.x="Gene.Name", by.y="Gene", all.x=T), y=Missense.sum.mut.gene.Count, by.x="Gene.Name", by.y="Gene", all.x=T) 

Obs.Exp.Pred$Gene.Name <- paste( "'" , Obs.Exp.Pred$Gene.Name , "'" , sep="" )


###########################
# P value table
###########################
print( "Step 8: P value table" )

if( opt$pli.source == "exac" ) {

	pli.score <- read.table("/media/Dazs/youngji/Resource/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt", header=T, sep="\t")
	pli.score$Gene.Name <- paste( "'" , pli.score$Gene.Name , "'" , sep="" )
	constraint.matrix <- pli.score[ , c("Gene.Name", "pLI", "mis_z")]
	anno.final.output <- merge( x=Obs.Exp.Pred , y=constraint.matrix, by="Gene.Name", all.x=T)
	
} else if( opt$pli.source == "gnomad" ) {

	pli.score <- read.table("/media/Dazs/youngji/Resource/gnomad.v2.1.1.lof_metrics.by_gene.txt", header=T, sep="\t")
	pli.score$gene <- paste( "'" , pli.score$gene , "'" , sep="" )
	temp.constraint.matrix <- pli.score[ , c("gene", "pLI", "mis_z")]
	constraint.matrix <- temp.constraint.matrix %>% group_by(gene) %>% filter(mis_z == max(mis_z)) %>% filter(pLI == max(pLI)) 
	anno.final.output <- merge( x=Obs.Exp.Pred , y=constraint.matrix, by.x="Gene.Name", by.y="gene", all.x=T)

}

uniq.anno.final.output <- uniquecombs(anno.final.output)

sorted.anno.final.output <- uniq.anno.final.output[ order( uniq.anno.final.output$P.value ) , ]

FINAL.output <- sorted.anno.final.output[ , c("Gene.Name", "LOF.Observed.Variant.Number", "FrameShift.Observed.Variant.Number", "Missense.Observed.Variant.Number", "Observed.Variant.Number", "LOF.mutation_rate", "FrameShift.mutation_rate", "Missense.mutation_rate", "Expected.Gene.Mutation", "Expected.Variant.Number", "fitted.value", "Enrichment", "P.value", "LOF.Predicted.Variants", "Missense.Predicted.Variants", "pLI", "mis_z") ]

dir.create(opt$output.dir, showWarnings=F)
setwd(opt$output.dir)
sink('analysis-output.log', append=TRUE)
print(opt)
print( paste( "X_linked.adjustment_factor:", X_linked.adjustment_factor ) )
print( paste( "ratio of frameshift to nonsense:", adj.factor ) )
sink()

print( "Top genes" )
head(FINAL.output)
write.csv( FINAL.output , "Final.Result.Pvalue.csv" , row.names=F , quote=F )

###########################
# qq plot
###########################
print( "Step 9: qq plot" )

burden.test.plot.input <- sorted.anno.final.output[ , c("Gene.Name", "P.value") ]
colnames(burden.test.plot.input) <- c("Gene", "P_DOM")
tempFile <- "Final.Result.burden.txt"
write.table( burden.test.plot.input, tempFile, sep="\t", quote=F, row.names=F )

system(paste("Rscript /media/Dazs/youngji/Scripts/trapd_base_QQ.R", "--prefix", "Final.Result"))
invisible(file.remove( tempFile ))