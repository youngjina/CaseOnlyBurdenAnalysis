# Usage: Rscript CaseOnlyBurdenAnalysis_Missense.R --help

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

# Count the number of observed variants per gene
observed.variant.gene.pair <- observed.genotypes[ , c("Variant.ID", "Gene.Name") ]
uniq.observed.variant.gene.pair <- observed.variant.gene.pair[ !duplicated( observed.variant.gene.pair$Variant.ID ) , ]
observed.variant.per.gene <- plyr::count( uniq.observed.variant.gene.pair, "Gene.Name" )
observed.variant.per.gene$Gene.Name <- gsub( "'", "", observed.variant.per.gene$Gene.Name )

out <- strsplit( as.character( observed.genotypes$Variant.ID ) , "-" )
chrom <- data.frame(t(sapply(out, `[`)))[,1]
Xlinked.genes <- as.character( gsub( "'", "", unique( observed.genotypes[ which( chrom == "X" ) , "Gene.Name" ] ) ))

# X-linked genes
Xlinked.observed.variant.per.gene <- observed.variant.per.gene %>% filter( Gene.Name %in% Xlinked.genes )

# Autosome genes
autosome.observed.variant.per.gene <- observed.variant.per.gene %>% filter( !( Gene.Name %in% Xlinked.genes ) )


##################################################################
# Mutation rate
##################################################################
print( "Step 2: Mutation rate" )

mut.data <- fread( opt$mutation.rate.file , header=T , sep="\t" )

uniq.variant.mut.data <- mut.data[ !duplicated( mut.data$Variant.ID ) , ]
uniq.variant.mut.gene <- uniq.variant.mut.data[ , c( "Gene", "Variant.ID", "Mut.Rate" )  ]
uniq.variant.mut.gene$Gene <- gsub( "'", "", uniq.variant.mut.gene$Gene )
sum.mut.gene.Count <- uniq.variant.mut.gene %>% group_by(Gene) %>% tally()
colnames(sum.mut.gene.Count) <- c( "Gene", "Predicted.Variants" )

# Aggregate mutation rate per position into mutation rate per gene by sum
sum.mut.gene <- uniq.variant.mut.gene %>% group_by(Gene) %>% tally(Mut.Rate)
expected.mutation_rate <- sum.mut.gene$n / mutation_ratio.adjustment_factor
names(expected.mutation_rate) <- sum.mut.gene$Gene

expected.mutation_rate <- rownames_to_column(as.data.frame(expected.mutation_rate))
colnames(expected.mutation_rate) <- c("Gene", "Expected.Gene.Mutation")


##################################################################
# Expected number of variants
##################################################################
print( "Step 3: Expected variants" )

observed.autosome.Xlinked <- rbind( autosome.observed.variant.per.gene , Xlinked.observed.variant.per.gene )

# Merged data
input.model <- merge( x = observed.autosome.Xlinked , y = expected.mutation_rate , by.x = "Gene.Name" , by.y = "Gene" , all.y = T )
input.model[ is.na( input.model ) ] <- 0

Observed.Variant.Number <- input.model$freq
Expected.Gene.Mutation <- input.model$Expected.Gene.Mutation
Total.Variant.Number <- sum(Observed.Variant.Number)

model <- lm( Observed.Variant.Number ~ poly( Expected.Gene.Mutation , 2 ) )
fitted.value <- fitted( model )

Expected.Variant.Number <- c()
for( i in 1 : nrow( input.model ) ) {
	Expected.Variant.Number <- append( Expected.Variant.Number , Total.Variant.Number * ( fitted.value[ i ] / sum( fitted.value ) ) )
}

Enrichment <- Observed.Variant.Number / Expected.Variant.Number


###########################
# Significance
###########################
print( "Step 4: Significance" )

P.value <- c()
for( i in 1 : nrow( input.model ) ) { 

	suppressWarnings(
	
	P.value <- append( P.value , 
	( ppois( q = Observed.Variant.Number[i] - 1 , 
	lambda = Expected.Variant.Number[i] ,
	lower.tail = FALSE ) ) )
	
	)
	
}
names( P.value ) <- input.model$Gene.Name

exp.pval <- ( rank( P.value , ties.method="first" )+.5 ) / ( length( P.value )+1 )

Obs.Exp.matrix <- cbind( input.model, fitted.value, Expected.Variant.Number, Enrichment, P.value )
colnames(Obs.Exp.matrix) <- c("Gene.Name", "Observed.Variant.Number", "Expected.Gene.Mutation", "fitted.value", "Expected.Variant.Number", "Enrichment", "P.value")
Obs.Exp.Pred <- merge(x=Obs.Exp.matrix, y=sum.mut.gene.Count, by.x="Gene.Name", by.y="Gene", all.x=T)

Obs.Exp.Pred$Gene.Name <- paste( "'" , Obs.Exp.Pred$Gene.Name , "'" , sep="" )


###########################
# P value table
###########################
print( "Step 5: P value table" )

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
FINAL.output <- sorted.anno.final.output[ , c("Gene.Name", "Observed.Variant.Number", "Expected.Gene.Mutation", "Expected.Variant.Number", "fitted.value", "Enrichment", "P.value", "Predicted.Variants", "pLI", "mis_z") ]

dir.create(opt$output.dir, showWarnings=F)
setwd(opt$output.dir)
sink('analysis-output.log', append=TRUE)
print(opt)
print( paste( "X_linked.adjustment_factor:", X_linked.adjustment_factor ) )
sink()

print("Top genes")
head(FINAL.output)
write.csv( FINAL.output , "Final.Result.Pvalue.csv" , row.names=F , quote=F )


###########################
# qq plot
###########################
print( "Step 6: QQ plot" )

burden.test.plot.input <- sorted.anno.final.output[ , c("Gene.Name", "P.value") ]
colnames(burden.test.plot.input) <- c("Gene", "P_DOM")
tempFile <- "Final.Result.burden.txt"
write.table( burden.test.plot.input, tempFile, sep="\t", quote=F, row.names=F )

system(paste("Rscript /media/Dazs/youngji/Scripts/trapd_base_QQ.R", "--prefix", "Final.Result"))
invisible(file.remove( tempFile ))