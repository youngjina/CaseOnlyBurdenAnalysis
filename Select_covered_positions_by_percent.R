args = commandArgs(trailingOnly=TRUE)
jobDir <- args[1]
CovSum_outFile <- args[2]
Filtered_outFile <- args[3]
n.case <- as.numeric( args[4] )
n.percent <- as.numeric( args[5] )
pheno <- args[6]

print( "Parameters" )
print( args )

setwd( jobDir )
print( getwd() )

if( pheno == "case" ) {
	selectCol = c(4:8)
} else if( pheno == "control" ) {
	selectCol = c(9:13)
}
print( selectCol )
inFile <- Sys.glob( "*_site.summary.csv" )
print( inFile )
input.data <- read.csv( inFile , head=T )
print( dim( input.data ) )

coverage.sum <- apply( input.data[ , selectCol ] , 1 , sum )
coverage.percent <- coverage.sum / n.case *100
merged.data <- cbind( input.data , coverage.sum , coverage.percent )
print( head( merged.data ) )
write.csv( merged.data , CovSum_outFile , row.names=F , quote=F )

covered.data <- subset( merged.data , coverage.percent >= n.percent )
print( head( covered.data ) )
print( nrow( covered.data ) )
write.csv( covered.data , Filtered_outFile , row.names=F , quote=F )

print( nrow( covered.data ) / nrow( input.data ) * 100 )
