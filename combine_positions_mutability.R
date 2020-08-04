args = commandArgs(trailingOnly=TRUE)
jobDir = args[1]
posFile = args[2]
mutFile = args[3]
outFile = args[4]
setwd( jobDir )

pos.data <- read.table( posFile , header=F , sep="\t" )
mut.data <- read.table( mutFile , header=F , sep="\t" )

colnames(pos.data) <- c( "Chr.Pos" , "Gene" )
colnames(mut.data) <- c( "Variant.ID" , "Mut.Rate" , "Chr.Pos" )
merged.data <- merge( x = mut.data , y = pos.data , by = "Chr.Pos" , all.x = T )
write.table( merged.data , outFile, row.names=F, quote=F, sep="\t" )