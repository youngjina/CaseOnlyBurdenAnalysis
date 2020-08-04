import os, sys

jobDir = sys.argv[1]
inFile = sys.argv[2]
mutationType = sys.argv[3]

if not os.path.exists( jobDir ):	os.makedirs( jobDir )
os.chdir( jobDir )

with open( inFile , "r" ) as f:
	for line in f:
		chrom = line.strip().split('-')[0]
		outFile = "%s%s%s%s%s" % ( "chr", chrom, ".", mutationType, "_cadd_anno_variant.txt" )
		if os.path.exists(outFile):
			append_write = 'a'
		else:
			append_write = 'w'
		fout = open( outFile , append_write )
		fout.write( line )		
