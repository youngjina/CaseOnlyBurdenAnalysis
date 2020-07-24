import os, re, sys

inDir = sys.argv[1]
in_postfix = sys.argv[2]
outDir = sys.argv[3]

out_postfix = re.sub( "sites_annotated" , "mutation_rate" , in_postfix )
if not os.path.exists( outDir ):	os.makedirs( outDir )

os.chdir( inDir )

try:
	chromosomeNumbers = range(1,23)
	chromosomes = [str(chromosome) for chromosome in chromosomeNumbers]
	chromosomes.append("X")
	for chromosome in chromosomes:
		print( "# chr" + str( chromosome ) )
		inFile = "%s%s_%s" % ( "chr" , chromosome , in_postfix )
		outFile = "%s%s%s_%s" % ( outDir , "chr" , chromosome , out_postfix )
		fout = open( outFile , 'w' )
		with open( inFile , 'r' ) as fin:
			header_line = next( fin )
			for line in fin:
				temp = line.strip().split('\t')
				chrom = temp[0]
				pos = temp[1]
				ref = temp[2]
				alt = temp[3]
				mut_rate = temp[5]
				output = "%s-%s-%s-%s\t%s\n" % ( chrom, pos , ref , alt , mut_rate )
				fout.write( output )	

except (BrokenPipeError, IOError):
    pass
