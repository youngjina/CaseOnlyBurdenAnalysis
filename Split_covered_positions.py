import os, sys

jobDir = sys.argv[1]
inFile = sys.argv[2]
prefix = sys.argv[3]

os.chdir( jobDir )

with open( inFile , "r" ) as f:
	header = next( f )
	for line in f:
		chrom = line.strip().split('-')[0]
		outFile = "%s%s%s%s%s" % ( "chr", chrom, "_Covered_90_percent_Positions_", prefix , "_min-coverage_10_chr_pos.txt" )
		if os.path.exists(outFile):
			append_write = 'a'
		else:
			append_write = 'w'
		fout = open( outFile , append_write )
		fout.write( line )		

os.system( "mkdir Covered_chr_pos" )
os.system( "mv chr*_chr_pos.txt ./Covered_chr_pos" )