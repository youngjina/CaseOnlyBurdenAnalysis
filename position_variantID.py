import os, sys

jobDir = sys.argv[1]
inFile = sys.argv[2]
outFile = sys.argv[3]

os.chdir( jobDir )
fout = open( outFile , 'w' )
with open( inFile , 'r' ) as f:
	header_line = next( f )
	for line in f:
		temp = line.strip().split(',')
		gene = temp[0]
		chrom= temp[1]
		pos  = temp[2]
		output = "%s-%s\t%s\n" % ( chrom , pos , gene )
		fout.write( output )
