import os, sys

jobDir = sys.argv[1]
inFile = sys.argv[2]
outFile = sys.argv[3]

os.chdir( jobDir )
fout = open( outFile , 'w' )
with open( inFile , 'r' ) as f:
	for line in f:
		temp = line.strip().split('\t')
		varID = temp[0]
		mut_rate = temp[1]
		chr_pos = '-'.join( varID.split('-')[:2] )
		output = "%s\t%s\t%s\n" % ( varID , mut_rate , chr_pos )
		fout.write( output )
