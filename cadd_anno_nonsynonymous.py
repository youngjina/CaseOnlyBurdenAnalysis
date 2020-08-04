import gzip, os, sys

jobDir = sys.argv[1]
inFile = sys.argv[2]
outFile = sys.argv[3]

os.chdir( jobDir )
fout = open( outFile , 'w' )
with gzip.open( inFile , 'r' ) as fin:        
	for line in fin:
		if line.startswith("#") == 0:
			temp = line.strip().split("\t")
			chrom = temp[0]
			pos = temp[1]
			ref = temp[2]
			alt = temp[3]
			cons = temp[7]
			if cons == "NON_SYNONYMOUS":
				variantID  = "%s-%s-%s-%s" % ( chrom , pos , ref , alt )
				output = "%s\t%s\n" % ( variantID , cons )				
				fout.write( output )	
