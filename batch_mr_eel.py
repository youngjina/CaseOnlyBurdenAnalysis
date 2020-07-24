import os, re, sys, subprocess

inDir = sys.argv[1]
in_postfix = sys.argv[2]
outDir = sys.argv[3]

out_postfix = re.sub( "variantID" , "sites_annotated" , in_postfix )
if not os.path.exists( outDir ):	os.makedirs( outDir )

try:
	chromosomeNumbers = range(1,23)
	chromosomes = [str(chromosome) for chromosome in chromosomeNumbers]
	chromosomes.append("X")
	for chromosome in chromosomes:
		print( chromosome )
		cmd = "%s%s%s%s_%s%s%s%s%s_%s" % ( "perl /media/Dazs/youngji/Resource/mr-eel/mr_eel.pl --in " , inDir , "chr" , chromosome , in_postfix , " --seq --rates /media/Dazs/youngji/Resource/mr-eel/ERV_7bp_rates.txt --ref /media/Dazs/youngji/Resource/human_g1k_v37.fasta > " , outDir , "chr" , chromosome , out_postfix )
		print( cmd )
		process = subprocess.Popen( cmd, shell=True )
		process.wait()
	
except (BrokenPipeError, IOError):
    pass

#sys.stderr.close()