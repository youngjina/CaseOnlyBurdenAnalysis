import os, sys, subprocess

jobDir = sys.argv[1]
postfix = sys.argv[2]

os.chdir( jobDir )

try:
	chromosomeNumbers = range(1,23)
	chromosomes = [str(chromosome) for chromosome in chromosomeNumbers]
	chromosomes.append("X")
	for chromosome in chromosomes:
		cmd = "%s%s_%s" % ( "vi -c \"%s/-/\t/g\" -c \"wq\" chr" , chromosome , postfix )
		process = subprocess.Popen( cmd, shell=True )
		process.wait()	

except (BrokenPipeError, IOError):
    pass