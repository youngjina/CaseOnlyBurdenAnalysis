### Write out positions in a CADD file that are known to be well-covered.
import os, sys

### Get well-covered positions for a chromosome:
### Input: well-covered positions file-name.
### Output: list of well-covered positions.

def getWellCoveredPositions(coveredPositionsFileName):

	print( "Getting covered positions..." )

	coveredPositions = set()
	with open(coveredPositionsFileName, "r" ) as f:
		for line in f:
			covered_position = line.strip()
			coveredPositions.add(covered_position)

	return coveredPositions


### Read through CADD file and write out the positions (in one chromosome) that are well-covered.
### Input: CADD filename, output filename, well-covered positions set.
### Output: nothing.

def writeCoveredPositionsFromCADDfile(caddFilename, outputFileName, coveredPositions):

	print( "Reading from CADD..." )

	with open(caddFilename , "r") as caddFile, open(outputFileName, "w") as out:
		for cadd_line in caddFile:
			cadd_temp = cadd_line.strip().split("-")
			cadd_chrom_pos = "-".join(cadd_temp[:2])
			cadd_variantID = cadd_line.strip()

			if cadd_chrom_pos in coveredPositions:
				out.write(cadd_variantID + "\n")
				

def doAllChromosomes(cohort, mutationType):

	chromosomeNumbers = range(1,23)
	chromosomes = [str(chromosome) for chromosome in chromosomeNumbers]
	chromosomes.append("X")
	for chromosome in chromosomes:

		print( "Working on chromosomes {}...".format(chromosome) )

		coveredPositionsFileName = "/media/Jerry/youngji/Coverage_Analysis/Coverage_Summary_{cohort}_min-coverage_10/Covered_chr_pos/chr{chr}_Covered_90_percent_Positions_{cohort}_min-coverage_10_chr_pos.txt".format(cohort = cohort, chr = chromosome)

		caddFilename = "/media/Dazs/youngji/Resource/absent_{mutationType}/absent_{mutationType}.chr{chr}.txt".format(chr = chromosome, mutationType = mutationType)

		outputDirectory = "/media/Jerry/youngji/Coverage_Analysis/Coverage_Summary_{cohort}_min-coverage_10/Absent_gnomAD_{mutationType}_variantID/".format(mutationType = mutationType, cohort = cohort)

		if not os.path.exists(outputDirectory):
			os.makedirs(outputDirectory)

		outputFileName = outputDirectory + "chr{chr}_Absent_gnomAD_{mutationType}_Covered_90_percent_Positions_{cohort}_min-coverage_10_variantID.txt".format(chr = chromosome, mutationType = mutationType, cohort = cohort)

		coveredPositions = getWellCoveredPositions(coveredPositionsFileName)
		writeCoveredPositionsFromCADDfile(caddFilename, outputFileName, coveredPositions)


if __name__ == "__main__":
	
	cohort = sys.argv[1]
	mutationType = sys.argv[2]
	doAllChromosomes(cohort, mutationType)