import os, sys, subprocess

def doMutationType(cohort, mutationType):
	job_dir_1 = "/media/Jerry/youngji/Coverage_Analysis/Coverage_Summary_{cohort}_min-coverage_10/".format(cohort = cohort)
	print( job_dir_1 )
	os.chdir( job_dir_1 )
	print( "cmd_1" )
		
	cmd_2 = "cp -r Absent_gnomAD_{mutationType}_variantID Absent_gnomAD_{mutationType}_input_anno".format(mutationType = mutationType)
	p_2 = subprocess.Popen( cmd_2, shell=True)
	p_2.wait()
	print( "cmd_2" )

	cmd_3 = "python input_mr_eel.py /media/Jerry/youngji/Coverage_Analysis/Coverage_Summary_{cohort}_min-coverage_10/Absent_gnomAD_{mutationType}_input_anno/ Absent_gnomAD_{mutationType}_Covered_90_percent_Positions_{cohort}_min-coverage_10_variantID.txt".format(mutationType = mutationType, cohort = cohort)
	p_3 = subprocess.Popen( cmd_3, shell=True)
	p_3.wait()
	print( "cmd_3" )

	cmd_4 = "python batch_mr_eel.py /media/Jerry/youngji/Coverage_Analysis/Coverage_Summary_{cohort}_min-coverage_10/Absent_gnomAD_{mutationType}_input_anno/ Absent_gnomAD_{mutationType}_Covered_90_percent_Positions_{cohort}_min-coverage_10_variantID.txt /media/Jerry/youngji/Coverage_Analysis/Coverage_Summary_{cohort}_min-coverage_10/Absent_gnomAD_{mutationType}_output_anno/".format(mutationType = mutationType, cohort = cohort)
	p_4 = subprocess.Popen( cmd_4, shell=True)
	p_4.wait()
	print( "cmd_4" )

	cmd_5 = "python summary_mr_eel.py /media/Jerry/youngji/Coverage_Analysis/Coverage_Summary_{cohort}_min-coverage_10/Absent_gnomAD_{mutationType}_output_anno/ Absent_gnomAD_{mutationType}_Covered_90_percent_Positions_{cohort}_min-coverage_10_sites_annotated.txt /media/Jerry/youngji/Coverage_Analysis/Coverage_Summary_{cohort}_min-coverage_10/Absent_gnomAD_{mutationType}_anno_mut_rate/".format(mutationType = mutationType, cohort = cohort)
	p_5 = subprocess.Popen( cmd_5, shell=True)
	p_5.wait()
	print( "cmd_5" )

	job_dir_2 = "/media/Jerry/youngji/Coverage_Analysis/Coverage_Summary_{cohort}_min-coverage_10/Absent_gnomAD_{mutationType}_anno_mut_rate/".format(mutationType = mutationType, cohort = cohort)
	print( job_dir_2 )
	os.chdir( job_dir_2 )
	print( "cmd_6" )

	cmd_7 = "cat *mutation_rate.txt | uniq > ../Absent_gnomAD_{mutationType}_{cohort}_Carlson_mutability_90p_covered.txt".format(mutationType = mutationType, cohort = cohort)
	p_7 = subprocess.Popen( cmd_7, shell=True)
	p_7.wait()
	print( "cmd_7" )

	cmd_8 = "python mutation_rate_variantID.py /media/Jerry/youngji/Coverage_Analysis/Coverage_Summary_{cohort}_min-coverage_10/ Absent_gnomAD_{mutationType}_{cohort}_Carlson_mutability_90p_covered.txt Absent_gnomAD_{mutationType}_{cohort}_variantID_Carlson_mutability_90p_covered.txt".format(mutationType = mutationType, cohort = cohort)
	p_8 = subprocess.Popen( cmd_8, shell=True)
	p_8.wait()
	print( "cmd_8" )

	cmd_9 = "Rscript combine_positions_mutability.R /media/Jerry/youngji/Coverage_Analysis/Coverage_Summary_{cohort}_min-coverage_10/ Covered_90_percent_Positions_Gene_variantID_{cohort}_min-coverage_10_site.summary.txt Absent_gnomAD_{mutationType}_{cohort}_variantID_Carlson_mutability_90p_covered.txt Absent_gnomAD_{mutationType}_{cohort}_merged_Carlson_MutRate_Gene_variantID_min-coverage_10_site.summary.txt".format(mutationType = mutationType, cohort = cohort)
	p_9 = subprocess.Popen( cmd_9, shell=True)
	p_9.wait()
	print( "cmd_9" )


if __name__ == "__main__":
	
	cohort = sys.argv[1]
	mutationType = sys.argv[2]
	print("In function")
	doMutationType(cohort, mutationType)
