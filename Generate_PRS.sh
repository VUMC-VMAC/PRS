#!/bin/bash

#fail on error
set -e

############################# Verify inputs ############################

#define usage
display_usage() {
    printf "Polygenic Risk Score Calculation Script

Usage:
Generate_PRS.sh -i GENOTYPES_STEM -s SUMMARY_STATISTICS -f OUTPUT_FOLDER -t OUTPUT_TAGS -o OUTPUT_STEM -p PVALUE_THRESHOLDS -r R2_THRESHOLD -w CLUMPING_WINDOW -b GENOME_BUILD -a -l -m MEMORY_LIMIT

GENOTYPES_STEM = the full path and stem to the genotypes for individuals to calculate the risk scores in (must be in plink binary format)
SUMMARY_STATISTICS = the full path and file name for the summary statistics to calculate the risk scores from (must have the columns SNP, A1, BETA/OR, and P, but can have other columns as well); note that you can specify multiple summary statistics files here separated by a comma (no space); if the -l flag is used, the results are assumed to be from logistic regression and having an odds ratio (OR) column; otherwise, results are assumed to be from a regular linear model and having the BETA column.
OUTPUT_FOLDER = the folder in which to calculate the PRS
OUTPUT_TAGS = a tag for each summary statistics file supplied (comma separated with no space); will be in the variable name in the outputted file with all the PRS
OUTPUT_STEM = name for the output file of all the PRS calculated
PVALUE_THRESHOLDS = comma separated p value thresholds for variants to be included in PRS calculation; default is plink's default of 0.0001 (optional)
R2_THRESHOLD = r-squared threshold for clumping; default is plink's default of 0.5 (optional)
CLUMPING_WINDOW = window in which to conduct the clumping; default is plink's default of 250 (optional)
GENOME_BUILD = the genome build for the input genotypes; only used if the -a flag is given; default is b38 (optional)
MEMORY_LIMIT = memory limit for plink calls in MB; especially helpful when running on badger rather than on slurm (optional)

-a indicates to exclude the APOE region from PRS calculation
-l indicates that the summary stats are from logistic regression, in which case the weight will be pulled from the OR column

-h will show this usage
"
        }

#set defaults (will be written over if the argument was supplied)
pvalues=0.0001
r2thresh=0.5
window=250
apoe_exclude=no
genome_build=b38
memory_limit=""
logistic=no

#parse arguments
while getopts 'i:s:f:t:o:p:r:w:b:m:alh' flag; do
  case "${flag}" in
    i) genotypes="${OPTARG}" ;;
    s) sumstats="${OPTARG}" ;;
    f) output_folder="${OPTARG}" ;;
    t) output_tags="${OPTARG}" ;;
    o) output="${OPTARG}" ;;
    p) pvalues="${OPTARG}" ;;
    r) r2thresh="${OPTARG}" ;;
    w) window="${OPTARG}" ;;
    b) genome_build="${OPTARG}" ;;
    m) memory_limit="${OPTARG}";;
    a) apoe_exclude=yes ;;
    l) logistic=yes ;;
    h) display_usage ; exit ;;
    \?|*) display_usage
       exit 1;;
  esac
done

#print out inputs to the log
printf "Polygenic Risk Score Calculation Script

Input genotypes: $genotypes
Summary statistics: $sumstats
Output folder: $output_folder
PRS output labels: $output_tags
Output stem: $output
P-value thresholds: $pvalues
R-squared threshold: $r2thresh
Window size: $window
"

if [ "$apoe_exclude" = "yes" ];
then 
    printf "PRS will be calculated with and without the APOE region (using coordinates from genome build ${genome_build}).\n"
fi

if [ "$logistic" = "no" ];
then 
    printf "Because the -l flag was not supplied, the summary statistics are assumed to be linear regression results with the test statistic to be used assumed to be a BETA.\n"
elif [ "$logistic" = "yes" ];
    printf "Because the -l flag was supplied, the summary statistics will be treated as logistic regression with the test statistic to be used assumed to be OR.\n"
fi

if [ "$memory_limit" ];
then 
    printf "Memory limit for plink calls: $memory_limit \n"
    memory_limit=$( echo "--memory $memory_limit" )
fi

#make sure the script is being run from the scripts directory
if [ ! -f Determine_overlapping_SNPs.R ];
then 
    printf "The working directory does not contain the scripts! Please rerun from within the scripts directory for PRS calculation!\n"
    exit 1
fi

#get genotypes_stem (sans path) using shell magic
genotypes_stem=${genotypes##*/}

#split summary stats and output tags into arrays for easy parsing
sumstats_array=(${sumstats//,/ })
output_tag_array=(${output_tags//,/ })

############################# Start actually doing things ############################

for i in ${!sumstats_array[@]};
do
    sumstats_current=${sumstats_array[i]}
    output_tag_current=${output_tag_array[i]}

    printf "\nStep 1: Determining overlapping variants between genotypes and $output_tag_current summary stats

current summary stats: $sumstats_current
current output tag: $output_tag_current\n"

    #check to make sure required columns are present
    if [[ $( head -n1 $sumstats_current | grep -o "SNP" ) ]] &&  [[ $( head -n1 $sumstats_current | grep -o "A1" ) ]] && [[ $( head -n1 $sumstats_current | grep -o "P" ) ]]; 
    then
	if [ "$logistic" = "no" ] && [[ $( head -n1 $sumstats_current | grep -o "BETA" ) ]]; 
	then 
	    printf "All necessary columns are present in the summary stats.\n" ; 
	elif [ "$logistic" = "yes"] && [[ $( head -n1 $sumstats_current | grep -o "OR" ) ]];
	    printf "All necessary columns are present in the summary stats.\n" ; 
	else
	    printf "One or more of the necessary columns (SNP, A1, BETA/OR, P) not present in the summary stats! Please confirm the columns are present and are named correctly and retry.\n"
        exit 1 ;
	fi
    else 
	printf "One or more of the necessary columns (SNP, A1, BETA/OR, P) not present in the summary stats! Please confirm the columns are present and are named correctly and retry.\n" 
	exit 1 ; 
    fi

    #get overlapping, non-palindromic variants
    Rscript Determine_overlapping_SNPs.R $sumstats_current $genotypes ${output_folder}/$output_tag_current $logistic
    #update summary stats variable to refer to the one updated for allele flips
    sumstats_current=${output_folder}/${output_tag_current}_summary_stats_updated.txt

    #Extract overlapping SNPs from the genotype data
    plink --bfile $genotypes --allow-no-sex --extract ${output_folder}/${output_tag_current}_overlapping_SNPs.txt --make-bed --out ${output_folder}/${genotypes_stem}_${output_tag_current} ${memory_limit} > /dev/null 2>&1 
    genotypes_new=${output_folder}/${genotypes_stem}_${output_tag_current}
    
    #check to see if the R script identified variants to 
    if [ -f "${output_folder}/${output_tag_current}_SNPs_to_flip.txt" ];
    then
	printf "Flipping variants...\n\n"
	#flip the variants identified by script to flip
	plink --bfile $genotypes_new --flip ${output_folder}/${output_tag_current}_SNPs_to_flip.txt --make-bed --out ${genotypes_new}_flip ${memory_limit} > /dev/null 2>&1
	genotypes_new=${genotypes_new}_flip 
    fi

    printf "Step 2: Performing LD clumping for $output_tag_current \n"
    #get the largest p value 
    maxp=$( echo $pvalues | sed 's/,/\n/g' | sort -n | tail -n1 )
    #Perform LD clumping, with the p value threshold of the largest p value
    plink --bfile ${genotypes_new} --allow-no-sex --clump ${output_folder}/${output_tag_current}_summary_stats_updated.txt --clump-p1 $maxp --clump-r2 $r2thresh --clump-kb $window --out ${genotypes_new} ${memory_limit} > /dev/null 2>&1
    
    #Create the input file for the score calculation
    Rscript Generate_score_input_file.R ${genotypes_new}.clumped $sumstats_current $pvalues

    #create range file, with one line for each p value threshold
    > ${output_folder}/${output_tag_current}_pvalue_range.txt
    for p in $( echo $pvalues | sed 's/,/ /g' ); do echo "Pval_$p 0 $p" >> ${output_folder}/${output_tag_current}_pvalue_range.txt ; done

    printf "Step 3: Calculating PRS for ${output_tag_current}\n"
    #Calculate PRS
    plink --bfile $genotypes_new --allow-no-sex --score ${genotypes_new}_score_input.txt --q-score-range ${output_folder}/${output_tag_current}_pvalue_range.txt ${genotypes_new}_score_input.txt 1 4 --out ${genotypes_new}_PRS ${memory_limit} > /dev/null 2>&1
    PRS_stem=${genotypes_new}_PRS

    #report the number of variants that were skipped, if any were skipped
    skipped_var=$( grep "lines skipped in --score" ${genotypes_new}_PRS.log | grep -o -E '[0-9]+' | head -n1 )
    if [ ! -z "$skipped_var" ]; then printf "$skipped_var variants skipped from PRS calculation. See ${genotypes_new}_PRS.log for more details.\n"; fi

    #if the no-apoe flag was specified, then recalculate the PRS without APOE
    if [ "$apoe_exclude" == "yes" ];
    then
	printf "Calculating PRS for ${output_tag_current} without the APOE region\n"
	#Calculate PRS
	plink --bfile $genotypes_new --allow-no-sex --exclude range APOE_region_${genome_build}.txt --score ${genotypes_new}_score_input.txt --q-score-range ${output_folder}/${output_tag_current}_pvalue_range.txt ${genotypes_new}_score_input.txt 1 4 --out ${genotypes_new}_noAPOE_PRS ${memory_limit} > /dev/null 2>&1
    fi
        
    #check for variants that failed to be incorporated into the score
    if [ -f  ${PRS_stem}.nopred ];
    then
	printf "Some variants had mismatched alleles for the PRS. Attempting to flip the strands and rebuild the score...\n"
    
	#Pull out list of SNPs with allele code mismatches
	awk '{print $2}' ${PRS_stem}.nopred > ${PRS_stem}.flipsnps
	
	#Flip strand for list of SNPs
	plink --bfile $genotypes_new --allow-no-sex --flip ${PRS_stem}.flipsnps --make-bed --out ${genotypes_new}_flipsnps ${memory_limit} > /dev/null 2>&1

	#Create .profile file
	plink --bfile ${genotypes_new}_flipsnps --allow-no-sex --score ${genotypes_new}_score_input.txt --q-score-range ${output_folder}/${output_tag_current}_pvalue_range.txt ${genotypes_new}_score_input.txt 1 4 --out $PRS_stem ${memory_limit} > /dev/null 2>&1
	#report the number of variants that were skipped, if any were skipped
	skipped_var=$( grep "lines skipped in --score" ${genotypes_new}_PRS.log | grep -o -E '[0-9]+' | head -n1 )
	if [ ! -z "$skipped_var" ]; then printf "$skipped_var variants skipped from PRS calculation. See ${genotypes_new}_PRS.log for more details.\n"; fi


	#if the no-apoe flag was specified, then recalculate the PRS without APOE
	if [ "$apoe_exclude" == "yes" ];
	then
	    printf "Recalculating the PRS excluding the APOE region.\n"
            #Calculate PRS
            plink --bfile ${genotypes_new}_flipsnps --allow-no-sex --exclude range APOE_region_${genome_build}.txt --score ${genotypes_new}_score_input.txt --q-score-range ${output_folder}/${output_tag_current}_pvalue_range.txt ${genotypes_new}_score_input.txt 1 4 --out ${genotypes_new}_noAPOE_PRS ${memory_limit} > /dev/null 2>&1
	fi


    fi

##Clean-up
#rm ${genotypes_new}*.bed ${genotypes_new}*.bim ${genotypes_new}*.fam ${genotypes_new}.clumped
rm ${genotypes_new}*.bed

done

printf "PRS generation complete! Find the PRS in the ${output_folder}/*.profile files.\n"

# printf "Step 4: Combining PRS into one file\n"

# #combine all PRS into one file
# ##but first echo the inputs
# printf "
# Inputs for the combo file:
# Output folder: $output_folder
# Genotypes stem: $genotypes_stem
# Output tags: $output_tags
# P values: $pvalues
# Output name: $output
# Exclude APOE? $apoe_exclude
# "
# Rscript Compile_PRS.R $output_folder $genotypes_stem $output_tags $pvalues $output $apoe_exclude

# printf "Step 5: Generating correlation plot of all generated scores\n"

# Rscript Create_PRS_corplot.R ${output_folder}/${output}.txt
