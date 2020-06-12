#!/bin/bash

#fail on error
set -e

#assumptions at the beginning of this script:
# the summary stats are already in shape to be able to figure out overlap
# will still need to determine the overlapping variants

############################# Verify inputs ############################

#parse arguments
while getopts 'i:s:f:t:o:p:r:w:h' flag; do
  case "${flag}" in
    i) genotypes="${OPTARG}" ;;
    s) sumstats="${OPTARG}" ;;
    f) output_folder="${OPTARG}" ;;
    t) output_tags="${OPTARG}" ;;
    o) output="${OPTARG}" ;;
    p) pvalues="${OPTARG}" ;;
    r) r2thresh="${OPTARG}" ;;
    w) window="${OPTARG}" ;;
    h) display_usage ; exit ;;
    \?|*) display_usage
       exit 1;;
  esac
done

#define usage

#test validity of inputs and run location

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

#get genotypes_stem (sans path) using shell magic
genotypes_stem=${genotypes##*/}

#split summary stats and output tags into arrays for easy parsing
IFS=',' read -r -a sumstats_array <<< "$sumstats"
IFS=',' read -r -a output_tag_array <<< "$output_tags"

############################# Start actually doing things ############################

for i in ${!sumstats_array[@]};
do
sumstats_current=${sumstats_array[i]}
output_tag_current=${output_tag_array[i]}

    printf "\nStep 1: Determining overlapping variants between genotypes and $output_tag_current summary stats

current summary stats: $sumstats_current
current output tag: $output_tag_current\n"

    #get overlapping, non-palindromic variants
    Rscript Determine_overlapping_SNPs.R $sumstats_current $genotypes ${output_folder}/$output_tag_current

    #Extract overlapping SNPs from the genotype data
    plink --bfile $genotypes --allow-no-sex --extract ${output_folder}/${output_tag_current}_overlapping_SNPs.txt --make-bed --out ${output_folder}/${genotypes_stem}_${output_tag_current} > /dev/null 2>&1 
    genotypes_new=${output_folder}/${genotypes_stem}_${output_tag_current}

    printf "Step 2: Performing LD clumping for Soutput_tag_current \n"
    #Perform LD clumping
    plink --bfile ${genotypes_new} --allow-no-sex --clump ${output_folder}/${output_tag_current}_summary_stats_updated.txt --clump-p1 1 --clump-r2 $r2thresh --clump-kb $window --out ${genotypes_new}  > /dev/null 2>&1
    
    #Create the input file for the score calculation
    Rscript Generate_score_input_file.R ${genotypes_new}.clumped $sumstats_current $pvalues

    #create range file, with one line for each p value threshold
    for p in $( echo $pvalues | sed 's/,/ /g' ); do echo "Pval_$p 0 $p" >> ${output_folder}/${output_tag_current}_pvalue_range.txt ; done

    printf "Step 3: Calculating PRS for ${output_tag_current}\n"
    #Calculate PRS
    plink --bfile $genotypes_new --allow-no-sex --score ${genotypes_new}_score_input.txt --q-score-range ${output_folder}/${output_tag_current}_pvalue_range.txt ${genotypes_new}_score_input.txt 1 3 --out ${genotypes_new}_PRS  > /dev/null 2>&1
    PRS_stem=${genotypes_new}_PRS
    
    #report the number of variants that were skipped, if any were skipped
    skipped_var=$( grep "lines skipped in --score" ${genotypes_new}_PRS.log | grep -o -E '[0-9]+' | head -n1 )
    if [ ! -z "$skipped_var" ]; then printf "$skipped_var variants skipped from PRS calculation. See ${genotypes_new}_PRS.log for more details.\n"; fi
    
    #check for variants that failed to be incorporated into the score
    if [ -f  ${PRS_stem}.nopred ];
    then
	printf "Some variants had mismatched alleles for the PRS. Attempting to flip the strands and rebuild the score...\n"
    
	#Pull out list of SNPs with allele code mismatches
	awk '{print $2}' ${PRS_stem}.nopred > ${PRS_stem}.flipsnps
	
	#Flip strand for list of SNPs
	plink --bfile $genotypes_new --allow-no-sex --flip ${PRS_stem}.flipsnps --make-bed --out ${genotypes_new}_flipsnps  > /dev/null 2>&1

	#Create .profile file
	plink --bfile ${genotypes_new}_flipsnps --allow-no-sex --score ${genotypes_new}_score_input.txt --q-score-range ${output_folder}/${output_tag_current}_pvalue_range.txt ${genotypes_new}_score_input.txt 1 3 --out $PRS_stem  > /dev/null 2>&1
	#report the number of variants that were skipped, if any were skipped
	skipped_var=$( grep "lines skipped in --score" ${genotypes_new}_PRS.log | grep -o -E '[0-9]+' | head -n1 )
	if [ ! -z "$skipped_var" ]; then printf "$skipped_var variants skipped from PRS calculation. See ${genotypes_new}_PRS.log for more details.\n"; fi
    fi

#Clean-up
rm ${genotypes_new}*.bed ${genotypes_new}*.bim ${genotypes_new}*.fam ${genotypes_new}.clumped

done

printf "Step 4: Combining PRS into one file\n"

#combine all PRS into one file
Rscript Compile_PRS.R $output_folder $genotypes_stem $output_tags $pvalues $output

printf "Step5: Generating correlation plot of all generated scores\n"

Rscript Create_PRS_corplot.R ${output_folder}/${output}.txt
