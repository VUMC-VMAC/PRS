#!/bin/bash

#fail on error
set -e

#assumptions at the beginning of this script:
# the summary stats are already in shape to be able to figure out overlap
# will still need to determine the overlapping variants

############################# Verify inputs ############################

#parse arguments

#test validity of inputs and run location

#print out inputs to the log
printf "Input genotypes: $genotypes
Summary statistics: $sumstats
Output stem: $output
P-value threshold: $pval
R-squared threshold: $r2thresh
Window size: $window
"

############################# Start actually doing things ############################

#Extract overlapping SNPs from the genotype data
Rscript Determing_overlapping_variants.sh $sumstats $genetics $output
plink --bfile $genetics --allow-no-sex --extract ${output}_overlapping_SNPs.txt --make-bed --out $output

#Perform LD clumping
plink --bfile $output --allow-no-sex --clump ${output}_summary_stats_updated.txt --clump-p1 1 --clump-r2 $r2thresh --clump-kb $window --out $output

#Create the input file for the score calculation
Rscript Generate_score_input_file.R $output $sumstats

#Calculate PRS
plink --bfile $output --allow-no-sex --score ${output}_score_input.txt --out ${output}_PRS
PRS_stem=${output}_PRS

#check for variants that failed to be incorporated into the score
if [ -f  ${PRS_stem}.nopred ];
then
    printf "\n\nSome variants had mismatched alleles for the PRS. Attempting to flip the strands and rebuild the score.\n"
    
    #Pull out list of SNPs with allele code mismatches
    awk '{print $2}' ${PRS_stem}.nopred > ${PRS_stem}.flipsnps

    #Flip strand for list of SNPs
    plink --bfile $output --allow-no-sex --flip ${PRS_stem}.flipsnps --make-bed --out ${output}_flipsnps

    #Create .profile file
    plink --bfile ${output}_flipsnps --allow-no-sex --score ${PRS_stem}.score --out ${output}_flipsnps_PRS

    PRS_stem=${output}_flipsnps_PRS
fi

