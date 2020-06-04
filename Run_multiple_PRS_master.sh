#!/bin/bash
#Generates Polygenic Risk Scores
#Derek Archer, August 13, 2019
#Modified from Logan's July 2019 script

top_dir=`pwd`


#Set variables
CLUMP1=.01    #Significance threshold for index SNPs
CLUMP2=.5     #Secondary significance threshold for clumped SNP
CLUMP3=.25    #LD threshold for clumping
CLUMP4=200    #Physical distance threshold for clumping

while read GENETICS SUMSTATS OVERLAPSCRIPT OUTPUTNAME APOEINCLUSION ; do 
    
    #Run R script which 1) removes ambiguous SNPs and 2) generates list of overlapping SNPs
    
    
    if [ $APOEINCLUSION == 1 ] ; then 
	
	#Extract overlapping SNPs from the genotype data
	Rscript /nfs/clarklc/PRS/Scripts/$OVERLAPSCRIPT $SUMSTATS $GENETICS $OUTPUTNAME
        echo "APOEINCLUSION is 1: including APOE in analysis"
        plink --bfile $GENETICS --allow-no-sex --extract ${OUTPUTNAME}_overlapping_SNPs.txt --make-bed --out $OUTPUTNAME

	

    else

	echo "APOEINCLUSION is 0: adding output \"_noAPOE tag\" and removing APOE from analysis"
	
	OUTPUTNAME=${OUTPUTNAME}_noAPOE
	Rscript /nfs/clarklc/PRS/Scripts/$OVERLAPSCRIPT $SUMSTATS $GENETICS $OUTPUTNAME
	plink --bfile $GENETICS --allow-no-sex --exclude range /nfs/clarklc/PRS/Scripts/APOErange.txt --extract ${OUTPUTNAME}_overlapping_SNPs.txt --make-bed --out $OUTPUTNAME

    fi 
    
      

echo $GENETICS
echo $SUMSTATS
echo $OUTPUTNAME

        
	#Perform LD clumping
	#NOTE: If just rerunning analyses using different clumping criteria, start script here and comment out previous steps

        plink --bfile $OUTPUTNAME --allow-no-sex --clump ${OUTPUTNAME}_summary_stats_updated.txt --clump-p1 $CLUMP1 --clump-p2 $CLUMP2 --clump-r2 $CLUMP3 --clump-kb $CLUMP4 --out ${OUTPUTNAME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

	#Pull list of SNPs from clumped file

        awk '{print $3}' ${OUTPUTNAME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.clumped > ${OUTPUTNAME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.SNPlist

	#Merge with updated summary statistics file and generate PRS input file

        Rscript /nfs/clarklc/PRS/Scripts/Generate_score_input_file.R $OUTPUTNAME $CLUMP1 $CLUMP2 $CLUMP3 $CLUMP4

	#Create .profile file

        plink --bfile $OUTPUTNAME --allow-no-sex --score ${OUTPUTNAME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTPUTNAME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

	if [ -f  ${OUTPUTNAME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred ];
	
	then

	    printf "\n\nSome variants had mismatched alleles for the PRS. Attempting to flip the strands and rebuild the score.\n"

	    #Pull out list of SNPs with allele code mismatches

            awk '{print $2}' ${OUTPUTNAME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred > ${OUTPUTNAME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps

	    #Flip strand for list of SNPs

            plink --bfile $OUTPUTNAME --allow-no-sex --flip ${OUTPUTNAME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps --make-bed --out ${OUTPUTNAME}_flipsnps

	    #Create .profile file

            plink --bfile ${OUTPUTNAME}_flipsnps --allow-no-sex --score ${OUTPUTNAME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTPUTNAME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}_flipsnps
       fi

done < $top_dir/PRS_driver.txt

###
#Clean up

#rm *.bim
#rm *.bed
#rm *.fam
