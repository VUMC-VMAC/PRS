# Polygenic Risk Score Pipeline

This repository contains the scripts to run the CNT polygenic risk score pipeline. These are contained within the CNT genomic processing container (current version: CNT_genomic_processing_v2.5.simg). 

For each summary statistics file supplied and for each p value threshold, these steps are completed:
* Subset to overlapping variants between input genotypes and summary statistics
* Perform LD clumping (default thresholds: p=0.01, r2=0.5, window = 250kb)
* Calculate the PRS using the standard averaging method (for missing data, the average frequency is used)
* Combine all PRS into one file

Before running the pipeline, the input genotypes of the samples to calculate the PRS for must be in plink format and the summary statistics must be formatted to have SNP, A1, BETA or OR, and P columns. They can have additional columns which will be ignored.

To run the pipeline, use a command similar to the following:

```
singularity exec --containall --bind /scratch/:/scratch CNT_genomic_processing_v1.simg \
	    /bin/bash -c "cd /scripts/PRS ; \
	    bash Generate_PRS.sh 
	    	 -i /scratch/mahone1/PRS/A4_NHW_imputed_final \
	    	 -s /scratch/mahone1/PRS/summary_stats/Jansenetal_summary_stats_b38.txt,/scratch/mahone1/PRS/summary_stats/Kunkle_all_summary_stats_cleaned_b38.txt \
		 -f /scratch/mahone1/PRS/A4/ \
		 -t Jansen_PRS,Kunkle_PRS \
		 -o A4_ADrisk_PRS \
		 -p 0.01 -r 0.25 -w 500 -a -b b37"
```

## Flags:
Define inputs: 
* -i specifies the input genotypes.
* -s specifies the summary statistics for calculating the score. Notice that multiple summary stats files can be supplied, separated by a comma.
* -l specifies that the input is from logistic regression and will assume the column from which to build the score will be OR rather than BETA. 

Define outputs:
* -f specifies the output folder in which results should be saved.
* -t specifies the tag or label for each of the specified summary statistics. This will be the name of the column of each score in the output file.
* -o specifies the output label for the final file.

Define score calculation thresholds:

[The following 3 options (p, r, and w) can be left unset. If they are not set, plink clumping defaults will be used (p=0.0001, r2=0.5, window=250).]
* -p specifies the p-value threshold to use for clumping. More than 1 p-value threshold can be supplied here, separated by a comma (no space).
* -r specifies the r2 threshold for clumping. 
* -w specifies the clumping window.

APOE inclusion/exclusion:
* -a specifies whether to calculate with and without APOE. 
* -b specifies the input genome build, only relevant if calculating without APOE since this just defines the base pair window for APOE exclusion. If this is not set and APOE is to be excluded, the build is assumed to be build 38.

* -m specifies a memory limit for all plink commands (in MB). This can be left unspecified which will result in no memory limit being applied. 
