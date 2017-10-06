# mutation_rainstorm
We implemented a multi-patient variant of the rainfall plot we have named a rainstorm plot. There is a fairly complicated set of pre-processing steps to help mitigate a variety of issues that are encountered when trying to plot more than one patient rainfall plot on one image. Our calculation also adjusts for local variation in mutation rate that can be observed at the level of the entire cohort.

## rainstorm.R 
This script takes a MAF file containing genome-wide mutations from many cancer genomes and determines the background mutation rate using the whole cohort. Then, one chromosome at a time, a patient-by-patient calculation similar to the rainfall plot calculation is used to infer the distance between each mutation and mutations in other genomes in the same cohort. 

```
usage: ./rainstorm.R [-h] [--input_maf INPUT_MAF]
                     [--output_base_name OUTPUT_BASE_NAME] [--cpu_num CPU_NUM]
                     [--genome_fai GENOME_FAI] [--plot PLOT]
                     [--max_mut MAX_MUT] [--off_by OFF_BY]
                     [--calc_background CALC_BACKGROUND]

Calculate rainstorm intermutation distance values for all mutations in a large
set of cancer genomes

optional arguments:
  -h, --help            show this help message and exit
  --input_maf INPUT_MAF, --m INPUT_MAF
                        MAF file containing mutation calls from many patient
                        genomes
  --output_base_name OUTPUT_BASE_NAME, --o OUTPUT_BASE_NAME
                        specify a base file name prefix for all outputs
  --cpu_num CPU_NUM, -c CPU_NUM
                        set to number of CPUs you would like to use to perform
                        calculation in parallel (consumes lots of RAM)
  --genome_fai GENOME_FAI, --g GENOME_FAI
                        provide the corresponding fasta index for the genome
                        you used. must match the chromosome naming style used
                        in your MAF!
  --plot PLOT, -p PLOT  ploduce rainstorm plot for each chromosome
  --max_mut MAX_MUT, -M MAX_MUT
                        genomes skipped if their total mutation load exceeds
                        this value
  --off_by OFF_BY, -k OFF_BY
                        take mean of the distance to the k closest mutations
                        to determine rainstorm distance value
  --calc_background CALC_BACKGROUND, -b CALC_BACKGROUND
                        if you have done this once for a cohort, you can
                        reload the result in future runs by setting this to 0

```
### Example
```
Rscript ./rainstorm.R --m ./cohort_mutations_merged.maf --o cohort_out --genome_fai ./hg19.ucsc.fa.fai --cpu_num 12
```

## rainstorm_plot.R 
A basic plotting script using ggplot2 to visualize chromosome-wide patterns or specific regions on chromosomes.

## rainstorm_peaks.R

### Example
```Rscript ./rainstorm_peaks.R --stringSplit mean_ --output_base fl_dlbcl_wavelet --input_maf ./DLBCL_FL.maf /home/rmorin/git/mutation_rainstorm/fl_dlbcl/FL_DLBCL_*mean*tsv```
