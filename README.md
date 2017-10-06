# mutation_rainstorm
We implemented a multi-patient variant of the rainfall plot we have named a rainstorm plot. There is a fairly complicated set of pre-processing steps to help mitigate a variety of issues that are encountered when trying to plot more than one patient rainfall plot on one image. Our calculation also adjusts for local variation in mutation rate that can be observed at the level of the entire cohort.

## rainstorm.R 
This script takes a MAF file containing genome-wide mutations from many cancer genomes and determines the background mutation rate using the whole cohort. Then, one chromosome at a time, a patient-by-patient calculation similar to the rainfall plot calculation is used to infer the distance between each mutation and mutations in other genomes in the same cohort. 

## rainstorm_plot.R 
A basic plotting script using ggplot2 to visualize chromosome-wide patterns or specific regions on chromosomes.

## rainstorm_peaks.R

### Example
```Rscript ./rainstorm_peaks.R --stringSplit mean_ --output_base fl_dlbcl_wavelet --input_maf ./DLBCL_FL.maf /home/rmorin/git/mutation_rainstorm/fl_dlbcl/FL_DLBCL_*mean*tsv```
