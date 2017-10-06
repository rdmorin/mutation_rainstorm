# mutation_rainstorm
We implemented a multi-patient variant of the rainfall plot we have named a rainstorm plot. There is a fairly complicated set of pre-processing steps to help mitigate a variety of issues that are encountered when trying to plot more than one patient rainfall plot on one image. Our calculation also adjusts for local variation in mutation rate that can be observed at the level of the entire cohort.

## Dependencies
The following R packages are needed by these tools:
```argparse```, ```GenomicRanges```, ```ggplot2```, ```maftools```, ```MassSpecWavelet```, ```parallel```

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
With genomes that have high mutation loads, we strongly recommend running this with as many CPU cores as possible.

## rainstorm_plot.R 
A basic plotting script using ggplot2 to visualize chromosome-wide patterns or specific regions on chromosomes.

```
usage: ./rainstorm_plot.R [-h] [--rainstorm_points RAINSTORM_POINTS]
                          [--xstart XSTART] [--xend XEND]

plot a genomic region with mutations and encode data

optional arguments:
  -h, --help            show this help message and exit
  --rainstorm_points RAINSTORM_POINTS, -r RAINSTORM_POINTS
                        Rainstorm output (points to plot) from a single
                        chromosome
  --xstart XSTART, -s XSTART
                        limit plot to points after this position
  --xend XEND, -e XEND  limit plot to points before this position
```
### Example
```
Rscript ./rainstorm_plot.R --rainstorm_points ./cohort_out_chr3.tsv
```
## rainstorm_peaks.R
Use the wavelet approach to identify regions enriched for mutations across multiple patients ("peaks").

```
usage: ./rainstorm_peaks.R [-h] [--stringSplit STRINGSPLIT]
                           [--input_maf INPUT_MAF]
                           [--output_base_file OUTPUT_BASE_FILE]
                           input_files [input_files ...]

wavelet searching argument

positional arguments:
  input_files           Input path and files

optional arguments:
  -h, --help            show this help message and exit
  --stringSplit STRINGSPLIT
                        characters before chr# or # in the input file names
  --input_maf INPUT_MAF
                        Input maf path and name
  --output_base_file OUTPUT_BASE_FILE
                        Output path and name base
```

### Example
```
Rscript ./rainstorm_peaks.R --stringSplit mean_ --output_base fl_dlbcl_wavelet \
 --input_maf ./DLBCL_FL.maf \
 /home/rmorin/git/mutation_rainstorm/fl_dlbcl/FL_DLBCL_*mean*tsv
```
