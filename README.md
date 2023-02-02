# Transposable Element Enrichment Estimator (T3E) 
## A tool for characterising the epigenetic profile of transposable elements using ChIP-seq data
![T3E Image](/figures/github_figure.jpg)
### Read about T3E in:
Almeida da Paz, M., Taher, L. T3E: a tool for characterising the epigenetic profile of transposable elements using ChIP-seq data. Mobile DNA 13, 29 (2022). https://doi.org/10.1186/s13100-022-00285-z

## Software requirements
T3E was developed for UNIX environments, written and tested with the following versions:
* `Python 3.8.5`
  * `sys 3.8.5 (depending on sys version: os, math, time)`
  * `argparse 1.1`
  * `numpy 1.19.0`
  * `random (default version = 2)`
* `Perl 5.30.0`
* `R 3.6.3`
  * `dplyr 1.0.9`
  * `ggplot2 3.3.6`
* `bedtools 2.27.1`
* `bedmap 2.4.37`
* `samtools 1.10`
* `git 2.25.1`
* `conda 4.13.0`
### External tool dependencies
T3E also requires external tools to run. They should be installed using a package manager or following their website instructions. 
#### Bedtools 
[Bedtools](https://bedtools.readthedocs.io/en/latest/) is a toolset for a wide-range of genomics analysis tasks. To verify if bedtools is installed (and its version), run this command:

    bedtools --help
    
If the usage instructions are printed to the terminal, it is installed. Check its version with the command (it is not recommended to use an older version than the one tested (v2.27.1), although older versions may work well):

    bedtools --version
    
#### Bedmap
[Bedmap](https://bedops.readthedocs.io/en/latest/content/reference/statistics/bedmap.html) is a program used to retrieve and process regions of interest in BED files. To verify if bedmap is installed (and its version), run this command:

    bedmap --help
    
If the usage instructions are printed to the terminal, it is installed. Check its version with the command (it is not recommended to use an older version than the one tested (v2.4.37), although older versions may work well):

    bedmap --version
    
#### Samtools
[Samtools](http://www.htslib.org/) is a suit of programs for interacting with high-throughput sequencing data. To verify if samtools is installed (and its version), run this command:

    samtools --help
    
If the usage instructions are printed to the terminal, it is installed. Check its version with the command (it is not recommended to use an older version than the one tested (v1.10), although older versions may work well):

    samtools --version
    
#### R with required packages
T3E runs [R](https://www.r-project.org/) script to filter out regions of extremely high signals (if asked to). The [dplyr](https://dplyr.tidyverse.org/) and [ggplot2](https://ggplot2.tidyverse.org/) are required. To verify if R is installed (and its version), run this command:

    R
    
If the usage instructions are printed to the terminal, it is installed. To check if the required packages are already installed, try to load them:

    library(dplyr)
    library(ggplot2)
    
If the packages load without any error, they are already installed. Otherwise, install them in R:

    install.packages("dplyr")
    install.packages("ggplot2")
    
## Installation
### Cloning T3E repository from GitHub
Before using [Git](https://git-scm.com/), make sure it is available. To verify if Git is installed, run this command:

    git
    
If the usage instructions are printed to the terminal, it is installed. To clone T3E repository into a new directory, download a .zip file from GitHub or run the command:

    git clone https://github.com/michelleapaz/T3E
   
### Installing with conda
We recommend using [conda](https://docs.conda.io/en/latest/) environment (with defined dependencies) to run T3E. To verify if conda is installed, run this command:

    conda
    
If the usage instructions are printed to the terminal, it is installed. Dependencies are defined in the file **`environment.yml`**. Create the environment using the command: 

    conda env create --file environment.yml
    
To see a list of all environments, run the command:

    conda env list
    
The recently created environment called **`t3e-env`** should appear in the list. Then, activate this environment using the command: 

    conda activate t3e-env
    
In case any dependency is not successfully installed, it can be installed, for example, using conda:

    conda install numpy
    
Or [pip](https://pypi.org/project/pip/):

    pip install numpy
    
After running T3E, do not forget to deactivate the t3e-env environment, running the command:

    conda deactivate
    
## Usage
We provide input files as examples and the repeat annotation files for human and mouse genomes used in our study. These files are available here: https://cloud.tugraz.at/index.php/s/Ec8HfnoasnMzcPS
<br />To run T3E the content of five folders should be considered:
1. **`./bam/`** - should contain the alignments (BAM files) for ChIP-seq samples and their corresponding input control. It is important that secondary reads are reported by the used mapper. Two files are provided as example (`test_sample.bam` and `test_control.bam`)

2. **`./references/`** - contains 
 * **`control_sample.txt`** file that contains the ChIP-seq input control and sample names (if more than one ChIP-seq sample, separate them using commas) separated by tab: <br />
 ```control	sample_1,sample_2,...sample_n``` 
 <br />For example: <br />
 ```test_control	test_sample```
 
 * **`parameters.txt`** file containing the following parameters: <br />
 
 | Arguments  | Explanation |
 | ------------- | ------------- |
 | species | hg38 (_Homo sapiens_) or mm10 (_Mus musculus_) |
 | iterations | number of iterations [Example: 100] |
 | alpha | level of significance to report enrichment [Example: 0.05] |
 | enrichment | log2FC threshold to report enrichment [Example: 1.0] |
 | filter | filter out regions of extremely high signals (0 for NO and 1 for YES) |
 
 Example: <br />
 
 ```
 species	hg38
 iterations	100
 alpha	0.05
 enrichment	1.0
 filter	0
 ```
 
 Like that, T3E considers the repeat annotation for _Homo sapiens_ (hg38), simulates 100 input libraries for each ChIP-seq sample, considers level of significance of 0.05 and log2FC of 1.0 for reporting enrichments, and it does not filter out regions of extremely high signals
 
 **It is recommended to filter out regions of extremely high signals for the mouse genome (filter = 1)!**

 * chromosome size files [**`hg38.genome`** (_Homo sapiens_) and **`mm10.genome`** (_Mus musculus_)]
 <br />Content of `hg38.genome` file (first 5 lines):
 
 ```
 chrom	size
 chr1	248956422
 chr2	242193529
 chr3	198295559
 chr4	190214555
 ```

 * **`path_dataset.csv`** file that is created by T3E and contains the path, library size and read length for each BAM file, separated by semicolon. Example:<br />
 
 ```
 ./T3E/bam/test_control.bam;697382;76
 ./T3E/bam/test_sample.bam;186923;76
 ```

3. **`./repeats/`** - contains transposable elements annotation [**`rmsk_hg38.bed`** (_Homo sapiens_) and **`rmsk_mm10.bed`** (_Mus musculus_)]. Custom repeat annotations for other species may be added and should follow the same format file (BED file) with information of TE individual copies. The order of the three first columns should be respected and must contain the chromosome, start coordinate on the chromosome and end coordinate on the chromosome. The fourth column should contain the information about the repeat (e.g. TE family/subfamily)
<br />Content of `rmsk_hg38.bed` file (first 5 lines):

```
chr1	11504	11675	L1MC5a
chr1	11677	11780	MER5B
chr1	15264	15355	MIR3
chr1	18906	19048	L2a
chr1	19971	20405	L3
```
4. **`./results/`** - contains the output files (one folder for each control and sample BAM files)
5. **`./scripts/`** - contains all Python, Perl and R scripts

The **`main.sh`** code uses the information contained in two files (`parameters.txt` and `control_sample.txt`) in `./references` folder, processes the datasets and runs T3E scripts **automatically**

### Run **main.sh** code:

    nohup bash main.sh > log_file.txt 2>&1 &

This command is all you have to run since you managed to configure everything until now. You can check the status of the run printing the end of the log file with the command: 

    tail log_file.txt

T3E also creates a log file (e.g. `log_test_sample.txt`) which can be checked in the same manner. It is also possible to run each script separately (**but it not necessary! T3E does it for you!**). The scripts consist in three steps:

#### Calculate input-based background probability distribution:
It estimates the probability of a read starting at an effective genomic position in the ChIP-seq input control experiment

    probabilities.py [-h] [--version] [--control <control_file>]
                     [--readlen <readlen>] [--species <species>]
                     [--outputfolder <outputfolder>]

| Arguments  | Explanation |
| ------------- | ------------- |
| -h, --help | shows help message and exits |
| --version | shows version message and exits |
| --control | ChIP-seq input control experiment [BED format] |
| --readlen | ChIP-seq input control experiment read length in base pairs [Example --readlen 76] |
| --species | hg38 (_Homo sapiens_) or mm10 (_Mus musculus_) [Example --species hg38] |
| --outputfolder | output folder path [Example: /probabilities] |

<br />Example of input BED file (`test_control.bed`) for `--control` parameter (first 5 lines):

```
chr1	10004	10080	NS500343:103:H72MMBGXY:1:21109:21707:18958
chr1	10016	10092	NS500343:103:H72MMBGXY:1:21109:21707:18958
chr1	10022	10098	NS500343:103:H72MMBGXY:1:21109:21707:18958
chr1	10028	10104	NS500343:103:H72MMBGXY:1:21109:21707:18958
chr1	10034	10110	NS500343:103:H72MMBGXY:1:21109:21707:18958
```
The file contains the chromosome (column 1), start (column 2) and end (column 3) coordinates on the chromosome and the read ID (column 4) for each read in the ChIP-seq input control library (note that all loci should be reported for multimappers)

<br />Example of command:

    python3 ./T3E/scripts/probabilities.py --control ./T3E/results/test_control/test_control.bed --readlen 76 --species hg38 --outputfolder ./T3E/results/test_control/probabilities
    
The output files are created in the specified folder. In total, each chromosome has one `.txt` output file containing the genomic position (column 1) and the corresponding cumulative probability (column 2). In the example above, 24 `.txt` files were generated. Example of `chr1_prob.txt` output file (first 5 lines):

```
10004 	 5.821325789847425e-10
10005 	 1.164265157969485e-09
10006 	 1.7463977369542275e-09
10007 	 2.32853031593897e-09
10008 	 2.910662894923712e-09
```

#### Construct the background probability distribution of read mappings:
It is the core script of T3E and it computes the background distribution of read mappings by randomly sampling read mappings based on the structure of the input control library

    t3e.py [-h] [--version] [--repeat <repeat_file>] [--sample <sample_file>] 
           [--readlen <readlen>] [--control <control_file>]
           [--controlcounts <control_counts>] [--probability <probability_folder>] 
           [--iter <iter>] [--species <species>]
           [--outputfolder <outputfolder>] [--outputprefix <outputprefix>]

| Arguments  | Explanation |
| ------------- | ------------- |
| -h, --help | shows help message and exits |
| --version | shows version message and exits |
| --repeat | transposable elements annotation [rmsk_hg38.bed (_Homo sapiens_) or rmsk_mm10.bed (_Mus musculus_)] |
| --sample | ChIP-seq sample experiment [BED format] |
| --readlen | ChIP-seq input control experiment read length in base pairs [Example --readlen 76] |
| --control | ChIP-seq input control experiment [BED format] |
| --controlcounts | ChIP-seq input control experiment counts [.txt format] |
| --probability | probability folder path [Example: /control/probability/] |
| --iter | number of iterations [Example: 100] |
| --species | hg38 (_Homo sapiens_) or mm10 (_Mus musculus_) [Example --species hg38] |
| --outputfolder | output folder path [Example: /results] |
| --outputprefix | prefix name of your analysis [Example: test_sample] |

<br />Example of `.txt` input file (`test_control_counts.txt`) for `--controlcounts` parameter (first 5 lines):

```
Alu	37.9599074240402
AluJb	8086.88216090022
AluJo	4291.85264145515
AluJr	4867.75105757094
AluJr4	1122.24456393148
```
The file contains the TE family/subfamily (column 1) and the corresponding read mapping counts (column 2) for the ChIP-seq input control

<br />Example of command:

    python3 ./T3E/scripts/t3e.py --repeat ./T3E/repeats/rmsk_hg38.bed --sample ./T3E/results/test_sample/test_sample.bed --readlen 76 --control ./T3E/results/test_control/test_control.bed --controlcounts ./T3E/results/test_control/test_control_counts.txt --probability ./T3E/results/test_control/probabilities --iter 100 --species hg38 --outputfolder ./T3E/results/test_sample/ --outputprefix test_sample > ./T3E/log_test_sample.txt
    
The output files are created in the specified folder. In the example, the background file is important for the next script for computing TE families/subfamilies enrichments. Example of `test_sample_background.txt` output file (first 5 lines):

```
iter1	Alu	8.782564266947237
iter1	AluJb	2219.1462339576756
iter1	AluJo	1164.8765166651572
iter1	AluJr	1304.9989875595957
iter1	AluJr4	313.5173699906347
```
The file contains the number of iteration (column 1), TE family/subfamily (column 2) and the corresponding read mapping counts for the simulated background (column 3)

#### Calculate TE families/subfamilies enrichments:
It computes ChIP-seq enrichment at TE families/subfamilies relative to a background

    enrichment.py [-h] [--version] [--background <background>]
                  [--signal <signal>] [--iter <iter>] [--alpha <alpha>]
                  [--enrichment <enrichment>] [--outputfolder <outputfolder>]
                  [--outputprefix <outputprefix>]
                  
| Arguments  | Explanation |
| ------------- | ------------- |
| -h, --help | shows help message and exits |
| --version | shows version message and exits |
| --background | background file created by T3E [Example: sample001_background.txt] |
| --signal | ChIP-seq sample experiment counts [.txt format] |
| --iter | number of iterations [Example: 100] |
| --alpha | level of significance to report enrichment [Example: 0.05] |
| --enrichment | log2FC threshold to report enrichment [Example: 1.0] |
| --outputfolder | output folder path [Example: /results] |
| --outputprefix | prefix name of your analysis [Example: test_sample] |

<br />Example of `.txt` input file (`test_sample_counts.txt`) for `--signal` parameter (first 5 lines):

```
Alu	14.2787735236149
AluJb	2088.6238435312
AluJo	1200.62110262026
AluJr	1333.68223686911
AluJr4	299.617299043992
```
The file contains the TE family/subfamily (column 1) and the corresponding read mapping counts for the ChIP-seq sample experiments (column 2)

<br />Example of command:

    python3 ./T3E/scripts/enrichment.py --background ./T3E/results/test_sample/test_sample_background.txt --signal ./T3E/results/test_sample/test_sample_counts.txt --iter 100 --alpha 0.05 --enrichment 1.0 --outputfolder ./T3E/results/test_sample/ --outputprefix test_sample
    
The output file contains all the enrichments and it is created in the specified folder. Example of `test_sample_enrichment.txt` output file (first 5 lines):

```
Alu	0.02	0.39804917554376207
AluJb	1.0	-0.07102969494637044
AluJo	0.09	0.046781617267046424
AluJr	0.34	0.014120002340119377
AluJr4	0.59	-0.015369887385721339
```

Note that the file contains the TE family/subfamily (column 1) and its corresponding P-value (column 2) and log2FC (column 3). But also, T3E prints the enriched TE families/subfamilies considering the chosen level of significance and log2FC thresholds:

```
AluYh7   0.04	1.2034230038490004
Charlie4 0.02	1.6156846715708846
DNA1_Mam 0.04	2.400492454970415
Eulor1   0.01	3.120769507990133
Eulor12  0.05	2.651369092553242
```
