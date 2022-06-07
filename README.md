# Transposable Element Enrichment Estimator (T3E) 
### A tool for characterising the epigenetic profile of transposable elements using ChIP-seq data
![T3E Image](/figures/github_figure.jpg)

## Requirements
T3E was developed for UNIX environments, written and tested with the following versions:
* Software
  * Python 3.8.5
    * Libraries: 
      * sys 3.8.5 (depending on sys version: os, math, time)
      * argparse 1.1
      * numpy 1.19.0
      * random (default version = 2)
  * Perl 5.30.0
  * R 3.6.3
  * bedtools 2.27.1
  * bedmap 2.4.37
  * samtools 1.10

## Instalation
Clone the T3E repository:

    git clone https://github.com/michelleapaz/T3E
    
## Usage
T3E contains the **main.sh** code and five folders:
1. **./bam/** - contains BAM files
2. **./references/** - contains 
 * **control_sample.txt**
 * **parameters.txt** 
 * chromosome size files [hg38.genome (_Homo sapiens_) and mm10.genome (_Mus musculus_)]
 * <em>path_dataset.csv file is created</em>
3. **./repeats/** - contains transposable elements annotation [rmsk_hg38.bed (_Homo sapiens_) and rmsk_mm10.bed (_Mus musculus_)]
4. **./results/** - contains the output files (one folder for each control and sample BAM files)
5. **./scripts/** - contains all Python, Perl and R scripts

The **main.sh** code uses the information contained in two files (**parameters.txt** and **control_sample.txt**) in **./references** folder
<br />
The **parameter.txt** file contains the following parameters:

| Arguments  | Explanation |
| ------------- | ------------- |
| species | hg38 (_Homo sapiens_) or mm10 (_Mus musculus_) |
| iterations | number of interations [Example: 100] |
| alpha | level of significance to report enrichment [Example: 0.05] |
| enrichment | log2FC threshold to report enrichment [Example: 1.0] |
| filter | filter out regions of extremely high signals (0 for NO and 1 for YES) |

The **control_sample.txt** file contains the ChIP-seq input control and sample names (if more than one ChIP-seq sample, separate them using commas) separated by tab:

```control sample_1,sample_2,...sample_n```

Directly, **main.sh** code processes your datasets and runs T3E

### Run **main.sh** code:

    nohup bash main.sh > log_file.txt 2>&1 &

### Or run each script separately:

#### Calculate input-based background probability distribution:

    probabilities.py [-h] [--version] [--control <control_file>]
                     [--readlen <readlen>] [--species <species>]
                     [--outputfolder <outputfolder>]

| Arguments  | Explanation |
| ------------- | ------------- |
| -h, --help | shows help message and exits |
| --version | shows version message and exits |
| --control | ChIP-seq input control experiment [BED format] |
| --readlen | ChIP-seq input control experiment read length in base pairs [Example --readlen 36] |
| --species | hg38 (_Homo sapiens_) or mm10 (_Mus musculus_) [Example --species hg38] |
| --outputfolder | output folder path [Example: /probabilities] |

#### T3E:

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
| --sample | ChIP-seq sample experiment [.dict format] |
| --readlen | ChIP-seq input control experiment read length in base pairs [Example --readlen 36] |
| --control | ChIP-seq input control experiment [BED format] |
| --controlcounts | ChIP-seq input control experiment counts [.txt format] |
| --probability | probability folder path [Example: /control/probability/] |
| --iter | number of interations [Example: 100] |
| --species | hg38 (_Homo sapiens_) or mm10 (_Mus musculus_) [Example --species hg38] |
| --outputfolder | output folder path [Example: /probabilities] |
| --outputprefix | prefix name of your analysis [Example: sample001] |

#### Calculate TE families/subfamilies enrichments:

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
| --iter | number of interations [Example: 100] |
| --alpha | level of significance to report enrichment [Example: 0.05] |
| --enrichment | log2FC threshold to report enrichment [Example: 1.0] |
| --outputfolder | output folder path [Example: /results] |
| --outputprefix | prefix name of your analysis [Example: sample001] |
