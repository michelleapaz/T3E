# Transposable Element Enrichment Estimator (T3E) 
A tool for characterising the epigenetic profile of transposable elements using ChIP-seq data
## Requirements
T3E was developed for UNIX environments, written and tested with the following versions:
* Software
  * Python 3.8.5
    * Libraries: 
      * sys 3.8.5
      * argparse 1.1
      * numpy 1.19.0
      * random (default version = 2)
      * os (depending on sys version)
      * math (depending on sys version)
      * time (depending on sys version)
  * Perl 5.30.0
  * R 3.6.3
  * bedtools 2.27.1
  * bedmap 2.4.37
  * samtools 1.10

## Instalation

## Usage
    nohup bash main.sh > log_file.txt 2>&1 &
Create **parameters** file and **control_sample** files

Calculate input-basd background probability distribution

    probabilities.py [-h] [--version] [--control <control_file>]
                     [--readlen <readlen>] [--species <species>]
                     [--outputfolder <outputfolder>]

| Arguments  | Explanation |
| ------------- | ------------- |
| -h, --help | Shows help message and exits |
| --version | Shows version message and exits |
| --control | ChIP-seq input control experiment [BED format] |
| --readlen | ChIP-seq input control experiment read length in base pairs [Example --readlen 36] |
| --species | hg38 (Homo sapiens) or mm10 (Mus musculus) [Example --species hg38] |
| --outputfolder | Output folder path [Example: /probabilities] |
