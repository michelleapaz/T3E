###########################################################################
# R (v3.6.3) script filter_highcov.R (v1.0) to filter artifact regions out
# Last update: 2022_06_07
# Author: Michelle Almeida da Paz
###########################################################################

#!/usr/bin/R

args <- commandArgs ( TRUE );
if ( length ( args ) != 4 ) {
	stop ( "[Error in R] Check the number of arguments!" );
}
percentile <- as.numeric(args[1]);
path <- args[2];
input <- args[3];
species <- args[4];

if ( ! file.exists(input) ) {
	stop( paste( input, "input does not exist!") );
}

library(dplyr)
library(ggplot2)

data <- read.delim(input, header = FALSE, sep = "\t")
data <- data %>% filter(V4 > 0)

if ( species == "hg19" || species == "hg38" ) {
	chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY");
} else if ( species == "mm10" ) {
	chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY");
} else { stop( paste( species, "species not valid!") ) };


threshold_all_chrs <- quantile(data$V4, probs = c(percentile), names = FALSE);

print(threshold_all_chrs);

for (chr in chrs) {
	file_name <- paste(chr, "window_counts.bed", sep= "_");
	file_path <- paste(path, file_name, sep="/");
	chr_windows <- filter(data, V1 == chr);
	threshold_this_chr <- quantile(chr_windows$V4, probs = c(percentile), names = FALSE);
	chr_windows <- chr_windows %>% filter(., V4 >= threshold_all_chrs | V4 >= threshold_this_chr);
	write.table(chr_windows, file=file_path, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE);
}

plot_path <- paste(path, "boxplot.png", sep= "/");
png(plot_path);
p <- ggplot(data, aes(x=V1, y=V4, group=V1)) + geom_boxplot();
p <- p + theme_bw();
p <- p + geom_hline(yintercept=threshold_all_chrs, linetype="dashed", color="red");
p <- p + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1));
p <- p + ylab("Counts in window100 and step50");
print(p);
dev.off();

