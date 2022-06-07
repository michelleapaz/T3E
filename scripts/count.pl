###########################################################################
# Perl (v5.30.0) script count.pl (v1.0) to count for TE families/subfamilies
# Last update: 2022_06_07
# Author: Michelle Almeida da Paz
###########################################################################

#!/usr/bin/perl
use strict;
use warnings;

unless ($#ARGV == 2){
	die "usage: perl $0 read_length bam intersect\n";
}

my $read_length = $ARGV[0];
my $bam = $ARGV[1];
my $intersect = $ARGV[2];
my $key;
my $sum_repeat;
my %frequency;
my $count_records = counting_records($bam);
my $count_repeat = counting_repeat($intersect);

foreach $key (sort keys %{$count_repeat})
{
	$sum_repeat = $$count_repeat{$key};
	print "$key\t$sum_repeat\n";
}

sub counting_records {
	my $file = $_[0];
	my @fields;
	my $read_id;
	
	open ( BAM, $file ) || die "Cannot open $file!\n";
	while ( <BAM> ){
		if (/^#/)
		{
			next;
		}
		chomp;
		@fields = split (/\t/);
		$read_id = $fields[3];
		$frequency{$read_id}++;
	}
	close ( BAM );
}

sub counting_repeat {
	my $file = $_[0];
	my @fields;
	my $repeats;
	my @repeat;
	my $reads;
	my @read;
	my $maps;
	my @map;
	my %repeat_counts;
	my %read_counts;
	my %count_reads_in_group;
	my $size_array;
	my @read_col;
	my ($perc_map_total, $perc_map);
	my ($repeat_id, $read_id);
	my $i;

	open ( READ, $file ) || die "Cannot open $file!\n";
	while ( <READ> ){
		if (/^#/)
		{
			next;
		}
		@fields = split (/\|/);
		$repeats = $fields[0];
		$reads = $fields[1];
		$maps = $fields[2];
		@repeat = split (/\t/, $repeats);
		$repeat_id = $repeat[3];
		
		if ($reads eq "")
		{
			next;
		}
		else
		{
			@read = split (/\;/, $reads);
			@map = split (/\;/, $maps);
			$size_array = scalar @read;
			for ($i=0; $i<=$#read; $i++)
			{
				@read_col = split (/\t/, $read[$i]);
				$read_id = $read_col[3];
				$perc_map = $map[$i]; 
				$read_counts{$read_id} += $perc_map;
				$repeat_counts{$repeat_id}{$read_id} += $perc_map;
			}
		}
	}
	close ( READ );

	foreach $repeat_id( sort keys %repeat_counts  ){
		foreach $read_id( sort keys %{$repeat_counts{$repeat_id}}  ){
			unless ( defined($read_counts{$read_id}) ){
				die "Error. No $read_id!\n";
			}
			$count_reads_in_group{$repeat_id} += $repeat_counts{$repeat_id}{$read_id}/($frequency{$read_id}*$read_length);
		}
	}

	return(\%count_reads_in_group);
}
