###########################################################################
# Perl (v5.30.0) script subsample.pl (v1.0) to Downsample ChIP-seq libraries
# Last update: 2022_06_07
# Author: Michelle Almeida da Paz
###########################################################################

#!/usr/bin/perl
use strict;
use warnings;

unless ($#ARGV == 1){
	die "usage: perl $0 bed nr_reads_to_subsample\n";
}

my $sam = $ARGV[0];
my $nr_reads_to_subsample = $ARGV[1];

my @result = subsampling($sam, $nr_reads_to_subsample);

sub subsampling {
	my $file = $_[0];
	my $n = $_[1];
	
	my $read_id;
	my $i;
	my @random_readids;
	my %read_freq;
	my $mapping;

	open ( READ, $file ) || die "Cannot open $file!\n";
	while ( <READ> ){
		chomp;
		if (/^#/)
		{
			next;
		}
		elsif ( /^([^\t]+\t[^\t]+\t[^\t]+)\t([^\t]+)/ ){
			$mapping = $1;
			$read_id = $2;
			$read_freq{$read_id}++;
		}
		else {
			die "$_: wrong format!\n";
		}
	}

	my @hash_keys = keys %read_freq;
	fisher_yates_shuffle(\@hash_keys);
	if ( $n > $#hash_keys ){
		$n = $#hash_keys + 1;
	}
	@random_readids = @hash_keys[0..($n-1)];
	%read_freq = ();
	foreach $read_id( @random_readids ){
		$read_freq{$read_id}++;
	}

	seek (READ, 0, 0);
	while ( <READ> ){
		chomp;
		if (/^#/)
		{
			next;
			
		}
		elsif ( /^([^\t]+\t[^\t]+\t[^\t]+)\t([^\t]+)/ ){
			$mapping = $1;
			$read_id = $2;
			if ( defined($read_freq{$read_id}) ){
				print "$mapping\t$read_id\n";
			}
		}
		else {
			die "$_: wrong format!\n";
		}
	}
	close ( READ );	

	return;
}

sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
}
