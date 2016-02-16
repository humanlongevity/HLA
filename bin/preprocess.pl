#!/usr/bin/env perl
# Author: Xie Chao
use strict;

my $Q = 20;
my $F = 5;
my $M = 100;

while(<>)
{
	my @fields = split(/\t/, $_);
	my $name = $fields[0];
	my $flag = $fields[1];
	my $seq = $fields[9];
	my $qual = $fields[10];

	# reverse complement if the read was flipped in the BAM
	if($flag & 16){
		$qual = reverse $qual;
		$seq = reverse $seq;
		$seq =~ tr/ATGC/TACG/;
	}
	my @qual = map { ord($_) - 33 } split('', $qual);

	# remove the N-prefix
	if($seq =~ m/^(N+)/){
		my $n = length $1;
		$seq = substr($seq, $n);
		@qual = @qual[$n..$#qual];
	}

	my @qual2 = @qual;

	# 3' trimming by quality
	my @trim;
	my $i = $#qual + 1;
	my $cum = 0;
	my $min = 100;
	my $which_min = $#qual + 100;
	while(my $q = pop @qual)
	{
		$i--;
		$cum += $q - $Q;
		unshift @trim, $cum;
		if($cum < $min)
		{
			$min = $cum;
			$which_min = $i;
		}
	}
	if($qual2[$which_min] < $Q){
		$seq = substr($seq, 0, $which_min);
		$qual = substr($qual, 0, $which_min);
		@qual = @qual2[0..($which_min-1)]
	}

	# min length after trimming
	next if length($seq) < $M;

	# max bad bases
	my $fail_count;
	for my $q(@qual){
		$fail_count++ if $q <= 3;
	}
	next if $fail_count > $F;

	# discard reads if contains N
	next unless $seq =~ m/^[ATGC]+$/;

	# name reads with /1 or /2
	if($flag & 128){
		$name .= '/2';
	}else{
		$name .= '/1';
	}

	print "\@$name\n$seq\n+\n$qual\n";
}

