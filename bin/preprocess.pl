#!/usr/bin/env perl
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
	next unless $seq =~ m/^[ATGC]+$/;
	my $qual = $fields[10];

	if($flag & 16){
		$qual = reverse $qual;
		$seq = reverse $seq;
		$seq =~ tr/ATGC/TACG/;
	}
	my @qual = map { ord($_) - 33 } split('', $qual);

	if($qual[$#qual] < $Q){
		my @qual2 = @qual;

		my @trim;
		my $i = $#qual + 1;
		my $cum = 0;
		my $min = 100;
		my $which_min;
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
		if($which_min){
			$seq = substr($seq, 0, $which_min);
			$qual = substr($qual, 0, $which_min);
			@qual = @qual2[0..($which_min-1)]
		}else{
			$seq = '';
		}
	}
	next if length($seq) < $M;

	my $fail_count;
	for my $q(@qual){
		$fail_count++ if $q <= 3;
	}
	next if $fail_count > $F;

	if($flag & 128){
		$name .= '/2';
	}else{
		$name .= '/1';
	}

	print "\@$name\n$seq\n+\n$qual\n";
}

