#!/usr/bin/env perl
use strict;

my $ref;
my %align;
my %raw;
my %ids;
my $i = 0;
while(<>)
{
	$i++;
	next if m/^\S/;
	if(m/^ Prot/)
	{
		$i = 0;
		next;
	}
	if(m/^ (\S+?)\s+(.+)/)
	{
		my $id = $1;
		my $seq = $2;
		$raw{$id} .= $seq;
#		print $seq, "\n";
		if($i == 2)
		{
			$ref = $seq;
		}else
		{
			$seq =~ s/(-+)/substr($ref, $-[1], $+[1]-$-[1])/eg;
#			print $seq, "\n";
		}
		$align{$id} .= $seq;
		my $type = $id;
#		$type =~ s/:(\d+).*$/:\1/ unless $id =~ m/\D$/;
		$ids{$type}->{$id} = 1;
	}
}

my %size;
for my $id(keys %align)
{
	$align{$id} =~ s/\s*//g;
	my $size = 0 + $align{$id}=~s/([A-Z])/\1/g;
	$size{$id} = $size;
}

for my $type(sort keys %ids)
{
	my @ids = sort{$size{$b} <=> $size{$a}} sort keys %{$ids{$type}};
	my $longest = shift @ids;
	my $ref = $align{$longest};
	my $from = $-[1] if $ref =~ m/([A-Z])/;
	my $to = $-[1] if $ref =~ m/([A-Z])[^A-Z]*$/;
	for my $id(@ids)
	{
		my $seq = $align{$id};
		$seq =~ s/X//;
		my $f = $-[1] if $seq =~ m/([A-Z])/;
		my $t = $-[1] if $seq =~ m/([A-Z])[^A-Z]*$/;
		$f = $from if $f < $from;
		$t = $to if $t > $to;
		print STDERR "WARNING: $longest and $id are of different sequences:\n$ref\n====\n$seq\n\n" if substr($ref, $f, $t-$f+1) ne substr($seq, $f, $t-$f+1);

	}
	my $raw = $raw{$longest};
	$raw =~ s/\s//g;
	my $ndiff = $raw=~s/([^\*.-])/\1/g;
	print "$type\t$size{$longest}\t$from\t$to\t$ndiff\t$raw\n";
	print "$type\t$size{$longest}\t$from\t$to\t$ndiff\t$ref\n";
}
