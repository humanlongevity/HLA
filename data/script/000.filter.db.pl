#!/usr/bin/env perl
use strict;

$/ = "\n>";
my %full;
my %seq;
my %len;
while(<>)
{
	s/^\s*>//;
	s/\s*>\s*$//;
	next unless $_;
	if(m/^\S+\s+(\w+\*\S+) (\d+) bp\n(.+)/s)
	{
		my $type = $1;
		my $len = $2;
		my $seq = $3;
		$seq =~ s/\s//g;
		my $len2 = length $seq;
		die "$len != $len2\n$_" unless $len == $len2;
		$seq =~ s/X+$//;
		$len = length $seq;
		$type =~ s/(:\d+)\b.+/\1/ unless $type =~ m/[NS]$/;
		if($len > $len{$type})
		{
			if($len{$type} && index($seq, $seq{$type}) == -1)
			{
				print STDERR "this is different from the shorter version:\n$_\n\n$full{$type}\n\n\n\n";
			}
			$seq{$type} = $seq;
			$len{$type} = $len;
			$full{$type} = $_;
		}else
		{
			if(index($seq{$type}, $seq) == -1)
			{
				print STDERR "this is different from the shorter version:\n$_\n\n$full{$type}\n\n\n\n";
			}
		}
	}else
	{
		print STDERR $_;
	}
}

for my $t(sort keys %seq)
{
#	print ">$t-$len{$t}\n$seq{$t}\n";
	print ">$t\n$seq{$t}\n";
}
		

