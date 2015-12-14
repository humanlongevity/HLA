#!/usr/bin/env perl
use strict;

die "usage: $0 infile temp-prefix outfile\n" unless $#ARGV == 2;
my ($in, $out, $final) = @ARGV;
open(IN, $in) or die $!;
my %exons;
while(my $h = <IN>)
{
	my $seq = <IN>;
	my $e = $1 if $h =~ m/-(E\d+)-/;
	$exons{$e}->{$h} = $seq;
}
for my $e(keys %exons)
{
	open OUT, ">$out-$e.faa" or die $!;
	for my $h(keys %{$exons{$e}})
	{
		print OUT $h;
		print OUT $exons{$e}->{$h};
	}
}

#system("ls $out-*.faa | parallel muscle -sv -in {} -out {.}.msa");

open(IN, "cat $out-*.msa |") or die $!;
my $name;
my %msa;
while(<IN>)
{
	chomp;
	if(m/^>(\S+)/)
	{
		$name = $1;
	}else
	{
		$msa{$name} .= $_;
	}
}

open(OUT, ">$final") or die $!;
for my $n(sort keys %msa)
{
	#A*32:19N-E7-L15-pre-GC-suf-G-suf2-	
	my ($type, $exon, $len, $tpre, $pre, $tsuf, $suf, $tsuf2, $suf2) = split(/-/, $n);
	$len =~ s/L//;
	my $seq = $msa{$n};
	my $LEN = length $seq;
	print OUT "$type\t$exon\t$len\t$LEN\t$pre\t$suf\t$suf2\t$seq\n";
}
close OUT;


