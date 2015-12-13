#!/usr/bin/env perl
use strict;

die "usage: $0 msa daa\n" unless $#ARGV == 1;
my $msa_file = shift;
my $daa_file = shift;

print STDERR "processing MSA file\n";
my %seq;
my %start;
my %end;
open(IN, $msa_file) or die $!;
while(<IN>)
{
	chomp;
	my ($id, $size, $from, $to, $diff, $seq) = split(/\t/, $_);
	$seq =~ s/\./-/g;
	$seq =~ s/\*/./g;
	$seq{$id} = $seq;
	$start{$id} = $from;
	$end{$id} = $to;
}
print STDERR "found ", scalar(keys %seq), " HLA types\n";

print STDERR "processing DAA file\n";
open(IN, "diamond view -a '$daa_file' -f sam -o /dev/stdout |") or die $!;
my %match;
my %matched;
my %mpos;
my %mseq;
my %mlen;
my %nonspec;
while(<IN>)
{
	next if m/^@/;
	my @f = split(/\t/, $_);
	next unless $f[16] =~ m/:100$/;
	my $q = $f[0];
	my $seq = $f[9];
	my $len = length $seq;
	if($len > $mlen{$q})
	{
		$mlen{$q} = $len;
		$mseq{$q} = $seq;
		delete $match{$q} if $match{$q};
	}
	my $t = $f[2];
	if(not(exists $start{$t}))
	{
		$nonspec{$q}++;
	}else
	{
		$matched{$t}++;
		if(not $match{$q})
		{
			$match{$q} = $t;
			$mpos{$q} = $start{$t} + $f[3] - 1;
		}
	}
}
print STDERR scalar(keys %match), " reads matched to ", scalar(keys %matched), " HLA types\n";
print STDERR scalar(keys %nonspec), " reads matched to HLA types not in the MSA file\n";

print STDERR "translating matches to MSA\n";
for my $q(keys %match)
{
	my $t = $match{$q};
	my $upstream = substr($seq{$t}, 0, $mpos{$q});
	my $upinsert = $upstream=~s/-/-/g;
	my $start = $upinsert + $mpos{$q};
	if(substr($seq{$t}, $start) =~ m/^(-+)/)
	{
		$start += $+[1] + 1;
	}
	my @ref = split(//, substr($seq{$t}, $start));
	my @seq = split(//, $mseq{$q});
	my $seq;
	my $qseq;
	while(my $rs = shift @ref)
	{
		last unless @seq;
		if($rs eq '-')
		{
			$seq .= $rs;
			$qseq .= $rs;
		}else
		{
			my $qs = shift @seq;
			$seq .= $rs;
			$qseq .= $qs;
		}
	}
	print STDERR "DIFFERENT: $q on $t at $mpos{$q} (=> $start):\nQ: $qseq\nT: $seq\nQuery: $mseq{$q}\nTARG-up: $upstream\nTARG: $seq{$t}\n\n" if $seq ne $qseq;
	$mseq{$q} = $qseq;
	$mlen{$q} = length $qseq;
	$mpos{$q} = $start;
}
print STDERR "done translating\n";

print STDERR "building match matrix\n";
my $spec = 0;
my $nspec = 0;
for my $q(keys %match)
{
	my $flag = $nonspec{$q} ? 0 : 1;
	if($flag)
	{
		$spec++;
	}else
	{
		$nspec++;
	}
	for my $t(keys %matched)
	{
		my $tseq = substr($seq{$t}, $mpos{$q}, $mlen{$q});
		my $good;
		if($mseq{$q} eq $tseq)
		{
			$good = $mlen{$q};
		}else
		{
			$tseq =~ s/\.//g;
			my $tlen = length($tseq);
			if($tlen && $mseq{$q} =~ m/$tseq/)
			{
				$good = $tlen;
			}
		}
		print "$q\t$t\t$mpos{$q}\t$mlen{$q}\t$good\t$flag\t$mseq{$q}\n" if $good;
	}
}
print STDERR "wrote $spec specific reads and $nspec non-specific reads\n";
