#!/usr/bin/env perl
use strict;
use List::Util qw(min);

die "usage: $0 fastq msa daa\n" unless $#ARGV == 2;
my $fastq_file = shift;
my $msa_file = shift;
my $daa_file = shift;

print STDERR "processing FASTQ file\n";
my %qseq;
my %qlen;
open(IN, $fastq_file) or die $!;
while(my $h = <IN>)
{
	my $seq = <IN>;
	my $temp = <IN>;
	my $temp = <IN>;
	$h =~ s/^\@(\S+).+/\1/s;
	chomp $seq;
	$qseq{$h} = $seq;
	$qlen{$h} = length $seq;
}
print STDERR "found ", scalar(keys %qseq), " reads";

print STDERR "processing MSA file\n";
my %tseq;
my %tlen;
my %pre;
my %suf;
my %type;
my %gene;
my %is_gene;
open(IN, $msa_file) or die $!;
while(<IN>)
{
	chomp;
	my ($id, $type, $exon, $size, $full, $pre, $suf, $suf2, $seq) = split(/\t/, $_);
	$tseq{$id} = $seq;
	$type{$id} = $type;
	$pre{$id} = $pre;
	$suf{$id} = $suf;
	$tlen{$id} = $seq=~s/[^-]//g;
	my $g = $1 if $type =~ m/^(\S+)\*/;
	$g = 'ClassI' if $g =~ m/^[ABC]$/;
	$gene{"$g\t$exon"}->{$id} = 1;
	$is_gene{$id} = "$g\t$exon";
}
print STDERR "found ", scalar(keys %tseq), " HLA exons\n";

print STDERR "processing DAA file\n";
open(IN, "diamond view -a '$daa_file' -o /dev/stdout |") or die $!;
my %mlen;
my %match;
my %matched;
my %nonspec;
while(<IN>)
{
	my ($q, $target, $identity, $len, $mis, $gap, $qs, $qe, $ts, $te, $e, $score) = split(/\t/, $_);
	next unless $identity == 100;
	my $len = $te - $ts + 1;
	my $g = $is_gene{$target};
	if($len > $mlen{$g}->{$q})
	{
		$mlen{$g}->{$q} = $len;
		$match{$g}->{$q} = [$target, $qs, $qe, $ts, $te];
		print STDERR "\t$target\n" unless $matched{$target};
	}
	$matched{$target}++ if $len >= $mlen{$g}->{$q};
	$nonspec{$q}++ if not($tseq{$target}) && $len >= $mlen{$g}->{$q};
}
print STDERR "matched to ", scalar(keys %matched), " HLA exons\n";
print STDERR scalar(keys %nonspec), " reads matched to HLA types not in the MSA file\n";

print STDERR "translating matches to MSA\n";
my %done;
for my $g(keys %match)
{
	print STDERR "\t$g\n";
	for my $q(keys %{$match{$g}})
	{
		print STDERR "\t\t$q\n";
		my ($target, $qs, $qe, $ts, $te) = @{$match{$g}->{$q}};
		my $spe = $nonspec{$q} ? 0 : 1;
		print STDERR "\t\t\t$q($qlen{$q}):$qs-$qe vs $target($tlen{$target}):$ts-$te\n";

		my $qseq = $qseq{$q};
		if($qs > $qe)
		{
			$qseq =~ tr/ATGC/TACG/;
			$qseq = reverse $qseq;
			print STDERR "\t\t\t\tflipping $qs-$qe to";
			$qs = $qlen{$q} - $qs + 1;
			$qe = $qlen{$q} - $qe + 1;
			print STDERR " $qs-$qe\n";
		}
		my $mseq = translate($qseq, $qs, $qe);

		my ($pre1, $pre2, $suf1, $suf2);
		if($ts == 1)
		{
			my $p1 = $qs-2;
			my $p2 = $qs-3;
			$pre1 = substr($qseq, $p1, 1) if $p1 >= 0;
			$pre2 = substr($qseq, $p2, 2) if $p2 >= 0;
		}
		if($te == $tlen{$target})
		{
			$suf1 = substr($qseq, $qe, 1);
			$suf2 = substr($qseq, $qe, 2);
		}
		print STDERR "\t\t\t\tchecking prefix or suffix leftover on the query:\n";
		print STDERR "\t\t\t\t\tpre1=$pre1, pre2=$pre2, suf1=$suf1, suf2=$suf2\n";

		my $tag = "$g-$pre2-$pre1-$mseq-$suf1-$suf2";
		print STDERR "\t\t\t\ttag code: $tag\n";
		next if $done{$tag};
		$done{$tag} = 1;
		print STDERR "\t\t\t\ttag not processed yet\n";

		my $left = min(int(($qs-1)/3), $ts-1);
		my $right = min(int(($qlen{$q}-$qe)/3), $tlen{$target} - $te);

		my $start = 0;
		my $end = 0;
		my $N = 0;
		my $n = 0;
		my @aa = split(//, $tseq{$target});
		for my $a(@aa)
		{
			$n++;
			$N++ unless $a eq '-';
			$start = $n if $N == $ts && not($start);
			if($N == $te)
			{
				$end = $n;
				last;
			}
		}
		my $tlen = $end - $start + 1;
		$start--;
		my $tseq = substr($tseq{$target}, $start, $tlen);
		print STDERR "\t\t\t\tconverting $ts-$te position to MSA $start-$end\n";

				
		# the sequence check part should be removed later. Only for debuging purpose
		#my $tseq2 = $tseq;
		#$tseq2 =~ s/-//g;
		#die "DIFFERENT: $q ($qs-$qe: $mseq) with $target ($start-$end: $tseq):\n$qseq$q}\n$tseq{$target}\n\n" if $mseq ne $tseq2;

		for my $t(keys %{$gene{$g}})
		{
			print STDERR "\t\t\t\t\tcomparing with $t\n";
			next unless $matched{$t};
			print STDERR "\t\t\t\t\tsaw this before ($matched{$t} times)\n";
			my $comp = substr($tseq{$t}, $start, $tlen);
			next unless $comp eq $tseq;
			print STDERR "\t\t\t\t\tsequence matched\n";
			my $pre_len = length $pre{$t};
			my ($pre, $suf);
			print STDERR "\t\t\t\t\tchecking prefix match\n" if $pre_len;
			if($pre_len == 1)
			{
				next if $pre1 && $pre{$t} ne $pre1;
				$pre = $pre1;
			}elsif($pre_len == 2)
			{
				next if $pre2 && $pre{$t} ne $pre2;
				$pre = $pre2;
			}
			print STDERR "\t\t\t\t\tprefix matched ($pre{$t} vs $pre)\n" if $pre_len;
			my $suf_len = length $suf{$t};
			print STDERR "\t\t\t\t\tchecking suffix match\n" if $suf_len;
			if($suf_len == 1)
			{
				next if $suf1 && $suf{$t} ne $suf1;
				$suf = $suf1;
			}elsif($suf_len == 2)
			{
				next if $suf2 && $suf{$t} ne $suf2;
				$suf = $suf2;
			}
			print STDERR "\t\t\t\t\tsuffix matched ($suf{$t} vs $suf)\n" if $pre_len;
			print "$q\t$qlen{$q}:$qs-$qe\t$t\t$tlen{$t}:$ts-$te\t$type{$t}\t$is_gene{$t}\t$spe\t$left\t$right\t$start\t$end\n";
		}
	}
}

sub translate
{
	my ($seq, $from, $to) = @_;
	my %codon = (
  	TTT => "F", TTC => "F", TTA => "L", TTG => "L",
  	TCT => "S", TCC => "S", TCA => "S", TCG => "S",
  	TAT => "Y", TAC => "Y", TAA => "X", TAG => "X",
  	TGT => "C", TGC => "C", TGA => "X", TGG => "W",
  	CTT => "L", CTC => "L", CTA => "L", CTG => "L",
  	CCT => "P", CCC => "P", CCA => "P", CCG => "P",
  	CAT => "H", CAC => "H", CAA => "Q", CAG => "Q",
  	CGT => "R", CGC => "R", CGA => "R", CGG => "R",
  	ATT => "I", ATC => "I", ATA => "I", ATG => "M",
  	ACT => "T", ACC => "T", ACA => "T", ACG => "T",
  	AAT => "N", AAC => "N", AAA => "K", AAG => "K",
  	AGT => "S", AGC => "S", AGA => "R", AGG => "R",
  	GTT => "V", GTC => "V", GTA => "V", GTG => "V",
  	GCT => "A", GCC => "A", GCA => "A", GCG => "A",
  	GAT => "D", GAC => "D", GAA => "E", GAG => "E",
  	GGT => "G", GGC => "G", GGA => "G", GGG => "G",
	);
	my $seq2 = substr($seq, $from - 1, $to - $from + 1);
	if($from > $to)
	{
		$seq2 = substr($seq, $to - 1, $from - $to + 1);
		$seq2 =~ tr/ATGC/TACG/;
		$seq2 = reverse $seq2;
	}
	$seq2 =~ s/(...)/$codon{$1}/ge;
	return $seq2;
}
