#!/usr/bin/env perl
# Author: Xie Chao
use strict;
use List::Util qw(min);
#use Text::Levenshtein::XS qw/distance/;
use FindBin;
my $bin=$FindBin::Bin;
my $root = $bin;
$root =~ s|/[^/]+?$||;

die "usage: $0 fastq output cores\n" unless $#ARGV >= 1;
my $doDNA = 1 if $#ARGV > 1;
my $fastq_file = shift;
my $out_file = shift;
my $cores = shift;

print "No of cores in align.pl ->>>>>>>> $cores";

print STDERR "processing FASTQ file\n";
my %qseq;
my %qlen;
my $open_fastq = $fastq_file;
if($open_fastq =~ m/\.gz$/)
{
	$open_fastq = "gzip -dc $fastq_file |";
}
open(IN, $open_fastq) or die $!;
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
print STDERR "\tfound ", scalar(keys %qseq), " reads\n";

print STDERR "processing MSA file\n";
my %tseq;
my %tlen;
my %pre;
my %suf;
my %type;
my %gene;
my %is_gene;
open(IN, "$root/data/hla.tsv") or die $!;
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
print STDERR "\tfound ", scalar(keys %tseq), " HLA exons\n";

print STDERR "processing FASTQ file\n";
open(IN, "diamond blastx -t . --threads $cores --max-hits 20000 --index-mode 1 --seg no --min-score 10 --top 20 --index-chunks 1 --db $root/data/hla --query $fastq_file --outfmt tab --quiet --out /dev/stdout |") or die $!;
my %mLEN;
my %mlen;
my %match;
my %matched;
my %nonspec;
while(<IN>)
{
	my ($qu, $target, $identity, $len, $mis, $gap, $qs, $qe, $ts, $te, $e, $score) = split(/\t/, $_);
	next unless $identity == 100;
	my $q = "$qu==$qs-$qe";
	my $qtag = "$qu==".abs($qs - $qe);
	my $len = $te - $ts + 1;
	my $g = $is_gene{$target};
#	if($len > $mlen{$g}->{$q})
	if($len >= $mLEN{$qu})
	{
		$mLEN{$qu} = $len;
		$mlen{$g}->{$q} = $len;
		$match{$g}->{$q} = [$target, $qs, $qe, $ts, $te, $len];
	}
	$matched{$target}++ if $len >= $mlen{$g}->{$q};
	$nonspec{$qtag}++ if not($tseq{$target}) && $len >= $mlen{$g}->{$q};
}
print STDERR "\tmatched to ", scalar(keys %matched), " HLA exons\n";
print STDERR "\t", scalar(keys %nonspec), " reads matched to HLA types not in the MSA file\n";

my %dna;
my %dna_name;
open(IN, "$root/data/hla.fna") or die $!;
while(<IN>)
{
	my @a = split(/\t/, $_);
	if($matched{$a[0]}){
		$dna{$a[0]}->{$a[2]} = 1;
		$dna_name{$a[0]}->{$a[2]}->{$a[1]} = 1;
	}
}

open OUT, ">$out_file" or die $!;
open OUT2, ">$out_file.dna" or die $!;
print STDERR "translating matches to MSA\n";
my %done;
my %done2;
my %save;
my %save2;
for my $g(keys %match)
{
	for my $q(keys %{$match{$g}})
	{
		my ($target, $qs, $qe, $ts, $te, $mlen) = @{$match{$g}->{$q}};
		my $qu = $1 if $q =~ m/(.+)==.+?/;
		my $qtag = "$qu==".abs($qs - $qe);
		next if $mlen < $mLEN{$qu};
		my $spe = $nonspec{$qtag} ? 0 : 1;

		my $qpart = "$qlen{$qu}:$qs-$qe";
		my $qseq = $qseq{$qu};
		if($qs > $qe)
		{
			$qseq =~ tr/ATGC/TACG/;
			$qseq = reverse $qseq;
			$qs = $qlen{$qu} - $qs + 1;
			$qe = $qlen{$qu} - $qe + 1;
		}
		my $mseq = translate($qseq, $qs, $qe);
		my $qcode = substr($qseq, $qs-1, $qe-$qs+1);

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

		my $tag = "$g-$pre2-$pre1-$qcode-$suf1-$suf2";
		next if $done{$tag.$qu};
		$done{$tag.$qu} = 1;

		if($save{$tag})
		{
			for my $k(keys %{$save{$tag}})
			{
				print OUT "$qu\t$qpart\t$qlen{$qu}:$qs-$qe\t$k";
			}
			if($doDNA){
				for my $k(keys %{$save2{$tag}})
				{
					print OUT2 "$qu\t$k";
				}
			}
			next;
		}

		my $left = min(int(($qs-1)/3), $ts-1);
		my $right = min(int(($qlen{$qu}-$qe)/3), $tlen{$target} - $te);

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

				
		# the sequence check part should be removed later. Only for debuging purpose
		#my $tseq2 = $tseq;
		#$tseq2 =~ s/-//g;
		#die "DIFFERENT: $q ($qs-$qe: $mseq) with $target ($start-$end: $tseq):\n$qseq$q}\n$tseq{$target}\n\n" if $mseq ne $tseq2;

		for my $t(keys %{$gene{$g}})
		{
			next unless $matched{$t};
			next if $done2{$qu}->{$t};
			my $comp = substr($tseq{$t}, $start, $tlen);
			next unless $comp eq $tseq;
			my $pre_len = length $pre{$t};
			my ($pre, $suf);
			if($pre_len == 1)
			{
				next if $pre1 && $pre{$t} ne $pre1;
				$pre = $pre1;
			}elsif($pre_len == 2)
			{
				next if $pre2 && $pre{$t} ne $pre2;
				$pre = $pre2;
			}
			my $suf_len = length $suf{$t};
			if($suf_len == 1)
			{
				next if $suf1 && $suf{$t} ne $suf1;
				$suf = $suf1;
			}elsif($suf_len == 2)
			{
				next if $suf2 && $suf{$t} ne $suf2;
				$suf = $suf2;
			}
			$done2{$qu}->{$t} = 1;
			my $dist = 1000;

			my $prefix = substr($tseq{$t}, 0, $start);
			my $tstart = scalar($prefix =~ s/\w//g) * 3;
			for my $tcode(keys %{$dna{$t}})
			{
				my $tt = substr($tcode, $tstart, ($te-$ts+1)*3);
#				my $edit = distance($qcode, $tt);
				my $edit = $qcode eq $tt ? 0 : 1;
				$dist = $edit if $edit < $dist;
				if($doDNA)
				{
					if($edit == 0)
					{
						for my $name(keys %{$dna_name{$t}->{$tcode}})

						{
							print OUT2 "$qu\t$t\t$name\n";
							$save2{$tag}->{"$t\t$name\n"} = 1;
						}
					}
				}
			}
			print OUT "$qu\t$qpart\t$qlen{$qu}:$qs-$qe\t$t\t$tlen{$t}\t$ts\t$te\t$dist\t$type{$t}\t$is_gene{$t}\t$spe\t$left\t$right\t$start\t$end\n";
			$save{$tag}->{"$t\t$tlen{$t}\t$ts\t$te\t$dist\t$type{$t}\t$is_gene{$t}\t$spe\t$left\t$right\t$start\t$end\n"} = 1;
		}
	}
}
close OUT;
close OUT2;
unlink "$out_file.dna" unless $doDNA;


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
