#!/usr/bin/env perl
use strict;

die "usage: $0 nuc_align.msa out.frame_shift out.faa\n" unless $#ARGV == 2;
my ($msa_file, $out_shift, $out_faa) = @ARGV;

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

open IN, $msa_file or die $!;
my %seq;
my %shift;
while(<IN>)
{
	chomp;
	my ($id, $len, $start, $to, $ndiff, $seq) = split(/\t/, $_);
	my @exons = split(/\|/, $seq);
	my $e = 0;
	# frame shifted types are translated past stop codon, and the exon after 
	# the frame-shifted exon are translated as original frame. This is deliberate,
	# because blastx type alinger can still map the part after frame shift
	for my $exon(@exons)
	{
		my $shift;
		$e++;
		my $pre = $1 if $exon =~ s/^(\S{1,2}) //;
		my $suf = $1 if $exon =~ s/ (\S{1,2})$//;
		my $i = 0;
		while(1)
		{
			my $replaced = $exon=~s/ (\S\S\S)(\S)/ \1 \2/;
			last unless $replaced;
		}
		for my $codon(split(/\s/, $exon))
		{
			$i++;
			if(length($codon) % 3)
			{
				$shift = $i;
				last;
			}
		}
		$exon =~ s/\s//g;
		$exon =~ s/(...)/\1 /g;
		$exon =~ s/ $//;
		my $suf2 = $1 if $exon =~ s/ (\S{1,2})$//;

		my @codons = split(/ /, $exon);
		my $prot;
		for my $co(@codons)
		{
			last if $co =~ m/\*/;
			my $aa = $codon{$co};
			$prot .= $aa;
#			last if $aa eq 'X';
		}
		my $len = length($prot);
		if($len >= 10)
		{
			my $id = "$id-E$e-L$len-pre-$pre-suf-$suf-suf2-$suf2";
			$seq{$id} = $prot;
			$shift{$id}->{"E$e"} = $shift if $shift;
		}
	}
}

open OUT, ">$out_faa" or die $!;
open SHIFT, ">$out_shift" or die $!;
my %done;
for my $id(sort keys %seq)
{
	my $newid = $id;
	$newid =~ s/^(.+?\*\d+:\d+):.+?(\D*)-E/\1\2-E/;
	if($done{$newid})
	{
		if($seq{$done{$newid}} ne $seq{$id})
		{
			print ">$id\n$seq{$id}\n";
			$done{$id} = $id;
			print STDERR "sequence differ between reported $done{$newid} and $id:\n$seq{$done{$newid}}\n$seq{$id}\n\n" if $seq{$done{$newid}} ne $seq{$id};
		}
	}else
	{
		print OUT ">$newid\n$seq{$id}\n";
		$done{$newid} = $id;
		for my $e(keys %{$shift{$id}})
		{
			print SHIFT "$newid\t$e\t$shift{$id}->{$e}\n";
		}
	}
}

