#!/usr/bin/env perl
use strict;

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

my %seq;
while(<>)
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
		$e++;
		my $pre = $1 if $exon =~ s/^(\S{1,2}) //;
		my $suf = $1 if $exon =~ s/ (\S{1,2})$//;
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
#		print ">$id-E$e-P$pos-L$len-pre-$pre-suf-$suf-suf2-$suf2\n$prot\n" if $len >= 10;
		$seq{"$id-E$e-L$len-pre-$pre-suf-$suf-suf2-$suf2"} = $prot if $len >= 10;
	}
}

my %done;
for my $id(sort keys %seq)
{
	my $newid = $id;
	$newid =~ s/^(.+?\*\d+:\d+):.+?(\D*)-E/\1\2-E/;
	if($done{$newid})
	{
		print STDERR "sequence differ between reported $done{$newid} and $id:\n$seq{$done{$newid}}\n$seq{$id}\n\n" if $seq{$done{$newid}} ne $seq{$id};
	}else
	{
		print ">$newid\n$seq{$id}\n";
		$done{$newid} = $id;
	}
}

		




__END__
A*01:01:01:01	1098	0	1923	1555	ATG GCC GTC ATG GCG CCC CGA ACC CTC CTC CTG CTA CTC TCG GGG GCC CTG GCC CTG ACC CAG ACC TGG GCG G|GC TCC CAC TCC ATG AGG TAT TTC TTC ACA TCC GTG TCC CGG CCC GGC CGC GGG GAG CCC CGC TTC ATC GCC GTG GGC TAC GTG GAC GAC ACG CAG TTC GTG CGG TTC GAC AGC GAC GCC GCG AGC CAG AAG ATG GAG CCG CGG GCG CCG TGG ATA GAG CAG GAG GGG CCG GAG TAT TGG GAC CAG GAG ACA CGG AAT ATG AAG GCC CAC TCA CAG ACT GAC CGA GCG AAC CTG GGG ACC CTG CGC GGC TAC TAC AAC CAG AGC GAG GAC G|GT TCT CAC ACC ATC CAG ATA ATG TAT GGC TGC GAC GTG GGG CCG GAC GGG CGC TTC CTC CGC GGG TAC CGG CAG GAC GCC TAC GAC GGC AAG GAT TAC ATC GCC CTG AAC GAG GAC CTG CGC TCT TGG ACC GCG GCG GAC ATG GCA GCT CAG ATC ACC AAG CGC AAG TGG GAG GCG GTC CAT GCG GCG GAG CAG CGG AGA GTC TAC CTG GAG GGC CGG TGC GTG GAC GGG CTC CGC AGA TAC CTG GAG AAC GGG AAG GAG ACG CTG CAG CGC ACG G|AC CCC CCC AAG ACA CAT ATG ACC CAC CAC CCC ATC TCT GAC CAT GAG GCC ACC CTG AGG TGC TGG GCC CTG GGC TTC TAC CCT GCG GAG ATC ACA CTG ACC TGG CAG CGG GAT GGG GAG GAC CAG ACC CAG GAC ACG GAG CTC GTG GAG ACC AGG CCT GCA GGG GAT GGA ACC TTC CAG AAG TGG GCG GCT GTG GTG GTG CCT TCT GGA GAG GAG CAG AGA TAC ACC TGC CAT GTG CAG CAT GAG GGT CTG CCC AAG CCC CTC ACC CTG AGA TGG G|AG CTG TCT TCC CAG CCC ACC ATC CCC ATC GTG GGC ATC ATT GCT GGC CTG GTT CTC CTT GGA GCT GTG ATC ACT GGA GCT GTG GTC GCT GCC GTG ATG TGG AGG AGG AAG AGC TCA G|AT AGA AAA GGA GGG AGT TAC ACT CAG GCT GCA A|GC AGT GAC AGT GCC CAG GGC TCT GAT GTG TCT CTC ACA GCT TGT AAA G|TG TGA 
A*01:01:01:02N	1098	0	1923	457	ATG GCC GTC ATG GCG CCC CGA ACC CTC CTC CTG CTA CTC TCG GGG GCC CTG GCC CTG ACC CAG ACC TGG GCG G|GC TCC CAC TCC ATG AGG TAT TTC TTC ACA TCC GTG TCC CGG CCC GGC CGC GGG GAG CCC CGC TTC ATC GCC GTG GGC TAC GTG GAC GAC ACG CAG TTC GTG CGG TTC GAC AGC GAC GCC GCG AGC CAG AAG ATG GAG CCG CGG GCG CCG TGG ATA GAG CAG GAG GGG CCG GAG TAT TGG GAC CAG GAG ACA CGG AAT ATG AAG GCC CAC TCA CAG ACT GAC CGA GCG AAC CTG GGG ACC CTG CGC GGC TAC TAC AAC CAG AGC GAG GAC G|GT TCT CAC ACC ATC CAG ATA ATG TAT GGC TGC GAC GTG GGG CCG GAC GGG CGC TTC CTC CGC GGG TAC CGG CAG GAC GCC TAC GAC GGC AAG GAT TAC ATC GCC CTG AAC GAG GAC CTG CGC TCT TGG ACC GCG GCG GAC ATG GCA GCT CAG ATC ACC AAG CGC AAG TGG GAG GCG GTC CAT GCG GCG GAG CAG CGG AGA GTC TAC CTG GAG GGC CGG TGC GTG GAC GGG CTC CGC AGA TAC CTG GAG AAC GGG AAG GAG ACG CTG CAG CGC ACG G|AC CCC CCC AAG ACA CAT ATG ACC CAC CAC CCC ATC TCT GAC CAT GAG GCC ACC CTG AGG TGC TGG GCC CTG GGC TTC TAC CCT GCG GAG ATC ACA CTG ACC TGG CAG CGG GAT GGG GAG GAC CAG ACC CAG GAC ACG GAG CTC GTG GAG ACC AGG CCT GCA GGG GAT GGA ACC TTC CAG AAG TGG GCG GCT GTG GTG GTG CCT TCT GGA GAG GAG CAG AGA TAC ACC TGC CAT GTG CAG CAT GAG GGT CTG CCC AAG CCC CTC ACC CTG AGA TGG G|AG CTG TCT TCC CAG CCC ACC ATC CCC ATC GTG GGC ATC ATT GCT GGC CTG GTT CTC CTT GGA GCT GTG ATC ACT GGA GCT GTG GTC GCT GCC GTG ATG TGG AGG AGG AAG AGC TCA G|AT AGA AAA GGA GGG AGT TAC ACT CAG GCT GCA A|GC AGT GAC AGT GCC CAG GGC TCT GAT GTG TCT CTC ACA GCT TGT AAA G|TG TGA 
A*01:01:02	546	240	1159	458	*** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|GC TCC CAC TCC ATG AGG TAT TTC TTC ACA TCC GTG TCC CGG CCC GGC CGC GGG GAG CCC CGC TTC ATC GCC GTG GGC TAC GTG GAC GAC ACG CAG TTC GTG CGG TTC GAC AGC GAC GCC GCG AGC CAG AAG ATG GAG CCG CGG GCG CCG TGG ATA GAG CAG GAG GGG CCG GAG TAT TGG GAC CAG GAG ACA CGG AAT ATG AAG GCC CAC TCA CAG ACT GAC CGA GCG AAC CTG GGG ACC CTG CGC GGC TAC TAC AAC CAG AGC GAG GAC G|GT TCT CAC ACC ATC CAG ATA ATG TAT GGC TGC GAC GTG GGG CCG GAC GGG CGC TTC CTC CGC GGG TAC CGG CAG GAC GCC TAC GAC GGC AAG GAT TAC ATC GCC CTG AAC GAG GAC CTG CGC TCT TGG ACC GCG GCG GAC ATG GCA GCT CAG ATT ACC AAG CGC AAG TGG GAG GCG GTC CAT GCG GCG GAG CAG CGG AGA GTC TAC CTG GAG GGC CGG TGC GTG GAC GGG CTC CGC AGA TAC CTG GAG AAC GGG AAG GAG ACG CTG CAG CGC ACG G|** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|** *** *** *** *** *** *** *** *** *** *** *|** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|** *** 
A*01:01:03	895	0	1536	458	ATG GCC GTC ATG GCG CCC CGA ACC CTC CTC CTG CTA CTC TCG GGG GCC CTG GCC CTG ACC CAG ACC TGG GCG G|GC TCC CAC TCC ATG AGG TAT TTC TTC ACA TCC GTG TCC CGG CCC GGC CGC GGG GAG CCC CGC TTC ATC GCC GTG GGC TAC GTG GAC GAC ACG CAG TTC GTG CGG TTC GAC AGC GAC GCC GCG AGC CAG AAG ATG GAG CCG CGG GCG CCG TGG ATA GAG CAG GAG GGG CCG GAG TAT TGG GAC CAG GAG ACA CGG AAT ATG AAG GCC CAC TCA CAG ACT GAC CGA GCG AAC CTG GGG ACC CTG CGC GGC TAC TAC AAC CAG AGC GAG GAC G|GT TCT CAC ACC ATC CAG ATA ATG TAT GGC TGC GAC GTG GGG CCG GAC GGG CGC TTC CTC CGC GGG TAC CGG CAG GAC GCC TAC GAC GGC AAG GAT TAC ATC GCC CTG AAC GAG GAC CTG CGC TCT TGG ACC GCG GCG GAC ATG GCA GCT CAG ATC ACC AAG CGC AAG TGG GAG GCG GTC CAT GCG GCG GAG CAG CGG AGA GTC TAC CTG GAG GGC CGG TGC GTG GAC GGG CTC CGC AGA TAC CTG GAG AAC GGG AAG GAG ACG CTG CAG CGC ACT G|AC CCC CCC AAG ACA CAT ATG ACC CAC CAC CCC ATC TCT GAC CAT GAG GCC ACC CTG AGG TGC TGG GCC CTG GGC TTC TAC CCT GCG GAG ATC ACA CTG ACC TGG CAG CGG GAT GGG GAG GAC CAG ACC CAG GAC ACG GAG CTC GTG GAG ACC AGG CCT GCA GGG GAT GGA ACC TTC CAG AAG TGG GCG GCT GTG GTG GTG CCT TCT GGA GAG GAG CAG AGA TAC ACC TGC CAT GTG CAG CAT GAG GGT CTG CCC AAG CCC CTC ACC CTG AGA TGG G|** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|** *** *** *** *** *** *** *** *** *** *** *|** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|** *** 
A*01:01:04	546	240	1159	458	*** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|GC TCC CAC TCC ATG AGG TAT TTC TTC ACA TCC GTG TCC CGG CCC GGC CGC GGG GAG CCC CGC TTC ATC GCC GTG GGC TAC GTG GAC GAC ACG CAG TTC GTG CGG TTC GAC AGC GAC GCC GCG AGC CAG AAG ATG GAG CCG CGG GCG CCG TGG ATA GAG CAG GAG GGG CCG GAG TAT TGG GAC CAG GAG ACA CGG AAT ATG AAG GCC CAC TCA CAG ACT GAC CGA GCG AAT CTG GGG ACC CTG CGC GGC TAC TAC AAC CAG AGC GAG GAC G|GT TCT CAC ACC ATC CAG ATA ATG TAT GGC TGC GAC GTG GGG CCG GAC GGG CGC TTC CTC CGC GGG TAC CGG CAG GAC GCC TAC GAC GGC AAG GAT TAC ATC GCC CTG AAC GAG GAC CTG CGC TCT TGG ACC GCG GCG GAC ATG GCA GCT CAG ATC ACC AAG CGC AAG TGG GAG GCG GTC CAT GCG GCG GAG CAG CGG AGA GTC TAC CTG GAG GGC CGG TGC GTG GAC GGG CTC CGC AGA TAC CTG GAG AAC GGG AAG GAG ACG CTG CAG CGC ACG G|** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|** *** *** *** *** *** *** *** *** *** *** *|** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|** *** 
A*01:01:05	822	240	1536	458	*** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|GC TCC CAC TCC ATG AGG TAT TTC TTC ACA TCC GTG TCC CGG CCC GGC CGC GGG GAG CCC CGC TTC ATC GCG GTG GGC TAC GTG GAC GAC ACG CAG TTC GTG CGG TTC GAC AGC GAC GCC GCG AGC CAG AAG ATG GAG CCG CGG GCG CCG TGG ATA GAG CAG GAG GGG CCG GAG TAT TGG GAC CAG GAG ACA CGG AAT ATG AAG GCC CAC TCA CAG ACT GAC CGA GCG AAC CTG GGG ACC CTG CGC GGC TAC TAC AAC CAG AGC GAG GAC G|GT TCT CAC ACC ATC CAG ATA ATG TAT GGC TGC GAC GTG GGG CCG GAC GGG CGC TTC CTC CGC GGG TAC CGG CAG GAC GCC TAC GAC GGC AAG GAT TAC ATC GCC CTG AAC GAG GAC CTG CGC TCT TGG ACC GCG GCG GAC ATG GCA GCT CAG ATC ACC AAG CGC AAG TGG GAG GCG GTC CAT GCG GCG GAG CAG CGG AGA GTC TAC CTG GAG GGC CGG TGC GTG GAC GGG CTC CGC AGA TAC CTG GAG AAC GGG AAG GAG ACG CTG CAG CGC ACG G|AC CCC CCC AAG ACA CAT ATG ACC CAC CAC CCC ATC TCT GAC CAT GAG GCC ACC CTG AGG TGC TGG GCC CTG GGC TTC TAC CCT GCG GAG ATC ACA CTG ACC TGG CAG CGG GAT GGG GAG GAC CAG ACC CAG GAC ACG GAG CTC GTG GAG ACC AGG CCT GCA GGG GAT GGA ACC TTC CAG AAG TGG GCG GCT GTG GTG GTG CCT TCT GGA GAG GAG CAG AGA TAC ACC TGC CAT GTG CAG CAT GAG GGT CTG CCC AAG CCC CTC ACC CTG AGA TGG G|** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|** *** *** *** *** *** *** *** *** *** *** *|** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|** *** 
A*01:01:06	546	240	1159	458	*** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|GC TCC CAC TCC ATG AGG TAT TTC TTC ACA TCC GTG TCC CGG CCC GGC CGC GGG GAG CCC CGC TTC ATC GCC GTG GGC TAC GTG GAC GAC ACG CAG TTC GTG CGG TTC GAC AGC GAC GCC GCG AGC CAG AAG ATG GAG CCG CGG GCG CCG TGG ATA GAG CAG GAG GGG CCG GAG TAT TGG GAC CAG GAG ACA CGG AAT ATG AAG GCC CAC TCA CAG ACT GAC CGA GCG AAC CTG GGG ACC CTG CGC GGC TAC TAC AAC CAG AGC GAG GAC G|GT TCT CAC ACC ATC CAG ATA ATG TAT GGC TGC GAC GTG GGG CCG GAC GGG CGC TTC CTC CGC GGG TAC CGG CAG GAC GCC TAC GAC GGC AAG GAT TAC ATC GCC CTG AAT GAG GAC CTG CGC TCT TGG ACC GCG GCG GAC ATG GCA GCT CAG ATC ACC AAG CGC AAG TGG GAG GCG GTC CAT GCG GCG GAG CAG CGG AGA GTC TAC CTG GAG GGC CGG TGC GTG GAC GGG CTC CGC AGA TAC CTG GAG AAC GGG AAG GAG ACG CTG CAG CGC ACG G|** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|** *** *** *** *** *** *** *** *** *** *** *|** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|** *** 
A*01:01:07	546	240	1159	458	*** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|GC TCC CAC TCC ATG AGG TAT TTC TTC ACA TCC GTG TCC CGG CCC GGC CGC GGG GAG CCC CGC TTC ATC GCC GTG GGC TAC GTG GAC GAC ACG CAG TTC GTG CGG TTC GAC AGC GAC GCC GCG AGC CAG AAG ATG GAG CCG CGG GCG CCG TGG ATA GAG CAG GAG GGG CCG GAA TAT TGG GAC CAG GAG ACA CGG AAT ATG AAG GCC CAC TCA CAG ACT GAC CGA GCG AAC CTG GGG ACC CTG CGC GGC TAC TAC AAC CAG AGC GAG GAC G|GT TCT CAC ACC ATC CAG ATA ATG TAT GGC TGC GAC GTG GGG CCG GAC GGG CGC TTC CTC CGC GGG TAC CGG CAG GAC GCC TAC GAC GGC AAG GAT TAC ATC GCC CTG AAC GAG GAC CTG CGC TCT TGG ACC GCG GCG GAC ATG GCA GCT CAG ATC ACC AAG CGC AAG TGG GAG GCG GTC CAT GCG GCG GAG CAG CGG AGA GTC TAC CTG GAG GGC CGG TGC GTG GAC GGG CTC CGC AGA TAC CTG GAG AAC GGG AAG GAG ACG CTG CAG CGC ACG G|** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|** *** *** *** *** *** *** *** *** *** *** *|** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|** *** 
A*01:01:08	546	240	1159	458	*** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|GC TCC CAC TCC ATG AGG TAT TTC TTC ACA TCC GTG TCC CGG CCC GGC CGC GGG GAG CCC CGC TTC ATC GCC GTG GGC TAC GTG GAC GAC ACG CAG TTC GTG CGG TTC GAC AGC GAC GCC GCG AGC CAG AAG ATG GAG CCG CGG GCG CCG TGG ATA GAG CAG GAG GGG CCG GAG TAT TGG GAC CAG GAG ACA CGG AAT ATG AAG GCC CAC TCA CAG ACT GAC CGA GCG AAC CTG GGG ACC CTG CGC GGC TAC TAC AAC CAG AGC GAG GAC G|GT TCT CAC ACC ATC CAG ATA ATG TAT GGC TGC GAC GTG GGG CCG GAC GGG CGC TTC CTC CGC GGG TAC CGG CAG GAC GCC TAC GAC GGC AAG GAT TAC ATC GCC CTG AAC GAG GAC CTG CGC TCT TGG ACC GCG GCG GAC ATG GCA GCT CAG ATC ACC AAG CGC AAG TGG GAG GCG GTC CAT GCG GCG GAG CAG CGG AGA GTC TAC CTG GAG GGC CGG TGC GTG GAC GGG CTC CGC AGA TAC TTG GAG AAC GGG AAG GAG ACG CTG CAG CGC ACG G|** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|** *** *** *** *** *** *** *** *** *** *** *|** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|** *** 
A*01:01:09	546	240	1159	458	*** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|GC TCC CAC TCC ATG AGG TAT TTC TTC ACA TCC GTG TCC CGG CCC GGG CGC GGG GAG CCC CGC TTC ATC GCC GTG GGC TAC GTG GAC GAC ACG CAG TTC GTG CGG TTC GAC AGC GAC GCC GCG AGC CAG AAG ATG GAG CCG CGG GCG CCG TGG ATA GAG CAG GAG GGG CCG GAG TAT TGG GAC CAG GAG ACA CGG AAT ATG AAG GCC CAC TCA CAG ACT GAC CGA GCG AAC CTG GGG ACC CTG CGC GGC TAC TAC AAC CAG AGC GAG GAC G|GT TCT CAC ACC ATC CAG ATA ATG TAT GGC TGC GAC GTG GGG CCG GAC GGG CGC TTC CTC CGC GGG TAC CGG CAG GAC GCC TAC GAC GGC AAG GAT TAC ATC GCC CTG AAC GAG GAC CTG CGC TCT TGG ACC GCG GCG GAC ATG GCA GCT CAG ATC ACC AAG CGC AAG TGG GAG GCG GTC CAT GCG GCG GAG CAG CGG AGA GTC TAC CTG GAG GGC CGG TGC GTG GAC GGG CTC CGC AGA TAC CTG GAG AAC GGG AAG GAG ACG CTG CAG CGC ACG G|** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|** *** *** *** *** *** *** *** *** *** *** *|** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *|** *** 
