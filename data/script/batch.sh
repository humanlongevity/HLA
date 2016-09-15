ls ../raw/alignments/*_nuc.txt | parallel './script/000.alignment.nuc.pl {} > nuc/{/}'
rm nuc/[ABC]_*
ls nuc | parallel ./script/005.exon.pl nuc/{} exon/{.}.shift exon/{.}.faa dna/{.}.fna
ls exon/*.faa | perl -pe 's|.+/||; s/_nuc.faa//' | parallel ./script/006.split.pl exon/{}_nuc.faa temp/{} align/{}.tsv

cat align/DQB.tsv align/DRB.tsv align/DPB.tsv align/ClassI.tsv dqb2/DQB2.tsv > hla.tsv
cat exon/*.shift > hla.shift
cat exon/*.faa dqb2/DQB2.faa > hla.faa
#diamond index hla.faa to hla.dmnd based on the DIAMOND version you are using
