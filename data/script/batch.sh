ls ../raw/alignments/*_nuc.txt | parallel './script/000.alignment.nuc.pl {} > nuc/{/}'
rm nuc/[ABC]_*
ls nuc | parallel ./script/005.exon.pl nuc/{} exon/{.}.shift exon/{.}.faa
ls exon/*.faa | perl -pe 's|.+/||; s/_nuc.faa//' | parallel ./script/006.split.pl exon/{}_nuc.faa temp/{} align/{}.tsv
#./script/006.split.pl nuc/DPB_nuc.txt temp/DPB align/DPB.tsv
#./script/006.split.pl nuc/DQB_nuc.txt temp/DQB align/DQB.tsv
#./script/006.split.pl nuc/DRB_nuc.txt temp/DRB align/DRB.tsv
#./script/006.split.pl nuc/classI_nuc.txt temp/classI align/classI.tsv
