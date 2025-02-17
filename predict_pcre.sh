
##============== step 1. compute the gene counts of 6- to 10-mer motifs in test and control gene sets
## 6- and 7-mer motifs were provided as example data, users need to run 8- to 10-mer motifs
perl predict_pcre.pl --detect_kmer_motifs 6 ./input/  reg_control.fa reg_BraA01g023140.3.5C_targ.fa ../output_pCRE/kmers/reg_BraA01g023140.3.5C_targ.pcre6  &
perl predict_pcre.pl --detect_kmer_motifs 7 ./input/  reg_control.fa reg_BraA01g023140.3.5C_targ.fa ../output_pCRE/kmers/reg_BraA01g023140.3.5C_targ.pcre7  &
perl predict_pcre.pl --detect_kmer_motifs 8 ./input/  reg_control.fa reg_BraA01g023140.3.5C_targ.fa ../output_pCRE/kmers/reg_BraA01g023140.3.5C_targ.pcre8  &
perl predict_pcre.pl --detect_kmer_motifs 9 ./input/  reg_control.fa reg_BraA01g023140.3.5C_targ.fa ../output_pCRE/kmers/reg_BraA01g023140.3.5C_targ.pcre9  &
perl predict_pcre.pl --detect_kmer_motifs 10 ./input/ reg_control.fa reg_BraA01g023140.3.5C_targ.fa ../output_pCRE/kmers/reg_BraA01g023140.3.5C_targ.pcre10 &

##============== step 2. use R script (predict_pcre.fisher.R) to do fisher.test and adjust p-values
##  6: minimum k-mer,  10: maximun k-mer
##  In the R script, I output motifs with adjusted p-value < 0.05, 
##    users can reset the cutoff
R --no-restore --slave --vanilla --args output_pCRE/kmers/reg_BraA01g023140.3.5C_targ.pcre   6  10   <  predict_pcre.fisher.R 

##============== step 3. combine 6- to 10-mer results
less output_pCRE/kmers/reg_BraA01g023140.3.5C_targ.pcre6.p > output_pCRE/annot/reg_BraA01g023140.3.5C_targ.p
less output_pCRE/kmers/reg_BraA01g023140.3.5C_targ.pcre7.p | grep -v -P '^id'  >> output_pCRE/annot/reg_BraA01g023140.3.5C_targ.p
less output_pCRE/kmers/reg_BraA01g023140.3.5C_targ.pcre8.p | grep -v -P '^id'  >> output_pCRE/annot/reg_BraA01g023140.3.5C_targ.p
less output_pCRE/kmers/reg_BraA01g023140.3.5C_targ.pcre9.p | grep -v -P '^id'  >> output_pCRE/annot/reg_BraA01g023140.3.5C_targ.p
less output_pCRE/kmers/reg_BraA01g023140.3.5C_targ.pcre10.p | grep -v -P '^id'  >> output_pCRE/annot/reg_BraA01g023140.3.5C_targ.p

##============== step 4. transform dna seq of motifs to meme format and run tomtom to annotate k-mer motifs
mkdir output_pCRE/annot/reg_BraA01g023140.3.5C_targ

##  users need to install the software 'MEMEsuite'
## input/reg_motifsDB.list is a combination of list files under 'MEMEsuit/motif_database'
cat ~/software/MEMEsuite/motif_databases/list.arabd  \
	~/software/MEMEsuite/motif_databases/list.cisbp2_ath  \
	~/software/MEMEsuite/motif_databases/list.jaspar_plants  \
 > input/reg_motifsDB.list

## provid the absolute directory of 'MEMEsuit/motif_database'
perl predict_pcre.pl --annotate_pCRE    \
	output_pCRE/annot/  reg_BraA01g023140.3.5C_targ  \
    input/reg_motifsDB.list  \
	~/software/MEMEsuite/motif_databases/
