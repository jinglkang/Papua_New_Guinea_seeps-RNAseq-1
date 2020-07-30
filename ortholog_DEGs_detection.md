Ortholgous gene detection by OrthoFinder && generate reads number matrix for DEGs detection  
========================================  
## before OrthoFinder
### all_contigs.second_pass.fa in setp 7 of DRAP as the de novo transcriptome for each species 
mkdir orthologue    
ls  
mv all_contigs.second_pass.fa Acura_tra_nuc.fa    
mv all_contigs.second_pass.fa Daru_tra_nuc.fa    
mv all_contigs.second_pass.fa ./Ocomp_tra_nuc.fa    
mv all_contigs.second_pass.fa ./Padel_tra_nuc.fa    
mv all_contigs.second_pass.fa ./Pmol_tra_nuc.fa    
  
### TransDecoder to obtain the protein sequences of ORFs  
nohup TransDecoder.LongOrfs -t Acura_tra_nuc.fa -O Acura_orf >Acura_transdecoder.process 2>&1 &    
nohup TransDecoder.LongOrfs -t Daru_tra_nuc.fa -O Daru_orf >Daru_transdecoder.process 2>&1 &    
nohup TransDecoder.LongOrfs -t Ocomp_tra_nuc.fa -O Ocomp_orf >Ocomp_transdecoder.process 2>&1 &    
nohup TransDecoder.LongOrfs -t Padel_tra_nuc.fa -O Padel_orf >Padel_transdecoder.process 2>&1 &    
nohup TransDecoder.LongOrfs -t Pmol_tra_nuc.fa -O Pmol_orf >Pmol_transdecoder.process 2>&1 &    
  
  
### put all the pep sequences of ORFs together in the same directory    
mkdir orthofinder_input_pep    
mv apoly_swath_Proteome.fasta orthofinder_input_pep/Apoly.fasta  
cp Acura_orf/longest_orfs.pep orthofinder_input_pep/Acura.fasta  
cp Daru_orf/longest_orfs.pep orthofinder_input_pep/Daru.fasta  
cp Ocomp_orf/longest_orfs.pep orthofinder_input_pep/Ocomp.fasta  
cp Padel_orf/longest_orfs.pep orthofinder_input_pep/Padel.fasta
cp Pmol_orf/longest_orfs.pep orthofinder_input_pep/Pmol.fasta
  
### copy these fasta file to new directory and change the gene name  
mkdir all_input # put all protein ORFs in all_input  
less Acura.fasta|perl -alne '$name=Acura;if (/>/) {$i++;$name1=$name."_"."$i";print ">$name1"}else{print}' >Acura-1.fasta  
mv Acura-1.fasta Acura.fasta  
less Apoly.fasta|perl -alne '$name=Apoly;if (/>/) {$i++;$name1=$name."_"."$i";print ">$name1"}else{print}' >Apoly-1.fasta  
mv Apoly-1.fasta Apoly.fasta  
less Daru.fasta|perl -alne '$name=Darua;if (/>/) {$i++;$name1=$name."_"."$i";print ">$name1"}else{print}' >Daru-1.fasta  
mv Daru-1.fasta Daru.fasta  
less Ocomp.fasta|perl -alne '$name=Padel;if (/>/) {$i++;$name1=$name."_"."$i";print ">$name1"}else{print}' >Ocomp-1.fasta  
mv Ocomp-1.fasta Ocomp.fasta  
less Padel.fasta|perl -alne '$name=Padel;if (/>/) {$i++;$name1=$name."_"."$i";print ">$name1"}else{print}' >Padel-1.fasta  
mv Padel-1.fasta Padel.fasta  
less Pmol.fasta|perl -alne '$name=Pmol;if (/>/) {$i++;$name1=$name."_"."$i";print ">$name1"}else{print}' >Pmol-1.fasta  
mv Pmol-1.fasta Pmol.fasta  

#################################  
## run OrthoFinder (-f: the input fasta file directory)  
#################################  
nohup orthofinder -f all_input -a 30 >orthofinder-process 2>&1 &  # the result will in "all_input" directory  
  
### select the orthologous group that contain transcripts of all species  
Orthogroups.GeneCount.tsv: in the "Orthogroups" directory of OrthoFinder results, blastp the selected groups to swiss-prot database respectively  
nohup perl blastp_uni.pl Orthogroups/Orthogroups.GeneCount.tsv >blast_process 2>&1 &  # the blastp result will be in blastp_result  
  

#################################  
## generate reads number matrix for each species  
#################################  

apoly: apoly_swath_Proteome.fasta (pep);   
apoly_primary_transcriptome_v1_shortHeader.fasta（nuc）  
### put the corresponding nucleotide sequences per species together (in orthofinder_input_nuc)  
mkdir orthofinder_input_nuc  
less ../Acura_orf/longest_orfs.cds|perl -alne '$name=Acura;if (/>/) {$i++;$name1=$name."_"."$i";print ">$name1"}else{print}' >Acura_nuc.fa  
less ../Daru_orf/longest_orfs.cds|perl -alne '$name=Daru;if (/>/) {$i++;$name1=$name."_"."$i";print ">$name1"}else{print}' >Daru_nuc.fa  
less ../Ocomp_orf/longest_orfs.cds|perl -alne '$name=Ocomp;if (/>/) {$i++;$name1=$name."_"."$i";print ">$name1"}else{print}' >Ocomp_nuc.fa  
less ../Padel_orf/longest_orfs.cds|perl -alne '$name=Padel;if (/>/) {$i++;$name1=$name."_"."$i";print ">$name1"}else{print}' >Padel_nuc.fa  
less ../Pmol_orf/longest_orfs.cds|perl -alne '$name=Pmol;if (/>/) {$i++;$name1=$name."_"."$i";print ">$name1"}else{print}' >Pmol_nuc.fa  
less apoly_primary_transcriptome_v1_shortHeader.fasta|perl -alne '$name=Apoly;if (/>/) {$i++;$name1=$name."_"."$i";print ">$name1"}else{print}' >Apoly_nuc.fa  
  
### RSEM caluculate reads number/per gene  
### build index  
nohup rsem-prepare-reference --bowtie2 Acura_nuc.fa Acura --bowtie2 &  
nohup rsem-prepare-reference --bowtie2 Apoly_nuc.fa Apoly --bowtie2 &  
nohup rsem-prepare-reference --bowtie2 Daru_nuc.fa Daru --bowtie2 &  
nohup rsem-prepare-reference --bowtie2 Ocomp_nuc.fa Ocomp --bowtie2 &  
nohup rsem-prepare-reference --bowtie2 Padel_nuc.fa Padel --bowtie2 &  
nohup rsem-prepare-reference --bowtie2 Pmol_nuc.fa Pmol --bowtie2 &  
  
#### Acura  
mkdir Acura_RSEM_output  
cd Acura_RSEM_output/  
vi RSEM_Acura.sh  
for Acura in ../Acura*_1.fastq.gz; do name=${Acura/_1.fastq.gz};name=${name##*/};name2=${Acura/_1/_2};rsem-calculate-expression -p 12 --bowtie2 --paired-end ${Acura} ${name2} ../orthofinder_input_nuc/Acura ${name};done  
nohup bash -x RSEM_Acura.sh  >RSEM-process 2>&1 &  

#### Apoly  
mkdir Apoly_RSEM_output  
vi RSEM_Apoly.sh  
for Apoly in ../Apoly*_1.fastq.gz; do name=${Apoly/_1.fastq.gz};name=${name##*/};name2=${Apoly/_1/_2};rsem-calculate-expression -p 16 --bowtie2 --paired-end ${Apoly} ${name2} ../orthofinder_input_nuc/Apoly ${name};done  
nohup bash -x RSEM_Apoly.sh  >RSEM-process 2>&1 &  

#### Daru  
mkdir Daru_RSEM_output  
cd Daru_RSEM_output  
vi RSEM_Daru.sh  
for Daru in ../Daru*_1.fastq.gz; do name=${Daru/_1.fastq.gz};name=${name##*/};name2=${Daru/_1/_2};rsem-calculate-expression -p 12 --bowtie2 --paired-end ${Daru} ${name2} ../orthofinder_input_nuc/Daru ${name};done  
nohup bash -x RSEM_Daru.sh  >RSEM-process 2>&1 &  

#### Ocomp  
mkdir Ocomp_RSEM_output  
cd Ocomp_RSEM_output/  
vi RSEM_Ocomp.sh  
for Ocomp in ../Ocomp*_1.fastq.gz; do name=${Ocomp/_1.fastq.gz};name=${name##*/};name2=${Ocomp/_1/_2};rsem-calculate-expression -p 12 --bowtie2 --paired-end ${Ocomp} ${name2} ../orthofinder_input_nuc/Ocomp ${name};done  
nohup bash -x RSEM_Ocomp.sh  >RSEM-process 2>&1 &  

#### Padel  
cd Padel_RSEM_output/  
vi RSEM_Padel.sh  
for Padel in ../Padel*_1.fastq.gz; do name=${Padel/_1.fastq.gz};name=${name##*/};name2=${Padel/_1/_2};rsem-calculate-expression -p 12 --bowtie2 --paired-end ${Padel} ${name2} ../orthofinder_input_nuc/Padel ${name};done  
nohup bash -x RSEM_Padel.sh  >RSEM-process 2>&1 &  

#### Pmol  
mkdir Pmol_RSEM_output  
cd Pmol_RSEM_output/  
vi RSEM_Pmol.sh  
for Pmol in ../Pmol*_1.fastq.gz; do name=${Pmol/_1.fastq.gz};name=${name##*/};name2=${Pmol/_1/_2};rsem-calculate-expression -p 12 --bowtie2 --paired-end ${Pmol} ${name2} ../orthofinder_input_nuc/Pmol ${name};done  
  
##########################################################  
## generate the initial reads number matrix of each species to select the better transcript with more reads number  
##########################################################  
  
#### Acura  
perl -e '@files=<*.genes.results>;print "../merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
sh merge.sh >total_Acura.gene.matrix  

#### Apoly  
perl -e '@files=<*.genes.results>;print "../merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
sh merge.sh >total_Apoly.gene.matrix  

#### Daru  
perl -e '@files=<*.genes.results>;print "../merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
sh merge.sh >total_Daru.gene.matrix  

#### Ocomp  
perl -e '@files=<*.genes.results>;print "../merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
sh merge.sh >total_Ocomp.gene.matrix  

#### Padel  
perl -e '@files=<*.genes.results>;print "../merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
sh merge.sh >total_Padel.gene.matrix  

#### Pmol  
perl -e '@files=<*.genes.results>;print "../merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
sh merge.sh >total_Pmol.gene.matrix  
  
## copy these reads number to ./six_fishes_reads_num  
mkdir six_fishes_reads_num  
cp total_Acura.gene.matrix total_Apoly.gene.matrix total_Dura.gene.matrix total_Ocomp.gene.matrix total_Padel.gene.matrix total_Pmol.gene.matrix > six_fishes_reads_num  

###########################################################################   
## obtain the representative sequence per species in each orthologous group   
###########################################################################  
nohup perl get_best_blast_orthogroup.pl -blast_result=blastp_result -reads_matrix=six_fishes_reads_num &  
ls -l | grep "^-" | wc -l  
14634 # finally there are 14634 qualified orthologous group  
  
## obtain the new reference per species by concatenated the final orthologous gene set  
./get_sequences_ref.pl -input=final_blast_orth_group -nuc=\~/Desktop/PapueNewGuinea-new/orthologue/orthofinder_input_nuc -output=final_reference  

########################################################################################  
## mapping again: map against the new references (concatenated by the orthologous genes)  
## RSEM estimate the reads number once again  
########################################################################################  
## index references in ./final_reference
nohup rsem-prepare-reference --bowtie2 Acura.fa Acura --bowtie2 &  
nohup rsem-prepare-reference --bowtie2 Apoly.fa Apoly --bowtie2 &  
nohup rsem-prepare-reference --bowtie2 Daru.fa Daru --bowtie2 &  
nohup rsem-prepare-reference --bowtie2 Ocomp.fa Ocomp --bowtie2 &  
nohup rsem-prepare-reference --bowtie2 Padel.fa Padel --bowtie2 &  
nohup rsem-prepare-reference --bowtie2 Pmol.fa Pmol --bowtie2 &  
  
#### Acura  
mkdir Acura_final_RSEM_output  
cd Acura_final_RSEM_output  
vi RSEM_Acura.sh  
for Acura in ../Acura*_1.fastq.gz; do name=${Acura/_1.fastq.gz};name=${name##*/};name2=${Acura/_1/_2};rsem-calculate-expression -p 12 --bowtie2 --paired-end ${Acura} ${name2} ../final_reference_3/Acura ${name};done  
nohup bash -x RSEM_Acura.sh  >RSEM-process 2>&1 &  

#### Apoly  
mkdir Apoly_final_RSEM_output  
cd Apoly_final_RSEM_output  
vi RSEM_Apoly.sh  
for Apoly in ../Apoly*_1.fastq.gz; do name=${Apoly/_1.fastq.gz};name=${name##*/};name2=${Apoly/_1/_2};rsem-calculate-expression -p 16 --bowtie2 --paired-end ${Apoly} ${name2} ../final_reference_3/Apoly ${name};done  
nohup bash -x RSEM_Apoly.sh  >RSEM-process 2>&1 &  

#### Daru  
mkdir Daru_final_RSEM_output  
cd Daru_final_RSEM_output  
vi RSEM_Daru.sh  
for Daru in ../Daru*_1.fastq.gz; do name=${Daru/_1.fastq.gz};name=${name##*/};name2=${Daru/_1/_2};rsem-calculate-expression -p 12 --bowtie2 --paired-end ${Daru} ${name2} ../final_reference_3/Daru ${name};done  
nohup bash -x RSEM_Daru.sh  >RSEM-process 2>&1 &  

#### Ocomp  
mkdir Ocomp_final_RSEM_output  
cd Ocomp_final_RSEM_output  
vi RSEM_Ocomp.sh  
for Ocomp in ../Ocomp*_1.fastq.gz; do name=${Ocomp/_1.fastq.gz};name=${name##*/};name2=${Ocomp/_1/_2};rsem-calculate-expression -p 12 --bowtie2 --paired-end ${Ocomp} ${name2} ../final_reference_3/Ocomp ${name};done  
nohup bash -x RSEM_Ocomp.sh  >RSEM-process 2>&1 &  

#### Padel  
mkdir Padel_final_RSEM_output  
cd Padel_final_RSEM_output  
vi RSEM_Padel.sh  
for Padel in ../Padel*_1.fastq.gz; do name=${Padel/_1.fastq.gz};name=${name##*/};name2=${Padel/_1/_2};rsem-calculate-expression -p 12 --bowtie2 --paired-end ${Padel} ${name2} ../final_reference_3/Padel ${name};done  
nohup bash -x RSEM_Padel.sh  >RSEM-process 2>&1 &  

#### Pmol  
mkdir Pmol_final_RSEM_output_3  
cd Pmol_final_RSEM_output_3  
vi RSEM_Pmol.sh  
for Pmol in ../Pmol*_1.fastq.gz; do name=${Pmol/_1.fastq.gz};name=${name##*/};name2=${Pmol/_1/_2};rsem-calculate-expression -p 12 --bowtie2 --paired-end ${Pmol} ${name2} ../final_reference_3/Pmol ${name};done  
nohup bash -x RSEM_Pmol.sh  >RSEM-process 2>&1 &  

#######################################################  
## put the reads number of all samples in the same file  
#######################################################  

### put all *.genes.results to the same directory  
mkdir final_RSEM_genes_results  
cp ../Acura_final_RSEM_output/*.genes.results ./  
cp ../Apoly_final_RSEM_output/*.genes.results ./  
cp ../Daru_final_RSEM_output/*.genes.results ./  
cp ../Ocomp_final_RSEM_output/*.genes.results ./  
cp ../Padel_final_RSEM_output/*.genes.results ./  
cp ../Pmol_final_RSEM_output/*.genes.results ./  

### merge_RSEM_frag_counts_single_table.pl to get the integrated reads number matrix  
perl -e '@files=<*.genes.results>;print "./merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
sh merge.sh >total_species.gene.matrix  
  
  
#######################################################  
## DEGs detection in DESeq2, generate sample information 
#######################################################  
less coldata.txt|perl -alne 's/\s+$//g;print "$_\tbatch" if /^\s+/;$name=$F[1]."_".$F[2];print "$_\t$name" if ! /^\s+/' >coldata.txt  


### generate reads number matrix, the reads number per transcript per species should be integer  
#### Acura  
perl -e '@files=<*.genes.results>;print "./merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
vi merge_RSEM_frag_counts_single_table.pl ## int()  
sh merge.sh >Acura.reads_num.matrix  

## merge_RSEM_frag_counts_single_table.pl to generate reads number matrix  
#### Apoly  
perl -e '@files=<*.genes.results>;print "./merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
sh merge.sh >Apoly.reads_num.matrix  

#### Daru  
perl -e '@files=<*.genes.results>;print "./merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
sh merge.sh >Daru.reads_num.matrix  

#### Ocomp  
perl -e '@files=<*.genes.results>;print "./merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
sh merge.sh >Ocomp.reads_num.matrix  

#### Padel  
perl -e '@files=<*.genes.results>;print "./merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
sh merge.sh >Padel.reads_num.matrix  
  
  
######################################  
## GABAa receptor genes  
######################################  
### Phylogeny tree: GABAa receptor genes in all species  
### filter out the GABAa receptor genes sequences of all species  
mkdir GABAa_genes  
cd GABAa_genes  
cat Acura.fasta Apoly.fasta Daru.fasta Ocomp.fasta Padel.fasta Pmol.fasta >all.fasta  
vi temp3.pl  

open fil2, "all.fasta";  
while (<fil2>) {  
	
        chomp;  
	if (/>/) {  
		s/>//;  
		$seq_name=$_;  
	} else {  
		$hash{$seq_name}.=$_;  
	}  
}   

open fil1, "GABAa_genes.txt";  from GABAa_genes.txt to get the GABAa receptor genes sequences of all six species  
while (<fil1>) {  

	chomp;  
	@a=split;  
	$file_name=$a[0].".blastp_result";  
	open fil, "../final_blast_orth_group_4/$file_name";  
	while (<fil>) {  
		chomp;  
		@b=split;  
		$Apoly_id=$b[1] if /Apoly/;  
		$Acura_id=$b[1] if /Acura/;  
		$Ocomp_id=$b[1] if /Ocomp/;  
		$Daru_id=$b[1] if /Daru/;  
		$Padel_id=$b[1] if /Padel/;  
		$Pmol_id=$b[1] if /Pmol/;  
	}  
	print ">$a[1]_Apoly\n$hash{$Apoly_id}\n";  
	print ">$a[1]_Acura\n$hash{$Acura_id}\n";  
	print ">$a[1]_Ocomp\n$hash{$Ocomp_id}\n";  
	print ">$a[1]_Daru\n$hash{$Daru_id}\n";  
	print ">$a[1]_Padel\n$hash{$Padel_id}\n";  
	print ">$a[1]_Pmol\n$hash{$Pmol_id}\n";  
}  
perl temp3.pl >all_species_GABAa_genes.fas  
### clustal alignment  
clustalo -i all_species_GABAa_genes.fas -t Protein -o all_species_GABAa_genes_align.fas --outfmt=fa  
vi temp2.pl  
open fil, "$ARGV[0]";  
while (<fil>) {  

        chomp;  
        if (/>/) {  
                s/>//;  
                $name=$_;  
        } else {  
                $seq{$name}.=$_;  
        }  
}  
foreach $key (sort keys %seq) {

        print ">$key\n$seq{$key}\n";  
}  
perl temp2.pl all_species_GABAa_genes_align.fas >all_species_GABAa_genes_align-final.fas  
### fasta2phy.pl tranform fasta format to phy format  
perl fasta2phy.pl all_species_GABAa_genes_align-final.fas >all_species_GABAa_genes_align-final.phy  
nohup raxmlHPC -f a -x 12345 -p 12345 -# 1000 -m PROTGAMMAAUTO -s all_species_GABAa_genes_align-final.phy -n all_species_GABAa_genes_align-final >raxml.process 2>&1 &  
