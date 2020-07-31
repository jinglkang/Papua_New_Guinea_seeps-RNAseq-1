Ortholgous gene detection by OrthoFinder && generate reads number matrix for DEGs detection  
========================================  
#####################  
## 1. before OrthoFinder
#####################  
### TransDecoder to detect the protein sequences of ORFs in the de novo transcriptome for each species
### _A. polyacanthus_ has a reference genome (needn't process in this step)
#### all_contigs.second_pass.fa in setp 7 of DRAP as the de novo transcriptome  
#### put these _de novo_ transcriptomes (Acura_tra_nuc.fa, Daru_tra_nuc.fa, Ocomp_tra_nuc.fa, Padel_tra_nuc.fa, Pmol_tra_nuc.fa) in ./orthologue
#### the fastq files of all samples in ./
mkdir orthologue    
cd orthologue  
nohup TransDecoder.LongOrfs -t Acura_tra_nuc.fa -O Acura_orf >Acura_transdecoder.process 2>&1 &    
nohup TransDecoder.LongOrfs -t Daru_tra_nuc.fa -O Daru_orf >Daru_transdecoder.process 2>&1 &    
nohup TransDecoder.LongOrfs -t Ocomp_tra_nuc.fa -O Ocomp_orf >Ocomp_transdecoder.process 2>&1 &    
nohup TransDecoder.LongOrfs -t Padel_tra_nuc.fa -O Padel_orf >Padel_transdecoder.process 2>&1 &    
nohup TransDecoder.LongOrfs -t Pmol_tra_nuc.fa -O Pmol_orf >Pmol_transdecoder.process 2>&1 &    
#### put all the predict pep sequences of ORFs together in the same directory from the TransDecoder result directory   
mkdir orthofinder_input_pep    
cp Acura_orf/longest_orfs.pep orthofinder_input_pep/Acura.fasta  
cp Daru_orf/longest_orfs.pep orthofinder_input_pep/Daru.fasta  
cp Ocomp_orf/longest_orfs.pep orthofinder_input_pep/Ocomp.fasta  
cp Padel_orf/longest_orfs.pep orthofinder_input_pep/Padel.fasta  
cp Pmol_orf/longest_orfs.pep orthofinder_input_pep/Pmol.fasta  
#### also put _A. polyacanthus_ protein sequence into ./orthofinder_input_pep
cp apoly_swath_Proteome.fasta orthofinder_input_pep/Apoly.fasta  
  
### copy these fasta file to ./all_input and change the gene name in each species  
mkdir all_input  
cd all_input  
less ../orthofinder_input_pep/Acura.fasta|perl -alne '$name=Acura;if (/>/) {$i++;$name1=$name."\_"."$i";print ">$name1"}else{print}' >Acura.fasta  
less ../orthofinder_input_pep/Apoly.fasta|perl -alne '$name=Apoly;if (/>/) {$i++;$name1=$name."\_"."$i";print ">$name1"}else{print}' >Apoly.fasta  
less ../orthofinder_input_pep/Daru.fasta|perl -alne '$name=Darua;if (/>/) {$i++;$name1=$name."\_"."$i";print ">$name1"}else{print}' >Daru.fasta  
less ../orthofinder_input_pep/Ocomp.fasta|perl -alne '$name=Ocomp;if (/>/) {$i++;$name1=$name."\_"."$i";print ">$name1"}else{print}' >Ocomp.fasta  
less ../orthofinder_input_pep/Padel.fasta|perl -alne '$name=Padel;if (/>/) {$i++;$name1=$name."\_"."$i";print ">$name1"}else{print}' >Padel.fasta  
less ../orthofinder_input_pep/Pmol.fasta|perl -alne '$name=Pmol;if (/>/) {$i++;$name1=$name."\_"."$i";print ">$name1"}else{print}' >Pmol.fasta  

#################################  
## 2. run OrthoFinder (-f: the input fasta file directory)  
#################################  
nohup orthofinder -f all_input -a 30 >orthofinder-process 2>&1 &  # the result will in "all_input" directory  
### select the orthologous group that contain transcripts of all species  
#### Orthogroups.GeneCount.tsv: in the "Orthogroups" directory of OrthoFinder results, blastp the selected groups to swiss-prot database respectively  
nohup perl blastp_uni.pl Orthogroups/Orthogroups.GeneCount.tsv >blast_process 2>&1 &  # the blastp result will be in ./blastp_result  
  

###########################################################################################  
## 3. Seltect the most representative transcript for each species in each orthologous group  
###########################################################################################  

########################################################################################### 
### 3.1 generate reads number matrix for each species based on nucleotide sequences of predict ORFs in the _de novo_ transcriptome  
###########################################################################################   
apoly: apoly_swath_Proteome.fasta (pep);   
apoly_primary_transcriptome_v1_shortHeader.fasta（nuc）  
#### put the corresponding nucleotide sequences per species together (in orthofinder_input_nuc)  
mkdir orthofinder_input_nuc  
cd orthofinder_input_nuc
less ../Acura_orf/longest_orfs.cds|perl -alne '$name=Acura;if (/>/) {$i++;$name1=$name."\_"."$i";print ">$name1"}else{print}' >Acura_nuc.fa  
less ../Daru_orf/longest_orfs.cds|perl -alne '$name=Daru;if (/>/) {$i++;$name1=$name."\_"."$i";print ">$name1"}else{print}' >Daru_nuc.fa  
less ../Ocomp_orf/longest_orfs.cds|perl -alne '$name=Ocomp;if (/>/) {$i++;$name1=$name."\_"."$i";print ">$name1"}else{print}' >Ocomp_nuc.fa  
less ../Padel_orf/longest_orfs.cds|perl -alne '$name=Padel;if (/>/) {$i++;$name1=$name."\_"."$i";print ">$name1"}else{print}' >Padel_nuc.fa  
less ../Pmol_orf/longest_orfs.cds|perl -alne '$name=Pmol;if (/>/) {$i++;$name1=$name."\_"."$i";print ">$name1"}else{print}' >Pmol_nuc.fa  
less apoly_primary_transcriptome_v1_shortHeader.fasta|perl -alne '$name=Apoly;if (/>/) {$i++;$name1=$name."\_"."$i";print ">$name1"}else{print}' >Apoly_nuc.fa  
#### RSEM caluculate reads number/per gene  
#### build index in ./orthofinder_input_nuc  
nohup rsem-prepare-reference --bowtie2 Acura_nuc.fa Acura --bowtie2 &  
nohup rsem-prepare-reference --bowtie2 Apoly_nuc.fa Apoly --bowtie2 &  
nohup rsem-prepare-reference --bowtie2 Daru_nuc.fa Daru --bowtie2 &  
nohup rsem-prepare-reference --bowtie2 Ocomp_nuc.fa Ocomp --bowtie2 &  
nohup rsem-prepare-reference --bowtie2 Padel_nuc.fa Padel --bowtie2 &  
nohup rsem-prepare-reference --bowtie2 Pmol_nuc.fa Pmol --bowtie2 &  
##########################################################  
### start to generate the initial reads number matrix of each species to select the better transcript with more reads number in each RSEM result direcotry  
#### merge_RSEM_frag_counts_single_table.pl to generate reads number matrix
##########################################################  
#### _A. curacao_  
mkdir Acura_RSEM_output  
cd Acura_RSEM_output/  
vi RSEM_Acura.sh  
for Acura in ../Acura*_1.fastq.gz; do name=${Acura/_1.fastq.gz};name=${name##*/};name2=${Acura/_1/_2};rsem-calculate-expression -p 12 --bowtie2 --paired-end ${Acura} ${name2} ../orthofinder_input_nuc/Acura ${name};done  
nohup bash -x RSEM_Acura.sh  >RSEM-process 2>&1 &  
perl -e '@files=<*.genes.results>;print "../merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
sh merge.sh >total_Acura.gene.matrix  

#### _A. polyacanthus_  
mkdir Apoly_RSEM_output  
Apoly_RSEM_output  
vi RSEM_Apoly.sh  
for Apoly in ../Apoly*_1.fastq.gz; do name=${Apoly/_1.fastq.gz};name=${name##*/};name2=${Apoly/_1/_2};rsem-calculate-expression -p 16 --bowtie2 --paired-end ${Apoly} ${name2} ../orthofinder_input_nuc/Apoly ${name};done  
nohup bash -x RSEM_Apoly.sh  >RSEM-process 2>&1 &  
perl -e '@files=<*.genes.results>;print "../merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
sh merge.sh >total_Apoly.gene.matrix  

#### _D. aruanus_  
mkdir Daru_RSEM_output  
cd Daru_RSEM_output  
vi RSEM_Daru.sh  
for Daru in ../Daru*_1.fastq.gz; do name=${Daru/_1.fastq.gz};name=${name##*/};name2=${Daru/_1/_2};rsem-calculate-expression -p 12 --bowtie2 --paired-end ${Daru} ${name2} ../orthofinder_input_nuc/Daru ${name};done  
nohup bash -x RSEM_Daru.sh  >RSEM-process 2>&1 &  
perl -e '@files=<*.genes.results>;print "../merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
sh merge.sh >total_Dura.gene.matrix  

#### _A. compressus_  
mkdir Ocomp_RSEM_output  
cd Ocomp_RSEM_output/  
vi RSEM_Ocomp.sh  
for Ocomp in ../Ocomp*_1.fastq.gz; do name=${Ocomp/_1.fastq.gz};name=${name##*/};name2=${Ocomp/_1/_2};rsem-calculate-expression -p 12 --bowtie2 --paired-end ${Ocomp} ${name2} ../orthofinder_input_nuc/Ocomp ${name};done  
nohup bash -x RSEM_Ocomp.sh  >RSEM-process 2>&1 &  
perl -e '@files=<*.genes.results>;print "../merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
sh merge.sh >total_Ocomp.gene.matrix  

#### _P. adelus_  
mkdir Padel_RSEM_output/  
cd Padel_RSEM_output/  
vi RSEM_Padel.sh  
for Padel in ../Padel*_1.fastq.gz; do name=${Padel/_1.fastq.gz};name=${name##*/};name2=${Padel/_1/_2};rsem-calculate-expression -p 12 --bowtie2 --paired-end ${Padel} ${name2} ../orthofinder_input_nuc/Padel ${name};done  
nohup bash -x RSEM_Padel.sh  >RSEM-process 2>&1 &  
perl -e '@files=<*.genes.results>;print "../merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
sh merge.sh >total_Padel.gene.matrix  

#### _P. molluscensis_  
mkdir Pmol_RSEM_output  
cd Pmol_RSEM_output/  
vi RSEM_Pmol.sh  
for Pmol in ../Pmol*_1.fastq.gz; do name=${Pmol/_1.fastq.gz};name=${name##*/};name2=${Pmol/_1/_2};rsem-calculate-expression -p 12 --bowtie2 --paired-end ${Pmol} ${name2} ../orthofinder_input_nuc/Pmol ${name};done  
perl -e '@files=<*.genes.results>;print "../merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
sh merge.sh >total_Pmol.gene.matrix  
  
## copy these reads number matrixes to ./six_fishes_reads_num  
mkdir six_fishes_reads_num  
cp Acura_RSEM_output/total_Acura.gene.matrix Apoly_RSEM_output/total_Apoly.gene.matrix Dura_RSEM_output/total_Dura.gene.matrix Ocomp_RSEM_output/total_Ocomp.gene.matrix Padel_RSEM_output/total_Padel.gene.matrix Pmol_RSEM_output/total_Pmol.gene.matrix > ./six_fishes_reads_num  

###########################################################################   
## 3.2 Obtain the representative sequence per species in each orthologous group   
###########################################################################  
#### finally there are 14634 qualified orthologous group, which would be in ./final_blast_orth_group  
nohup perl get_best_blast_orthogroup.pl -blast_result=blastp_result -reads_matrix=six_fishes_reads_num &  
ls -l | grep "^-" | wc -l  
14634  
#### obtain the new reference per species by concatenated the final orthologous gene set  
./get_sequences_ref.pl -input=./final_blast_orth_group -nuc=./orthofinder_input_nuc -output=final_reference  

########################################################################################  
## 4. mapping again: map against the new references (concatenated by the orthologous genes)  
## RSEM estimate the reads number once again  
########################################################################################  
### index references in ./final_reference
cd ./final_reference  
nohup rsem-prepare-reference --bowtie2 Acura.fa Acura --bowtie2 &  
nohup rsem-prepare-reference --bowtie2 Apoly.fa Apoly --bowtie2 &  
nohup rsem-prepare-reference --bowtie2 Daru.fa Daru --bowtie2 &  
nohup rsem-prepare-reference --bowtie2 Ocomp.fa Ocomp --bowtie2 &  
nohup rsem-prepare-reference --bowtie2 Padel.fa Padel --bowtie2 &  
nohup rsem-prepare-reference --bowtie2 Pmol.fa Pmol --bowtie2 &  
#### _A. curacao_  
mkdir Acura_final_RSEM_output  
cd Acura_final_RSEM_output  
vi RSEM_Acura.sh  
for Acura in ../Acura*_1.fastq.gz; do name=${Acura/_1.fastq.gz};name=${name##*/};name2=${Acura/_1/_2};rsem-calculate-expression -p 12 --bowtie2 --paired-end ${Acura} ${name2} ../final_reference_3/Acura ${name};done  
nohup bash -x RSEM_Acura.sh  >RSEM-process 2>&1 &  
perl -e '@files=<*.genes.results>;print "./merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
vi merge_RSEM_frag_counts_single_table.pl ## int()  
sh merge.sh >Acura.reads_num.matrix  

#### _A. polyacanthus_  
mkdir Apoly_final_RSEM_output  
cd Apoly_final_RSEM_output  
vi RSEM_Apoly.sh  
for Apoly in ../Apoly*_1.fastq.gz; do name=${Apoly/_1.fastq.gz};name=${name##*/};name2=${Apoly/_1/_2};rsem-calculate-expression -p 16 --bowtie2 --paired-end ${Apoly} ${name2} ../final_reference_3/Apoly ${name};done  
nohup bash -x RSEM_Apoly.sh  >RSEM-process 2>&1 &  
perl -e '@files=<*.genes.results>;print "./merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
sh merge.sh >Apoly.reads_num.matrix  

#### _D. aruanus_  
mkdir Daru_final_RSEM_output  
cd Daru_final_RSEM_output  
vi RSEM_Daru.sh  
for Daru in ../Daru*_1.fastq.gz; do name=${Daru/_1.fastq.gz};name=${name##*/};name2=${Daru/_1/_2};rsem-calculate-expression -p 12 --bowtie2 --paired-end ${Daru} ${name2} ../final_reference_3/Daru ${name};done  
nohup bash -x RSEM_Daru.sh  >RSEM-process 2>&1 &  
perl -e '@files=<*.genes.results>;print "./merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
sh merge.sh >Daru.reads_num.matrix  

#### _A. compressus_  
mkdir Ocomp_final_RSEM_output  
cd Ocomp_final_RSEM_output  
vi RSEM_Ocomp.sh  
for Ocomp in ../Ocomp*_1.fastq.gz; do name=${Ocomp/_1.fastq.gz};name=${name##*/};name2=${Ocomp/_1/_2};rsem-calculate-expression -p 12 --bowtie2 --paired-end ${Ocomp} ${name2} ../final_reference_3/Ocomp ${name};done  
nohup bash -x RSEM_Ocomp.sh  >RSEM-process 2>&1 &  
perl -e '@files=<*.genes.results>;print "./merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
sh merge.sh >Ocomp.reads_num.matrix  

#### _P. adelus_  
mkdir Padel_final_RSEM_output  
cd Padel_final_RSEM_output  
vi RSEM_Padel.sh  
for Padel in ../Padel*_1.fastq.gz; do name=${Padel/_1.fastq.gz};name=${name##*/};name2=${Padel/_1/_2};rsem-calculate-expression -p 12 --bowtie2 --paired-end ${Padel} ${name2} ../final_reference_3/Padel ${name};done  
nohup bash -x RSEM_Padel.sh  >RSEM-process 2>&1 &  
perl -e '@files=<*.genes.results>;print "./merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
sh merge.sh >Padel.reads_num.matrix  

#### _P. molluscensis_  
mkdir Pmol_final_RSEM_output_3  
cd Pmol_final_RSEM_output_3  
vi RSEM_Pmol.sh  
for Pmol in ../Pmol*_1.fastq.gz; do name=${Pmol/_1.fastq.gz};name=${name##*/};name2=${Pmol/_1/_2};rsem-calculate-expression -p 12 --bowtie2 --paired-end ${Pmol} ${name2} ../final_reference_3/Pmol ${name};done  
nohup bash -x RSEM_Pmol.sh  >RSEM-process 2>&1 &  

#######################################################  
### put the reads number of all samples into the same tab file  
#######################################################  
#### put all *.genes.results to the same directory  
mkdir final_RSEM_genes_results  
cd inal_RSEM_genes_results  
cp ../Acura_final_RSEM_output/*.genes.results ./  
cp ../Apoly_final_RSEM_output/*.genes.results ./  
cp ../Daru_final_RSEM_output/*.genes.results ./  
cp ../Ocomp_final_RSEM_output/*.genes.results ./  
cp ../Padel_final_RSEM_output/*.genes.results ./  
cp ../Pmol_final_RSEM_output/*.genes.results ./  
#### merge_RSEM_frag_counts_single_table.pl to get the integrated reads number matrix  
perl -e '@files=<*.genes.results>;print "./merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "}'>merge.sh  
sh merge.sh >total_species.gene.matrix  
  
#######################################################  
## these reads number matrixes could be used in DEGs detection by DESeq2   
#######################################################  

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
