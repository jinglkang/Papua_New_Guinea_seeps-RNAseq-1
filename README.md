# Papue_New_Guinea-RNA-seq

## 1. after the running of OrthoFinder, select Orthogroup that contain at least one reads per species; and blastp to swiss-prot database, Orthogroups/Orthogroups.GeneCount.tsv is OrthoFinder result file; the qualified orthogroup will be in "./orthogroup_needed_sequences"; the blastp result will be in "blastp_result".  
perl blastp_uni.pl Orthogroups/Orthogroups.GeneCount.tsv  
## 2. select the most representative transcript for each species in each orthologous group. "./six_fishes_reads_num" includes reads number matrix of all species, which are used to select the trancript with one mapped reads. the final orthologous group will be in "./final_blast_orth_group".  
perl get_best_blast_orthogroup.pl -blast_result=blastp_result -reads_matrix=six_fishes_reads_num  
## 3. obtain the new reference per species by concatenated the final orthologous gene set. The input is the ortholog group in "final_blast_orth_group" and the nucleotide sequence of de novo sequences of all species (the fasta files are in "~/Desktop/PapueNewGuinea-new/orthologue/orthofinder_input_nuc"), it will concatenate transcript according to the orthologous group of each species, the references of each species will be in "./final_reference".    
./get_sequences_ref.pl -input=final_blast_orth_group -nuc=~/Desktop/PapueNewGuinea-new/orthologue/orthofinder_input_nuc -output=final_reference
