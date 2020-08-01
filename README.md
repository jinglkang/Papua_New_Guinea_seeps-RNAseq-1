# Papue_New_Guinea-RNA-seq    
  
## 1. Quality Control    
1-Quality_Control.md    
  
## 2. _de novo_    
2-denovo.md    
  
## 3. Orthologous gene detection    
3-ortholog_DEGs_detection.md  
### 3.1 blastp_uni.pl  
#### after the running of OrthoFinder, select Orthogroup that contain at least one reads per species; and blastp to swiss-prot database, Orthogroups/Orthogroups.GeneCount.tsv is OrthoFinder result file; the qualified orthogroup will be in "./orthogroup_needed_sequences"; the blastp result will be in "blastp_result".    
Example: perl blastp_uni.pl Orthogroups/Orthogroups.GeneCount.tsv    
  
### 3.2 perl get_best_blast_orthogroup.pl  
#### select the most representative transcript for each species in each orthologous group. "./six_fishes_reads_num" includes reads number matrix of all species, which are used to select the trancript with one mapped reads. the final orthologous group will be in "./final_blast_orth_group".    
Example: perl get_best_blast_orthogroup.pl -blast_result=blastp_result -reads_matrix=six_fishes_reads_num    
  
### 3.3 get_sequences_ref.pl  
#### obtain the new reference per species by concatenated the final orthologous gene set. The input is the ortholog group in "final_blast_orth_group" and the nucleotide sequence of de novo sequences of all species (the fasta files are in "./orthofinder_input_nuc"), it will concatenate transcript according to the orthologous group of each species, the references of each species will be in "./final_reference".      
Example:./get_sequences_ref.pl -input=final_blast_orth_group -nuc=./orthofinder_input_nuc -output=final_reference    
  
### 3.4 merge_RSEM_frag_counts_single_table.pl  
#### revised merge_RSEM_frag_counts_single_table.pl from RSEM, remove the suffix of the sample name and get the integer reads number  
Example: ./merge_RSEM_frag_counts_single_table.pl RSEM_result_file      
  
## 4. Assessment  
4-assessment.md  
  
## 5. SNPs calling  
5-SNPs_calling.md  
#### perl gatk_rna.pl: --fasta <reference fasta file> --tmp <buffer memory space of GATK> --out <result_directory>  
Example: perl gatk_rna.pl --fasta ./Acura.fa --bam . --tmp tmp --output Acura_gatk >Acura_gatk.process    
#### transform vcf file to SNP genetype file, which could be used in BayeScan    
#### need provide the information of sample in which population    
Example: perl vcf_bayescan.pl --vcf Acomp.all.snp.final.passed-1.vcf --pop_def sample_def.txt --id_correlation_num  loci-numb-id.bayescan >Acomp_bayecan.input    

## 6. EVE model analysis
6-EVE_model_analysis.md
#### get_sequence_of_orthogroup.pl: get the protein sequences of each transcript from the ORFs fasta files of each species and do the alignment by clustalo  
#### -output "direcotry that store orthologous groups (fasta files contain the sequences of each transcript per species)"  
Example: perl get_sequence_of_orthogroup.pl -fasta=all.fasta -blast_result=final_blast_orth_group -output=final_orth_sequences  

#### perl glocks_process.pl: -input "alignment of fasta files" -Gblocks_output "Gblocks output directory" -long_seq "Gblocks output direcotry that each fasta with more than 50 codons"  
Exaple: perl glocks_process.pl -input final_orth_sequences -gblocks_output Gblocks_output -long_seq Gblocks_output_long_than_50  

#### fasta2phy.pl: transformed fasta file to phylip format  
perl fasta2phy.pl concatenate.fas >concatenate.phy  
