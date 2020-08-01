EVE model analysis  
==================  
#### notice: must have "./results" in the running directory of EVE model  
mkdir results  

#########################################################  
## 1 Phylogenetic tree construction for EVE model analysis
#########################################################  
### 1.1 preproccess in each orthogroup before Phylogeny  
#### 1.1.1 get_sequence_of_orthogroup.pl: get the protein sequences of each transcript from the ORFs fasta files of each species and do the alignment by clustalo  
##### -fasta "ORFs sequences in all sepcies before orthologous gene detecion"  
##### -blast_result "the result directory that store your final orthologous groups"  
##### -output "direcotry that store orthologous groups (fasta files contain the sequences of each transcript per species)"  
cat Acura_nuc.fa Apoly.fa Daru.fa Ocomp.fa Padel.fa Pmol.fa > all.fasta  
nohup perl get_sequence_of_orthogroup.pl -fasta=all.fasta -blast_result=final_blast_orth_group -output=final_orth_sequences >get_sequences.process 2>&1 &  
#### 1.1.2 remove the gaps within the alignment by Gblocks, and keep the orthologous group with more than 50 codons  
#### perl glocks_process.pl: -input "alignment of fasta files" -Gblocks_output "Gblocks output directory" -long_seq "Gblocks output direcotry that each fasta with more than 50 codons"  
nohup perl glocks_process.pl -input final_orth_sequences -gblocks_output Gblocks_output -long_seq Gblocks_output_long_than_50 >gblocks.process 2>&1 &  
#### 1.1.3 concatenate all sequences in ./Gblocks_output_long_than_50 according orthogroup id by species  
vi temp1.pl  
##### concatenated all orthogroups (condons > 50)  
@fas=<Gblocks_output_long_than_50/*.fas>;  
foreach $fas (@fas) {  

	open fil, $fas;  
	while (<fil>) {  
		chomp;  
		if (/>/) {  
			s/>//;  
			($spe)=$_=~/(.*)\_\d+/;  
		} else {  
			$hash{$spe}.=$_;  
		}  
	}  
}  
foreach $key (sort keys %hash) {  

	print ">$key\n$hash{$key}\n";  
}  
perl temp1.pl >concatenate.fas  
  
### 1.2 RAxML-NG to construct the maximum likelihood (ML) tree  
mkdir RAxML-NG  
mv concatenate.fas RAxML-NG  
#### transformed fasta file to phylip  
perl fasta2phy.pl concatenate.fas >concatenate.phy  
#### first check that the MSA can actually be read and doesn't contain sites with only undetermined characters or sequences with undetermined characters or duplicate taxon names  
#### model: LG+G4  
raxml-ng --parse --msa concatenate.phy --model LG+G4 --prefix T1    
#### Alignment can be successfully read by RAxML-NG  
Gaps: 0.00 %  
Invariant sites: 83.61 %  
* Estimated memory requirements                : 392 MB  
  
* Recommended number of threads / MPI processes: 51  
  
Please note that numbers given above are rough estimates only.  
Actual memory consumption and parallel performance on your system may differ!  
  
Alignment can be successfully read by RAxML-NG.  
  
#### start to estimate best ML tree  
raxml-ng --bootstrap --msa T1.raxml.rba --model LG+G4 --prefix T2 --threads 32 --seed 2 --tree pars{50},rand{50} --bs-cutoff 0.01 >raxml-ng.process 2>&1 &  
raxml-ng --bsconverge --bs-trees T2.raxml.bootstraps --prefix T3 --seed 2 --threads 2 --bs-cutoff 0.01  
Performing bootstrap convergence assessment using autoMRE criterion  
    trees        avg WRF       avg WRF in %       # perms: wrf <= 1.00 %     converged?  
      50          0.000              0.000                         1000        YES  
Bootstopping test converged after 50 trees  
  
##### get the best tree
raxml-ng --msa T1.raxml.rba --model LG+G4 --prefix T5 --threads 32 --seed 2 # T5.raxml.bestTree  
raxml-ng --support --tree T5.raxml.bestTree --bs-trees T2.raxml.bootstraps --prefix T6 --threads 32 --bs-metric tbe	# T6.raxml.support  
##### the T6.raxml.support would be the input tree in the EVE model analysis  
  
#############################  
## 2. start EVE model analysis  
#############################  
mkdir EVE_model_analysis  
cp RAxML-NG/T6.raxml.support EVE_model_analysis  
mv T6.raxml.support 1_Tree.newick # first line: species number (1 to N to represents the species from left to right); second line: phylogeny tree (Notice: remove bootstraps in 1_Tree.newick)  
vi 2_Nindivs.indiv # individual numbers (from 1 to N to represents the species, and the second line is the corresponding sample numbers per species)  
vi 3_sampleExpr.dat # which is transformed by normalized reads number matrix from DESeq2  

#### only samples from CO2 seep  
~/software/EVE_release/EVEmodel -S -n 14634 -t 1_Tree.newick -i 2_Nindivs.indiv -d 3_sampleExpr.dat -f \_CO2 -v 10 >eve.process  
287 genes have higher variance among than within lineages at 0.05 FDR.   
8 genes have higher variance within than among lineages at 0.05 FDR.  

#### only samples from Control site  
nohup ~/software/EVE_release/EVEmodel -S -n 14634 -t 1_Tree.newick -i 2_Nindivs.indiv -d 3_sampleExpr.dat -f \_Control -v 10 >eve.process  

#### All samples
~/software/EVE_release/EVEmodel -S -n 14634 -t 1_Tree.newick -i 2_Nindivs.indiv -d 3_sampleExpr.dat -f _total -v 10  


#### all sample excluding _A. polyacanthus_ from mid  
$~/software/EVE_release/EVEmodel -S -n 14634 -t 1_Tree.newick -i 2_Nindivs.indiv -d 3_sampleExpr_exclu_mid.dat -f _total_exclu_mid -v 10  
