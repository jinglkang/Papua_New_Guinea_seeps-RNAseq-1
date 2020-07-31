GATK SNPs calling    
================  
## The concatenated sequences of orthologous genes as the references in the SNPs calling    
### gatk_rna.pl to do the SNPs calling and initial filtering    
### Acura.fa, Daru.fa, Ocomp.fa, Padel.fa, Pmol.fa, Apoly.fa in ./final_reference    
## _A. curacao_    
mkdir Acura_gatk  # the alignment bam file of Acura should be moved here    
cd Acura_gatk    
### run perl gatk_rna.pl, the resulted final vcf file (Acura.all.snp.final.vcf) would be in ./Acura_gatk    
nohup perl gatk_rna.pl --fasta ./Acura.fa --bam . --tmp tmp --output Acura_gatk >Acura_gatk.process 2>&1 &  
cd ./Acura_gatk
vi temp1.pl
#### vcf SNPs file：keep loci with "PASS" flag, maximun heterozygosity <0.8 and minor allele frequency >=0.1    
use List::Util qw/max min/;  
open fil, "$ARGV[0]";  
while (<fil>) {  

        %hash=();  
        %het=();  
        $het=0;  
        @allel_freq=();  
        @het_ratio=();  
        chomp;  
        if (/label/i) {  
                print "$_\n";  
                $count = $_ =~ s/Label//g;  
        }elsif (/^#/) {  
                print "$_\n";  
                } elsif (! /^#/) {  
                if (/pass/i) {  
                        @a=split;  
                        for ($i = -$count; $i <= -1; $i++) {  
                                @c=split /\:/, $a[$i];  
                                @d=split /\//, $c[0];  
                                if ($d[0] != $d[1]) {  
                                        $het{$c[0]}++;  
                                }  
                                foreach $d (@d) {  
                                        $hash{$d}++;  
                                }  
                        }  
                        foreach $key1 (keys %het) {  
                                $het_ratio=$het{$key1}/$count;  
                                push @het_ratio, $het_ratio;  
                        }  
                        foreach $key (keys %hash) {  
                                $allel_freq=$hash{$key}/(2*$count);  
                                push @allel_freq, $allel_freq;  
                        }  
                        $MAF=min(@allel_freq);  
                        $het_max=max(@het_ratio);  
                        @hash=keys %hash;  
                        $num=@hash;  
                        #       print "$num\n";  
                        print "$_\n" if $num>1 && $het_max<0.8 && $MAF >=0.1;  
                        #                        print "$num\t$het_max\t$MAF\n";  
                }  
        }  
}  
perl temp1.pl Acura.all.snp.final.vcf >Acura.all.snp.final.passed-1.vcf  
grep -v '^#' Acura.all.snp.final.passed-1.vcf|wc -l    
26 # 26 SNPs in Acura  

############
### Bayescan
############
#### transform vcf file to genetype file, which would be the input of bayescan  
#### sample_def.txt is the infromation where the sample come from, the output file (loci-numb-id.bayescan) tells the correslation bettwen the SNPs loci and the id that used in the Bayescan  
#### change the information of population in the sample_def.txt to set the population num  
perl vcf_bayescan.pl --vcf Acura.all.snp.final.passed-1.vcf --pop_def sample_def.txt --id_correlation_num loci-numb-id.bayescan >Acura_bayecan.input  
BayeScan -snp Acura_bayecan.input  
    
#############
### admixture
#############
plink --vcf Acura.all.snp.final.passed.vcf --recode --out my_data --allow-extra-chr 0 --make-bed  
for K in 1 2 3 4 5 6 7 8 9 10; do admixture --cv my_data.bed $K | tee log${K}.out; done  
grep -h CV log*.out  
CV error (K=10): 0.89330  
CV error (K=1): 0.36177  
CV error (K=2): 0.49057  
CV error (K=3): 0.63131  
CV error (K=4): 0.61580  
CV error (K=5): 0.87449  
CV error (K=6): 0.90621  
CV error (K=7): 0.94013  
CV error (K=8): 0.78433  
CV error (K=9): 0.68934  
  
## _A. polyacanthus_  
mkdir Apoly_gatk  
cd Apoly_gatk    
nohup perl gatk_rna.pl --fasta ./Apoly.fa --bam . --tmp tmp --output Apoly_gatk >Apoly_gatk.process 2>&1 &  
cd ./Apoly_gatk  
#### vcf SNPs file：keep loci with "PASS" flag, maximun heterozygosity <=0.8 and minor allele frequency >=0.1    
perl temp1.pl Apoly.all.snp.final.vcf >Apoly.all.snp.final.passed-1.vcf     
grep -v '^#' Apoly.all.snp.final.passed-1.vcf|wc -l    
1035 # 1035 SNPs in Apoly  
  
############
### Bayescan
############
perl vcf_bayescan.pl --vcf Apoly.all.snp.final.passed-1.vcf --pop_def sample_def.txt --id_correlation_num loci-numb-id.bayescan >Apoly_bayecan.input  
BayeScan -snp Apoly_bayecan.input  
  
#############
### admixture
#############
plink --vcf Apoly.all.snp.final.passed-1.vcf --recode --out my_data --allow-extra-chr 0 --make-bed  
for K in 1 2 3 4 5 6 7 8 9 10; do admixture --cv my_data.bed $K | tee log${K}.out; done  
grep -h CV log*.out  
CV error (K=10): 1.27078  
CV error (K=1): 0.44976  
CV error (K=2): 0.59192  
CV error (K=3): 0.69454  
CV error (K=4): 0.80680  
CV error (K=5): 1.04558  
CV error (K=6): 1.11638  
CV error (K=7): 1.32218  
CV error (K=8): 1.51571  
CV error (K=9): 1.34219  
  


############################################  
############################################  
## _A. compressus_    
mkdir ./Acomp_gatk  
nohup perl gatk_rna.pl --fasta ./Ocomp.fa --bam . --tmp tmp --output Acomp_gatk >Acomp_gatk.process 2>&1 &  
cd ./Acomp_gatk  
#### vcf SNPs file：keep loci with "PASS" flag, maximun heterozygosity <=0.8 and minor allele frequency >=0.1    
perl temp1.pl Acomp.all.snp.final.vcf >Acomp.all.snp.final.passed-1.vcf     
grep -v '^#' Acomp.all.snp.final.passed-1.vcf|wc -l    
452 # 452 SNPs in Acomp  

############
### Bayescan
############
perl vcf_bayescan.pl --vcf Acomp.all.snp.final.passed-1.vcf --pop_def sample_def.txt --id_correlation_num loci-numb-id.bayescan >Acomp_bayecan.input  
BayeScan -snp Acomp_bayecan.input  

#############
### admixture
#############
plink --vcf Acomp.all.snp.final.passed-1.vcf --recode --out my_data --allow-extra-chr 0 --make-bed  

# admixture  
grep -h CV log*.out  
CV error (K=10): 1.91282  
CV error (K=1): 0.49162  
CV error (K=2): 0.61411  
CV error (K=3): 0.75024  
CV error (K=4): 1.08110  
CV error (K=5): 1.30350  
CV error (K=6): 1.61921  
CV error (K=7): 1.50137  
CV error (K=8): 1.38740  
CV error (K=9): 1.01894  
  

## _D. aruanus_    
mkdir Daru_gatk  
nohup perl gatk_rna.pl --fasta ./Daru.fa --bam . --tmp tmp --output Daru_gatk >Daru_gatk.process 2>&1 &  
cd ./Daru_gatk  
#### vcf SNPs file：keep loci with "PASS" flag, maximun heterozygosity <=0.8 and minor allele frequency >=0.1    
perl temp1.pl Daru.all.snp.final.vcf >Daru.all.snp.final.passed-1.vcf     
grep -v '^#' Daru.all.snp.final.passed-1.vcf|wc -l    
2105 # 2105 SNPs in Daru  

############
### Bayescan
############
perl vcf_bayescan.pl --vcf Daru.all.snp.final-1.vcf --pop_def pop_def.txt --id_correlation_num Daru_loci-numb-id.bayescan >Daru_bayecan.input  
BayeScan -snp Daru_bayecan.input  

#############
### admixture
#############
plink --vcf Daru.all.snp.final.passed-1.vcf --recode --out my_data --allow-extra-chr 0 --make-bed  
grep -h CV log*.out  
CV error (K=10): 1.91880  
CV error (K=1): 0.48868  
CV error (K=2): 0.65945  
CV error (K=3): 0.85918  
CV error (K=4): 1.05459  
CV error (K=5): 1.29056  
CV error (K=6): 1.46555  
CV error (K=7): 1.68151  
CV error (K=8): 1.48244  
CV error (K=9): 2.15463  
  
## _P. adelus_  
mkdir Padel_gatk  
nohup perl gatk_rna.pl --fasta ./Padel.fa --bam . --tmp tmp --output Padel_gatk >Padel_gatk.process 2>&1 &  
cd ./Padel_gatk  
#### vcf SNPs file：keep loci with "PASS" flag, maximun heterozygosity <=0.8 and minor allele frequency >=0.1    
perl temp1.pl Padel.all.snp.final.vcf >Padel.all.snp.final.passed-1.vcf     
grep -v '^#' Padel.all.snp.final.passed-1.vcf|wc -l    
52 # 52 SNPs in Padel  

############
### Bayescan
############
perl vcf_bayescan.pl --vcf Padel.all.snp.final-1.vcf --pop_def pop_def.txt --id_correlation_num Pmol_loci-numb-id.bayescan >Padel_bayecan.input  
BayeScan -snp Padel_bayecan.input  

#############
### admixture
#############
plink --vcf Padel.all.snp.final.passed-1.vcf --recode --out my_data --allow-extra-chr 0 --make-bed  
grep -h CV log*.out  
CV error (K=10): 0.73828  
CV error (K=1): 0.39519  
CV error (K=2): 0.57515  
CV error (K=3): 0.60623  
CV error (K=4): 0.57874  
CV error (K=5): 0.78310  
CV error (K=6): 0.63551  
CV error (K=7): 0.76925  
CV error (K=8): 0.89937  
CV error (K=9): 0.76462  
  
## _P. molluscensis_  
mkdir Pmol_gatk  
nohup perl gatk_rna.pl --fasta ./Pmol.fa --bam . --tmp tmp --output Pmol_gatk >Pmol_gatk.process 2>&1 &  
cd ./Pmol_gatk  
#### vcf SNPs file：keep loci with "PASS" flag, maximun heterozygosity <=0.8 and minor allele frequency >=0.1    
perl temp1.pl Pmol.all.snp.final.vcf >Pmol.all.snp.final.passed-1.vcf     
grep -v '^#' Pmol.all.snp.final.passed-1.vcf|wc -l    
63 # 63 SNPs in Pmol  

############
### Bayescan
############
perl vcf_bayescan.pl --vcf Pmol.all.snp.final-1.vcf --pop_def pop_def.txt --id_correlation_num Pmol_loci-numb-id.bayescan >Pmol_bayecan.input  
BayeScan -snp Pmol_bayecan.input  

#############
### admixture
#############
plink --vcf Pmol.all.snp.final.passed-1.vcf --recode --out my_data --allow-extra-chr 0 --make-bed  
grep -h CV log*.out  
CV error (K=10): 0.83266  
CV error (K=1): 0.43347  
CV error (K=2): 0.65058  
CV error (K=3): 0.80360  
CV error (K=4): 0.93647  
CV error (K=5): 0.82158  
CV error (K=6): 0.80991  
CV error (K=7): 0.74335  
CV error (K=8): 0.97828  
CV error (K=9): 0.79244  
