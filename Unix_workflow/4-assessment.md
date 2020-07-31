Assessment: _de novo_ assemblies, ORFs and Orthlogous gene sets of all species   
======================================
###################  
## 1 Busco assessment   
###################  
mkdir Busco_assessment  
cd Busco_assessment  
### 1.1 on _de novo_ assemblies of all six species (./orthologue)  
#### Acura_tra_nuc.fa, Daru_tra_nuc.fa, Ocomp_tra_nuc.fa, Padel_tra_nuc.fa, Pmol_tra_nuc.fa
mkdir denovo_assembl  
cd denovo_assembl
nohup python ~/software/Busco/scripts/run_BUSCO.py -m tran -i Acura_tra_nuc.fa -o Acura -l ~/software/Busco/lineage/actinopterygii_odb9 -c 12 -t Acura-tmp >Busco_acura.process 2>&1 &  
nohup python ~/software/Busco/scripts/run_BUSCO.py -m tran -i Daru_tra_nuc.fa -o Daru -l ~/software/Busco/lineage/actinopterygii_odb9 -c 12 -t Daru-tmp >Busco_daru.process 2>&1 &  
nohup python ~/software/Busco/scripts/run_BUSCO.py -m tran -i Ocomp_tra_nuc.fa -o Ocomp -l ~/software/Busco/lineage/actinopterygii_odb9 -c 12 -t Ocomp-tmp >Busco_ocomp.process 2>&1 &  
nohup python ~/software/Busco/scripts/run_BUSCO.py -m tran -i Padel_tra_nuc.fa -o Padel -l ~/software/Busco/lineage/actinopterygii_odb9 -c 12 -t Padel-tmp >Busco_Padel.process 2>&1 &  
nohup python ~/software/Busco/scripts/run_BUSCO.py -m tran -i Pmol_tra_nuc.fa -o Pmol -l ~/software/Busco/lineage/actinopterygii_odb9 -c 12 -t Pmol-tmp >Busco_Pmol.process 2>&1 &  

### 1.2 On the ORFs after transdecoder  
#### Acura_nuc.fa, Daru_nuc.fa, Ocomp_nuc.fa, Padel_nuc.fa, Pmol_nuc.fa
mkdir ORF_nuc_transdecoder  
cd ORF_nuc_transdecoder  
nohup python ~/software/Busco/scripts/run_BUSCO.py -m tran -i Acura_nuc.fa -o Acura -l ~/software/Busco/lineage/actinopterygii_odb9 -c 12 -t Acura-tmp >Busco_acura.process 2>&1 &  
nohup python ~/software/Busco/scripts/run_BUSCO.py -m tran -i Daru_nuc.fa -o Daru -l ~/software/Busco/lineage/actinopterygii_odb9 -c 12 -t Daru-tmp >Busco_daru.process 2>&1 &  
nohup python ~/software/Busco/scripts/run_BUSCO.py -m tran -i Ocomp_nuc.fa -o Ocomp -l ~/software/Busco/lineage/actinopterygii_odb9 -c 12 -t Ocomp-tmp >Busco_ocomp.process 2>&1 &  
nohup python ~/software/Busco/scripts/run_BUSCO.py -m tran -i Padel_nuc.fa -o Padel -l ~/software/Busco/lineage/actinopterygii_odb9 -c 12 -t Padel-tmp >Busco_padel.process 2>&1 &  
nohup python ~/software/Busco/scripts/run_BUSCO.py -m tran -i Pmol_nuc.fa -o Pmol -l ~/software/Busco/lineage/actinopterygii_odb9 -c 12 -t Pmol-tmp >Busco_pmol.process 2>&1 &  
  
### 1.3 On Orthologous gene sets  
#### Acura.fa, Daru.fa, Ocomp.fa, Padel.fa, Pmol.fa, Apoly.fa in ./final_reference
cd ./final_reference
nohup python ~/software/Busco/scripts/run_BUSCO.py -m tran -i Acura.fa -o Acura -l ~/software/Busco/lineage/actinopterygii_odb9 -c 12 -t Acura-tmp >Busco_acura.process 2>&1 &  
nohup python ~/software/Busco/scripts/run_BUSCO.py -m tran -i Daru.fa -o Daru -l ~/software/Busco/lineage/actinopterygii_odb9 -c 12 -t Daru-tmp >Busco_daru.process 2>&1 &  
nohup python ~/software/Busco/scripts/run_BUSCO.py -m tran -i Ocomp.fa -o Ocomp -l ~/software/Busco/lineage/actinopterygii_odb9 -c 12 -t Ocomp-tmp >Busco_ocomp.process 2>&1 &  
nohup python ~/software/Busco/scripts/run_BUSCO.py -m tran -i Padel.fa -o Padel -l ~/software/Busco/lineage/actinopterygii_odb9 -c 12 -t Padel-tmp >Busco_padel.process 2>&1 &  
nohup python ~/software/Busco/scripts/run_BUSCO.py -m tran -i Pmol.fa -o Pmol -l ~/software/Busco/lineage/actinopterygii_odb9 -c 12 -t Pmol-tmp >Busco_pmol.process 2>&1 &  
nohup python ~/software/Busco/scripts/run_BUSCO.py -m tran -i Apoly.fa -o Apoly -l ~/software/Busco/lineage/actinopterygii_odb9 -c 12 -t Apoly-tmp >Busco_apoly.process 2>&1 &  
  
#######################  
## 2 transrate assessment
#######################  
  
### 2.1 on _de novo_ assemblies of all six species (in the Drap_trinity/e-rmbt_editing of each species)  
nohup transrate --assembly all_contigs.second_pass.fa --left Acura_total_1.fastq.gz --right Acura_total_2.fastq.gz --threads 12 >Acura_transrate.process 2>1 &  
nohup transrate --assembly all_contigs.second_pass.fa --left Ocomp_total_1.fastq.gz --right Ocomp_total_2.fastq.gz --threads 12 >Ocomp_transrate.process 2>1 &  
nohup transrate --assembly all_contigs.second_pass.fa --left Daru_total_1.fastq.gz --right Daru_total_2.fastq.gz --threads 12 >Daru_transrate.process 2>1 &  
nohup transrate --assembly all_contigs.second_pass.fa --left Padel_total_1.fastq.gz --right Padel_total_2.fastq.gz --threads 12 >Padel_transrate.process 2>1 &  
nohup transrate --assembly all_contigs.second_pass.fa --left Pmol_total_1.fastq.gz --right Pmol_total_2.fastq.gz --threads 12 >Pmol_transrate.process 2>1 &  

### 2.2 On the ORF after transdecoder  
#### Acura_nuc.fa, Daru_nuc.fa, Ocomp_nuc.fa, Padel_nuc.fa, Pmol_nuc.fa  
nohup transrate --assembly Acura_nuc.fa --left Acura_total_1.fastq.gz --right Acura_total_2.fastq.gz --threads 12 >Acura_transrate.process 2>1 &  
nohup transrate --assembly Daru_nuc.fa --left Ocomp_total_1.fastq.gz --right Ocomp_total_2.fastq.gz --threads 12 >Ocomp_transrate.process 2>1 &  
nohup transrate --assembly Ocomp_nuc.fa --left Daru_total_1.fastq.gz --right Daru_total_2.fastq.gz --threads 12 >Daru_transrate.process 2>1 &  
nohup transrate --assembly Padel_nuc.fa --left Padel_total_1.fastq.gz --right Padel_total_2.fastq.gz --threads 12 >Padel_transrate.process 2>1 &  
nohup transrate --assembly Pmol_nuc.fa --left Pmol_total_1.fastq.gz --right Pmol_total_2.fastq.gz --threads 12 >Pmol_transrate.process 2>1 & 
  
### 2.3 On Orthologous gene set  
#### Acura.fa, Daru.fa, Ocomp.fa, Padel.fa, Pmol.fa, Apoly.fa in ./final_reference
nohup transrate --assembly Acura.fa --left Acura_total_1.fastq.gz --right Acura_total_2.fastq.gz --threads 12 >Acura_transrate.process 2>1 &  
nohup transrate --assembly Daru.fa --left Ocomp_total_1.fastq.gz --right Ocomp_total_2.fastq.gz --threads 12 >Ocomp_transrate.process 2>1 &  
nohup transrate --assembly Ocomp.fa --left Daru_total_1.fastq.gz --right Daru_total_2.fastq.gz --threads 12 >Daru_transrate.process 2>1 &  
nohup transrate --assembly Padel.fa --left Padel_total_1.fastq.gz --right Padel_total_2.fastq.gz --threads 12 >Padel_transrate.process 2>1 &  
nohup transrate --assembly Pmol.fa --left Pmol_total_1.fastq.gz --right Pmol_total_2.fastq.gz --threads 12 >Pmol_transrate.process 2>1 &  
