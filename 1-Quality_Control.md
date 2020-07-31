Quality Control 
===============  
#######################
## 1. FastQc (first time)    
#######################
### run FastQC
### fastq files of all samples are in the currently directory  
fastqc \*.fastq -o ./fastqc --extract -t 32  
### extract the Overrepresented sequences into Overrepresented_sequences.txt   
find -name "fastqc_data.txt" >11  
less 11|perl -ne 'chomp;print "$_\t"'|perl -alne 'print "cat $_ > all.samples.fastqc_data.txt"' > 1.sh  
sh 1.sh  
less all.samples.fastqc_data.txt|perl -alne 'print if /^[ATCG]{8,}/'|cut -f 1|sort -u > Overrepresented_sequences.txt  
### generate TruSeq2-PE-new.fa into the adapter file of Trimmomatic, add overrepresent sequences into TruSeq2-PE-new.fa  
less Overrepresented_sequences.txt |perl -alne '$i++;$num=$i+2;print ">FlowCell$num\n$\_"' > Overrepresented_sequences.fa  
### cat the index file (TruSeq2-PE.fa) in the install directory of Trimmomatic and Overrepresented_sequences.fa together
### and put TruSeq2-PE-new.fa and trimmomatic-0.30.jar in the current directory ./
cat TruSeq2-PE.fa Overrepresented_sequences.fa > TruSeq2-PE-new.fa  

#################################  
## 2. Trimmomatic  
#################################  
mkdir paired  
mkdir unpaired  
ll *.gz|perl -alne '@a=split /\./,$F[-1];print"$a[0]"'|perl -alne 's/\_R[1|2]//;print'|sort -u >sample.name  
vi temp2.pl  
open fil, "sample.name";  
while (<fil>) {  

        chomp;  
        $forward="_R1.fastq.gz";  
        $reverse="_R2.fastq.gz";  
        $seq1=$_.$forward;  
        $seq2=$_.$reverse;  
        $seq1_paired=$_.".paired".$forward;  
        $seq2_paired=$_.".paired".$reverse;  
        $seq1_unpaired=$_.".unpaired".$forward;  
        $seq2_unpaired=$_.".unpaired".$forward;  
        print "java -jar trimmomatic-0.39.jar PE $seq1 $seq2 ./paired/$seq1_paired ./unpaired/$seq1_unpaired ./paired/$seq2_paired ./unpaired/$seq2_unpaired ILLUMINACLIP:TruSeq2-PE-new.fa:2:30:10 LEADING:4 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:40 -threads 32\n";  
}  
perl temp2.pl >1.sh  
# split 1.sh into 4 files  
split -l 33 1.sh trim  # "trim" as the prefix  
vi temp3.pl  
open fil, "$ARGV[0]";  
while (<fil>) {  

        chomp;  
        `$_`;  
}  
nohup perl temp3.pl trimaa &  
nohup perl temp3.pl trimab &  
nohup perl temp3.pl trimac &  
nohup perl temp3.pl trimad &  
#### extract Trimmomatic process results  
vi temp5.pl  
open fil,"nohup.out";  
while (<fil>) {  

        chomp;  
        if (/(\D+\d+_\D\d+)_R1\.fastq.gz/){  
                $name=$1;  
                print "$name\t";  
        }  
        if (/Input/){  
                @num=$_=~/\d+/g;  
                $per1=$num[2].".".$num[3];  
                $per2=$num[5].".".$num[6];  
                $per3=$num[8].".".$num[9];  
                $per4=$num[11].".".$num[12];  
                print "$num[0]\t$num[1]\t$per1\t$num[4]\t$per2\t$num[7]\t$per3\t$num[10]\t$per4\n";  
        }  
}  
perl temp5.pl > Trimmomatic.result  

######################  
## 3. FastQc (Second time)  
######################  
cd paired
mkdir fastqc  
fastqc *.fastq -o ./fastqc --extract -t 32  
  
### compare the differences between fastqc result files (summary.txt) of two FastQc process  
#### ./paired (FastQc second time)
cd ./paired
find -name "summary.txt"|perl -ne 'chomp;print "$_\t"'|perl -alne 'print "cat $_ > all.samples.fastqc.txt"' >1.sh  
sh 1.sh  
perl temp1.pl >all.samples.second.fastqc.result  
#### ./ (FastQc first time)
find -name "summary.txt"|perl -ne 'chomp;print "$_\t"'|perl -alne 'print "cat $_ > all.samples.fastqc.txt"' >1.sh  
sh 1.sh  
perl temp1.pl >all.samples.first.fastqc.result  
  
#######################################  
## 4. Kraken to move potential contaminates
#######################################  

### construct Kraken library (archaea, bacteria, fungi, viral)  in Kraken install directory 
mkdir library  
cd library  
#### download archaea, bacteria, fungi, viral  
kraken2-build --download-library archaea --db archaea  
kraken2-build --download-library bacteria --db bacteria  
kraken2-build --download-library fungi --db fungi  
kraken2-build --download-library viral --db viral  
### run Kraken
mkdir kraken  
cd kraken  
mkdir reports   ## kraken reports  
mkdir total_output ## kraken output file  
#### working in ./paired  
for file in *.paired_R1.fastq.gz; do name=${file/.paired_R1.fastq.gz};kraken2 --db \~/software/kraken2/library --paired --threads 12 --gzip-compressed --unclassified-out kraken/${name}#.fastq ${name}.paired_R1.fastq.gz ${name}.paired_R2.fastq.gz --report kraken/reports/${name}.kraken_report --use-name --confidence 0.7 > kraken/total_output/${name}.out;done  
  
## kraken result  
## get the dropped reads number and the drapped reads ratio in total reads  
cd ./kraken/reports  
vi temp1.pl  
@reports=<*report>;  
foreach $report(@reports){  

	open fil,"$report";  
	while (<fil>){  
		chomp;  
		$name=$report=~/(.*)\.kraken_report/;  
		next if /unclassified/;  
		next if /0\.00/;  
		if (/Archaea/) {  
			print "$_\t$name\n";  
		}  
		if (/Viruses/) {  
			print "$_\t$name\n";  
		}  
		if (/Bacteria$/) {  
			print "$_\t$name\n";  
		}  
		if (/fungi$/i) {  
			print "$_\t$name\n";  
		}  
	}  
}  
perl temp1.pl >kraken.contaminate.ratio  
