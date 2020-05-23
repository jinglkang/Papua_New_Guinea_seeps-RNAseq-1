# GATK work flow
# the reference fasta must in the current directory
#!/usr/bin/perl -w
use File::Basename;
use Getopt::Long;
my $opt = GetOptions( 'fasta:s', \$fasta,
        'bam:s', \$bam_dir,
    'output:s', \$output);
my $help;
if (!($opt && $fasta && $output && $bam_dir) || $help) {
    print STDERR "\nExample usage:\n";
    print STDERR "\n$0 -fasta=\"reference sequence\" -bam=\"bam files directory\" -output=\"output directory\" \n\n";
    exit;
}
unless (-e $output) {
        `mkdir $output`;
}
my ($spec)=basename($fasta)=~/(.*)\.fa/;
my @bam=<$bam_dir/*.transcript.bam>;
my $PICARD = "/home/Kang/software/picard/picard.jar";
my $GATK = "/home/Kang/software/gatk3.8/GenomeAnalysisTK.jar";
my $GATKREF = $fasta;
my $DIC = $spec.".dict";
my $fai=$fasta.".fai";

unless (-e $DIC && -e $fai) {
        `java -jar $PICARD CreateSequenceDictionary R=$GATKREF O=$DIC`; # build index of fasta file
        `samtools faidx $fasta`;
}

foreach (@bam) {
        my ($name)=basename($_)=~/(.*)\.transcript\.bam/;
        my $input_1=$_; # sortsam input
        my $bam_sort=$name.".sorted.bam";
        my $bam_sort_add=$name.".sort.add.bam"; # add read group
        my $bam_sort_add_dup=$name.".sort.add.dup.bam"; # mark duplicated reads
        my $dump="name".".mdup.metrics";
        my $resort=$name.".resorted.bam"; # resorted
        my $split=$name.".split.bam"; # SplitNCigarReads
        my $vcf=$name.".vcf"; my $vcf_idx=$vcf.".idx"; # GATK Haplotype
        my $snp_vcf=$name.".snp.vcf"; my $snp_vcf_idx=$snp_vcf.".idx";
        if (-e $snp_vcf) {
                print "congrat, $name already has snp file! let's start the next one.\n";
                `mv $vcf $vcf_idx $snp_vcf $snp_vcf_idx $output`;
                `rm $bam_sort $bam_sort_add $bam_sort_add_dup $resort`;
                next;
        } else {
                `java -jar $PICARD SortSam INPUT=$input_1 OUTPUT=$bam_sort SO=coordinate`; # sortsam
                if (-e $bam_sort) {
                        `java -jar $PICARD AddOrReplaceReadGroups INPUT=$bam_sort OUTPUT=$bam_sort_add SO=coordinate RGID=Label RGLB=Label RGPL=illumina RGPU=Label RGSM=Label`; # add read group
                } else {
                        die "There is no $bam_sort, error occurs in the first step: SortSam\n";
                }
                if (-e $bam_sort_add) {
                        `java -jar $PICARD MarkDuplicates INPUT=$bam_sort_add OUTPUT=$bam_sort_add_dup M=$dump CREATE_INDEX=true`; # mark duplicated reads
                } else {
                        die "This is no $bam_sort_add, error occurs in the second step: AddOrReplaceReadGroups\n";
                }
                if (-e $bam_sort_add_dup) {
                        `java -jar $PICARD ReorderSam INPUT=$bam_sort_add_dup OUTPUT=$resort SEQUENCE_DICTIONARY= $DIC CREATE_INDEX=TRUE`; # resorted
                } else {
                        die "There is no $bam_sort_add_dup, error occurs in the third step: MarkDuplicates\n";
                }
                if (-e $resort) {
                        `java -jar $GATK -T SplitNCigarReads -R $GATKREF -I $resort -o $split -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 --filter_reads_with_N_cigar`; # SplitNCigarReads
                } else {
                        die "There is no $resort, error occurs in the fourth step: ReorderSam\n";
                }
                if (-e $split) {
                        `java -jar $GATK -T HaplotypeCaller -R $GATKREF -I $split -o $vcf`; # GATK Haplotype
                } else {
                        die "There is no $split, error occurs in the fifth step: SplitNCigarReads\n";
                }
                if (-e $vcf) {
                        `java -jar $GATK -T SelectVariants -selectType SNP -R $GATKREF -selectType MNP -V $vcf -o $snp_vcf`; # SelectVariants
                } else {
                        die "There is no $vcf, error occurs in the sixth step: HaplotypeCaller\n";
                }
                if (-e $snp_vcf) {
                        `rm $bam_sort $bam_sort_add $bam_sort_add_dup $resort`;
                        `mv $vcf $vcf_idx $snp_vcf $snp_vcf_idx $output`;
                } else {
                        die "There is no $snp_vcf, error occurs in the seventh sep: SelectVariants\n";
                }
        }
}
