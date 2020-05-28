# GATK work flow
# the reference fasta must in the current directory
#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
my $opt = GetOptions( 'fasta:s', \ my $fasta,
    'bam:s', \ my $bam_dir,
    'output:s', \ my $output,
    'tmp:s', \my $tmp);

if (!($opt && $fasta && $output &&$tmp)) {
    print STDERR "\nExample usage:\n";
    print STDERR "\n$0 --fasta \"reference sequence\" --bam \"bam files directory \(don't add the last \/\)\" --output \"output directory\" --tmp \"java tmp directory\"\n\n";
    exit;
}

unless (-e $output) {
    `mkdir $output`;
}
unless (-e $tmp) {
    `mkdir $tmp`;
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
    if (-e "./$output/$snp_vcf") {
        print STDERR "congrat, $name already has snp file! let's start the next one.\n";
        next;
    } else {
        print STDERR "\nIt starts the snp calling in $name\n";
        &gatk_prepare_bam($input_1, $bam_sort, $bam_sort_add, $bam_sort_add_dup, $resort, $split, $dump, $vcf);
        &gatk_generate_vcf($vcf, $snp_vcf, $name, $bam_sort, $bam_sort_add, $bam_sort_add_dup, $resort);
    }
    print STDERR "\nsnp calling in $name is done, now is next one\n";
}

# combine all snp files: CombineVariants
my $snp_vcf_1=$spec.".1.snp.vcf";
&merge_snp_vcf();

# get loci that appears in all samples
my $all_snp_1=$spec.".all.snp.1.vcf";
my $all_snp_2=$spec.".all.snp.2.vcf";
# get the headers all the loci with flag "set=Intersection"
`grep -e "^#" -e "section" $all_snp_1 > $all_snp_2`; 

# filter snp
&filter_snp();



sub gatk_prepare_bam {
    my ($input_1, $bam_sort, $bam_sort_add, $bam_sort_add_dup, $resort, $split, $dump, $vcf)=@_;
    `java -jar -Djava.io.tmpdir=$tmp $PICARD SortSam INPUT=$input_1 OUTPUT=$bam_sort SO=coordinate`; # sortsam
    if (-e $bam_sort) {
        `java -jar -Djava.io.tmpdir=$tmp $PICARD AddOrReplaceReadGroups INPUT=$bam_sort OUTPUT=$bam_sort_add SO=coordinate RGID=Label RGLB=Label RGPL=illumina RGPU=Label RGSM=Label`; # add read group
    } else {
        die "There is no $bam_sort, error occurs in the first step: SortSam\n";
    }
    if (-e $bam_sort_add) {
        `java -jar -Djava.io.tmpdir=$tmp $PICARD MarkDuplicates INPUT=$bam_sort_add OUTPUT=$bam_sort_add_dup M=$dump CREATE_INDEX=true`; # mark duplicated reads
    } else {
        die "This is no $bam_sort_add, error occurs in the second step: AddOrReplaceReadGroups\n";
    }
    if (-e $bam_sort_add_dup) {
        `java -jar -Djava.io.tmpdir=$tmp $PICARD ReorderSam INPUT=$bam_sort_add_dup OUTPUT=$resort SEQUENCE_DICTIONARY= $DIC CREATE_INDEX=TRUE`; # resorted
    } else {
        die "There is no $bam_sort_add_dup, error occurs in the third step: MarkDuplicates\n";
    }
    if (-e $resort) {
        `java -jar -Djava.io.tmpdir=$tmp $GATK -T SplitNCigarReads -R $GATKREF -I $resort -o $split -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 --filter_reads_with_N_cigar`; # SplitNCigarReads
    } else {
        die "There is no $resort, error occurs in the fourth step: ReorderSam\n";
    }
    if (-e $split) {
        `java -jar -Djava.io.tmpdir=$tmp $GATK -T HaplotypeCaller -R $GATKREF -I $split -o $vcf --annotation MappingQualityZero`; # GATK Haplotype
    } else {
        die "There is no $split, error occurs in the fifth step: SplitNCigarReads\n";
    }
}

sub gatk_generate_vcf {
    my ($vcf, $snp_vcf, $name, $bam_sort, $bam_sort_add, $bam_sort_add_dup, $resort)=@_;
    if (-e $vcf) {
        `java -jar -Djava.io.tmpdir=$tmp $GATK -T SelectVariants -selectType SNP -R $GATKREF -selectType MNP -V $vcf -o $snp_vcf`; # SelectVariants
    } else {
        die "There is no $vcf, error occurs in the sixth step: HaplotypeCaller\n";
    }
    if (-e $snp_vcf) {
        open SNP, "$snp_vcf" or die "can not open $snp_vcf\n";
        my $snp_vcf_1=$name.".1.snp.vcf";
        open SNP1, ">$snp_vcf_1" or die "can not open $snp_vcf_1\n";
        # filter the low quality genotype (QG<30.0)
        while (<SNP>) {
            chomp;
            my @info=();
            my @info1=();
            if (/^#/) {
            print SNP1 "$_\n";
        } else {
            @info=split;
            @info1=split /\:/, $info[-1];
            if ($info1[-2] && $info1[-2]>=30) {
                print SNP1 "$_\n";
                }
            }
        }
        `cp $snp_vcf_1 $output/`;
        `rm $bam_sort $bam_sort_add $bam_sort_add_dup $resort`;
        } else {
            die "There is no $snp_vcf, error occurs in the seventh sep: SelectVariants\n";
        }
}

sub merge_snp_vcf {
    my @snp=<$output/*.snp.vcf>;
    my $all_snp_1=$spec.".all.snp.1.vcf";
    my $merge="java -jar -Djava.io.tmpdir=$tmp $GATK -T CombineVariants -R $GATKREF";
    foreach (@snp) {
        (my $sample)=basename($_)=~/(.*).snp.vcf/;
        $merge.=" --variant:$sample $_";
    }
    $merge.=" -o $all_snp_1 -genotypeMergeOptions UNIQUIFY";
    `$merge`;
    unless (-e $all_snp_1) {
        die "There is no $all_snp_1\n"
    }
}

sub filter_snp {
    # filter SNPs
    my $all_snp_3=$spec.".all.snp.3.vcf";
    my $filter="java -jar -Djava.io.tmpdir=$tmp $GATK -T VariantFiltration -R $GATKREF -V $all_snp_2 -o $all_snp_3 ";
    # in 35 bp window shouldn't have more than 3 SNPs. shows "SnpCluster" if failed
    $filter .= '--clusterSize 3 --clusterWindowSize 35 ';
    # QD: variant confidence/quality by depth; should more than 5
    $filter .= '--filterName "low coverage" --filterExpression "QD < 5.0" ';
    # DP: read depth. should more 10
    $filter .= '--filterName "no reads" --filterExpression "DP < 10" ';
    # ReadPosRankSum: z-score from wilcoxon rank sum test of Alt vs. ref read position bias; should more than -8.0
    $filter .= '--filterName "failed RPRS" --filterExpression "ReadPosRankSum < -8.0" ';
    # MQRankSum: z-score from wilcoxon rank sum test of Alt vs. Ref read mapping quality; should more than -12.5
    $filter .= '--filterName "failed MQRS" --filterExpression "MQRankSum < -12.5" ';
    # MQ: "RMS (root mean square) of the mapping quality of reads across all samples"; should more than 40
    $filter .= '--filterName "failed MQ" --filterExpression "MQ < 40.0" ';
    # FS: "phred-scaled p-value using Fisher's exact test to detect strand bias"; should more than 60
    $filter .= '--filterName "failed FS" --filterExpression "FS > 60.0" ';
    # MQ0: "Total Mapping Quality Zero Reads"; MP0 should not more than 4 and variants for which reads with 0 mapping quality constitute 
    # should not more than 10%
    $filter .= '--filterName "failed 10% reads MQ0" --filterExpression "MQ0 >= 4 && MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)"';
    `$filter`;
    if (-e $all_snp_3) {
        print STDERR "congrat, GATK finish success\n";
        my $all_snp_4=$spec.".all.snp.final.vcf";
        `mv $all_snp_3 $output/$all_snp_4`;
    } else {
        die "The last step crashed: there is no $all_snp_3\n";
    }
}
