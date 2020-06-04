use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $opt = GetOptions( 'vcf:s', \ my $vcf,
    'pop_def:s', \ my $pop_def,
    'loci:s', \my $loci);

if (!($vcf)) {
    print STDERR "\nExample usage:\n";
    print STDERR "\n$0 --vcf \"vcf file\" --pop_def \"which sample in which population\" --loci \"loci file\"\n\n";
    exit;
}
my @loci=();
my %genetype=();
my @sample=();
my %hash;
my @samp_pop;
my $sample_num;
open VCF, "$vcf" or die "cannot open $vcf\n";
while (<VCF>) {
        chomp;
        if (/^##GATKCommandLine/) { # get the sample in the merged vcf
            my $variant="";
            ($variant)=$_=~/variant=(\[.*\])\s+out\=/;
            my @variant=split /\,/, $variant if $variant;
            my $i=0;
            foreach my $var (@variant) {
                $i++;
                (my $sample)=$var=~/name=(.*)\.1\ssource=/;
                $sample=~s/^\s+//;
                $sample=~s/\s+$//;
                $hash{$i}=$sample;
                push @sample, $sample;
            }
        }
        if (/^#/) {
            next;
        } elsif (/PASS/) {  # get the snp with the flag "PASS", and get the genetype of each sample in all loci
            my @a=split;
            my $loci=$a[0]."_".$a[1];
            push @loci, $loci;
            $sample_num=scalar(keys %hash);
            my $j=0;
            for (my $i = -$sample_num; $i <= -1; $i++) {
                $j++;
                my @b=split /\:/, $a[$i];
                my @c=split /\//, $b[0];
                my $allel1=$c[0]+1;
                my $allel2=$c[1]+1;
                $genetype{$hash{$j}}->{$loci}="0"."$allel1"."0"."$allel2";
            }
        }
}

my @loci_need;
open LOCI, "$loci" or die "cannot open $loci\n";
while (<LOCI>) {
    chomp;
    my @a=split;
    push @loci_need, $a[0];
}

my $snp_num=@loci_need;
my ($spe)=$vcf=~/(.*?)\./;
print "# $snp_num snps in $spe\n";

my $header;
foreach my $loci (@loci_need) { # print the header of output (all loci)
    $header.=$loci."\t";
}
$header=~s/\s+$//;

my %pop=();
open POP, "$pop_def" or die "cannot open $pop_def\n"; # get the sample in each population
while (<POP>) {
    chomp;
    s/^\s+//;
    s/\s+$//;
    my @a=split;
    $pop{$a[1]}.=$a[0]."\t";
}

my $pop_num=scalar(keys %pop);

print "$sample_num\t$pop_num\t$snp_num\t2\t2\n";

my $i;
foreach my $pop (sort keys %pop) {
    $i++;
    print "pop$i\n";
}

print "pop\tind\t$header\n";
$i=0;
foreach my $pop (sort keys %pop) { # print the genetype of loci in samples according population
    $i++;
    $pop{$pop}=~s/\s+$//;
    $pop{$pop}=~s/^\s+//;
    @samp_pop=split /\t/, $pop{$pop};
    &write_genetype_sample();
    @samp_pop=();
}

sub write_genetype_sample {
    foreach my $sample (@samp_pop) {
        my $total_genetype="";
        foreach my $loci (@loci_need) {
            my $genetype=$genetype{$sample}->{$loci};
            $total_genetype.=$genetype."\t";
        }
        $total_genetype=~s/\s+$//;
        print "$i\t$sample\t$total_genetype\n";
    }
}
