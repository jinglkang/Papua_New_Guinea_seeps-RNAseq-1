use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $opt = GetOptions( 'vcf:s', \ my $vcf,
    'pop_def:s', \ my $pop_def);

if (!($vcf)) {
    print STDERR "\nExample usage:\n";
    print STDERR "\n$0 --vcf \"vcf file\" --pop_def \"which sample in which population\" \n\n";
    exit;
}
my @loci=();
my %genetype=();
my @sample=();
my %hash;
my @samp_pop;
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
            my $sample_num=scalar(keys %hash);
            my $j=0;
            for (my $i = -$sample_num; $i <= -1; $i++) {
                $j++;
                my @b=split /\:/, $a[$i];
                $genetype{$hash{$j}}->{$loci}=$b[0];
            }
        }
}

my $header;
foreach my $loci (@loci) { # print the header of output (all loci)
    $header.=$loci."\t";
}
$header=~s/\s+$//;
print "\t$header\n";

my %pop=();
open POP, "$pop_def" or die "cannot open $pop_def\n"; # get the sample in each population
while (<POP>) {
    chomp;
    s/^\s+//;
    s/\s+$//;
    my @a=split;
    $pop{$a[1]}.=$a[0]."\t";
}

foreach my $pop (sort keys %pop) { # print the genetype of loci in samples according population
    $pop{$pop}=~s/\s+$//;
    $pop{$pop}=~s/^\s+//;
    @samp_pop=split /\t/, $pop{$pop};
    &write_genetype_sample();
    @samp_pop=();
}

sub write_genetype_sample {
    foreach my $sample (@samp_pop) {
        my $total_genetype="";
        foreach my $loci (@loci) {
            my $genetype=$genetype{$sample}->{$loci};
            $total_genetype.=$genetype."\t";
        }
        $total_genetype=~s/\s+$//;
        print "$sample\t$total_genetype\n";
    }
}
