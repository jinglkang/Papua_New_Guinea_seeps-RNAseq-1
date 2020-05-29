use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $opt = GetOptions( 'vcf:s', \ my $vcf);

if (!($vcf)) {
    print STDERR "\nExample usage:\n";
    print STDERR "\n$0 --vcf \"vcf file\" \n\n";
    exit;
}
my @loci=();
my %genetype=();
my @sample=();
my %hash;
open VCF, "$vcf" or die "cannot open $vcf\n";
while (<VCF>) {
        chomp;
        if (/^##GATKCommandLine/) {
            my $variant="";
            ($variant)=$_=~/variant=(\[.*\])\s+out\=/;
            my @variant=split /\,/, $variant if $variant;
            my $i=0;
            foreach my $var (@variant) {
                $i++;
                (my $sample)=$var=~/name=(.*)\ssource=/;
                $sample=~s/^\s+//;
                $sample=~s/\s+$//;
                $hash{$i}=$sample;
                push @sample, $sample;
            }
        }
        if (/^#/) {
            next;
        } elsif (/PASS/) {
            my @a=split;
            my $loci=$a[0]."_".$a[1];
            push @loci, $loci;
            my $j=0;
            for (my $i = -19; $i <= -1; $i++) {
                $j++;
                my @b=split /\:/, $a[$i];
                $genetype{$hash{$j}}->{$loci}=$b[0];
            }
        }
}

my $header;
foreach my $loci (@loci) {
    $header.=$loci."\t";
}
$header=~s/\s+$//;
print "\t$header\n";
foreach my $sample (@sample) {
    my $total_genetype="";
    foreach my $loci (@loci) {
        my $genetype=$genetype{$sample}->{$loci};
        $total_genetype.=$genetype."\t";
    }
    $total_genetype=~s/\s+$//;
    print "$sample\t$total_genetype\n";
}
