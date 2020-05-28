# transform vcf file to input genetype file in bayescan
#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $opt = GetOptions( 'vcf:s', \ my $vcf,
    'pop_def:s', \ my $pop_def
    );

if (!($vcf && $pop_def)) {
    print STDERR "\nExample usage:\n";
    print STDERR "\n$0 --vcf \"vcf file\" --pop_def \"the sample correlation with population\"\n\n";
    exit;
}

my %hash=();
my @loci=();
my %genetype=();
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

my %pop=();
open POP, "$pop_def";
while (<POP>) {
        chomp;
        s/^\s+//;
        s/\s+$//;
        my @a=split;
        $pop{$a[1]}.=$a[0]."\t";
}

foreach my $pop (keys %pop) {
        $pop{$pop}=~s/\s+$//;
        print "[population]= $pop\t$pop{$pop}\n";
        my @sample=split /\t/, $pop{$pop};
        my $gene_num=2*scalar(@sample);
        foreach my $loci (@loci) {
                my %hash2=();
                my $info;
                foreach my $sample (@sample) {
                        my $gtype;
                        $gtype=$genetype{$sample}->{$loci};
                        if ($hash2{$gtype} && $hash2{$gtype}>=1) {
                                $hash2{$gtype}++;
                        } elsif (!$hash2{$gtype}) {
                                $hash2{$gtype}=1;
                        }
                }

                my %allel_num=();
                foreach my $key (sort keys %hash2) {
                        my @a=split /\//, $key;
                        for (my $i = 0; $i < 2; $i++) {
                                if ($allel_num{$a[$i]}) {
                                        $allel_num{$a[$i]}+=$hash2{$key};
                                } elsif (!$allel_num{$a[$i]}) {
                                        $allel_num{$a[$i]}=$hash2{$key};
                                }
                        }
                }

                my $allel_total_type=scalar(keys %allel_num);
                foreach my $key (sort {$a<=>$b} keys %allel_num) {
                        $info.=$allel_num{$key}."\t";
                }
                $info=~s/\s+$//;
                print "$loci $gene_num $allel_total_type $info\n";
        }
}
