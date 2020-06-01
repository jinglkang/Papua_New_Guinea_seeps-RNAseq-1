# transform vcf file to input genetype file in bayescan
#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $opt = GetOptions( 'vcf:s', \ my $vcf,
    'pop_def:s', \ my $pop_def,
    'id_correlation_num:s', \my $cor_tab
    );

if (!($vcf && $pop_def && $cor_tab)) {
    print STDERR "\nExample usage:\n";
    print STDERR "\n$0 --vcf \"vcf file\" --pop_def \"the sample correlation with population\" --id_correlation_num \"correlation table output\"\n\n";
    exit;
}

my %hash=();
my @loci=();
my %genetype=();
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
            }
        }
        if (/^#/) {
            next;
        } elsif (/PASS/) {
            my @a=split;
            my $loci=$a[0]."_".$a[1];
            push @loci, $loci;
            my $j=0;
            my $sample_num=scalar(keys %hash);
            for (my $i = -$sample_num; $i <= -1; $i++) {
                $j++;
                my @b=split /\:/, $a[$i];
                $genetype{$hash{$j}}->{$loci}=$b[0];
            }
        }
}

my $loci_num=@loci;
print "[loci]=$loci_num\n\n";

my %pop=();
my @sample=();
my %allel_type=();
my %allel_num=();
open POP, "$pop_def" or die "cannot open $pop_def\n";
while (<POP>) {
    chomp;
    s/^\s+//;
    s/\s+$//;
    my @a=split;
    push @sample, $a[0];
    $pop{$a[1]}.=$a[0]."\t";
}

my $pop_num=scalar(keys %pop);
print "[populations]=$pop_num\n\n";

# get the total genetype of each loci in all individuals
&get_allel_from_all_sample();

my @samp_pop;
my $gene_num;
my $i=0;
foreach my $pop (keys %pop) {
    $i++;
    $pop{$pop}=~s/\s+$//;
    $pop{$pop}=~s/^\s+//;
    print "[pop]=$i\n";
    @samp_pop=split /\t/, $pop{$pop};
    $gene_num=2*@samp_pop;
    &write_info_according_loci();
}


sub get_allel_from_all_sample {
    my @gtype;
    # get all allel types of all sample in each loci
    foreach my $loci (@loci) {
        my %hash;
        foreach my $sample (@sample) {
            my $gtype;
            $gtype=$genetype{$sample}->{$loci};
            my @allel=split /\//, $gtype;
            foreach my $allel (@allel) {
                if ($hash{$allel} && $hash{$allel} >1) {
                    $hash{$allel}++;
                } elsif (!$hash{$allel}) {$hash{$allel}=1}
            }
        }
        my @allel_type=sort {$a<=>$b} keys %hash;
        @{$allel_type{$loci}}=@allel_type;
    }
}

sub write_info_according_loci {
    my $i=0;
    foreach my $loci (@loci) {
        $i++;
        my %hash2=();
        foreach my $sample (@samp_pop) {
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
            for (my $i = 0; $i < @a; $i++) {
                if ($allel_num{$a[$i]}) {
                    $allel_num{$a[$i]}+=$hash2{$key};
                } elsif (!$allel_num{$a[$i]}) {
                    $allel_num{$a[$i]}=$hash2{$key};
                }
            }
        }

        my $allel_total_type=@{$allel_type{$loci}};
        my @allel=@{$allel_type{$loci}};
        my $info="";
        foreach my $allel (@allel) {
            if ($allel_num{$allel}) {
                $info.=$allel_num{$allel}." ";
                } elsif (! $allel_num{$allel}) {
                    $allel_num{$allel}=0;
                    $info.=$allel_num{$allel}." ";
                }
        }
        $info=~s/\s+$//;
        print "$i $gene_num $allel_total_type $info\n";
        open COR, ">>$cor_tab";
        print COR "$loci $i\n";
    }
}
