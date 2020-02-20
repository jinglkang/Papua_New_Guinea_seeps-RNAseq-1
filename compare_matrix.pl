#!/usr/bin/perl -w
use File::Basename;
use Getopt::Long;
# compare the difference of reads number between using all transcript and using selected transcripts
my $opt = GetOptions( 'total_matrix:s', \$total_matrix,
        'matrixs:s', \$matrix);
if (!($total_matrix && $matrix)) {
        print STDERR "\nExample usage:\n";
        print STDERR "\n$0 -total_matrix=\"total reads number matrix\"";
        print STDERR "   -matrix = directory of a list of matrixs\n\n";
        exit;
}
my @matrix = <$matrix/*.matrix>;
my (@ind, @org, @in);
my (%hash, %num1, %num2);
open "final_total", "$total_matrix";
while (<final_total>) {
        chomp;
        if (/^\s+/) {
                s/^\s+//;
                s/\s+$//;
                @ind=split;
        } else {
                my @a=split;
                $hash{$a[0]}=$a[1];
                my $org=$a[0];
                push @org,$org;
                for (my $i = 2; $i < @a; $i++) {
                        $num1{$org}->{$ind[$i-2]}=$a[$i];
                }
        }
}
foreach my $file (@matrix) {
        open "fil", "$file";
        while (<fil>) {
                chomp;
                if (/^\s+/) {
                        s/^\s+//;
                        s/\s+$//;
                        @in=split;
                } else {
                        my @a=split;
                        for (my $i = 1; $i < @a; $i++) {
                                $num2{$a[0]}->{$in[$i-1]}=$a[$i];
                        }
                }
        }
}
print "\t";
foreach my $ind (@ind) {
        print "$ind\t";
}
print "\n";
foreach my $org (sort @org) {
        print "$org\t$hash{$org}\t";
        foreach my $ind (@ind) {
                my $read_num1=$num1{$org}->{$ind};
                my $read_num2=$num2{$org}->{$ind};
                print "$read_num1($read_num2)\t";
        }
        print "\n";
}
