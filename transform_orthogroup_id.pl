#!/usr/bin/perl -w
# transform gene id to the common orthogroup id
use File::Basename;
use Getopt::Long;
my ($blast_result);
my $opt = GetOptions( 'input:s', \$matrix,
        'fasta:s', \$fasta,
        'output:s', \$output);
my $help;
if (!($opt && $matrix && $fasta) || $help) {
    print STDERR "\nExample usage:\n";
    print STDERR "\n$0 -input=\"reads matrix\" -fasta=\"fasta file\" -output=\"new matrix output\"\n\n";
    exit;
}
my %hash;
open "fil1", "$fasta";
while (<fil1>) {
        chomp;
        if (/>/) {
                s/>//;
                my @a=split;
                $hash{$a[-1]}=$a[0];
        }
}
open "fil2", "$matrix";
open "fil3", ">$output";
while (<fil2>) {
        chomp;
        print fil3 "$_\n" if /^\s+/;
        my @a=split;
        my $info="";
        if (defined $hash{$a[0]}){
                print fil3 "$hash{$a[0]}\t";
                for (my $i = 1; $i <= @a-1; $i++) {
                        $info.="$a[$i]\t";
                }
                $info=~s/\s+$//;
                print fil3 "$info\n";
        }
}
