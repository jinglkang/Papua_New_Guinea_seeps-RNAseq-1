#!/usr/bin/perl -w
use File::Basename;
use Getopt::Long;
my ($blast_result);
my $opt = GetOptions( 'fasta:s', \$fasta,
    'blast_result:s', \$blast_result,
    'output:s', \$output);
my $help;
if (!($opt && $fasta && $blast_result) || $help) {
    print STDERR "\nExample usage:\n";
    print STDERR "\n$0 -fasta=\"transcript sequence\" -blast_result=\"final blast result directory\" -output=\"output directory\" \n\n";
    print STDERR "Option:\n";
    exit;
}
mkdir $output;

my %seq;
# build hash between the sequence id and sequence
open "fasta", "$fasta" or die "can not open $fasta";
while (<fasta>) {
        chomp;
        if (/>/) {
                s/>//;
                $seq_name=$_;
        } else {
                $seq{$seq_name}.=$_;
        }
}
my @blast;
@blast=<$blast_result/*.blastp_result>;
foreach my $blast (@blast) {
        my ($orth_id)=basename($blast)=~/(.*)\.blastp_result/;
        my $orth_fas=$orth_id.".fas"; # construct new fasta file which contains the orhogroup sequence
        open "orth_fas", ">$output/$orth_fas";
        open "blast", "$blast";
        while (<blast>) {
                chomp;
                my @a=split;
                print orth_fas ">$a[1]\n$seq{$a[1]}\n" if $seq{$a[1]} or die "there is no hash %seq";
        }
        close orth_fas;
        my $new_fas1=$orth_id."-1.fas";
        system("clustalo -i $output\/$orth_fas -t Protein -o $output\/$new_fas1 --outfmt=fa");
        my $new_fas2=$orth_id."-2.fas";
        open "new_fas1", "$output/$new_fas1";
        open "new_fas2", ">$output/$new_fas2";
        my %new_seq=();
        while (<new_fas1>) {
                chomp;
                if (/>/) {
                        s/>//;
                        $a=$_;
                } else {
                        $new_seq{$a}.=$_;
                }
        }
        foreach my $key (sort keys %new_seq) {
                print new_fas2 ">$key\n$new_seq{$key}\n" if $new_seq{$key} or die "$new_fas1 hash is unvalid";
        }
        system("rm $output\/$orth_fas") if $orth_fas or die "cannot rm $orth_fas";
        system("rm $output\/$new_fas1") if $new_fas1 or die "cannot rm $new_fas1";
        system("mv $output\/$new_fas2 $output\/$orth_fas");

}
