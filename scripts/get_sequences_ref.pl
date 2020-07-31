#!/usr/bin/perl -w
use File::Basename;
use Getopt::Long;
my $opt = GetOptions( 'input:s', \$blast_result,
        'nuc:s', \$nuc,
        'output:s', \$output);
if (!($opt && $blast_result && $nuc)) {
    print STDERR "\nExample usage:\n";
    print STDERR "\n$0 -input=\"blast_result diretory\" /
    -nuc=\"nucleotide sequeces directory\" /
    -output=\"final reference sequences direcotry\"\n\n";
    exit;
}
my @fas=<$nuc/*fa>;
my %seq;
my ($spe,$tran);
foreach my $fa (@fas) {
        open "fil", "$fa";
        my $spec=basename($fa);
        (my $spe)=$spec=~/(.*)\_nuc\.fa/;
        while (<fil>) {
                chomp;
                if (/>/) {
                        s/>//;
                        $tran=$_;
                } else {
                        $seq{$spe}->{$tran}.=$_ if defined $tran;
                }
        }
}
my @blasts=<$blast_result/*result>;
mkdir $output;
foreach my $blast (@blasts) {
        my $file=basename($blast);
        (my $orth)=$file=~/(.*)\.blastp\_result/;
        open "fil", "$blast";
        while (<fil>) {
                chomp;
                my @a=split;
                (my $spec)=$a[1]=~/(.*)\_\d+/;
                $a[1]=~s/\s+$//;
                my $ref=$spec.".fa";
                open $spec, ">>$output/$ref";
                my $info=$orth."\t".$a[0]."\t".$a[1];
                print $spec ">$info\n$seq{$spec}->{$a[1]}\n" if defined $seq{$spec}->{$a[1]};
        }
}
