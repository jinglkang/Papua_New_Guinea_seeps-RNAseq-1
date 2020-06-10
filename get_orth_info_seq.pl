#!/usr/bin/perl -w
use File::Basename;
use Getopt::Long;
my $opt = GetOptions( 'input:s', \$blast_result,
        'nuc:s', \$nuc);
if (!($opt && $blast_result && $nuc)) {
    print STDERR "\nExample usage:\n";
    print STDERR "\n$0 -input=\"blast_result diretory\" /
    -nuc=\"nucleotide sequeces directory\"\n\n";
    exit;
}
my @fas=<$nuc/*fasta>;
my %seq;
my ($spe,$tran);
foreach my $fa (@fas) {
        open "fil", "$fa";
        my $spec=basename($fa);
        (my $spe)=$spec=~/(.*)\.fasta/;
        while (<fil>) {
                chomp;
                if (/>/) {
                        s/>//;
                        my @a=split;
                        $tran=$a[0];
                } else {
                        $seq{$spe}->{$tran}.=$_ if defined $tran;
                }
        }
}
my @blasts=<$blast_result/*result>;

my $j=0;
foreach my $blast (@blasts) {
    $j++;
        my $file=basename($blast);
        (my $orth)=$file=~/(.*)\.blastp\_result/;
        open "fil", "$blast";
        my $i=0;
        while (<fil>) {
            $i++;
                chomp;
                my @a=split;
                my $uni_seq=$seq{"uniprot"}->{$a[0]};
                (my $spec)=$a[1]=~/(.*)\_\d+/;
                $a[1]=~s/\s+$//;
                my $info=$orth."\t".$a[0]."\t".$a[1];
                print "$j\t$a[0]\t$uni_seq\t$a[1]\t$seq{$spec}->{$a[1]}\t$a[-3]\t$a[-2]\t$a[-1]\n" if defined $seq{$spec}->{$a[1]} && $i==1;
                print "\t\t\t$a[1]\t$seq{$spec}->{$a[1]}\t$a[-3]\t$a[-2]\t$a[-1]\n" if defined $seq{$spec}->{$a[1]} && $i>1;
        }
}
