# output the Gblocks result files and the fasta files with long sequence (>50 condons);
#!/usr/bin/perl -w
use File::Basename;
use Getopt::Long;
my ($blast_result);
my $opt = GetOptions( 'input:s', \$input,
    'gblocks_output:s', \$output1,
    'long_seq:s', \$output2);
my $help;
if (!($opt && $input && $output1 && $output2) || $help) {
    print STDERR "\nExample usage:\n";
    print STDERR "\n$0 -input=\"alignment files directory\" -gblocks_output=\"Gblocks output directory\" -long_seq=\"fasta files directory from Gblocks with longer sequences\" \n\n";
    print STDERR "Option:\n";
    exit;
}

mkdir $output1;
mkdir $output2;
my @fas=<$input/*.fas>;
foreach my $fas (@fas) {
        system("Gblocks $fas -t=p");
        my ($name)=basename($fas)=~/(.*)\.fas/;
        my $output=$name.".fas-gb";
        if (-z $output) { # if the output file is empty, delete it; else move it to $output1
                system("rm $input/$output");
        } else {
                system("mv $input/$output $output1");
        }
        system("rm $input/*htm");
}
my @gb=<$output1/*gb>;
foreach my $gb (@gb) {
        open GB, "$gb" or die "can not open $gb";
        my ($fas)=basename($gb)=~/(.*)-gb/;
        my %seq=();
        while (<GB>) {
                chomp;
                if (/>/) {
                        s/>//;
                        $seq_id=$_;
                } else {
                        s/\s+//g;
                        $seq{$seq_id}.=$_;
                }
        }
        my $length=0;
        foreach my $key (sort keys %seq) {
                $length=length($seq{$key});
                if ($length > 0) {
                        open FAS, ">>$output1/$fas" or die "cannot open $output1/$fas";
                        print FAS ">$key\n$seq{$key}\n";
                }
        }
        system("rm $gb");
        if ($length >=50) {
                system("cp $output1/$fas $output2");
        }
}
