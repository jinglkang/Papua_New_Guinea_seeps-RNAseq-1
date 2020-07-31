#!/usr/bin/perl -w
open orth_gro, "$ARGV[0]"; # open orthogroup table (Orthogroups.GeneCount.tsv)
open orth_gro_need, ">orthogroup_needed";
while (<orth_gro>) {
        chomp;
        next if /Acura/; # ignore the header
        my @a=split;
        my %num;
        for (my $i = 1; $i <= @a-1; $i++) { # if this orthogroup have one species without transcript, then remove it
                if ($a[$i]==0) {
                        $num{$a[0]}=0;
                        last;
                } else {
                        $num{$a[0]}=1;
                }
        }
        print orth_gro_need "$_\n" if $num{$a[0]}==1;
}
my $dir1="orthogroup_needed_sequences"; # creat a directory to save sequences of the orthogroup_needed
mkdir $dir1 or die "can not creat $dir1";
open orth_gro_need, "orthogroup_needed";
my ($count1,$count2);
while (<orth_gro_need>) {
        chomp;
        $count1++;
        my @a=split;
        my $fa=$a[0].".fa";
        print "now copy $count1 files to $dir1\n";
        `cp ./Orthogroup_Sequences/$fa $dir1`;
}
my $dir2="blastp_result"; # creat a directory to save blast results
mkdir $dir2 or die "can not creat $dir2";
my @fa=<$dir1/*.fa>;
my $uni_pro="~/Desktop/Annotation_database/swiss-prot/uniprot-filtered-reviewed_yes.fasta";
foreach my $fa(@fa) {
        $count2++;
        my ($fa_name)=$fa=~/\/(.*)\.fa/;
        my $result_blast=$fa_name.".blastp_result";
        print "now is blastp the $count2 file to uniprot";
        `blastp -query $fa -db $uni_pro -outfmt 6 -num_threads 32 -out $dir2/$result_blast`;
}
