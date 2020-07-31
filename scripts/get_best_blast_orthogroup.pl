#!/usr/bin/perl -w
use File::Basename;
use Getopt::Long;
my ($blast_result);
my $opt = GetOptions( 'blast_result:s', \$blast_result,
    'reads_matrix:s', \$matrix);
my $help;
if (!($opt && $blast_result && $matrix) || $help) {
    print STDERR "\nExample usage:\n";
    print STDERR "\n$0 -blast_result=\"blast result directory\" -reads_matrix=\"reads matrix directory\"\n\n";
    print STDERR "Option:\n";
    exit;
}
my @matrix=<$matrix/*.matrix>; # get the reads number from reads number matrix directory
my %numb;
foreach my $mat (@matrix) {
    open "fil", "$mat";
    (my $fish)=basename($mat)=~/total_(.*)\.gene/;
    while (<fil>) {
        chomp;
        my @a=split;
        my $read_num=0;
        next if /^\s+/;
        for (my $i = 1; $i < @a; $i++) {
            $read_num=$read_num+$a[$i];
        }
        $numb{$fish}->{$a[0]}=$read_num;
    }
}
my @blast = <$blast_result/*result>; # get blast result from blast result directory
foreach my $blast(@blast) {
    open Blast_result, "$blast" or die "cannot open $blast"; # open blast result file
    open Blast_top10, ">blast_top10.txt" or die "cannot open blast_top10.txt"; # save the top10 hit per transcript to this file
    my %num=();
    while (<Blast_result>) {
        chomp;
        my @a=split;
        my $trans=$a[0];
        $num{$trans}++; # save the top10 hit per transcript
        if ($num{$trans} < 11) {
                print Blast_top10 "$_\n";
        }
    }
    open Blast_top10, "blast_top10.txt" or die "cannot open blast_top10.txt"; # open the top10 result file
    open Best_pro_tra, ">blast_pro_tra.txt" or die "cannot open blast_pro_tra.txt"; # save the best hit (protein to species)
    my %info=();
    my @trans=();
    while (<Blast_top10>) {
        chomp;
        my @a=split;
        my $pro=$a[1]; # protein id
        my $tran=$a[0]; # transcript id
        push @trans, $tran;
        (my $spec)=$tran=~/(.*)_/; # species name
        my $eval=$a[-2]; # e-value of each transcript
        my $info1;
        for (my $i = 2; $i < @a; $i++) {
                $info1.=$a[$i]."\t";
        }
        $info1=~s/\t$//;
        my $info2=$a[1]."\t".$a[0]."\t".$info1; # get transcript info of per protein
        # select the best blast result if a transcript were blasted to different position of the same protein
        if ((defined $info{$pro}->{$spec}->{$tran} && $info{$pro}->{$spec}->{$tran}->{eval} >= $a[-2]) || (!defined $info{$pro}->{$spec}->{$tran})) {
                    $info{$pro}->{$spec}->{$tran}={
                        eval => $a[-2],
                        info => $info2,
                        read_um => $numb{$spec}->{$tran}
                    };
            }
    }
    my %count;
    my @tran=grep { ++$count{ $_ } < 2; } @trans;
    my @specs=("Apoly","Padel","Daru","Acura","Ocomp","Pmol"); # get the best hit of each protein of these species;
    foreach my $key (keys %info) {
        foreach my $spec (@specs) {
            foreach my $tran (@tran) {
                if (defined $info{$key}->{$spec}->{$tran} && defined $numb{$spec}->{$tran}) {
                    print Best_pro_tra "$info{$key}->{$spec}->{$tran}->{info}\t$numb{$spec}->{$tran}\n";
                }
            }
        }
    }
    open Best_pro_tra, "blast_pro_tra.txt" or die "cannot open blast_best_pro_tra.txt"; # obtain proteins that has a best hit to all species
    open Most_reads, ">blast_pro_tra_spec_most_reads.txt" or die "cannot open blast_best_pro_tra_5_spec.txt";
    my %info1=();
    my %hash=();
    while (<Best_pro_tra>) {
        chomp;
        my @a=split;
        my $pro=$a[0];
        next if $a[-3] > 1E-25; # filter hit with high e-value
        (my $spec)=$a[1]=~/(.*)_\d+/;
        if ((defined $info1{$pro}->{$spec}->{num} && $info1{$pro}->{$spec}->{num} <= $a[-1]) || (!defined $info1{$pro}->{$spec}->{num})) {
            $info1{$pro}->{$spec}={
                info => $_,
                num => $a[-1]
            };
        }
    }
    foreach my $key (keys %info1) {  # sort blast hits according proteins
        foreach my $spec (@specs) {
            print Most_reads "$info1{$key}->{$spec}->{info}\n" if defined $info1{$key}->{$spec};
        }
    }
    open Most_reads, "blast_pro_tra_spec_most_reads.txt" or die "cannot open blast_pro_tra_spec_most_reads.txt"; # obtain proteins that has a hit to all species
    open spec_6, ">blast_best_pro_tra_6_spec.txt" or die "cannot open blast_best_pro_tra_6_spec.txt";
    my %info1_1=();
    my %hash1=();
    while (<Most_reads>) {
        chomp;
        my @a=split;
        my $pro=$a[0];
        (my $spec)=$a[1]=~/(.*)_\d+/;
        my $eval=$a[-3];
        my $score=$a[-2];
        if ($info1_1{$pro}) {
            $hash1{$pro}++;
            my $info2 = "$info1_1{$pro}->{info}"."\n"."$_";
            my $eval2=$info1_1{$pro}->{eval}+$eval;
            $info1_1{$pro}={
                num_spec => $hash1{$pro},
                info => $info2,
                eval => $eval2
            }
        }
        else {$info1_1{$pro}={
                info => $_,
                eval => $eval
                };
            }
    }
    foreach my $key (keys %info1_1) {  # select if the protein has a hit to all six species
        if (exists $info1_1{$key}->{num_spec} && $info1_1{$key}->{num_spec}==5 && $info1_1{$key}->{eval} <=1E-30) {
            print spec_6 "$info1_1{$key}->{info}\n";
        }
    }
    open spec_6, "blast_best_pro_tra_6_spec.txt" or die "cannot open blast_best_pro_tra_6_spec.txt"; # get the best one if we get more than one final protein
    open fil1, ">temp" or die "cannot open temp";
    my %info2=();
    while (<spec_6>) {
        chomp;
        my @a=split;
        if ($info2{$a[0]}) {
            my $info="$info2{$a[0]}->{infor}"."\t"."$_"."\t";
            my $score=$info2{$a[0]}->{score} + $a[-2];
            $info2{$a[0]}={
                score => $score,
                infor => $info
            };
        }
        else{$info2{$a[0]}={
            score => $a[-2],
            infor => $_
            };
        }
    }
    foreach my $key (keys %info2) {
        my $value=$info2{$key}->{infor};
        my $score=$info2{$key}->{score};
        print fil1 "$value\t$score\n";
    }
    my %score=();
    my %info3=();
    my $name="";
    open fil1, "temp";
    while (<fil1>) {
        chomp;
        @a=split;
        next if /^\s+$/;
        s/^\s+//;
        ($name)=$a[0]=~/(sp)\|.*/;
        if ($score{$name}) {
            if ($score{$name} < $a[-1]) {
                $info3{$name}=$_;
                $score{$name}=$a[-1];
            }
        } else {
            $score{$name}=$a[-1];
            $info3{$name}=$_;}
    }
    my $final="";
    $final=$info3{$name};
    $final=~s/\t\d+$// if defined $final;
    $final=~s/^\s+// if defined $final;
    $final=~s/\tsp/\nsp/g if defined $final;
    $orth=basename($blast);
    my $final_blast_orth_group="final_blast_orth_group";
    unless (-f $final_blast_orth_group) { # if there is no directory of $final_blast_orth_group, mkdir it
        mkdir $final_blast_orth_group;
    }
    open orth_gr, ">$final_blast_orth_group/$orth";
    print orth_gr "$final\n" if defined $final;
    `rm blast_top10.txt blast_pro_tra.txt blast_pro_tra_spec_most_reads.txt blast_best_pro_tra_6_spec.txt temp`; # delete all of the temp files
}
my @files=<final_blast_orth_group/*result>;
foreach my $file (@files) {
    if (-z $file) {
        `rm $file`; # delete the empty files
    }
}
