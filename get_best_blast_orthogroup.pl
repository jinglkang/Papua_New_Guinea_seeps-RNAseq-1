open Blast_result, "$ARGV[0]"; # open blast result file
open Blast_top10, ">blast_top10.txt"; # save the top10 hit per transcript to this file
while (<Blast_result>) {
        chomp;
        @a=split;
        $trans=$a[0];
        $num{$trans}++; # save the top10 hit per transcript
        if ($num{$trans} < 11) {
                print Blast_top10 "$_\n";
        }
}
open Blast_top10, "blast_top10.txt"; # open the top10 result file
open Best_pro_tra, ">blast_best_pro_tra.txt"; # save the best hit (protein to species)
while (<Blast_top10>) {
        chomp;
        @a=split;
        $pro=$a[1]; # protein id
        $tran=$a[0]; # transcript id
        ($spec)=$tran=~/(.*)_/; # species name
        $eval=$a[-1]; # e-value of each transcript
        for ($i = 2; $i < @a; $i++) {
                $info1.=$a[$i]."\t";
        }
    $info1=~s/$\t//;
    $info2=$a[1]."\t".$a[0]."\t".$info1;
    $info1="";
        if ($info{$pro}->{$spec}) {
                if ($info{$pro}->{$spec}->{eval} < $eval) { # if the new e-value is lower, save the new one
                        $info{$pro}->{$spec}={
                                tran => $tran,
                                eval => $eval,
                                info => $info2
                        };
                }
        } else {
                        $info{$pro}->{$spec}={
                        tran => $tran,
                        eval => $eval,
                        info => $info2
                };
        }
}
@specs=("Apoly","Padel","Daru","Acura","Ocomp"); # get the best hit: protein to species;
foreach $key (keys %info) {
        foreach $spec (@specs) {
                if ($info{$key}->{$spec}) {
                        print Best_pro_tra "$info{$key}->{$spec}->{info}\n";
                }
        }
}
open Best_pro_tra, "blast_best_pro_tra.txt"; # esure protein has hit to each species
open spec_5, ">blast_best_pro_tra_5_spec.txt";
while (<Best_pro_tra>) {
        chomp;
        @a=split;
        $pro=$a[0];
        ($spec)=$a[1]=~/(.*)_\d+/;
        $eval=$a[-2];
        $score=$a[-1];
        if ($info1{$pro}) {
                $hash{$pro}++;
                $eval = $eval + $info1{$pro}->{eval};
                $info2 = "$info1{$pro}->{info}"."\n"."$_";
                $info1{$pro}={
                        num_spec => $hash{$pro},
                        eval => $eval,
                        info => $info2
                }
        }
        else {$info1{$pro}={
                specs => $spec,
                eval => $eval,
                score => $score,
                info => $_
        };}
}
foreach $key (keys %info1) {
        if ($info1{$key}->{num_spec}==4) {
                $ave_eval=($info1{$key}->{eval})/5;
                if ($ave_eval < 1E-30) {
                        print spec_5 "$info1{$key}->{info}\n";
                        }
        }
}
open spec_5, "blast_best_pro_tra_5_spec.txt";
open fil1, ">temp";
while (<spec_5>) {
    chomp;
    @a=split;
    if ($info{$a[0]}) {
        $info="$info{$a[0]}->{infor}"."\t"."$_"."\t";
        $score=$info{$a[0]}->{score} + $a[-1];
        $info{$a[0]}={
            score => $score,
            infor => $info
        };
    }
    else{$info{$a[0]}={
        score => $a[-1],
        infor => $_
    };}
}
foreach $key (keys %info) {
    $value=$info{$key}->{infor};
    $score=$info{$key}->{score};
    print fil1 "$value\t$score\n";
}
open fil1, "temp";
while (<fil1>) {
    chomp;
    @a=split;
    next if /^\s+$/;
    s/^\s+//;
        ($name)=$a[0]=~/(sp)\|.*/;
    if ($score{$name}) {
        if ($score{$name} < $a[-1]) {
            $infor{$name}=$_;
        $score{$name}=$a[-1];
        }
    }
    else {
        $score{$name}=$a[-1];
        $infor{$name}=$_;}
}
$final=$infor{$name};
$final=~s/\t\d+$//;
$final=~s/^\s+//;
$final=~s/\tsp/\nsp/g;
print "$final\n";
