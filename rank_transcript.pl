#!/usr/bin/perl -w
open "fil", "OG0000001.blastp_result";
my (%hash1, %hash2, %hash3, %infor1, %infor2);
while (<fil>) {
        chomp;
        my @a=split;
        $hash1{$a[0]}++;
        $hash2{$a[1]}++;
        my $spec=$a[0]=~/(.*)\_\d+/;
        $hash3{$spec}++;
        if ((defined $infor1{$a[0]}->{$a[1]}->{score}>=$a[-2]) || (!defined $infor1{$a[0]}->{$a[1]}->{score})) {
                $infor1{$a[0]}->{$a[1]}={
                        infor1 => $_,
                        score => $a[-2]
        };
        my $info;
        for (my $i = 2; $i < @a; $i++) {
                $info.=$a[$i]."\t";
        }
        $info=~s/\s+$//;
        $info=$a[1]."\t".$a[0]."\t".$info;
        $infor2{$a[1]}->{$spec}->{$a[0]}={
                infor2 => $info,
                score => $a[-2]
                };
        }
} 

foreach my $pro (keys %hash2) {
        foreach my $spec (keys %hash3) {
                foreach my $tran (sort {($infor2{$pro}->{$spec}->{$a}->{score}||0) <=> ($infor2{$pro}->{$spec}->{$b}->{score}||0)} keys %infor1){
                        if (defined  $infor2{$pro}->{$spec}->{$tran}->{infor2}) {
                                $i++;
                                print "$infor2{$pro}->{$tran}->{infor2}\t$i\n";
                        }
                }
                $i=0;
        }
}
