open fil, "Trimmomatic.result-remove-unneed.txt";
while (<fil>) {
        chomp;
        @a=split;
        ($spe, $num)=$a[0]=~/(\w+?)(\d+)\_\w+/;
        $hash1{$spe}++;
        $hash2{$num}++;
        push @spe, $spe if $hash1{$spe}==1;
        push @num, $num if $hash2{$num}==1;
        if ($info{$spe}->{$num}) {
                $tot_read_num=$info{$spe}->{$num}->{'tot_read_num'}+$a[1];
                $sur_read_num=$info{$spe}->{$num}->{'sur_read_num'}+$a[2];
                $dro_read_num=$info{$spe}->{$num}->{'dro_read_num'}+$a[-2];
                $dro_ratio=($dro_read_num/$tot_read_num)*100;
                $info{$spe}->{$num}={
                                tot_read_num => $tot_read_num,
                                sur_read_num => $sur_read_num,
                                dro_read_num => $dro_read_num,
                                dro_ratio => $dro_ratio
                        };
                } else {
                        $info{$spe}->{$num}={
                                tot_read_num => $a[1],
                                sur_read_num => $a[2],
                                dro_read_num => $a[-2],
                                dro_ratio => ($a[-2]/$a[1])*100
                        };
                }
}

foreach $spe (sort {$a cmp $b} @spe) {
        foreach $num (sort {$a <=> $b} @num) {
                if ($info{$spe}->{$num}) {
                        $id=$spe.$num;
                        $tot_read_num=$info{$spe}->{$num}->{'tot_read_num'};
                        $sur_read_num=$info{$spe}->{$num}->{'sur_read_num'};
                        $dro_read_num=$info{$spe}->{$num}->{'dro_read_num'};
                        $dro_ratio=$info{$spe}->{$num}->{'dro_ratio'};
                        $dro_ratio=sprintf "%.2f", $dro_ratio;
                        print "$id\t$tot_read_num\t$sur_read_num\t$dro_read_num\t$dro_ratio\n";
                }
        }
}
