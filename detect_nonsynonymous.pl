# identified the synonymous or nonsynonymous snp by the vcf file
#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my $opt = GetOptions ('vcf:s' => \my $vcf,
        'ref:s' => \my $ref
        );
if (! ($vcf && $ref)) {
        print STDERR "\nExample usage:\n";
    print STDERR "\n$0 --vcf \"the vcf file\" --ref \"the reference fasta file\" \n\n";
    exit;
}

my %ref;
my $orf_id="";

open REF, "$ref" or die "cannot open $ref\n";
while (<REF>) {
        chomp;
        #       my $orf_id="";
        if (/>/) {
                s/>//;
                my @a=split;
                $orf_id=$a[0];
        } else {
                $ref{$orf_id}.=$_;
        }
}

my $new_vcf="new."."$vcf";
`grep -v "^##" $vcf|grep -i "pass" > $new_vcf`;
open VCF, "$new_vcf" or die "cannot open $new_vcf\n";
while (<VCF>) {
        chomp;
        my @info=split;
        my $orf_id=$info[0];
        my $snp_pos=$info[1];
        my $ref_nuc=$info[3];
        my $snp_nuc;
        if ($info[4]=~/\,/) {
                next;
        } else {
                $snp_nuc=$info[4];
        }
        # confirm the location of snp in the coden
    my $coden_pos=&confirm_snp_location_coden($snp_pos);
    # get the old coden and the new coden in the corresponding position
    my ($old_coden, $new_coden)=&get_coden_nuc_with_snp($snp_pos, $coden_pos, $orf_id, $snp_nuc);
    # translate the olde coden and the new coden
    my $old_coden_pep=&translate_nuc_to_protein($old_coden);
    my $new_coden_pep=&translate_nuc_to_protein($new_coden);
    my $synonymous;
    if ($old_coden_pep eq $new_coden_pep) {
        $synonymous="synonymous";
    } else {
        $synonymous="nonsynonymous";
    }
    print "$orf_id\t$snp_pos\t$ref_nuc\t$snp_nuc\t$old_coden\t$new_coden\t$old_coden_pep\t$new_coden_pep\t$synonymous\n";
}
`rm $new_vcf`;

# confirm the position of snp in the coden
sub confirm_snp_location_coden {
        my ($snp_pos)=@_;
        my $remainder = ($snp_pos) % 3;
    my $codon_pos; #position in codon (1st, 2nd, or 3rd)
    if ($remainder == 0){
        $codon_pos = 3;
      }elsif ($remainder == 1){
        $codon_pos = 1;
      }elsif ($remainder == 2){
        $codon_pos = 2;
      }else{
        die "Impossible codon position\n";
      }
    return $codon_pos;
}
sub get_coden_nuc_with_snp {
        my ($snp_pos, $coden_pos, $orf_id, $snp_nuc)=@_;
        $snp_pos=$snp_pos-1;
        my ($pos1, $pos2, $pos3, $nt1, $nt2, $nt3, $old_coden, $new_coden);
        die "there is no $orf_id in hash" if ! $ref{$orf_id};
        if($coden_pos == 3){ # snp is the end nuc of the coden
        $pos1 = $snp_pos-2;
        $pos2 = $snp_pos-1;
        $pos3 = $snp_pos;
        $nt1 = substr($ref{$orf_id}, $pos1 , 1 );
        $nt2 = substr($ref{$orf_id}, $pos2 , 1 );
        $nt3 = substr($ref{$orf_id}, $pos3 , 1 );
        $old_coden = ($nt1.$nt2.$nt3);
        $new_coden = ($nt1.$nt2.$snp_nuc);
      }elsif($coden_pos == 1){ # snp is the first nuc of the coden
        $pos1 = $snp_pos;
        $pos2 = $snp_pos+1;
        $pos3 = $snp_pos+2;
        $nt1 = substr($ref{$orf_id}, $pos1 , 1 );
        $nt2 = substr($ref{$orf_id}, $pos2 , 1 );
        $nt3 = substr($ref{$orf_id}, $pos3 , 1 );
        $old_coden = ($nt1.$nt2.$nt3);
        $new_coden = ($snp_nuc.$nt2.$nt3);
      }elsif($coden_pos == 2){ # snp is the second nuc of the coden
        $pos1 = $snp_pos-1;
        $pos2 = $snp_pos;
        $pos3 = $snp_pos+1;
        $nt1 = substr($ref{$orf_id}, $pos1 , 1 );
        $nt2 = substr($ref{$orf_id}, $pos2 , 1 );
        $nt3 = substr($ref{$orf_id}, $pos3 , 1 );
        $old_coden = ($nt1.$nt2.$nt3);
        $new_coden = ($nt1.$snp_nuc.$nt3);
      }else{
        die "The remainder is not a multiple of 1/3\n";
      }
      return ($old_coden, $new_coden);
}
sub translate_nuc_to_protein {
        my ($codon)=@_;
        my $aa;
        if($codon eq "GCT" or $codon eq "GCC" or $codon eq "GCA" or $codon eq "GCG"){
        $aa = "A";
      }elsif($codon eq "CGT" or $codon eq "CGC" or $codon eq "CGA" or $codon eq "CGG" or $codon eq "AGA" or $codon eq "AGG"){
        $aa = "R";
      }elsif($codon eq "AAT" or $codon eq "AAC"){
        $aa = "N";
      }elsif($codon eq "GAT" or $codon eq "GAC"){
        $aa = "D";
      }elsif($codon eq "TGT" or $codon eq "TGC"){
        $aa = "C";
      }elsif($codon eq "CAA" or $codon eq "CAG"){
        $aa = "Q";
      }elsif($codon eq "GAA" or $codon eq "GAG"){
        $aa = "E";
      }elsif($codon eq "GGT" or $codon eq "GGC" or $codon eq "GGA" or $codon eq "GGG"){
        $aa = "G";
      }elsif($codon eq "CAT" or $codon eq "CAC"){
        $aa = "H";
      }elsif($codon eq "ATT" or $codon eq "ATC" or $codon eq "ATA"){
        $aa = "I";
      }elsif($codon eq "TTA" or $codon eq "TTG" or $codon eq "CTT" or $codon eq "CTC" or $codon eq "CTA" or $codon eq "CTG"){
        $aa = "L";
      }elsif($codon eq "AAA" or $codon eq "AAG"){
        $aa = "K";
      }elsif($codon eq "ATG"){
        $aa = "M";
      }elsif($codon eq "TTT" or $codon eq "TTC"){
        $aa = "F";
      }elsif($codon eq "CCT" or $codon eq "CCC" or $codon eq "CCA" or $codon eq "CCG"){
        $aa = "P";
      }elsif($codon eq "TCT" or $codon eq "TCC" or $codon eq "TCA" or $codon eq "TCG" or $codon eq "AGT" or $codon eq "AGC"){
        $aa = "S";
      }elsif($codon eq "ACT" or $codon eq "ACC" or $codon eq "ACA" or $codon eq "ACG"){
        $aa = "T";
      }elsif($codon eq "TGG"){
        $aa = "W";
      }elsif($codon eq "TAT" or $codon eq "TAC"){
        $aa = "Y";
      }elsif($codon eq "GTT" or $codon eq "GTC" or $codon eq "GTA" or $codon eq "GTG"){
        $aa = "V";
      }elsif($codon eq "TAA" or $codon eq "TGA" or $codon eq "TAG"){
        $aa = "STOP";
      }else{
        print STDERR "The reference codon cannot be translated\n";
        $aa = "X";
      }
      return $aa;
}
