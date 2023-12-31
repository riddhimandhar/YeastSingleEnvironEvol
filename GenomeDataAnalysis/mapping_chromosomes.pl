#!/usr/bin/perl -w 
#
# Mapping chromosome id and sequence 

$flag=0; 
open(FP,"yeast_ref_genome/S288C_reference_genome_R64-2-1_20150113/S288C_reference_sequence_R64-2-1_20150113.fsa") or die; 
open(WR,">../results/S288C_chromosome_id_mapping.txt") or die; 
while($fp=<FP>)
{
   chomp($fp);  
   if($fp=~/^\>/) 
   {
      @arr=split(/\s+/,$fp); 
      @qw=split(/\>/,$arr[0]); 
      $id=$qw[1]; 
      chop($arr[5]); 
      @we=split(/\=/,$arr[5]); 
      print "$id $we[1]\n"; 
      if($we[1] eq "mitochondrion") { $we[1]="M"; } 
      $chr="Chr".$we[1]; 
      undef(@qw); 
      undef(@we); 
      undef(@arr); 
      if($flag==0) 
      {
	$seq=""; 
	$flag=1; 
      }
      else
      {
	$len=length($seq); 
        print WR "$previd\t$prevchr\t$len\t$seq\n"; 	
	undef $len; 
	##print "$seq\n"; 
	$seq=""; 
      }
      $prevchr=$chr; 
      $previd=$id; 
      undef $chr; 
      undef $id; 
      next;  
   }    
   chop($fp); 
   $seq.=$fp; 
}
close(FP); 

$len=length($seq); 
print WR "$previd\t$prevchr\t$len\t$seq\n"; 	
undef $len; 
undef $seq; 
undef $previd; 
undef $prevchr; 

close(WR); 

