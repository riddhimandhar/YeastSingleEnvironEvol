#!/usr/bin/perl -w 
#
# Identify CNVs from normalized log2 coverage ratios 
# Value of 1 = no change in copy number 
# Coverage ratio in sliding window 
# V2 

use warnings;

use Statistics::R;

$R=Statistics::R->new();
$R->start_sharedR;

open(FP,"../results/CNV_analysis/Normalized_log2_coverage_ratio.txt") or die; 
open(WR1,">../results/CNV_analysis/Regions_N1_log2covRatio.txt") or die; 
open(WR2,">../results/CNV_analysis/Regions_N2_log2covRatio.txt") or die; 
open(WR3,">../results/CNV_analysis/Regions_N3_log2covRatio.txt") or die; 

$cn=0; 
while($fp=<FP>) 
{
   chomp($fp); 
   @sp=split(/\s+/,$fp); 
   $covlist[$cn][0]=$sp[0];
   $covlist[$cn][1]=$sp[1];
   $covlist[$cn][2]=$sp[2];
   $covlist[$cn][3]=$sp[3];
   $covlist[$cn][4]=$sp[4];
   $cn++; 
}
close(FP); 

$pp=0; 
$pp1=0; $pp2=0; $pp3=0; 
$prevchr = "ChrI"; 

$k=1;
$pos=$k; 
$retk=0; 
open(TM1,">TMP1.txt") or die; 
open(TM2,">TMP2.txt") or die; 
open(TM3,">TMP3.txt") or die; 

while($k<$cn)
{
   if($pp==1000) 
   {
      close(TM1);
      close(TM2);
      close(TM3);

      if($pp1>500) 
      {
	#$R->set('ARR1',\@list1);       
        $R->run(qq'ARR1=read.table("TMP1.txt")');       
	$R->run(qq'mc1=mean(ARR1\$V1)');
        $mc1=$R->get('mc1'); 	
      }
      else { $mc1="NA";  } 

      if($pp2>500) 
      {
	#$R->set('ARR2',\@list2);       
        $R->run(qq'ARR2=read.table("TMP2.txt")');       
	$R->run(qq'mc2=mean(ARR2\$V1)');
        $mc2=$R->get('mc2'); 	
      }
      else { $mc2="NA";  } 

      if($pp3>500) 
      {
	#$R->set('ARR3',\@list3);       
        $R->run(qq'ARR3=read.table("TMP3.txt")');       
	$R->run(qq'mc3=mean(ARR3\$V1)');
        $mc3=$R->get('mc3'); 	
      }
      else { $mc3="NA"; } 

      $R->run(q'rm(ARR1)');
      $R->run(q'rm(ARR2)');
      $R->run(q'rm(ARR3)');
      $R->run(q'rm(mc1)');
      $R->run(q'rm(mc2)');
      $R->run(q'rm(mc3)');

      ## Cut-off 0.0025 
      print "$prevchr $pos $mc1 $mc2 $mc3\n"; 
      
      if($mc1 ne "NA" && ($mc1>=1.5 || $mc1<=0.2)) 
      {
	 print WR1 "$prevchr $covlist[$k][1] $mc1\n"; 
      }
      if($mc2 ne "NA" && ($mc2>=1.5 || $mc2<=0.2)) 
      {
	 print WR2 "$prevchr $covlist[$k][1] $mc2\n"; 
      }
      if($mc3 ne "NA" && ($mc3>=1.5 || $mc3<=0.2)) 
      {
	 print WR3 "$prevchr $covlist[$k][1] $mc3\n"; 
      }

      if($retk==0)
      {
	 $k+=100; 
      }
      else
      {
         $k=$retk;
      }
      if($k>=$cn) { last; } 

      $retk=0; 
      $prevchr=$covlist[$k][0]; 
      $pos=$covlist[$k][1]; 

      $pp=0; 
      $pp1=0; 
      $pp2=0; 
      $pp3=0; 
      undef(@lpos1); 
      undef(@lpos2); 
      undef(@lpos3); 
      #undef(@list1); 
      #undef(@list2); 
      #undef(@list3); 
      open(TM1,">TMP1.txt") or die; 
      open(TM2,">TMP2.txt") or die; 
      open(TM3,">TMP3.txt") or die; 
   }

   if($covlist[$k+$pp][0] ne $prevchr) 
   { 
     $retk=$k+$pp;
     $pp=1000; 
     $pp1=0; $pp2=0; $pp3=0; 
     #undef(@list1);
     #undef(@list2);
     #undef(@list3);
     #close(TM1);
     #close(TM2);
     #close(TM3);
     redo; 
   } 

   if($covlist[$k+$pp][2] ne "NA") 
   {
     $lpos1[$pp1]=$pos;
     print TM1 "$covlist[$k+$pp][2]\n"; 
     #$list1[$pp1]=$covlist[$pos][2];  
     $pp1++; 
   }

   if($covlist[$k+$pp][3] ne "NA") 
   {
     $lpos2[$pp2]=$pos;
     #$list2[$pp2]=$covlist[$pos][3];  
     print TM2 "$covlist[$k+$pp][3]\n"; 
     $pp2++; 
   }

   if($covlist[$k+$pp][4] ne "NA") 
   {
     $lpos3[$pp3]=$pos;
     #$list3[$pp3]=$covlist[$k][4];  
     print TM3 "$covlist[$k+$pp][4]\n"; 
     $pp3++; 
   }

   $prevchr=$covlist[$k+$pp][0]; 
   $pp++; 
   $pos++; 

   if($k+$pp>=$cn) { $pp=1000; redo; } 
   #if($pp==10000) { last; } 
}
close(WR3); 
close(WR2); 
close(WR1); 
close(FP); 
close(TM1);
close(TM2);
close(TM3);

undef(@covlist); 

$R->stopR; 
