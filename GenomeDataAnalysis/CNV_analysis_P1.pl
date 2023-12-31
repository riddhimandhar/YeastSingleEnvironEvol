#!/usr/bin/perl -w 
#
# Check for Copy Number Variations 
# Calculating normalized log-ratios 

use warnings;

use Statistics::R;


$file[0]="R0-1_coverage.txt"; 
$file[1]="R0-2_coverage.txt"; 
$file[2]="R0-3_coverage.txt"; 
$file[3]="R100-1_coverage.txt"; 
$file[4]="R100-2_coverage.txt"; 
$file[5]="R100-3_coverage.txt"; 

print "Calculating geometric mean ...\n"; 

for($i=0;$i<6;$i++)
{
   $pp=0; 

   open(FP,"../results/coverage/$file[$i]") or die; 
   while($fp=<FP>)
   {
      chomp($fp); 
      @arr=split(/\s+/,$fp); 
      if($i==0) 
      {
	$pos[$pp][0]=$arr[0];
	$pos[$pp][1]=$arr[1];
      }
      $covlist[$pp][$i]=$arr[2];         

      undef(@arr); 
      $pp++; 
      ##if($pp==100000) { last; } 
   }
   close(FP); 
}

undef(@file); 

print "Calculating ratios ...\n"; 

open(WR,">Ratios_columns.txt") or die; 
for($i=0;$i<$pp;$i++)
{
   $val=1; 
   for($j=0;$j<6;$j++)
   {
      $val*=$covlist[$i][$j];
   }
   $val=($val**(1/6)); 

   for($j=0;$j<6;$j++)
   {
     if($val!=0) 
     {
        $rat=sprintf("%0.2f",$covlist[$i][$j]/$val); 
     }
     else { $rat="NA"; } 

     print WR "$rat\t"; 
     undef $rat; 
   }
   print WR "\n"; 
   undef $val; 
}
close(WR); 

$R=Statistics::R->new();
$R->start_sharedR;

$R->run(qq'tab=read.table(\"Ratios_columns.txt\")');
$R->run(qq'med1=median(na.omit(tab\$V1))'); 
$R->run(qq'med2=median(na.omit(tab\$V2))'); 
$R->run(qq'med3=median(na.omit(tab\$V3))'); 
$R->run(qq'med4=median(na.omit(tab\$V4))'); 
$R->run(qq'med5=median(na.omit(tab\$V5))'); 
$R->run(qq'med6=median(na.omit(tab\$V6))'); 

$med[0]=$R->get('med1');
$med[1]=$R->get('med2');
$med[2]=$R->get('med3');
$med[3]=$R->get('med4');
$med[4]=$R->get('med5');
$med[5]=$R->get('med6');

$R->run(q'rm(med6)');
$R->run(q'rm(med5)');
$R->run(q'rm(med4)');
$R->run(q'rm(med3)');
$R->run(q'rm(med2)');
$R->run(q'rm(med1)');
$R->run(q'rm(tab)');

$R->stopR; 

##system "rm Ratios_columns.txt"; 

print "Writing Normalized log2 coverage ratios ...\n"; 

open(WR,">../results/CNV_analysis/Normalized_coverage.txt") or die; 
open(WR2,">../results/CNV_analysis/Normalized_log2_coverage_ratio.txt") or die; 

for($i=0;$i<$pp;$i++)
{
   print WR "$pos[$i][0]\t$pos[$i][1]\t"; 
   print WR2 "$pos[$i][0]\t$pos[$i][1]\t"; 

   for($j=0;$j<6;$j++)
   {
     $rat=sprintf("%0.2f",$covlist[$i][$j]/$med[$j]); 
     print WR "$rat\t"; 
     undef $rat; 
   }
   print WR "\n"; 

   for($j=3;$j<6;$j++)
   {
     $rat1=sprintf("%0.2f",$covlist[$i][$j]/$med[$j]); 
     $rat0=sprintf("%0.2f",$covlist[$i][$j-3]/$med[$j-3]); 

     if($rat0!=0 && $rat0 ne "NA") 
     {
        $crat=sprintf("%0.2f",log(1+$rat1/$rat0)/log(2)); 
     }
     else { $crat="NA"; } 
     print WR2 "$crat\t"; 
     undef $crat; 
   }
   print WR2 "\n"; 
}

close(WR2); 
close(WR); 

undef(@med); 

undef(@pos);
undef(@covlist); 
