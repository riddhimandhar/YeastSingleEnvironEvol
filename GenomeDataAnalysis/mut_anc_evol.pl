#!/usr/bin/perl -w 
#

$file[0]="../results/coverage/R0-1_coverage.txt";
$file[1]="../results/coverage/R0-2_coverage.txt";
$file[2]="../results/coverage/R0-3_coverage.txt";
$file[3]="../results/coverage/R100-1_coverage.txt";
$file[4]="../results/coverage/R100-2_coverage.txt";
$file[5]="../results/coverage/R100-3_coverage.txt";

for($i=0;$i<6;$i++)
{
   open(FP,$file[$i]) or die; 
   while($fp=<FP>)
   {
      chomp($fp); 
      @arr=split(/\s+/,$fp);
      $id=$arr[0]."*$arr[1]"; 
      if(exists $covlist{$id}) 
      {
        $covlist{$id}[1][$i]=$arr[2]; 
      }
      else 
      {
	for($j=0;$j<6;$j++)
	{
	   $covlist{$id}[1][$j]=0; 
	}
        $covlist{$id}[1][$i]=$arr[2]; 
      }
      undef $id; 
      undef(@arr); 
   }
   close(FP); 
}
undef(@file); 


$file[0]="../results/FinalMutations/R0-1_finalMutations.txt";
$file[1]="../results/FinalMutations/R0-2_finalMutations.txt";
$file[2]="../results/FinalMutations/R0-3_finalMutations.txt";
$file[3]="../results/FinalMutations/R100-1_finalMutations.txt";
$file[4]="../results/FinalMutations/R100-2_finalMutations.txt";
$file[5]="../results/FinalMutations/R100-3_finalMutations.txt";

for($i=0;$i<6;$i++)
{
   open(FP,$file[$i]) or die; 
   while($fp=<FP>)
   {
      chomp($fp); 
      @arr=split(/\s+/,$fp);
      $id=$arr[0]."*$arr[1]*$arr[2]*$arr[3]"; 
      if(exists $list{$id}) 
      {
        $list{$id}[1][$i]=$arr[4]; 
        $list{$id}[2][$i]=$arr[5]; 
        $list{$id}[3][$i]=$arr[6]; 
      }
      else 
      {
	for($j=0;$j<6;$j++)
	{
	   $list{$id}[1][$j]=0; 
	   if(exists $covlist{$arr[0]."*$arr[1]"}) 
	   {
	     $list{$id}[2][$j]=$covlist{$arr[0]."*$arr[1]"}[1][$j]; 
	   }
	   else { $list{$id}[2][$j]=0; }
	   $list{$id}[3][$j]=0; 
	}
        $list{$id}[1][$i]=$arr[4]; 
        $list{$id}[2][$i]=$arr[5]; 
        $list{$id}[3][$i]=$arr[6]; 
      }
      undef $id; 
      undef(@arr); 
   }
   close(FP); 
}

undef(@file); 

undef %covlist; 


open(WR,">../results/FinalMutations/Mutation_Staistics_All.txt") or die; 
open(WR2,">../results/FinalMutations/Mutation_Staistics_AncOnly.txt") or die; 
open(WR3,">../results/FinalMutations/Mutation_Staistics_EvoOnly.txt") or die; 

foreach $key (sort keys %list)
{
   @arr=split(/\*/,$key); 

   $fl=0;

   if($list{$key}[1][3]==0 && $list{$key}[1][4]==0 && $list{$key}[1][5]==0) { $fl=1; } ## AncOnly
   if($list{$key}[1][0]==0 && $list{$key}[1][1]==0 && $list{$key}[1][2]==0) { $fl=2; } ## EvoOnly

   print WR "$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\t"; 
   if($fl==1) 
   { 
      print WR2 "$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\t"; 
   }
   if($fl==2) 
   { 
      print WR3 "$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\t"; 
   }
   for($i=0;$i<6;$i++)
   {
      print WR "$list{$key}[1][$i] $list{$key}[2][$i] $list{$key}[3][$i]  \t"; 
      if($fl==1) 
      { 
         print WR2 "$list{$key}[1][$i] $list{$key}[2][$i] $list{$key}[3][$i]  \t"; 
      }
      if($fl==2) 
      { 
         print WR3 "$list{$key}[1][$i] $list{$key}[2][$i] $list{$key}[3][$i]  \t"; 
      }
   }
   print WR "\n"; 
   if($fl==1) { print WR2 "\n"; }
   if($fl==2) { print WR3 "\n"; }
   undef(@arr);
}
undef %list; 

close(WR); 
