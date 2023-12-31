#!/usr/bin/perl -w 
#
# Check statistics regarding mutations 
#
# Mutations in common genes/common mutations 
# Transitions/Transversions Syn/Nonsyn


$file[0]="../results/Format_Circplot/N1_formatted.txt"; 
$file[1]="../results/Format_Circplot/N2_formatted.txt"; 
$file[2]="../results/Format_Circplot/N3_formatted.txt"; 

$frc=4.5; ## frequency cut-off 

open (WR,">../results/Basic_Statistics_Mutations.txt") or die; 

for($i=0;$i<3;$i++)
{
   $tot=0; 
   $syn=0; $nonsyn=0; $prom=0; $htot=0;  
   $miss=0; $mitoh=0; $mitoall=0;  
   $indel=0; $hindel=0;  $codindel=0; $hcodindel=0;  

   open(FP,$file[$i]) or die; 
   if($i==0) { $file[$i]="N1"; } 
   if($i==1) { $file[$i]="N2"; } 
   if($i==2) { $file[$i]="N3"; } 
   while($fp=<FP>)
   {
      chomp($fp); 
      if($fp=~/Chr\s+/) { next; } 
      @arr=split(/\s+/,$fp); 
      if($arr[4]!=0) { $tot++; } 
      if($arr[4]>=$frc) { $htot++; } 
      if($arr[2]!~/\-/ && $arr[3]!~/\-/ && $arr[4]>=$frc && $arr[5]=~/^Coding\;/) 
      {
	 @qw=split(/\;/,$arr[5]); 
         if($qw[2] eq "Syn") { $syn++; }
         if($qw[2] eq "Nonsyn") { $nonsyn++; }
         if($qw[2] eq "Nonsense") { $miss++; }
	 undef(@qw);  

      }

      if($arr[4]>=$frc) 
      {
	 if($arr[0] eq "ChrM") 
	 {
            $mitoh++; 
	 }

	 $lfl=0; 
         if($arr[5]=~/^Coding\;/) 
	 {
            $lfl=1;
	    @qw=split(/\;/,$arr[5]); 
	 }
	 elsif($arr[5]=~/^Coding\_/) 
	 {
            $lfl=1; 
	    @qw=split(/\./,$arr[5]); 
	 }

	 $id1=$arr[0]."*$arr[1]*$arr[2]*$arr[3]"; 

	 if(exists $mutlist1{$id1} && !exists $dup1{$file[$i]."*$id1"}) 
         {
	    $cnt=$mutlist1{$id1}[0]; 
	    $mutlist1{$id1}[1][$cnt]=$file[$i]; 
	    $mutlist1{$id1}[2][$cnt]=$arr[4]; 
            $mutlist1{$id1}[0]++; 
	    $dup1{$file[$i]."*$id1"}=1; 
	    undef $cnt; 
         }
	 elsif(!exists $mutlist1{$id1})
         {
	    $mutlist1{$id1}[1][0]=$file[$i]; 
	    $mutlist1{$id1}[2][0]=$arr[4]; 
	    $mutlist1{$id1}[3]=$arr[5]; 
	    $mutlist1{$id1}[0]=1; 
	    $dup1{$file[$i]."*$id1"}=1; 
         }

	 if($lfl==1)
	 {
           $id2=$qw[1]; 

	   if(exists $mutlist2{$id2} && !exists $dup2{$file[$i]."*$id2"}) 
           {
	     $cnt=$mutlist2{$id2}[0]; 
	     $mutlist2{$id2}[1][$cnt]=$file[$i]; 
	     $mutlist2{$id2}[2][$cnt]=$arr[4]; 
	     $mutlist2{$id2}[3][$cnt]="$arr[0];$arr[1];$arr[2];$arr[3]"; 
	     $mutlist2{$id2}[4][$cnt]=$arr[5]; 
             $mutlist2{$id2}[0]++; 
	     $dup2{$file[$i]."*$id2"}=1; 
	     undef $cnt; 
           }
	   elsif(!exists $mutlist2{$id2}) 
           {
	     $mutlist2{$id2}[1][0]=$file[$i]; 
	     $mutlist2{$id2}[2][0]=$arr[4]; 
	     $mutlist2{$id2}[3][0]="$arr[0];$arr[1];$arr[2];$arr[3]"; 
	     $mutlist2{$id2}[4][0]=$arr[5]; 
             $mutlist2{$id2}[0]=1; 
	     $dup2{$file[$i]."*$id2"}=1; 
           }

	   undef $id2; 
	 }
	 undef $id1;  
      }

      if(($arr[2]=~/\-/ || $arr[3]=~/\-/) && ($arr[4]!=0))
      {
	$indel++; 
	if($arr[4]>=$frc)
        {
	   $hindel++;  	
	   if($arr[5]=~/^Coding\_/) 
	   {
	      $hcodindel++; 
	   }
	}
	if($arr[5]=~/^Coding\_/) 
	{
	   $codindel++; 
	}
      }
      if($arr[0] eq "ChrM" && $arr[4]!=0) 
      {
	 $mitoall++; 
      } 
      elsif($arr[4]>=$frc && $arr[5]=~/Prom_/) { $prom++; }
      undef(@arr); 
   }

   print WR "$file[$i] Total: $tot HighFreqTot: $htot Nonsyn: $nonsyn Syn: $syn Prom: $prom Nonsense (not included in missense): $miss\nMitoAll: $mitoall  MitoHigh: $mitoh\nTotIndel: $indel HighFreqIndel: $hindel CodIndel: $codindel HighFreqCodIndel: $hcodindel\n\n"; 
}

undef %dup1; 
undef %dup2; 
undef(@file); 

open(WR,">../results/Common_Unique_mutations_Evoonly.txt") or die; 
foreach $key (sort keys %mutlist1)
{
   @qw=split(/\*/,$key);
   print WR "$qw[0] $qw[1] $qw[2] $qw[3]\t$mutlist1{$key}[0]\t"; 
   for($i=0;$i<$mutlist1{$key}[0];$i++)
   {
      print WR "$mutlist1{$key}[1][$i] $mutlist1{$key}[2][$i],"; 
   } 
   print WR " $mutlist1{$key}[3]\n"; 
   undef(@qw);    
}
close(WR); 


open(WR2,">../results/Common_Unique_genes_mutated_Evoonly.txt") or die; 
foreach $key (sort keys %mutlist2)
{
   print WR2 "$key\t$mutlist2{$key}[0]\t"; 
   for($i=0;$i<$mutlist2{$key}[0];$i++)
   {
      print WR2 "$mutlist2{$key}[1][$i] $mutlist2{$key}[2][$i] $mutlist2{$key}[3][$i] $mutlist2{$key}[4][$i], "; 
   } 
   print WR2 "\n"; 
}
close(WR2); 

undef %mutlist1; 
undef %mutlist2; 
