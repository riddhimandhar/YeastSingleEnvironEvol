#!/usr/bin/perl -w 
#
# Calculate coverage across the genome 
#

$file[0]="R0-1_mutations_consPE.txt"; 
$file[1]="R0-2_mutations_consPE.txt"; 
$file[2]="R0-3_mutations_consPE.txt"; 
$file[3]="R100-1_mutations_consPE.txt"; 
$file[4]="R100-2_mutations_consPE.txt"; 
$file[5]="R100-3_mutations_consPE.txt"; 

open(GR,">../results/coverage/Average_coverage.txt") or die; 

for($gi=0;$gi<6;$gi++)
{
   print "$file[$gi]\n"; 
   open(FP,"../results/S288C_chromosome_id_mapping.txt") or die; 
   while($fp=<FP>) 
   {
      chomp($fp); 
      @arr=split(/\s+/,$fp); 
      for($i=1;$i<=$arr[2];$i++)
      {
         $coverage{$arr[1]}[$i]=0;  
      }
      undef(@arr); 
   }
   close(FP); 

   @we=split(/\_/,$file[$gi]); 
   print GR "$we[0]\n"; 
   $wrf=$we[0]."_coverage.txt"; 
   undef(@we); 

   open(FP,"../results/Mutations_consolidatedPE/$file[$gi]") or die; 
   while($fp=<FP>)
   {
      chomp($fp); 
      @arr=split(/\t/,$fp); 
      $chr=$arr[1]; 
      @st=split(/\,/,$arr[2]); 
      @en=split(/\,/,$arr[3]); 
      $nn=@arr; 

      $no=0;
      if($nn==7)
      {
        @qw=split(/\,/,$arr[6]); 
        $no=@qw; 
      }
      undef $nn; 

      $cn=0; $hq=0;  
      for($k=0;$k<$no;$k++)
      {
         @we=split(/\s+/,$qw[$k]); 
	 ##print "$we[0] $we[1] $we[2] $we[4]\n"; 
	 if($we[3] eq "NA" || $we[3]>=30 || $we[4]==2) { $hq++; } ## 
	 if($we[4]==2) { $cn++; } 
	 undef(@we); 
      }
      if($no<=5) ## || $cn/$no==1) ## Filtering on number of mutations  
      #if($hq<=4) ## || $cn/$no==1) 
      {
	 if(!exists $dupid{$arr[1]."*$arr[2]*$arr[3]"})
	 {
	   for($m=0;$m<2;$m++)
	   {
             for($k=$st[$m];$k<=$en[$m];$k++) 
	     {
               $coverage{$chr}[$k]++; 
	     }
           }	  
	   $dupid{$arr[1]."*$arr[2]*$arr[3]"}=1; 
	 } 
      }
      undef $no; 
      undef(@qw); 
      undef(@en); 
      undef(@st); 
      undef $chr;  
      undef(@arr); 
   }
   close(FP); 
   undef %dupid; 

   open(WR,">../results/coverage/$wrf") or die; 
   open(FP,"../results/S288C_chromosome_id_mapping.txt") or die; 
   $nuclcov=0; $nnc=0; 
   while($fp=<FP>) 
   {
      chomp($fp); 
      @arr=split(/\s+/,$fp); 
      $avcov=0; 
      for($i=1;$i<=$arr[2];$i++)
      {
         print WR "$arr[1]\t$i\t$coverage{$arr[1]}[$i]\n";  
	 $avcov+=$coverage{$arr[1]}[$i]; 
	 if($arr[1] ne "ChrM") 
	 {
	    $nuclcov+=$coverage{$arr[1]}[$i]; 
	 }
      }
      if($arr[1] ne "ChrM") 
      {
         $nnc+=$arr[2]; 
      }
      $avcov=sprintf("%0.2f",$avcov/$arr[2]); 
      print GR "  $arr[1]\t$avcov\n"; 
      undef(@arr); 
   }
   $nuclcov=sprintf("%0.2f",$nuclcov/$nnc); 
   print GR "Avg Nucl Cov: $nuclcov\n"; 
   close(FP); 
   close(WR); 
}

undef(@file); 

close(GR); 
