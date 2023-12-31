#!/usr/bin/perl -w 
#
# Find mutations in all lines 

$file[0]="R0-1_mutations_consPE.txt"; 
$file[1]="R0-2_mutations_consPE.txt"; 
$file[2]="R0-3_mutations_consPE.txt"; 
$file[3]="R100-1_mutations_consPE.txt"; 
$file[4]="R100-2_mutations_consPE.txt"; 
$file[5]="R100-3_mutations_consPE.txt"; 

for($gi=0;$gi<6;$gi++) 
{
   print "$file[$gi]\n"; 
   @qw=split(/\_/,$file[$gi]); 
   $cov=$qw[0]."_coverage.txt";
   $wrf=$qw[0]."_finalMutations.txt"; 
   undef(@qw); 

   open(FP,"../results/coverage/$cov") or die; 
   while($fp=<FP>) 
   {
      chomp($fp); 
      @arr=split(/\s+/,$fp); 
      $covlist{$arr[0]."*$arr[1]"}=$arr[2];  
      undef(@arr); 
   }
   close(FP); 

   open(FP,"../results/Mutations_consolidatedPE/$file[$gi]") or die; 
   open(WR,">../results/FinalMutations/$wrf") or die; 
   while($fp=<FP>) 
   {
      chomp($fp); 
      @arr=split(/\t/,$fp); 
      $chr=$arr[1]; 
      $nn=@arr; 

      $no=0; 
      if($nn==7) 
      {
        @qw=split(/\,/,$arr[6]); 
        $no=@qw; 
      }
      else
      {
	undef(@arr);
	undef $chr; 
	undef $nn; 
	next; 
      }
      undef $nn; 

      $hq=0; 
      for($k=0;$k<$no;$k++)
      {
         @we=split(/\s+/,$qw[$k]); 
	 ##print "$we[0] $we[1] $we[2] $we[4]\n"; 
	 if($we[3] eq "NA" || $we[3]>=30 || $we[4]==2) { $hq++; } 
	 undef(@we); 
      }

      if($no<=5) ## Filtering on the number of mutations 
      {
      for($i=0;$i<$no;$i++)
      {
        @se=split(/\s+/,$qw[$i]); 
	if(($se[3] eq "NA" && $se[4]==2) || ($se[3] ne "NA" && $se[3]>=30) || ($se[3] ne "NA" && $se[3]>=20 && $se[4]==2)) 
	{
	   $id=$chr."*$se[0]*$se[1]*$se[2]"; 
	   if(exists $mutlist{$id}) 
	   {
	      $cnt=$mutlist{$id}[0]; 
	      $lfl=0;
	      for($k=0;$k<$cnt;$k++)
	      { 
                if((($mutlist{$id}[1][$k] eq $arr[2]) && ($mutlist{$id}[2][$k] eq $arr[3])))
		{
		   $lfl=1;
		}
	      } 
	      if($lfl==0)
	      {
		$mutlist{$id}[1][$cnt]=$arr[2];
		$mutlist{$id}[2][$cnt]=$arr[3];
	        $mutlist{$id}[0]++; 
	      }
	   }
	   else
	   {
              $mutlist{$id}[1][0]=$arr[2];
	      $mutlist{$id}[2][0]=$arr[3];
	      $mutlist{$id}[0]=1; 
	   }
	   undef $id; 
	}
	undef(@se); 
      }
      }
      undef $chr; 
      undef $no; 
      undef(@qw); 
      undef(@arr); 
   }
   close(FP); 

   $tot=0; 
   $hfr=0; $maxf=0;  
   foreach $key (sort keys %mutlist) 
   {
      if(($mutlist{$key}[0]>=10 && $key!~/\-/) || ($mutlist{$key}[0]>=10 && $key=~/\-/))
      {
	 $tot++; 
	 @qw=split(/\*/,$key); 
	 $id=$qw[0]."*$qw[1]"; 
	 if(exists $covlist{$id} && $covlist{$id}!=0) 
	 {
	    $freq=sprintf("%0.2f",$mutlist{$key}[0]/$covlist{$id}*100); 
	    if($freq>=5) { $hfr++; } 
	    if($freq>$maxf) { $maxf=$freq; } 
	 }
	 else
	 {
            $covlist{$id}="NA"; 
            $freq="NA"; 
	 }
	 print WR "$qw[0]\t$qw[1]\t$qw[2]\t$qw[3]\t$mutlist{$key}[0]\t$covlist{$id}\t$freq\n"; 
	 undef $freq; 
	 undef(@qw); 
      }
   }

   print "TOT: $tot   High_Freq: $hfr  Max Mut freq: $maxf\n"; 
   close(WR); 
   undef $wrf; 

   undef %mutlist; 
   undef %covlist; 
}

undef(@file); 
