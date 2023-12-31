#!/usr/bin/perl -w 
#
# Annotate genome information 
#

## Dubious list 
#

$file[0]="../results/FinalMutations/Mutation_Staistics_All.txt"; 
$file[1]="../results/FinalMutations/Mutation_Staistics_AncOnly.txt"; 
$file[2]="../results/FinalMutations/Mutation_Staistics_EvoOnly.txt"; 

$cnt=0; 
open(FP,"../datasets/Yeast_genome_annot_func.txt") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   @arr=split(/\t+/,$fp); 
   $gene=$arr[1]; 
   $chr=$arr[2]; 
   $st=$arr[3]; 
   $en=$arr[4]; 
   $strand=$arr[5]; 
   $intron=$arr[6]; 
   $seq=$arr[7]; 
   $func=$arr[8]; 

   if(exists $genelist{$chr}) 
   {
      $cnt=$genelist{$chr}[0]; 
      for($i=0;$i<$cnt;$i++)
      {
	 if($st<$genelist{$chr}[2][$i] && $en<$genelist{$chr}[3][$i]) 
	 {
            last;
	 }
      }
      for($k=$cnt;$k>$i;$k--)
      {
	 $genelist{$chr}[1][$k]=$genelist{$chr}[1][$k-1]; 
	 $genelist{$chr}[2][$k]=$genelist{$chr}[2][$k-1]; 
	 $genelist{$chr}[3][$k]=$genelist{$chr}[3][$k-1]; 
	 $genelist{$chr}[4][$k]=$genelist{$chr}[4][$k-1]; 
	 $genelist{$chr}[5][$k]=$genelist{$chr}[5][$k-1]; 
	 $genelist{$chr}[6][$k]=$genelist{$chr}[6][$k-1]; 
	 $genelist{$chr}[7][$k]=$genelist{$chr}[7][$k-1]; 
      }
      $genelist{$chr}[1][$i]=$gene; 
      $genelist{$chr}[2][$i]=$st; 
      $genelist{$chr}[3][$i]=$en; 
      $genelist{$chr}[4][$i]=$strand; 
      $genelist{$chr}[5][$i]=$intron; 
      $genelist{$chr}[6][$i]=$seq; 
      $genelist{$chr}[7][$i]=$func; 
      $genelist{$chr}[0]++; 
   }
   else
   {
      $genelist{$chr}[1][0]=$gene; 
      $genelist{$chr}[2][0]=$st; 
      $genelist{$chr}[3][0]=$en; 
      $genelist{$chr}[4][0]=$strand; 
      $genelist{$chr}[5][0]=$intron; 
      $genelist{$chr}[6][0]=$seq; 
      $genelist{$chr}[7][0]=$func; 
      $genelist{$chr}[0]=1; 
   }

   undef $gene; 
   undef $chr; 
   undef $st; 
   undef $en; 
   undef $strand; 
   undef $intron; 
   undef $seq; 
   undef $func; 
   undef(@arr); 
}
close(FP); 


open(WR,">../datasets/Yeast_genome_annot_func_sorted.txt") or die; 
foreach $key (sort keys %genelist) 
{
   for($k=0;$k<$genelist{$key}[0];$k++)
   {
      print WR "$genelist{$key}[1][$k]\t$key\t"; 
      for($m=2;$m<=7;$m++) 
      {
	  print WR "$genelist{$key}[$m][$k]\t"; 
      }
      print WR "\n"; 
   }
}
close(WR); 

for($mk=0;$mk<3;$mk++)
{
  @fle=split(/\//,$file[$mk]); 
  print "$fle[3]\n"; 
  @qq=split(/\_/,$fle[3]); 
  $wrf="../Mutation_annotations/Mutation_annotation_$qq[2]"; 
  undef(@qq);
  undef(@fle); 

  open(WR,">$wrf") or die; 
  undef $wrf; 

  open(FP,$file[$mk]) or die; 
  while($fp=<FP>)
  {
     chomp($fp); 
     @arr=split(/\s+/,$fp); 
     $chr=$arr[0];
     $pos=$arr[1]; 
     $mut=$arr[2]; 
     $anc=$arr[3]; 

     $flag=""; 
     $prevpos=0; 
     $prevgn="";
     $prevstr="";
     if(exists $genelist{$chr}) 
     {
        for($i=0;$i<$genelist{$chr}[0];$i++)
	{
	   $gn=$genelist{$chr}[1][$i]; 

           $strdir=$genelist{$chr}[4][$i]; 

           if($pos>=$genelist{$chr}[2][$i] && $pos<=$genelist{$chr}[3][$i])
	   {
	      $flag="Coding"; 
              if($genelist{$chr}[5][$i] ne "NA")
	      {
		 @qw=split(/\,/,$genelist{$chr}[5][$i]); 
		 $nn=@qw; 
		 for($k=0;$k<$nn;$k++)
		 {
		    @we=split(/\-/,$qw[$k]); 
		    ##print "$chr $pos $we[0] $we[1]\n"; 
	            if($pos>=$we[0] && $pos<=$we[1])
		    {
		       $flag="Intron.$gn"; 
		    }
		    undef(@we); 
		 }
		 undef $nn;
		 undef(@qw); 
	      }
	      if($flag eq "Coding" && $anc ne "-" && $mut ne "-") 
	      {
		 if($strdir eq "+") ## Gene in positive strand  
		 {
	            $cc=0; 
	            for($k=$genelist{$chr}[2][$i];$k<$genelist{$chr}[3][$i];$k++)
                    { 
		       $lfl=0; 
                       if($genelist{$chr}[5][$i] ne "NA")
	               {
		         @qw=split(/\,/,$genelist{$chr}[5][$i]); 
		         $nn=@qw; 
		         for($j=0;$j<$nn;$j++)
		         {
		            @we=split(/\-/,$qw[$j]); 
			    if($k>=$we[0] && $k<=$we[1]) 
			    {
			       $lfl=1;
			       undef(@we);  
			       last; 
			    }
			    undef(@we); 
			 }
			 undef $nn; 
			 undef(@qw); 
		       }
                       if($lfl==1) { next; }  

		       $cc++; 
		       if($k==$pos) 
		       {
			  $base=substr($genelist{$chr}[6][$i],$cc-1,1); 
                          $res=$cc % 3; 
			  if($res==0) 
			  { 
		             $oricod=substr($genelist{$chr}[6][$i],$cc-3,3); 
			     $mutcod=substr($genelist{$chr}[6][$i],$cc-3,2).$mut; 
			     $nextcod=substr($genelist{$chr}[6][$i],$cc,3);
			  }
			  if($res==1) 
			  { 
		             $oricod=substr($genelist{$chr}[6][$i],$cc-1,3); 
			     $mutcod=$mut.substr($genelist{$chr}[6][$i],$cc,2); 
			     $nextcod=substr($genelist{$chr}[6][$i],$cc+2,3);
			  }
			  if($res==2) 
			  { 
		             $oricod=substr($genelist{$chr}[6][$i],$cc-2,3); 
			     $mutcod=substr($genelist{$chr}[6][$i],$cc-2,1).$mut.substr($genelist{$chr}[6][$i],$cc,1); 
			     $nextcod=substr($genelist{$chr}[6][$i],$cc+1,3);
			  }

			  if($res==0) { $aar=$cc/3; } 
			  else { $aar=($cc-$res)/3+1; }

			  $val=checkcodons($mutcod,$oricod); 
			  @lo=split(/\*/,$val); 
			  $ex1=$lo[0].";$lo[2]$aar$lo[1]";
			  $ex2=$oricod."->$mutcod"; 
			  $flag.=";$genelist{$chr}[1][$i];$ex1;$ex2"; 
			  undef $ex1;
			  undef $ex2; 
			  undef(@lo);
			  #print "$cc $chr $genelist{$chr}[1][$i] $strdir $mut $anc $base $res $mutcod $oricod $nextcod $pos $flag $aar $genelist{$chr}[2][$i] $genelist{$chr}[3][$i]\t"; 
			  #for($u=0;$u<length($genelist{$chr}[6][$i])-3;$u+=3)
			  #{
			  #  $ss=substr($genelist{$chr}[6][$i],$u,3); 
			  #  print "$ss "; 
			  #}
			  #print "\n"; 
		       }    
		    }
		 }
		 elsif($strdir eq "-") 
		 {
	            $cc=0; 
	            for($k=$genelist{$chr}[3][$i];$k>$genelist{$chr}[2][$i];$k--)
                    { 
		       $lfl=0; 
                       if($genelist{$chr}[5][$i] ne "NA")
	               {
		         @qw=split(/\,/,$genelist{$chr}[5][$i]); 
		         $nn=@qw; 
		         for($j=0;$j<$nn;$j++)
		         {
		            @we=split(/\-/,$qw[$j]); 
			    if($k>=$we[0] && $k<=$we[1]) 
			    {
			       $lfl=1;
			       undef(@we);  
			       last; 
			    }
			    undef(@we); 
			 }
			 undef $nn; 
			 undef(@qw); 
		       }
                       if($lfl==1) { next; }  

		       $cc++; 
		       if($k==$pos) 
		       {
			  $base=substr($genelist{$chr}[6][$i],$cc-1,1); 

			  if($anc eq "A") { $canc="T"; }
			  if($anc eq "T") { $canc="A"; }
			  if($anc eq "G") { $canc="C"; }
			  if($anc eq "C") { $canc="G"; }
			  if($mut eq "A") { $cmut="T"; }
			  if($mut eq "T") { $cmut="A"; }
			  if($mut eq "G") { $cmut="C"; }
			  if($mut eq "C") { $cmut="G"; }

                          $res=$cc % 3; 
			  if($res==0) 
			  { 
		             $oricod=substr($genelist{$chr}[6][$i],$cc-3,3); 
			     $mutcod=substr($genelist{$chr}[6][$i],$cc-3,2).$cmut; 
			     $nextcod=substr($genelist{$chr}[6][$i],$cc,3);
			  }
			  if($res==1) 
			  { 
		             $oricod=substr($genelist{$chr}[6][$i],$cc-1,3); 
			     $mutcod=$cmut.substr($genelist{$chr}[6][$i],$cc,2); 
			     $nextcod=substr($genelist{$chr}[6][$i],$cc+2,3);
			  }
			  if($res==2) 
			  { 
		             $oricod=substr($genelist{$chr}[6][$i],$cc-2,3); 
			     $mutcod=substr($genelist{$chr}[6][$i],$cc-2,1).$cmut.substr($genelist{$chr}[6][$i],$cc,1); 
			     $nextcod=substr($genelist{$chr}[6][$i],$cc+1,3);
			  }

			  if($res==0) { $aar=$cc/3; } 
			  else { $aar=($cc-$res)/3+1; }

			  $val=checkcodons($mutcod,$oricod); 
			  @lo=split(/\*/,$val); 
			  $ex1=$lo[0].";$lo[2]$aar$lo[1]";
			  $ex2=$oricod."->$mutcod"; 
			  $flag.=";$genelist{$chr}[1][$i];$ex1;$ex2"; 
			  undef $ex1;
			  undef $ex2; 
			  undef(@lo);
			  #print "$cc $chr $genelist{$chr}[1][$i] $strdir $cmut $canc $base $mut $anc $res $mutcod $oricod $nextcod $pos $flag $aar $genelist{$chr}[2][$i] $genelist{$chr}[3][$i]\t"; 
			  #for($u=0;$u<length($genelist{$chr}[6][$i])-3;$u+=3)
			  #{
			  #  $ss=substr($genelist{$chr}[6][$i],$u,3); 
			  #  print "$ss "; 
			  #}
			  #print "\n"; 
		       }    
		    }
		 }
	      }
	      elsif($flag eq "Coding")
	      {
		 $flag="Coding_Indel.$gn"; 
	      }
	      last; 
	   }
	   else
	   {
	      ## Intergenic noncoding mutations - check whether they fall in the promoter region of a gene 
	      #$flag="";
	      $intg=""; 
	      if($pos>$prevpos && $pos<$genelist{$chr}[2][$i])
	      {
                 if($genelist{$chr}[5][$i] eq "+" && $genelist{$chr}[2][$i]-$pos<=1000) { $intg="Prom_$gn,"; } 
                 if($prevstr eq "-" && $pos-$prevpos<=1000) { $intg="Prom_$prevgn,"; } 
	      }
	      $flag.=$intg; 
	      $prevpos=$genelist{$chr}[3][$i]; 
	      $prevgn=$genelist{$chr}[1][$i]; 
	      $prevstr=$genelist{$chr}[4][$i]; 
	   }
	   undef $gn; 
	   undef $strdir; 
	}
     }

     if($flag eq "") { $flag="X"; } 
     print WR "$fp $flag\n"; 
     undef $flag; 
     undef $chr; 
     undef $pos; 
     undef $mut; 
     undef $anc; 
     undef(@arr); 
  }
  close(WR); 
  close(FP); 
}

undef(@file);

undef(@genelist); 

sub checkcodons 
{
  local($lmut,$lori)=@_; 

  local ($lch,$oriaa,$mutaa); 
  
  $lcod{"TTT"}="F"; 
  $lcod{"TTC"}="F"; 
  $lcod{"TTA"}="L"; 
  $lcod{"TTG"}="L"; 
  $lcod{"TCT"}="S"; 
  $lcod{"TCC"}="S"; 
  $lcod{"TCA"}="S"; 
  $lcod{"TCG"}="S"; 
  $lcod{"TAT"}="T"; 
  $lcod{"TAC"}="T"; 
  $lcod{"TAA"}="STOP"; 
  $lcod{"TAG"}="STOP"; 
  $lcod{"TGA"}="STOP"; 
  $lcod{"TGT"}="C"; 
  $lcod{"TGC"}="C"; 
  $lcod{"TGG"}="W"; 
  $lcod{"CTT"}="L"; 
  $lcod{"CTC"}="L"; 
  $lcod{"CTA"}="L"; 
  $lcod{"CTG"}="L"; 
  $lcod{"CCT"}="P"; 
  $lcod{"CCC"}="P"; 
  $lcod{"CCA"}="P"; 
  $lcod{"CCG"}="P"; 
  $lcod{"CAT"}="H"; 
  $lcod{"CAC"}="H"; 
  $lcod{"CAA"}="Q"; 
  $lcod{"CAG"}="Q"; 
  $lcod{"CGT"}="R"; 
  $lcod{"CGC"}="R"; 
  $lcod{"CGA"}="R"; 
  $lcod{"CGG"}="R"; 
  $lcod{"ATT"}="I"; 
  $lcod{"ATC"}="I"; 
  $lcod{"ATA"}="I"; 
  $lcod{"ATG"}="M"; 
  $lcod{"ACT"}="T"; 
  $lcod{"ACC"}="T"; 
  $lcod{"ACA"}="T"; 
  $lcod{"ACG"}="T"; 
  $lcod{"AAT"}="N"; 
  $lcod{"AAC"}="N"; 
  $lcod{"AAA"}="K"; 
  $lcod{"AAG"}="K"; 
  $lcod{"AGT"}="S"; 
  $lcod{"AGC"}="S"; 
  $lcod{"AGA"}="R"; 
  $lcod{"AGG"}="R"; 
  $lcod{"GTT"}="V"; 
  $lcod{"GTC"}="V"; 
  $lcod{"GTA"}="V"; 
  $lcod{"GTG"}="V"; 
  $lcod{"GCT"}="A"; 
  $lcod{"GCC"}="A"; 
  $lcod{"GCA"}="A"; 
  $lcod{"GCG"}="A"; 
  $lcod{"GAT"}="D"; 
  $lcod{"GAC"}="D"; 
  $lcod{"GAA"}="E"; 
  $lcod{"GAG"}="E"; 
  $lcod{"GGT"}="G"; 
  $lcod{"GGC"}="G"; 
  $lcod{"GGA"}="G"; 
  $lcod{"GGG"}="G"; 

  if(exists $lcod{$lori}) { $oriaa=$lcod{$lori}; }
  if(exists $lcod{$lmut}) { $mutaa=$lcod{$lmut}; }

  if($oriaa eq $mutaa) 
  {
     $lch="Syn"; 
  }
  if($oriaa ne $mutaa) 
  {
     $lch="Nonsyn"; 
  }
  if($mutaa eq "STOP")
  {
     $lch="Nonsense"; 
  } 
  return $lch."*$mutaa*$oriaa"; 
}
