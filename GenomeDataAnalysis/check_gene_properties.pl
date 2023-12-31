#!/usr/bin/perl -w 
#
# Check details of genes which contain mutations only in the evolved lines
# Associate these genes with stress response, gene regulation 
# Functional association 
#

open(FP,"../datasets/GENE_Common_Names.txt") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   @arr=split(/\s+/,$fp); 
   $comgene{$arr[0]}=$arr[1]; 
   undef(@arr); 
}
close(FP); 

open(FP,"../datasets/all_regulators.tsv") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   if($fp=~/^Gene/) { next; }

   @arr=split(/\s+/,$fp); 
   $reglist{$arr[2]}=1; 
   undef(@arr); 
}
close(FP);

open(FP,"../datasets/paxdb_expr_check/paxdb_cerevisiae_abundance.txt") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   if($fp=~/string/) { next; } 
   @arr=split(/\s+/,$fp); 
   @qw=split(/\./,$arr[0]); 
   
   if(exists $comgene{$qw[1]}) { $qw[1]=$comgene{$qw[1]}; } 
   $paxlist{$qw[1]}=$arr[1]; 
   
   undef(@qw); 
   undef(@arr); 
}
close(FP); 


open(FP,"../datasets/SUMMARY_All_PROMOTERS.txt") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   if($fp=~/^Gene/) { next; }
   @arr=split(/\s+/,$fp); 
   $nslist{$arr[1]}[0]=$arr[5];
   $nslist{$arr[1]}[1]=$arr[6];
   undef(@arr); 
}
close(FP);


open(FP,"../datasets/genes_STRE.txt") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   @arr=split(/\s+/,$fp); 
   $strlist{$arr[0]}=1; 
   undef(@arr); 
}
close(FP); 

## Essential ORFs

open(FP,"../datasets/Yeast_essential_ORFs.txt") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   if($fp=~/^rec\_/ || $fp=~/\=+/) { next; } 
   @arr=split(/\s+/,$fp); 
   if(exists $comgene{$arr[1]}) { $arr[1]=$comgene{$arr[1]}; } 
   $essential{$arr[1]}=1; 
   undef(@arr); 
}
close(FP); 


## Our eLife mitochondria data 

open(FP,"../datasets/eLife2019_T4_vs_T1_data.txt") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   if($fp=~/^Gene/) { next; } 
   @arr=split(/\s+/,$fp); 
   $t4t1{$arr[0]}=$arr[2]; 
   undef(@arr); 
}
close(FP); 

open(FP,"../datasets/eLife2019_T2_vs_T1_data.txt") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   if($fp=~/^Gene/) { next; } 
   @arr=split(/\s+/,$fp); 
   $t2t1{$arr[0]}=$arr[2]; 
   undef(@arr); 
}
close(FP); 


## Gasch data 

open(FP,"../datasets/Gasch_data/Gasch_complete_dataset.txt") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   if($fp=~/^UID/) 
   { 
      @arr=split(/\t/,$fp);  
      $no=@arr; 
      for($i=3;$i<$no;$i++)
      {
	while($arr[$i]=~/\s+/) 
	{
           $arr[$i]=~s/\s+/\_/; 
	}
	$head_gs[1][$i-3]="GS_".$arr[$i]; 
	#print "$i $arr[$i]\n";
      }
      $head_gs[0]=$no-3; 
      undef $no;  
      undef(@arr);       
      next; 
   }

   @arr=split(/\t/,$fp); 
   $no=@arr; 

    
   @qw=split(/\s/,$arr[1]); 
   if($qw[1] eq "") 
   {
      if(exists $comgene{$qw[0]}) { $qw[1]=$comgene{$qw[0]}; }
      else { $qw[1]=$qw[0]; } 
   }
   for($i=3;$i<$no;$i++)
   {
      if($arr[$i] eq "") { $arr[$i]="NA"; } 
      $gasch{$qw[1]}[$i-3]=$arr[$i]; 
      #print "$qw[1] $arr[$i]\n"; 
   }
   undef(@qw); 

   #print "$no $head[0] $gasch{$arr[2]}[$head[0]-1]\n"; 
   undef $no; 
   undef(@arr); 
}
close(FP); 

open(FP,"../datasets/Causton_et_al_paper/Causton_processed_data.txt") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   chop($fp); 
   if($fp=~/^ORF/) 
   { 
      @arr=split(/\t/,$fp);  
      $no=@arr; 
      for($i=2;$i<$no;$i++)
      {
	while($arr[$i]=~/\s+/) 
	{
           $arr[$i]=~s/\s+/\_/; 
	}
	$head_caust[1][$i-2]="CT_".$arr[$i]; 
	#print "$i $arr[$i]\n";
      }
      $head_caust[0]=$no-2; 
      undef $no;  
      undef(@arr);       
      next; 
   }

   @arr=split(/\t/,$fp); 
   $no=@arr; 
    
   if(exists $comgene{$arr[1]}) { $arr[1]=$comgene{$arr[1]}; }

   for($i=2;$i<$no;$i++)
   {
      if($arr[$i] eq "") { $arr[$i]="NA"; } 
      if($arr[$i]=~/\#DIV/) { $arr[$i]="NA"; } 
      $causton{$arr[1]}[$i-2]=$arr[$i]; 
      #print "$arr[1] $arr[$i]\n"; 
   }

   ##print "$no $head_caust[0] $causton{$arr[2]}[$head_caust[0]-1]\n"; 
   undef $no; 
   undef(@arr); 
}
close(FP); 


# Global_GO_mapping.txt

open(FP,"../Mutation_annotations/Mutation_annotation_EvoOnly.txt") or die; 

open(WR,">../results/Analysis_EvoOnly_MutationGenes_properties.txt") or die; 

print WR "Chr\tPos\tMut\tAnc\tN1_Mut N1_nCov N1_Freq\tN2_Mut N2_nCov N2_Freq\tN3_Mut N3_nCov N3_Freq\tMut_Annot\tIsRegulator\tpaxdb_abund\tDM_YPD\tTATAbox\tSTRE_elem\tIsEssential\tT4vsT1_LFC\tT2vsT1_LFC\t"; 

for($i=0;$i<$head_gs[0];$i++)
{
   print WR "$head_gs[1][$i]\t"; 
}

for($i=0;$i<$head_caust[0];$i++)
{
   print WR "$head_caust[1][$i]\t"; 
}

print WR "\n"; 

while($fp=<FP>)
{
   chomp($fp); 
   @arr=split(/\t+/,$fp); 
   if($arr[10] ne " X") 
   {
      @q1=split(/\s+/,$arr[7]);
      @q2=split(/\s+/,$arr[8]);
      @q3=split(/\s+/,$arr[9]);
      print WR "$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\t$q1[2]\t$q2[2]\t$q3[2]\t$arr[10]\t"; 
      undef(@q1);
      undef(@q2);
      undef(@q3);

      $gnname="X"; 

      if($arr[10]=~/Coding\;/) 
      {
	@qw=split(/\;/,$arr[10]);
	$gnname=$qw[1]; 
	undef(@qw); 
      }
      if($arr[10]=~/Intron/ || $arr[10]=~/Coding_Indel/) 
      {
	@qw=split(/\./,$arr[10]);
	$gnname=$qw[1]; 
	undef(@qw); 
      }
      if($arr[10]=~/Prom\_/) 
      {
	chop($arr[10]); 
	@qw=split(/\_/,$arr[10]);
	$gnname=$qw[1]; 
	undef(@qw); 
      }

      if(exists $reglist{$gnname})
      {
	print WR "YES\t"; 
      }
      else
      {
	print WR "No\t"; 
      }

      if(exists $paxlist{$gnname}) 
      {
	print WR "$paxlist{$gnname}\t"; 
      }
      else
      {
	print WR "NA\t"; 
      }

      if(exists $nslist{$gnname})
      {
	print WR "$nslist{$gnname}[0]\t$nslist{$gnname}[1]\t"; 
      }
      else
      {
	print WR "NA\tNA\t"; 
      }

      if(exists $strlist{$gnname})
      {
	print WR "YES\t"; 
      }
      else
      {
	print WR "No\t"; 
      }

      if(exists $essential{$gnname}) 
      {
	 print WR "$essential{$gnname}\t"; 
      }
      else
      {
	 print WR "No\t"; 
      }

      if(exists $t4t1{$gnname}) 
      {
	 print WR "$t4t1{$gnname}\t"; 
      }
      else
      {
	 print WR "NA\t"; 
      }

      if(exists $t2t1{$gnname}) 
      {
	 print WR "$t2t1{$gnname}\t"; 
      }
      else
      {
	 print WR "NA\t"; 
      }

      if(exists $gasch{$gnname}) 
      {
	for($j=0;$j<$head_gs[0];$j++)
	{
	   print WR "$gasch{$gnname}[$j]\t"; 
	}
      }
      else
      {
	for($j=0;$j<$head_gs[0];$j++)
	{
           print WR "NA\t"; 
	}
      }

      if(exists $causton{$gnname}) 
      {
	for($j=0;$j<$head_caust[0];$j++)
	{
	   print WR "$causton{$gnname}[$j]\t"; 
	}
      }
      else
      {
	for($j=0;$j<$head_caust[0];$j++)
	{
           print WR "NA\t"; 
	}
      }

      print WR "\n"; 
      undef $gnname; 
   }

   undef(@arr); 
}
close(WR); 
close(FP); 

undef %essential; 
undef %t2t1; 
undef %t4t1; 
undef(@head_gs);
undef(@head_caust); 
undef %causton; 
undef %gasch; 
undef %reglist; 
undef %paxlist; 
undef %nslist; 
undef %comgene; 
