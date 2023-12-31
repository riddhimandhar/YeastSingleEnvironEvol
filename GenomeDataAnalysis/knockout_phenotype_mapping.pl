#!/usr/bin/perl -w 
#
# Check mutated gene knockout phenotypes from published datasets 

open(FP,"../datasets/GENE_Common_Names.txt") or die;
while($fp=<FP>)
{
   chomp($fp);
   @arr=split(/\s+/,$fp); 
   $comlist{$arr[1]}=$arr[0]; 
   undef(@arr); 
}
close(FP); 

## Datasets to be used 

open(FP,"../datasets/yeast_deletion_fitness_data/fitness_data/Hillenmeyer_hom.pub") or die;
while($fp=<FP>)
{
   chomp($fp); 
   if($fp=~/^Orf/) 
   {
      @arr=split(/\t/,$fp); 
      $no=@arr;
      for($i=1;$i<$no;$i++)
      {
         @qw=split(/\:/,$arr[$i]);
	 
	 while($qw[1]=~/\s+/) { $qw[1]=~s/\s+//; } 

	 if($qw[2] ne "" && $qw[7] ne "") 
	 {
	    $id=$qw[1].$qw[2]."_$qw[7]";
	 }
	 else
	 {
            $id=$qw[1]; 
	 }
	 $head[$i-1]="Hi_".$id; 
	 #print "$id\n"; 
	 undef $id; 
	 undef(@qw);
      }
      $chead=$no-1; 
      undef $no; 
      undef(@arr); 
      next; 
   }

   @arr=split(/\t/,$fp); 

   @qw=split(/\:/,$arr[0]); 
   $no=@arr; 
   for($i=1;$i<$no;$i++)
   {
      if($arr[$i] ne "NA") 
      {
        $hilist{$qw[0]}[$i-1]=(-$arr[$i]); 
      }
      else
      {
        $hilist{$qw[0]}[$i-1]=$arr[$i]; 
      }
   }
   $hic=$no-1; 
   undef $no; 
   undef(@qw); 

   undef(@arr); 
}
close(FP); 

open(FP,"../datasets/yeast_deletion_fitness_data/fitness_data/Stonge_data.txt") or die;
while($fp=<FP>)
{
   chomp($fp); 
   if($fp=~/^orf/) 
   {
      @arr=split(/\s+/,$fp); 
      $head[$chead]="Stonge_MMS1"; 
      $head[$chead+1]="Stonge_MMS2"; 
      $chead+=2; 
      undef(@arr); 
      next; 
   }
   @arr=split(/\s+/,$fp); 

   if($arr[2] eq "") { $arr[2]="NA"; }  
   if($arr[3] eq "") { $arr[3]="NA"; }  
   if($arr[2] ne "NA") { $arr[2]=(-$arr[2]); } 
   if($arr[3] ne "NA") { $arr[3]=(-$arr[3]); } 
   $stongelist{$arr[0]}[0]=$arr[2]; 
   $stongelist{$arr[0]}[1]=$arr[3]; 
   #print "$arr[0] $arr[2] $arr[3]\n"; 
   undef(@arr);    
}
close(FP); 


open(FP,"../datasets/yeast_deletion_fitness_data/fitness_data/Jonas_processed_data.txt") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   if($fp=~/^Gene/) 
   { 
      $head[$chead]="Jonas_Salt"; 
      $chead+=1; 
      next; 
   }
   @arr=split(/\t/,$fp); 
   $arr[0]=uc($arr[0]); 
   #print "$arr[0] $arr[1]\n"; 
   if($arr[1] eq "") { $arr[1]="NA"; } 
   $jonlist{$arr[0]}=$arr[1]; 
   undef(@arr); 
}
close(FP); 


open(FP,"../datasets/yeast_deletion_fitness_data/fitness_data/Tucker_data.txt") or die;
while($fp=<FP>)
{
   chomp($fp); 
   if($fp=~/^ORF/) 
   {
      @arr=split(/\t/,$fp); 
      $head[$chead]="Tucker_".$arr[3]; 
      $head[$chead+1]="Tucker_".$arr[4]; 
      $head[$chead+2]="Tucker_".$arr[5]; 
      $head[$chead+3]="Tucker_".$arr[6]; 
      $chead+=4; 
      undef(@arr);
      next; 
   }
   @arr=split(/\t/,$fp); 
   $arr[0]=uc($arr[0]); 
   if($arr[3] ne "") 
   {
     $tucklist{$arr[0]}[0]=(-$arr[3]); 
   }
   else
   {
     $tucklist{$arr[0]}[0]="NA"; 
   }

   if($arr[4] ne "") 
   {
     $tucklist{$arr[0]}[1]=(-$arr[4]); 
   }
   else
   {
     $tucklist{$arr[0]}[1]="NA"; 
   }

   if($arr[5] ne "") 
   {
     $tucklist{$arr[0]}[2]=(-$arr[5]); 
   }
   else
   {
     $tucklist{$arr[0]}[2]="NA"; 
   }

   if($arr[6] ne "") 
   {
     $tucklist{$arr[0]}[3]=(-$arr[6]); 
   }
   else
   {
     $tucklist{$arr[0]}[3]="NA"; 
   }
   undef(@arr); 
}
close(FP); 

$brc=0; 
open(FP,"../datasets/yeast_deletion_fitness_data/fitness_data/Brown_data.txt") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   chop($fp); 
   if($fp=~/^2004/) 
   {
      @arr=split(/\t/,$fp); 
      $no=@arr; 
      for($i=4;$i<$no;$i++)
      {
        $head[$chead]="Brown_".$arr[$i]; 
	$chead++;
      }
      undef $no; 
      undef(@arr); 
      next; 
   }

   @arr=split(/\t/,$fp); 
   $no=@arr; 

   for($i=4;$i<$no;$i++)
   {
      if($arr[$i] eq "") { $arr[$i]="NA"; } 
      $brownlist{$arr[0]}[$i-4]=$arr[$i]; 
   }
   if($no-4>$brc) { $brc=$no-4; } 
   for($j=$i;$j<$brc;$j++)
   {
      $brownlist{$arr[0]}[$j-4]="NA"; 
   }
   #print "$brc $no $arr[3]\n"; 
   undef $no; 
   undef(@arr); 
}
close(FP); 

open(FP,"../datasets/yeast_deletion_fitness_data/fitness_data/Baryshnikova_data.txt") or die;
while($fp=<FP>)
{
   chomp($fp); 
   if($fp=~/^Gene/)
   { 
     next; 
   } 
   @arr=split(/\s+/,$fp); 
   $barylist{$arr[0]}=$arr[1]; 
   undef(@arr); 
}
close(FP); 


open(FP,"../datasets/Need_based_upregulation_duplicate_genes/DeLuna_need_based_duplicates_YPD.txt") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   if($fp=~/^\[/ || $fp=~/^ID/) { next; } 
   @arr=split(/\s+/,$fp); 
   if(exists $comlist{$arr[1]}) 
   {
      $arr[1]=$comlist{$arr[1]}; 
   }
   if(exists $comlist{$arr[2]}) 
   {
      $arr[2]=$comlist{$arr[2]}; 
   }
   $duplist{$arr[1]}=1;
   $duplist{$arr[2]}=1;
   undef(@arr);
}
close(FP); 


open(FP,"../datasets/YGOB_WolfeLab.txt") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   if($fp=~/^Number/) { next; } 
   @arr=split(/\,/,$fp); 
   $ygoblist{$arr[2]}=1; 
   $ygoblist{$arr[4]}=1; 
   undef(@arr); 
}
close(FP); 


open(WR,">../results/Mutant_gene_knockout_phenotype.txt") or die; 

print WR "Gene\tComName\tN1_freq N2_freq N3_freq Mut_annot\tBaryshnikova_data ";

for($i=0;$i<$chead;$i++)
{
   print WR "$head[$i] "; 
}
print WR "NeedBasedDup YGOBdup\n"; 

open(FP,"../Mutation_annotations/Mutation_annotation_EvoOnly.txt") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   @arr=split(/\t+/,$fp); 
   
   $gnname="X"; 

   if($arr[10] ne " X") 
   {
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
   }

   if($gnname eq "X") { next; }

   $gnname2=$gnname; 
   if(exists $comlist{$gnname}) 
   {
      $gnname=$comlist{$gnname}; 
      print "$gnname $gnname2\n"; 
   }

   @q1=split(/\s+/,$arr[7]); 
   @q2=split(/\s+/,$arr[8]); 
   @q3=split(/\s+/,$arr[9]); 
   
   print WR "$gnname2\t$gnname\t$q1[2] $q2[2] $q3[2] $arr[10]\t";   

   undef(@q1);
   undef(@q2);
   undef(@q3);

   if(exists $barylist{$gnname})
   {
      print WR "$barylist{$gnname} "; 
   }
   else { print WR "NA "; } 

   if(exists $hilist{$gnname}) 
   {
     for($i=0;$i<$hic;$i++)
     { 
       print WR "$hilist{$gnname}[$i] ";
     }
   }
   else
   {
     for($i=0;$i<$hic;$i++)
     { 
       print WR "NA ";
     }
   }

   if(exists $stongelist{$gnname})
   {
      print WR "$stongelist{$gnname}[0] $stongelist{$gnname}[1] ";
   }
   else
   {
      print WR "NA NA "; 
   }

   if(exists $jonlist{$gnname})
   {
      print WR "$jonlist{$gnname} "; 
   }
   else
   {
      print WR "NA "; 
   }

   if(exists $tucklist{$gnname})
   {
      print WR "$tucklist{$gnname}[0] $tucklist{$gnname}[1] $tucklist{$gnname}[2] $tucklist{$gnname}[3] "; 
   }
   else
   {
      print WR "NA NA NA NA "; 
   }

   if(exists $brownlist{$gnname})
   {
     for($i=0;$i<$brc;$i++)
     { 
       print WR "$brownlist{$gnname}[$i] ";
     }
   }
   else
   {
     for($i=0;$i<$brc;$i++)
     { 
       print WR "NA ";
     }
   }

   if(exists $duplist{$gnname}) 
   {
      print WR "Yes "; 
   }
   else { print WR "No "; }

   if(exists $ygoblist{$gnname}) 
   {
      print WR "Yes "; 
   }
   else { print WR "No "; }

   print WR "\n"; 

   undef $gnname; 
   undef(@arr); 
}
close(WR); 
close(FP); 

undef %ygoblist; 
undef %duplist; 
undef %hilist; 
undef %stongelist; 
undef %jonlist;
undef %tucklist; 
undef %brownlist; 
undef %comlist;
