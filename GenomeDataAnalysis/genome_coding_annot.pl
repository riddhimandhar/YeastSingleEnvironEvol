#!/usr/bin/perl -w 
# Find yeast gene coding sequences, locations, functions 
#

open(FP,"yeast_ref_genome/S288C_reference_genome_R64-2-1_20150113/orf_coding_all_R64-2-1_20150113.fasta") or die; 
open(WR,">../datasets/Yeast_genome_annot_func.txt") or die; 

$flag=0; 
$seq=""; 
while($fp=<FP>)
{
   chomp($fp); 
   chop($fp); 
   if($fp=~/^\>/) 
   {
      if($fp=~/\sDubious\s/) { next; } 
      if($flag==1) 
      {
	 if($intron eq "") {  $intron="NA"; } 
         print WR "$gene\t$comgene\t$chr\t$st\t$en\t$sign\t$intron\t$seq\t$func\n"; 
	 undef $gene;
	 undef $comgene; 
	 undef $chr; 
	 undef $st;
	 undef $en; 
	 undef $sign; 
	 undef $intron;
	 undef $seq; 
	 undef $func; 
      }
      
      $flag=0;
      $seq="";   

      @fsp=split(/\"/,$fp); 
      @head=split(/\, /,$fsp[0]);

      @qw=split(/\s+/,$head[0]);
      @er=split(/\>/,$qw[0]); 
      $gene=$er[1]; 
      $comgene=$qw[1]; 
      undef(@er); 
      undef(@qw); 

      @qw=split(/\s+/,$head[1]); 
      if($qw[1] eq "Mito") { $qw[1]="M"; } 
      $chr=$qw[0]."$qw[1]"; 
      @er=split(/\,/,$qw[3]); 
      $nn=@er; 
      @we=split(/\-/,$er[0]); 
      $st=$we[0]; 
      $en=$we[1]; 
      $sign="+"; 
      if($st>$en) { $sign="-"; } 
      undef(@we); 
      $intron=""; 

      for($k=1;$k<$nn;$k++) 
      {
        @we=split(/\-/,$er[$k]); 
	if($sign eq "+") 
	{
          $p1=$en+1;
	  $p2=$we[0]-1;
	  $intron.="$p1-$p2,";
	  undef $p1; 
	  undef $p2; 
	  if($we[0]<$st) { $st=$we[0]; }
	  if($we[1]>$en) { $en=$we[1]; } 
	}
	if($sign eq "-") 
	{
	  $p1=$st+1; 
	  $p2=$we[1]-1; 
	  $intron.="$p1-$p2,";
	  undef $p1; 
	  undef $p2; 
	  if($we[0]>$st) { $st=$we[0]; }
	  if($we[1]<$en) { $en=$we[1]; } 
	}
        undef(@we); 
      }
      if($sign eq "-") 
      {
	$tmp=$st; 
	$st=$en; 
	$en=$tmp; 
	undef $tmp; 
      }
      undef $nn;
      undef(@er); 
      undef(@qw); 

      $func=$fsp[1];
      undef(@head);  
      undef(@fsp);
      $flag=1; 
      next;
   } 
   if($flag==1) 
   {
      $seq.=$fp; 
      #print "$flag $fp $seq\n"; 
   }
}
close(FP); 

if($flag==1) 
{
    if($intron eq "") {  $intron="NA"; } 

    print WR "$gene\t$comgene\t$chr\t$st\t$en\t$sign\t$intron\t$seq\t$func\n"; 
    undef $gene;
    undef $comgene; 
    undef $chr; 
    undef $st;
    undef $en; 
    undef $sign; 
    undef $intron; 
    undef $seq; 
    undef $func; 
}
      
$flag=0;

close(WR); 
