#!/usr/bin/perl -w 
#
# Process sam files generated after alignment 
#

use LWP::Simple;
use Parallel::ForkManager;

open(FP,"../results/S288C_chromosome_id_mapping.txt") or die; 
while($fp=<FP>) 
{
   chomp($fp); 
   @arr=split(/\s+/,$fp); 
   $chrlist{$arr[0]}[0]=$arr[1];  
   $chrlist{$arr[0]}[1]=$arr[3];  
   undef(@arr); 
}
close(FP); 


$file[0]="samOut_R0-1.sam";
$file[1]="samOut_R0-2.sam";
$file[2]="samOut_R0-3.sam";
$file[3]="samOut_R100-1.sam";
$file[4]="samOut_R100-2.sam";
$file[5]="samOut_R100-3.sam";

$pm = Parallel::ForkManager->new(20);

LINKS:

for($gi=0;$gi<6;$gi++)
{

$pm->start and next LINKS; # do the fork

@eq=split(/\_/,$file[$gi]); 
@re=split(/\./,$eq[1]); 

$wrf1="../results/MutationsFromAlignment_$re[0].txt"; 
$wrf2="../results/SamAlignment_$re[0].txt";

undef(@re); 
undef(@eq); 

print "$file[$gi]\n"; 

open(FP,"/home/shreya/SO_10693_rawdata/SO_10693_rawdata/newSam_output_files/$file[$gi]") or die; 
system "sudo sed -n \'\$=\' /home/shreya/SO_10693_rawdata/SO_10693_rawdata/newSam_output_files/$file[$gi] > TMP_$gi";

open(TM,"TMP_$gi") or die; 
while($tm=<TM>) 
{
   chomp($tm); 
   $totln=$tm; 
}

system "rm TMP_$gi"; 

$ln=0; 
open(WR1,">$wrf1") or die; 
open(WR2,">$wrf2") or die; 

while($fp=<FP>) 
{
   chomp($fp); 
   if($fp=~/^\@/) { next; } 
   $ln++; 
   @arr=split(/\s+/,$fp); 
   $chid=$arr[2]; 
   $chrseq=""; 
   #print "$arr[0]\n";

   if(exists $chrlist{$chid}) 
   {
      $chrseq=$chrlist{$chid}[1]; 
      $chid=$chrlist{$chid}[0]; 
   }

   $pos=$arr[3];
   $mapq=$arr[4]; 
   $rdseq=$arr[9]; 
   $rdqual=$arr[10]; 
   $orcig=$arr[5]; 

   $no=@arr; 

   $ndf=0; $dfs=""; $als=0; $alts=0; $mode="";  
   for($i=11;$i<$no;$i++) 
   {
      if($arr[$i]=~/NM/) 
      {
	 @qw=split(/\:/,$arr[$i]); 
         $ndf=$qw[2]; 
	 undef(@qw); 
      }
      if($arr[$i]=~/MD/) 
      {
	 @qw=split(/\:/,$arr[$i]); 
         $dfs=$qw[2];
	 undef(@qw);  
      }
      if($arr[$i]=~/AS/) 
      {
	 @qw=split(/\:/,$arr[$i]); 
         $als=$qw[2];
	 undef(@qw);  
      }
      if($arr[$i]=~/XS/) 
      {
	 @qw=split(/\:/,$arr[$i]); 
         $alts=$qw[2];
	 undef(@qw);  
      }
      if($arr[$i]=~/YT/) 
      {
	 @qw=split(/\:/,$arr[$i]); 
         $mode=$qw[2];
	 undef(@qw);  
      }
   }

   $cig=$orcig; 

   #print "$arr[0]\t$chid\t$pos\t$len\t$mapq\t$als\t$alts\t$mode\t$ndf\t$dfs\t$als\t$alts\n"; 

   if(($mode eq "CP")) ## && ($als ne $alts)) ## && ($ndf!=0))
   {
      if($cig=~/S/) 
      {
	 $nrdseq=$rdseq; 
         if($cig=~/S$/) 
	 {
            $cnt=0; 
	    for($i=length($cig)-2;$i>=0;$i--) 
            {
		$ft=substr($cig,$i,1);
		#print "$ft\n"; 
                if($ft!~/[0-9]/) { last; } 
		$cnt++; 
	    }
	    $cc=""; 
	    $cc=substr($cig,$i+1,$cnt); 
	    $cig=substr($cig,0,length($cig)-$cnt-1); 
	    $nrdseq=substr($rdseq,0,length($rdseq)-$cc); ## 
	    #print "HERE $cig $cnt $cc S\n"; 
            undef $cc; 
	    undef $cnt;
	 }
	 
	 @qw=split(/S/,$cig); 
	 $nn=@qw; 
	 if($nn==2)
	 {
           $nrdseq=substr($nrdseq,$qw[0],length($rdseq)-$qw[0]); 
	 }
         undef(@qw); 

	 #print "$cig\t$rdseq\n$nrdseq\n"; 
         $rdseq=$nrdseq; 
         undef $nrdseq; 	 
      }

      $len=length($rdseq);
      $substr=substr($chrseq,$pos-1,$len+5); 
      if($ndf==0)
      {
	 #print WR2 "$substr\n"; 
	 #print WR2 "$rdseq\n"; 
         $end=$pos+$len; 
         print WR1 "$arr[0]\t$chid\t$ndf\tAS:$als\tXS:$alts\t$pos\t$end\n"; 
      }
      elsif($ndf!=0) 
      {
         print WR1 "$arr[0]\t$chid\t$ndf\tAS:$als\tXS:$alts\t"; 
         print WR2 ">$arr[0]\t$chid\t$ndf\n"; 
         align($rdseq,$substr,$pos,$rdqual); 
      }
      undef $end; 
      undef $len; 
      undef $substr; 
   }

   undef $cig; 
   undef $orcig; 
   undef $chid; 
   undef $pos; 
   undef $mode; 
   undef $mapq; 
   undef $ndf; 
   undef $dfs; 
   undef $als; 
   undef $alts; 
   undef $rdseq; 
   undef $rdqual; 
   undef $chrseq; 

   undef(@arr); 

   $rat=$ln/$totln*100; 
   $rem=$rat/2; 
   $fr=sprintf("%0.0f",$rem); 
    
   if($rem-$fr==0) 
   {
      print " $file[$gi]\t$rat% complete\n"; 
   }
   undef $fr; 
   undef $rem; 
   undef $rat; 
}
close(FP); 
close(WR2); 
close(WR1); 
undef $wrf2; 
undef $wrf1; 

$pm->finish; # do the exit in the child process

}

$pm->wait_all_children;


undef %chrlist; 

sub align
{
   local($que,$ref,$lpos,$lql)=@_;
   local ($matsc,$mismatch,$gap,$cnmism,@lar1,@lar2,$li,$lj,$matsc,$mismatch,$gap,$sc1,$sc2,$sc3,$misfl,$cnmism,@las,@locsc,$lst1,$lst2,$lid,%track,$id1,$id2,$ll,$ls1,$ls2,$mp,$lq,$lg,$enpos,$tmp,$prev1,$prev2);

   $matsc=5;
   $mismatch=-2;
   $gapop=-20;

   #print WR2 "$substr\n"; 
   #print WR2 "$rdseq\n"; 

   for($li=1;$li<=length($que);$li++)
   {
      $lar1[$li]=substr($que,$li-1,1);
   }
   for($lj=1;$lj<=length($ref);$lj++)
   {
      $lar2[$lj]=substr($ref,$lj-1,1);
   }
   for($li=0;$li<=length($que);$li++)
   {
      for($lj=0;$lj<=length($ref);$lj++)
      {
         if($li==0) { $locsc[$li][$lj]=$gapop*$lj; next; }
         if($lj==0) { $locsc[$li][$lj]=$gapop*$li; next; }
         $misfl=0;
         if($lar1[$li] eq $lar2[$lj]) { $sc1=$locsc[$li-1][$lj-1]+$matsc; }
         else { $sc1=$locsc[$li-1][$lj-1]+$mismatch; $misfl=1; }
        
	 $sc2=$locsc[$li-1][$lj]+$gapop;
         $sc3=$locsc[$li][$lj-1]+$gapop;
      																		             ##print "$li $lj $lar1[$li] $lar2[$lj] $locsc[$li-1][$lj-1] $locsc[$li-1][$lj] $locsc[$li][$lj-1] $sc1 $sc2 $sc3\n";
        if($sc1>=$sc2 && $sc1>=$sc3)
        {
          $locsc[$li][$lj]=$sc1;
          $id1=$li-1;
          $id2=$lj-1;
          $track{$li."*$lj"}=$id1."*$id2*$misfl";
          undef $id1;
          undef $id2;
        }
        if($sc2>=$sc1 && $sc2>=$sc3)
        {
          $locsc[$li][$lj]=$sc2;
          $id1=$li-1;
          $id2=$lj;
	  $misfl=2;
          $track{$li."*$lj"}=$id1."*$id2*$misfl";
          undef $id1;
          undef $id2;
        }
        if($sc3>=$sc1 && $sc3>=$sc2)
        {
          $locsc[$li][$lj]=$sc3;
          $id1=$li;
          $id2=$lj-1;
	  $misfl=2;
          $track{$li."*$lj"}=$id1."*$id2*$misfl";
          undef $id1;
          undef $id2;
        }
      }
   }


   $lst1=length($que);
   $lst2=length($ref);

   #print "$lst1 $lst2 | \t"; 

   $max=0;  
   $xr=$lst1; $yr=$lst2; 

   for($li=$lst1;$li>=10;$li--)
   {
      for($lj=$lst2;$lj>=10;$lj--)
      {
        if($locsc[$li][$lj]>$max)
        {
           $max=$locsc[$li][$lj];
           $xr=$li; $yr=$lj;
        }
      }
  }
  $lst1=$xr; $lst2=$yr;

   ##print "MAX: $lst1 $lst2 $max\n"; 
  
   $qstr=""; $rstr="";

#   if($lst1<length($que))
#   {
#      for($li=length($que);$li>=$lst1+1;$li--) 
#      {
#	 $rstr.=substr($que,$li-1,1); 
#	 $qstr.="-";
#      }
#   }
#   if($lst2<length($ref))
#   {
#      for($li=length($ref);$li>=$lst2+1;$li--) 
#      {
#	 $rstr.="-";
#	 $qstr.=substr($ref,$li-1,1); 
#      }
#   }

   $cnmism=0;
   $prev1=-1; $prev2=-1; 
   while($lst1!=0 && $lst2!=0)
   {
     $lid=$lst1."*$lst2";
     if(exists $track{$lid})
     {
       if($lst1==$prev1) 
       {
          $qstr.="-";  
          $rstr.=$lar2[$lst2];  
       }
       elsif($lst2==$prev2) 
       {
          $rstr.="-";  
          $qstr.=$lar1[$lst1];  
       }
       elsif($lst1!=$prev1 && $lst2!=$prev2)
       {
         $qstr.=$lar1[$lst1]; 
         $rstr.=$lar2[$lst2]; 
       }
       $prev1=$lst1; 
       $prev2=$lst2; 

       @las=split(/\*/,$track{$lid});
       $lst1=$las[0];
       $lst2=$las[1];
       $cnmism+=$las[2];
       undef(@las);

     }
     undef $lid;
   }

   $tmp=$rstr; 
   while($tmp=~/\-/) 
   {
      $tmp=~s/\-//; 
   }

   $lg=length($tmp); 
   $enpos=$lpos+$lg; 
   print WR1 "$lpos\t$enpos\t"; 
   undef $tmp; 
   undef $lg;
   undef $enpos; 


   undef $que; 
   undef $ref; 
   undef %track;
   undef $lst1;
   undef $lst2;

   undef(@locsc);
   undef(@lar1);
   undef(@lar2);

   $qstr=reverse($qstr);
   $rstr=reverse($rstr); 
   print WR2 "$qstr\n$rstr\n"; 

   $ll=length($rstr); 
   $lc=0; $lr=0;  
   for($li=0;$li<$ll;$li++)
   {
      $ls1=substr($qstr,$li,1); 
      $ls2=substr($rstr,$li,1); 
      if($ls1 ne "-") {  $lc++; } 
      if($ls2 ne "-") {  $lr++; } 

      if($ls1 ne $ls2)
      {
	 $mp=$lpos+$lr-1; 
	 if($ls1 ne "-") 
	 {
	   $lq=checkqual($lql,$lc-1); 
	   ##if($lc>=length($lql)) { print "HERE $lc $ll\n"; } 
	 }
	 else { $lq="NA"; } 
	 print WR1 "$mp $ls1 $ls2 $lq,"; 
      }
   }
   print WR1 "\n"; 
   undef $ll; 
   undef $ls1; 
   undef $ls2; 

   undef $qstr; 
   undef $rstr; 
   return $cnmism;
}

sub checkqual
{
   local($quals,$llp)=@_;
   local($lval,$nlval);  

   $lval=substr($quals,$llp,1); 

   #print "$lval\t"; 
      
   if($lval eq "\!") { $nlval=0; } 
   if($lval eq "\"") { $nlval=1; } 
   if($lval eq "\#") { $nlval=2; } 
   if($lval eq "\$") { $nlval=3; } 
   if($lval eq "\%") { $nlval=4; } 
   if($lval eq "\&") { $nlval=5; } 
   if($lval eq "\'") { $nlval=6; } 
   if($lval eq "\(") { $nlval=7; } 
   if($lval eq "\)") { $nlval=8; } 
   if($lval eq "\*") { $nlval=9; } 
   if($lval eq "\+") { $nlval=10; } 
   if($lval eq "\,") { $nlval=11; } 
   if($lval eq "\-") { $nlval=12; } 
   if($lval eq "\.") { $nlval=13; } 
   if($lval eq "\/") { $nlval=14; } 
   if($lval eq "0") { $nlval=15; } 
   if($lval eq "1") { $nlval=16; } 
   if($lval eq "2") { $nlval=17; } 
   if($lval eq "3") { $nlval=18; } 
   if($lval eq "4") { $nlval=19; } 
   if($lval eq "5") { $nlval=20; } 
   if($lval eq "6") { $nlval=21; } 
   if($lval eq "7") { $nlval=22; } 
   if($lval eq "8") { $nlval=23; } 
   if($lval eq "9") { $nlval=24; } 
   if($lval eq "\:") { $nlval=25; } 
   if($lval eq "\;") { $nlval=26; } 
   if($lval eq "\<") { $nlval=27; } 
   if($lval eq "\=") { $nlval=28; } 
   if($lval eq "\>") { $nlval=29; } 
   if($lval eq "\?") { $nlval=30; } 
   if($lval eq "\@") { $nlval=31; } 
   if($lval eq "A") { $nlval=32; } 
   if($lval eq "B") { $nlval=33; } 
   if($lval eq "C") { $nlval=34; } 
   if($lval eq "D") { $nlval=35; } 
   if($lval eq "E") { $nlval=36; } 
   if($lval eq "F") { $nlval=37; } 
   if($lval eq "G") { $nlval=38; } 
   if($lval eq "H") { $nlval=39; } 
   if($lval eq "I") { $nlval=40; } 
   if($lval eq "J") { $nlval=41; } 
   if($lval eq "K") { $nlval=42; } 
   if($lval eq "L") { $nlval=43; } 
   if($lval eq "M") { $nlval=44; } 
   if($lval eq "N") { $nlval=45; } 
   if($lval eq "O") { $nlval=46; } 
   if($lval eq "P") { $nlval=47; } 
   if($lval eq "Q") { $nlval=48; } 
   if($lval eq "R") { $nlval=49; } 
   if($lval eq "S") { $nlval=50; } 
   if($lval eq "T") { $nlval=51; } 
   if($lval eq "U") { $nlval=52; } 
   if($lval eq "V") { $nlval=53; } 
   if($lval eq "W") { $nlval=54; } 
   if($lval eq "X") { $nlval=55; } 
   if($lval eq "Y") { $nlval=56; } 
   if($lval eq "Z") { $nlval=57; } 
   if($lval eq "\[") { $nlval=58; } 
   if($lval eq "\\") { $nlval=59; } 
   if($lval eq "\]") { $nlval=60; } 
   if($lval eq "\^") { $nlval=61; } 
   if($lval eq "\_") { $nlval=62; } 
   if($lval eq "\`") { $nlval=63; } 
   if($lval eq "a") { $nlval=64; } 
   if($lval eq "b") { $nlval=65; } 
   if($lval eq "c") { $nlval=66; } 
   if($lval eq "d") { $nlval=67; } 
   if($lval eq "e") { $nlval=68; } 
   if($lval eq "f") { $nlval=69; } 
   if($lval eq "g") { $nlval=70; } 
   if($lval eq "h") { $nlval=71; } 
   if($lval eq "i") { $nlval=72; } 
   if($lval eq "j") { $nlval=73; } 
   if($lval eq "k") { $nlval=74; } 
   if($lval eq "l") { $nlval=75; } 
   if($lval eq "m") { $nlval=76; } 
   if($lval eq "n") { $nlval=77; } 
   if($lval eq "o") { $nlval=78; } 
   if($lval eq "p") { $nlval=79; } 
   if($lval eq "q") { $nlval=80; } 
   if($lval eq "r") { $nlval=81; } 
   if($lval eq "s") { $nlval=82; } 
   if($lval eq "t") { $nlval=83; } 
   if($lval eq "u") { $nlval=84; } 
   if($lval eq "v") { $nlval=85; } 
   if($lval eq "w") { $nlval=86; } 
   if($lval eq "x") { $nlval=87; } 
   if($lval eq "y") { $nlval=88; } 
   if($lval eq "z") { $nlval=89; } 
   if($lval eq "z") { $nlval=90; } 
   if($lval eq "\{") { $nlval=91; } 
   if($lval eq "\|") { $nlval=92; } 
   if($lval eq "\}") { $nlval=93; } 
   if($lval eq "\~") { $nlval=94; } 

   #print "$lval\n"; 

   return $nlval; 
}
