#!/usr/bin/perl -w

## Calculation of lag time, r and carrying capacity 

use Statistics::R;
$R=Statistics::R->new();

$file[0]="ModelFit/C_Cas.txt";
$file[1]="ModelFit/C_FCN.txt";
$file[2]="ModelFit/C_H2O2_48h.txt";
$file[3]="ModelFit/C_H2O2_94h.txt";
$file[4]="ModelFit/C_NaCl_48h.txt";
$file[5]="ModelFit/C_SCD_40h.txt";
$file[6]="ModelFit/C_YPD_24h.txt";
$file[7]="ModelFit/F_Cas.txt";
$file[8]="ModelFit/F_FCN.txt";
$file[9]="ModelFit/F_H2O2_48h.txt"; ##
$file[10]="ModelFit/F_H2O2_91h.txt"; ##
$file[11]="ModelFit/F_NaCl.txt";
$file[12]="ModelFit/F_SCD_40h.txt";
$file[13]="ModelFit/F_YPD_24h.txt";
$file[14]="ModelFit/G_Cas.txt"; ##
$file[15]="ModelFit/G_FCN.txt";
$file[16]="ModelFit/G_H2O2.txt";
$file[17]="ModelFit/G_NaCl.txt";
$file[18]="ModelFit/G_SCD_40h.txt";
$file[19]="ModelFit/G_YPD_24h.txt";
$file[20]="ModelFit/Glycerol_growth_curve.txt"; 

for($i=0;$i<21;$i++)
{
  print "$file[$i]\n";

  @qq=split(/\//,$file[$i]); 
  $wr="Estimated_parameters/Estimates_"."$qq[1]"; 
  open(WR,">$wr") or die; 
  undef $wr; 
  undef(@qq);
  ## READ FILE
  open(FP,$file[$i]) or die;
  $ln=0; 
  while($fp=<FP>)
  {
    chomp($fp);
    @arr=split(/\s+/,$fp);
    #if($arr[0]!~/H06/) { undef(@arr); next; }
    if($fp=~/^Time/) 
    { 
       $nn=@arr; 
       for($k=1;$k<$nn;$k++)
       {
	 $head[$k]=$arr[$k]; 
       }
       undef(@arr); 
       next; 
    }
    for($k=0;$k<$nn;$k++)
    {
       $list[$k][$ln]=$arr[$k]; 
    }
    undef(@arr);
    $ln++; 
  }
  close(FP);

  ## ARRANGE DATA

  for($k=1;$k<$nn;$k++)
  { 
     $stp=20; 
     $cnt=0;

     for($li=0;$li<$ln;$li++) 
     {
         $xval[$cnt]=$list[0][$li];
         $yval[$cnt]=$list[$k][$li];
	 $cnt++; 
     }

     $bmax=0; $amax=0;
     $xmax=0; $ymax=0;
     $satmax=0; 

     $cc=0; 

       for($up=$stp;$up<$ln-$stp;$up++)
       {
          $sumx=0; $sumy=0;
          for($np=$up-$stp;$np<$up+$stp;$np++)
          {
             $sumx+=$xval[$np];
             $sumy+=$yval[$np];
          }
          $sumx/=($stp*2);
          $sumy/=($stp*2);

          $sum1=0; $sum2=0;
          for($np=$up-$stp;$np<$up+$stp;$np++)
          {
            $sum1+=(($xval[$np]-$sumx)*($yval[$np]-$sumy));
            $sum2+=(($xval[$np]-$sumx)**2);
          }
          $b=sprintf("%0.5f",$sum1/$sum2);
          $a=$sumy-($b*$sumx);
	  $fit[$cc][0]=$xval[$up];
	  $fit[$cc][1]=$yval[$up];
	  $fit[$cc][2]=$b;
	  $cc++; 
	  #print "$xval[$up] $a $b\n"; 
          if($b>$bmax) { $bmax=$b; $amax=$a; $xmax=$xval[$up]; $ymax=$yval[$up]; }
	  if($satmax<$yval[$up]) { $satmax=$yval[$up]; }
          undef $b;
          undef $a;
       }  

       #$xlag=(-$amax/$bmax); 
       
       $nc=0; 
       for($g=0;$g<$cc;$g++)
       {
	  if(abs($fit[$g][2]-$bmax)<=0.01)
	  {
	     $newx[$nc]=$fit[$g][0]; 
	     $newy[$nc]=$fit[$g][1]; 
	     $nc++;
	  }
       }

       ##print "$nc\n"; 

       $b=0; $a=0;

       $sumx=0; $sumy=0;
       for($np=0;$np<$nc;$np++)
       {
          $sumx+=$newx[$np];
          $sumy+=$newy[$np];
       }
       $sumx/=($nc);
       $sumy/=($nc);

       $sum1=0; $sum2=0;
       for($np=0;$np<$nc;$np++)
       {
          $sum1+=(($newx[$np]-$sumx)*($newy[$np]-$sumy));
          $sum2+=(($newx[$np]-$sumx)**2);
       }
       $b=sprintf("%0.5f",$sum1/$sum2);
       $a=$sumy-($b*$sumx);
       ##print "$xval[$up] $a $b\n"; 
         
       $xlag=sprintf("%0.2f",(-$a/$b)); 

       undef(@newy);
       undef(@newx);


       $nc=0; 
       $sumymax=0; 
       for($g=0;$g<$cc;$g++)
       {
	  if(abs($fit[$g][1]-$satmax)<=0.01)
	  {
	     $sumymax+=$fit[$g][1]; 
	     $nc++;
	  }
       }
       $sumymax=sprintf("%0.3f",$sumymax/$nc); 

       $tsat=sprintf("%0.2f",($sumymax-$a)/$b); 

       undef $a; 
       undef $b; 


       $R->start_sharedR; 

       $R->set('XAR',\@xval); 
       $R->set('YAR',\@yval); 
       
       $stod=0;
       for($u=0;$u<5;$u++)
       {
	  $stod+=$yval[$u]; 
       }
       $stod/=10; 

       $R->set('stod',$stod); 
       $R->set('satod',$sumymax); 

       $grsl=0.2; 
       if($file[$i]=~/SCD/ || $file[$i]=~/YPD/) { $grsl=0.4; } 
       #$grsl=$bmax*3; 
       $R->set('grsl',$grsl); 
       ##print "$head[$k] $stod $sumymax $grsl\n"; 

       $R->run(qq'm=nls(YAR ~ a/(1+((a-b)/b)*exp(-c*XAR)),start=list(a=satod,b=stod,c=grsl),control=list(warnOnly=T,minFactor=0.000001,maxiter=100))'); 

       $R->run(q'k=summary(m)[[10]][[1]]'); 
       $R->run(q'N0=summary(m)[[10]][[2]]'); 
       $R->run(q'r=summary(m)[[10]][[3]]'); 

       $kk=$R->get('k');
       $n0=$R->get('N0');
       $rr=$R->get('r');

       $R->run(q'rm(k)'); 
       $R->run(q'rm(N0)'); 
       $R->run(q'rm(r)'); 
       $R->run(q'rm(m)'); 
       $R->run(q'rm(XAR)');
       $R->run(q'rm(YAR)');
       $R->run(q'rm(stod)');
       $R->run(q'rm(satod)');
       $R->run(q'rm(grsl)');

       $R->stopR; 

       undef(@yval);
       undef(@xval);

       #print "$bmax $rr\n"; 
       print WR "$head[$k] Lagtime: $xlag(h) k: $kk N0: $n0 r: $rr  TimeToStationary: $tsat(h) MaxSlope: $bmax\n";

       undef $xlag; 
       undef $bmax; 
       undef $sumymax; 
       undef $kk;
       undef $n0; 
       undef $rr; 
       undef $tsat; 
       undef(@fit); 
       #exit();
  }
  close(WR); 

  ##print "$bmax $xmax\n"; 
  undef(@head); 
  undef(@list);  
}

undef(@file); 
