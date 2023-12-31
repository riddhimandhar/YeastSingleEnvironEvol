#!/usr/bin/perl -w 
#
## Merge alignment data and mutation data from PE read mapping 

$file[0]="MutationsFromAlignment_R0-1.txt"; 
$file[1]="MutationsFromAlignment_R0-2.txt"; 
$file[2]="MutationsFromAlignment_R0-3.txt"; 
$file[3]="MutationsFromAlignment_R100-1.txt"; 
$file[4]="MutationsFromAlignment_R100-2.txt"; 
$file[5]="MutationsFromAlignment_R100-3.txt"; 


for($gi=0;$gi<6;$gi++)
{
   print "$file[$gi]\n"; 

   @we=split(/\_/,$file[$gi]); 
   @sq=split(/\./,$we[1]); 
   $wrf1=$sq[0]."_mutations_consPE.txt"; 
   undef(@sq); 
   undef(@we); 

   open(WR,">../results/Mutations_consolidatedPE/$wrf1") or die; 
   open(FP,"../results/$file[$gi]") or die; 
   $ln=0; $flg=0; 
   while($fp=<FP>) 
   {
     $ln++; 
     chomp($fp); 
     @arr=split(/\t/,$fp); 
     $no=@arr; 

     $cc=($ln % 2); 
     if($cc==0) 
     {
        $temp[$cc][0]=$arr[5]; 
        $temp[$cc][1]=$arr[6]; 
	@qw1=split(/\:/,$arr[3]);
	@qw2=split(/\:/,$arr[4]);
        $temp[$cc][2]=$qw1[1]; 
        $temp[$cc][3]=$qw2[1]; 
	undef(@qw1);
	undef(@qw2); 

	if($no>7) 
	{
	   $temp[$cc][4]=$arr[7]; 
	}
	else
	{
	   $temp[$cc][4]="X"; 
	}

	if(($temp[0][2] eq $temp[0][3]) && ($temp[1][2] eq $temp[1][3]))
	{
	   undef(@temp); 
	   next; 
	}
        

	$ov1=0; $ov2=0; 
        if($temp[0][0]>=$temp[1][0] && $temp[0][0]<=$temp[1][1])
        {
            $ov1=$temp[0][0]; 
	}	
        if($temp[1][0]>=$temp[0][0] && $temp[1][0]<=$temp[0][1])
        {
            $ov1=$temp[1][0]; 
	}	
        if($temp[0][1]>=$temp[1][0] && $temp[0][1]<=$temp[1][1])
        {
            $ov2=$temp[0][1]; 
	}	
        if($temp[1][1]>=$temp[0][0] && $temp[1][1]<=$temp[0][1])
        {
            $ov2=$temp[1][1]; 
	}	

	if($ov1!=0 && $ov2!=0)
	{
	  $st=0; $en=0;  
        
	  if($temp[0][0]<=$temp[1][0]) { $st=$temp[0][0]; } 
	  elsif($temp[0][0]>$temp[1][0]) { $st=$temp[1][0]; } 

	  if($temp[0][1]>$temp[1][1]) { $en=$temp[0][1]; }
	  elsif($temp[0][1]<=$temp[1][1]) { $en=$temp[1][1]; }
	  $st.=",0";
	  $en.=",0";
	}
	else 
	{ 
          $st=$temp[0][0].",$temp[1][0]"; 
	  $en=$temp[0][1].",$temp[1][1]"; 
        }

	##print "$arr[0] $arr[1] R1: $temp[0][0] $temp[0][1]   R2: $temp[1][0]  $temp[1][1]\tST:$st EN:$en\tOV1:$ov1 OV2:$ov2\n"; 
	
	#if($temp[0][4] eq "X" && $temp[1][4] eq "X")
	#{
            $readlist{$arr[0]}[0]=0; 
	    $readlist{$arr[0]}[1]=$arr[1];  
	    $readlist{$arr[0]}[2]=$st;  
	    $readlist{$arr[0]}[3]=$en;  
	    $readlist{$arr[0]}[4]=$ov1;  
	    $readlist{$arr[0]}[5]=$ov2;  
	    #}

	@qw1=split(/\,/,$temp[0][4]); 
	$nn1=@qw1; 
	@qw2=split(/\,/,$temp[1][4]); 
	$nn2=@qw2; 

	if($temp[0][4] ne "NA" && $temp[0][4] ne "X") 
	{
	   for($i=0;$i<$nn1;$i++) 
	   {
		   ##print "$i $nn1 $qw1[$i]\n"; 
	      @er=split(/\s+/,$qw1[$i]); 
	      $lid=$er[0]."*$er[2]"; 
	      if(exists $lcdup{$lid})
	      { 
                 $pp=$lcdup{$lid}; 
		 @we=split(/\s+/,$qw1[$pp]); 
		 $we[1].=$er[1]; 
		 $we[3]=sprintf("%0.0f",($er[3]+$we[3])/2);
		 #print "$arr[0] $qw1[$pp] $nn1\n"; 
		 $qw1[$pp]="$we[0] $we[1] $we[2] $we[3]"; 
		 undef(@we);

		 #print "$arr[0] $qw1[$pp] $nn1\n"; 

		 for($k=$i;$k<$nn1-1;$k++)
		 {
	            $qw1[$k]=$qw1[$k+1]; 
		 } 
		 $nn1--; 
		 #print "$nn1\n"; 
		 undef $pp; 
		 if($i<$nn1) { redo; }
	      }
	      if($er[2] eq "-" && !exists $lcdup{$lid}) 
	      {
		 $lcdup{$lid}=$i; 
	      }
	      undef $lid; 
	      undef(@er);
           }	   
	   undef %lcdup; 

	   for($i=0;$i<$nn1;$i++) 
	   {
	      @er=split(/\s+/,$qw1[$i]); 
	      print "$arr[0] $er[0] $er[1] $er[2] $er[3]\n"; 
	      $id=$er[0]."*$er[1]*$er[2]"; 
	      if($er[0]>=$ov1 && $er[0]<=$ov2) 
	      {
	         $checkmut{$id}[0]=0; 
	         $checkmut{$id}[1]=$er[3]; 
	      }
	      else
	      {
		 if($nn1<=5) 
		 {
	           $checkmut{$id}[0]=1; 
	           $checkmut{$id}[1]=$er[3]; 
		 }
	      } 
	      undef $id; 
	      undef(@er); 
	   }
	}

	if($temp[1][4] ne "NA" && $temp[1][4] ne "X") 
	{
	   for($i=0;$i<$nn2;$i++) 
	   {
	      @er=split(/\s+/,$qw2[$i]); 
	      $lid=$er[0]."*$er[2]"; 
	      if(exists $lcdup{$lid})
	      { 
                 $pp=$lcdup{$lid}; 
		 @we=split(/\s+/,$qw2[$pp]); 
		 $we[1].=$er[1]; 
		 $we[3]=sprintf("%0.0f",($er[3]+$we[3])/2);
		 #print "$arr[0] $qw1[$pp] $nn1\n"; 
		 $qw2[$pp]="$we[0] $we[1] $we[2] $we[3]"; 
		 undef(@we);

		 #print "$arr[0] $qw1[$pp] $nn1\n"; 

		 for($k=$i;$k<$nn2-1;$k++)
		 {
	            $qw2[$k]=$qw2[$k+1]; 
		 } 
		 $nn2--; 
		 #print "$nn1\n"; 
		 undef $pp; 
		 if($i<$nn2) { redo; }
	      }
	      if($er[2] eq "-" && !exists $lcdup{$lid}) 
	      {
		 $lcdup{$lid}=$i; 
	      }
	      undef $lid; 
	      undef(@er);
           }	   
	   undef %lcdup; 

	   for($i=0;$i<$nn2;$i++) 
	   {
	      @er=split(/\s+/,$qw2[$i]); 
	      $id=$er[0]."*$er[1]*$er[2]"; 
              if(exists $checkmut{$id}) 
	      {
	         print ":$arr[0] $er[0] $er[1] $er[2] $er[3]\n"; 

		 if($er[3] ne "NA" && $checkmut{$id}[1] ne "NA") 
		 {
		    $checkmut{$id}[1]=sprintf("%0.0f",($checkmut{$id}[1]+$er[3])/2); 
		 }
		 else
		 {
		    $checkmut{$id}[1]="NA"; 
		 }
		 if($nn2<=5) 
		 {
		   $checkmut{$id}[0]=2; 
		 }
	      }
	      else
	      {
		 if($er[0]<$ov1 || $er[0]>$ov2)
		 {
	           if($nn2<=5) 
		   { 
	             $checkmut{$id}[0]=1; 
	             $checkmut{$id}[1]=$er[3]; 
		   }
		 }
	      }
	      undef $id; 
	      undef(@er); 
	   }
	}

	##print "$nn1 $nn2\n"; 
       
	undef $nn1; 
	undef(@qw1); 
	undef $nn2; 
	undef(@qw2); 

	foreach $key (keys %checkmut) 
	{
	   @we=split(/\*/,$key); 
	   if($checkmut{$key}[0]>=1)
	   {
	      if(exists $readlist{$arr[0]}) 
	      {
		  $cnt=$readlist{$arr[0]}[0]; 
		  $readlist{$arr[0]}[6][$cnt]=$we[0];  
		  $readlist{$arr[0]}[7][$cnt]=$we[1];  
		  $readlist{$arr[0]}[8][$cnt]=$we[2];  
		  $readlist{$arr[0]}[9][$cnt]=$checkmut{$key}[1];  
		  $readlist{$arr[0]}[10][$cnt]=$checkmut{$key}[0];  
                  $readlist{$arr[0]}[0]++; 
		  undef $cnt; 
	      }
	      #else
	      #{
	      #    $readlist{$arr[0]}[0]=1; 
#		  $readlist{$arr[0]}[1]=$arr[1];  
#		  $readlist{$arr[0]}[2]=$st;  
#		  $readlist{$arr[0]}[3]=$en;  
#		  $readlist{$arr[0]}[4]=$ov1;  
#		  $readlist{$arr[0]}[5]=$ov2;  
#		  $readlist{$arr[0]}[6][0]=$we[0];  
#		  $readlist{$arr[0]}[7][0]=$we[1];  
#		  $readlist{$arr[0]}[8][0]=$we[2];  
#		  $readlist{$arr[0]}[9][0]=$checkmut{$key}[1];  
#		  $readlist{$arr[0]}[10][0]=$checkmut{$key}[0];  
#	      }
	      #print WR "$arr[0]\t$checkmut{$key}[0]\t$we[0]\t$we[1]\t$we[2]\t$checkmut{$key}[1]\n"; 
	   }
	   undef(@we); 
	}
	undef %checkmut; 
	undef $key; 

	foreach $key (sort keys %readlist) 
	{
	   print WR "$key\t$readlist{$key}[1]\t$readlist{$key}[2]\t$readlist{$key}[3]\t$readlist{$key}[4]\t$readlist{$key}[5]\t"; 
	   for($k=0;$k<$readlist{$key}[0];$k++) 
	   {
	      print WR "$readlist{$key}[6][$k] $readlist{$key}[7][$k] $readlist{$key}[8][$k] $readlist{$key}[9][$k] $readlist{$key}[10][$k],"; 
	   }	
	   print WR "\n"; 
	}

	undef %readlist; 
	undef(@temp); 
     } 
     else
     {
	  $temp[$cc][0]=$arr[5]; 
	  $temp[$cc][1]=$arr[6]; 
	  @qw1=split(/\:/,$arr[3]);
	  @qw2=split(/\:/,$arr[4]);
          $temp[$cc][2]=$qw1[1]; 
          $temp[$cc][3]=$qw2[1]; 
	  undef(@qw1);
	  undef(@qw2); 

	if($no>7) 
	{
	     $temp[$cc][4]=$arr[7]; 
	}
	else
	{
	   $temp[$cc][4]="X"; 
	}
     }

     undef $no; 
     undef(@arr); 
   }
   close(WR); 
   close(FP); 

   undef $wrf1; 
}
undef(@file); 
