#!/usr/bin/perl

#******************** perl script for psi4 dma output to pull multipole moments
#******************* 
#*******************   IMPORTANT:  Units should all be A.U., we need to convert coordinates from Angstrom to Bohr
$conv=1.88973;

$ifile=$ARGV[0];
open ifile, $ifile , or die;

do { $_=<ifile> } until(/Multipole moments/);

@atomtype;@xyz;
$i=0;
while(<ifile>)
{ 
$_=<ifile>;
if (/Total/){last}
@array=split; push( @atomtype , $array[0] );
push(@xyz , [ $array[3] , $array[6] , $array[9] ] );
# make sure rank 4
$_=<ifile>;@array=split;my@tensor=();
$rank=$array[3]; if ( $rank !~/4/ ){ die "must have rank 4 multipoles\n" }
# rank 0
$_=<ifile>;@array=split;push(@tensor,[ $array[2] ] ) ;
# rank 1
$_=<ifile>;@array=split;push(@tensor,[ $array[5], $array[8], $array[11]] ) ;
#rank 2
$_=<ifile>;@array=split;$_=<ifile>;@array1=split;
push(@tensor,[ $array[5], $array[8], $array[11] , $array1[2] , $array1[5] ] ) ;
#rank 3
$_=<ifile>;@array=split;$_=<ifile>;@array1=split;$_=<ifile>;@array2=split;
push(@tensor,[ $array[5], $array[8], $array[11] , $array1[2] , $array1[5], $array1[8] , $array2[2] ] ) ;
#rank 4
$_=<ifile>;@array=split;$_=<ifile>;@array1=split;$_=<ifile>;@array2=split;
push(@tensor,[ $array[5], $array[8], $array[11] , $array1[2] , $array1[5], $array1[8] , $array2[2] , $array2[5], $array2[8] ] ) ;

   # now print, IMPORTANT, convert coordinates from angstrom to Bohr
printf("%3s%16.6f%16.6f%16.6f%7d\n",$atomtype[$i], $xyz[$i][0]*$conv, $xyz[$i][1]*$conv, $xyz[$i][2]*$conv, 4 ) ;
   # print tensor components for each rank
for(my$j=0;$j<5;$j++)
{
    $ncomp= 2 * $j + 1 ;
    for(my$k=0;$k<$ncomp;$k++)
{
    printf("%9.6f\t", $tensor[$j][$k] );
}
print"\n";
}
 $i++;   
}
