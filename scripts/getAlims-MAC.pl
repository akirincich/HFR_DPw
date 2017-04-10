#!/usr/bin/perl -w
use File::Basename;
use Time::Local;
use POSIX;
use strict;

######################
## Local Variables  ##
######################

my $filename =  "";
my @data;
my $nAntenna=3;
my $ldrb;
my $ldlb;
my $rdlb;
my $rdrb;
my $noise_floor;
my $i;
my $k;
#####################
## Main Processing ##
#####################

if($filename = $ARGV[0]) {
   unless((-s $filename) && (-B $filename)){
	    die("Need to specify valid file.\n");
   }
# read and print header to STDOUT    
   @data = &HeaderRead($filename);
#   &printHeader(@data);

# use header variables to setup spectra read

   my $nRangeCells=$data[12];
   my $nDopplerCells=$data[11];
   my $nCsaKind=$data[2];
   my $headerLength=$data[15];
   my $posBragg=$data[19];
   my $negBragg=$data[20];

# find first order region and noise floor from resource fork

   my ($ldlb,$rdlb,$ldrb,$rdrb,$noise_floor)=getAlims($filename,$nDopplerCells);

# print output file	
   for ($i=0;$i<$nRangeCells;$i++){
	print  "$i $ldlb->[$i] $rdlb->[$i] $ldrb->[$i] $rdrb->[$i] \n" if ($ldlb->[$i]);
    }
}
#################
## Subroutines ##
#################

sub HeaderRead{
#####################
## Local Variables ##
#####################

my $file = $_[0];
my $buff = "";
my $nCsaFileVersion = 0;
my $nDateTime = 0;
my $nV1Extent = 0;
my $nCsaKind = 0;
my $nV2Extent = 0;
my $nSiteCodeName = "";
my $nV3Extent = 0;
my $nAverMinutes = 0;
my $bDeleteCSWhenUsed = 0;
my $bOverrideHeaderInfo = 0;
my $fStartFreqMHz = 0;
my $fRepFreqHz = 0;
my $fBandwidthKHz = 0;
my $bSweepUp = 0;
my $nDopplerCells = 0;
my $nRangeCells = 0;
my $nFirstRangeCell = 0;
my $fRangeCellDistKm = 0;
my $nV4Extent = 0;
my $file_position= 0;
my $fCenterFreqMHZ=0;
my $fBraggFreq=0;
my $fBraggWavelength=0;
my $nForwardBragg=0;
my $nReceedingBragg=0; 
###########################
## HeaderRead Processing ##
###########################
	
if(open(FILE, "<:raw", $file) or die "can't open $file: $!"){
	read(FILE, $buff, 10);
	($nCsaFileVersion, $nDateTime, $nV1Extent) = unpack("n1N1N1", $buff);
	if(scalar($nV1Extent) > 0){
	    read(FILE, $buff, 6);
	    ($nCsaKind, $nV2Extent) = unpack("n1N1", $buff);
	    if(scalar($nV2Extent) > 0){
		read(FILE, $buff, 8);
		($nSiteCodeName, $nV3Extent) = unpack("a4N1", $buff);
		if(scalar($nV3Extent) > 0){
		    read(FILE, $buff, 12);
		    ($nAverMinutes, $bDeleteCSWhenUsed, $bOverrideHeaderInfo) 
		      = unpack("N1N1N1", $buff);
		    read(FILE, $buff, 4);
		    $fStartFreqMHz=unpack "f",  $buff;
		    read(FILE, $buff, 4);
		    $fRepFreqHz=unpack "f",  $buff;
		    read(FILE, $buff, 4);
		    $fBandwidthKHz=unpack "f",  $buff;
		    read(FILE, $buff, 16);
		    ($bSweepUp, $nDopplerCells, $nRangeCells,$nFirstRangeCell) 
		      = unpack("N1N1N1N1", $buff);
		    read(FILE, $buff, 4);
		    $fRangeCellDistKm=unpack "f",  $buff;
		    read(FILE, $buff, 4);
		    $nV4Extent=unpack "N1", $buff; 
		    $file_position=tell(FILE);
		    close(FILE);
# find center freq = start_freq + bandwidth/1000/2 if sweep down
# center freq = start_freq - bandwidth/1000/2 if sweep up
		    $fCenterFreqMHZ=$fStartFreqMHz-$fBandwidthKHz/1000/2 unless ($bSweepUp);
		    $fCenterFreqMHZ=$fStartFreqMHz+$fBandwidthKHz/1000/2 if ($bSweepUp);
# find Bragg Frequency = +- sqrt(g/piL)
		    $fBraggWavelength=299792458/(1000000*$fCenterFreqMHZ);
		    $fBraggFreq=sqrt(9.86/(3.14159265*$fBraggWavelength));
		    my $dc=$nDopplerCells/2-1;
		    my $bragg_offset=$fBraggFreq*$dc;
		    $nReceedingBragg=sprintf("%3.3d",round($dc-$bragg_offset));
		    $nForwardBragg=sprintf("%3.3d",round($dc+$bragg_offset));
		    return ($nCsaFileVersion, $nDateTime, $nCsaKind, 
			$nSiteCodeName, $nAverMinutes, $bDeleteCSWhenUsed, 
			$bOverrideHeaderInfo, $fStartFreqMHz, $fRepFreqHz, 
			$fBandwidthKHz, $bSweepUp, $nDopplerCells, 
			$nRangeCells, $nFirstRangeCell, $fRangeCellDistKm,
			$file_position,$fCenterFreqMHZ,$fBraggFreq,$fBraggWavelength,
			$nForwardBragg,$nReceedingBragg);
		}
		else{
		    close(FILE);
		    return ($nCsaFileVersion, $nDateTime, $nCsaKind, $nSiteCodeName);
		}
	    }
	    else{
		close(FILE);
		return ($nCsaFileVersion, $nDateTime, $nCsaKind);
	    }
	}
	else{
	    close(FILE);
	    return ($nCsaFileVersion, $nDateTime);
	}
}
else{
	return;
}
}

sub printHeader {
#################################################################################
# Script for reading the header of codar CSS files.				#
#										#
# USAGE:									#
# codar_css_HeaderRead.pl <filename>						#
#										#
# Typing in the script name and the file name runs the script.  The script will #
# will die if the file name is not valid or no file name was given.  The script #
# opens the file and reads in the binary header information.  This script will  #
# only retrieve information up to ver. 4 for codar css files.			#
#	Version 1 ( @_ length = 2)						#
#        	$_[0] = nCsaFileVersion						#
#        	$_[1] = nDateTime						#
#       Version 2(length = 3)							#
#		$_[2] = nCsaKind						#
#       Version 3(length = 4)							#
#		$_[3] = nSiteCodeName						#
#       Version 4(length = 15) 							#
#		$_[4] = nAverMinutes						#
#        	$_[5] = bDeleteCSWhenUsed					#
#        	$_[6] = bOverrideHeaderInfo					#
#        	$_[7] = fStartFreqMHz						#
#        	$_[8] = fRepFreqHz						#
#        	$_[9] = fBandwidthKHz						#
#        	$_[10] = bSweepUp						#
#        	$_[11] = nDopplerCells						#
#        	$_[12] = nRangeCells						#
#        	$_[13] = nFirstRangeCell					#
#        	$_[14] = fRangeCellDistKm					#
#        	$_[15] = file_position  					#
#        	$_[16] = fCenterFreqMHZ 					#
#        	$_[17] = fBraggFreqHZ   					#    
#        	$_[18] = fBraggWavelength   					#    
#        	$_[19] = nForwardBragg                                          #
#        	$_[20] = nReceedingBragg                                        #
# 										#
# 2005/07/28 - Timothy Harris							#
# 2005/12/08 - TC update unpack of float values for little endian machines      #
#################################################################################

my @data=@_;
print("File version number:\t\t", $data[0], "\n");
print("Date Created:\t\t\t", scalar(localtime($data[1] - 2082844800)), "\n");
print("Type of CrossSpectra Data:\t", $data[2], "\n");
print("Site Name:\t\t\t", $data[3], "\n");
print("Coverage Time:\t\t\t", $data[4], " min\n");
print("\'Cross Spectra #\' deleted:\t");
if($data[5]){
	print("Yes\n");
}
else{
	print("No\n");
}
print("CSPro override header:\t\t");
if($data[6]){
	print("Yes\n");
}
else{
	print("No\n");
}
print("Transmit Center Freq:\t\t", $data[16], " MHz\n");
print("Transmit Start Freq:\t\t", $data[7], " MHz\n");
print("Transmit Sweep Rate:\t\t", $data[8], " Hz\n");
print("Transmit Bandwidth:\t\t", $data[9], " kHz\n");
print("Sweep Direction:\t\t");
if($data[10]){
	print("Up\n");
}
else{
	print("Down\n");
}
print("Number of Doppler Cells:\t", $data[11], "\n");
print("Number of Range Cells:\t\t", $data[12], "\n");
print("Index of First Range Cell\t", $data[13], "\n");
print("Distance between Range Cells:\t", $data[14], " km\n");
print("Bragg Frequency:\t\t", $data[17], " Hz\n");
print("Bragg Wavelength:\t\t", $data[18], " m\n");
print("Bragg Lines (doppler cell):\t", "$data[20]\t$data[19]","\n");
}

sub round{
my($number) = shift;
return int($number + .5 * ($number <=> 0));
}

sub max {
my($max) = shift;
my $count=0;
my $index=0;

foreach (@_) {
	if ($max<$_){
	    $max = $_;
	    $index = $count;
	}
	$count++;
}
return ($max,$index);
} 

sub min {
my($min) = shift;
my($count,$index);
$count=0;

foreach (@_) {
	if ($min>$_){
	    $min = $_;
	    $index = $count;
	}
	$count++;
}
return ($min,$index);
} 

sub getAlims{
my ($filename,$nDopplerCells)=@_;
my $buff=0;
my $rsrcfile=$filename.".rsrc";
my $i=0;
my $j=0;
my @Ldlb=();
my @Ldrb=();
my @Rdlb=();
my @Rdrb=();
my @NF=();

system("cp $filename/..namedfork/rsrc $rsrcfile");

open(FILE,"<:raw",$rsrcfile);
seek(FILE,260,0);
	read(FILE,$buff,8);
	my ($tnam,$num_rc)=unpack "a4N",$buff;
	read(FILE,$buff,4);
	my $first_rc=unpack "f", $buff;	
	read(FILE,$buff,4);
	my $bearing=unpack "f", $buff;	
	read(FILE,$buff,16);
	my ($start_rc,$r1,$r2,$r3)=unpack "NNNN",$buff;
	for ($i=$start_rc;$i<=$num_rc;$i++){
	    read(FILE,$buff,16);
	    ($Ldlb[$j],$Rdlb[$j],$Ldrb[$j],$Rdrb[$j])=unpack "NNNN", $buff;
	    $Ldlb[$j]--;
	    $Rdlb[$j]--;
	    $Ldrb[$j]=$Ldrb[$j]+$nDopplerCells/2;
	    $Rdrb[$j]=$Rdrb[$j]+$nDopplerCells/2;	    
	    $j++;
	}

    read(FILE,$buff,36);
	my ($dum,$tnam1,$num_rc1,$start_rc1,$num_ant,$r11)=unpack "Na4NNNN",$buff;
	$j=0;
	for ($i=$start_rc1;$i<=$num_rc1;$i++){
	    read(FILE,$buff,2);
	    $NF[0][$j]=unpack "s", $buff;
	    $j++;
	}
	$j=0;
	for ($i=$start_rc1;$i<=$num_rc1;$i++){
	    read(FILE,$buff,2);
	    $NF[1][$j]=unpack "s", $buff;
	    $j++;
	}
	$j=0;
	for ($i=$start_rc1;$i<=$num_rc1;$i++){
	    read(FILE,$buff,2);
	    $NF[2][$j]=unpack "s", $buff;
	    $j++;
	}
#	for ($i=$start_rc;$i<=$num_rc;$i++){
#	    print "$i $nf1[$i] $nf2[$i] $nf3[$i]\n";
#	}

system("rm $rsrcfile");
return (\@Ldlb,\@Rdlb,\@Ldrb,\@Rdrb,\@NF);
}

sub calcSNR{
	my ($ldlb,$rdlb,$ldrb,$rdrb,$noise_floor,@ndata)=@_;
	my $i=0;
	my $lbsnr=0;
	my $rbsnr=0;
# start with left bragg	
    my $sigsum=0;
    my $sigcount=0;
# find average signal over first order limit    
	for ($i=$ldlb;$i<=$rdlb;$i++){
		$sigsum=$sigsum+$ndata[$i];
		$sigcount++;
	}
	my $sigavg=$sigsum/$sigcount;
	$lbsnr=($sigavg-$noise_floor);
# start with right bragg	
    $sigsum=0;
    $sigcount=0;
# find average signal over first order limit    
	for ($i=$ldrb;$i<=$rdrb;$i++){
		$sigsum=$sigsum+$ndata[$i];
		$sigcount++;
	}
	$sigavg=$sigsum/$sigcount;
	$rbsnr=($sigavg-$noise_floor);
	my ($snr,$in)=sprintf("%3.1f",max($rbsnr,$lbsnr));
	return $snr;	
}
