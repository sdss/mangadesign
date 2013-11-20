#!/usr/bin/perl -w
#
#writes out plateDefinition files
#APOGEE will add their inputs and priorities in later

use strict;
use PDL;
use PDL::Graphics::PGPLOT::Window;
use PDL::NiceSlice;

#inputs
my $inputdir = "inputs/manga/2013.02.x.manga/";

#read in plate centers for current run
my ($rap,$decp,$designID,$plateID,$locationID,$fld) = rcols 'MANGA-plates.dat',1,2,3,4,5,{perlcol=>[0]};

foreach my $i (0 .. $rap->nelem-1)
  {
    #definitions dir
    my $plo = $designID->at($i);
    my $desidst = substr $plo, 0,-2;
    my $defdir = "definitions/".$desidst."XX/";

    #check to see if directory exists
    #if not creat it
    if (-e $defdir){}else{mkdir $defdir;}

    my $name = $defdir."plateDefinition-".$plo.".par";
    my $inputname = "manga/2013.02.x.manga/plateInput-".$plo.".par";

    my $ra = $rap->at($i);
    my $dec = $decp->at($i);


    #check to see if definition file already exists
    #if it does ask if we want to overwrite it
    if (-e $name)
      {
	print "$name already exists. Do you wish to overwrite it (y/n)?\n";
	my $text = <STDIN>;
	chomp $text;
	if ($text =~ /y/){unlink $name;}else{next;}
      }

    print "writing $name\n";


    open OUT, "> $name";

    print OUT "\# Basic design info
designID $plo
platedesignversion v1
plateType APOGEE-MANGA
plate_lead MANGA
 
\# What are the pointings and offsets? 
raCen $ra
decCen $dec
 
\# Inputs and priority order of each. 
nInputs 1
plateInput $inputname

priority 1

";
}  
