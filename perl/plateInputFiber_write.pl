#!/usr/bin/perl -w
#
#writes out plateInputFiber yanny file
#run after plateInput_write.pl

use strict;
use PDL;
use efftickle;
require "yanny_sub.pl";

#define inputs directory
my $dir = $ENV{MANGA_PLATES}."/inputs/manga/2013.02.x.manga/";

#loop over each plateInput file
my $pattern = $dir."plateInput-*.par";
foreach my $name (glob $pattern)
  {
    print "Reading $name\n";
    my %par = read_yanny($name);
    
    #Want one line per fiber
    #Columns for mangaID, ifuDesign, ifuDesignSize, manga_target1/manga_target2 as for plateInput
    #provide empty fiber mags

    #define data structure
    my %instr1 = %{$par{STRUCT1}};#input struct1
    
    #loop over each IFU and create one entry per fiber
    my (@ra,@dec,@mangaID,@ifuDesign,@ifuDesignSize,@manga_target1,@manga_target2);
    my $ifu = 0;
    foreach my $nfib (@{$instr1{IFUDESIGNSIZE}})
      {
	foreach my $i (0 .. $nfib-1)
	  {
	    push(@ra,${$instr1{RA}}[$ifu]);
	    push(@dec,${$instr1{DEC}}[$ifu]);
	    push(@mangaID,${$instr1{MANGAID}}[$ifu]);
	    push(@ifuDesign,${$instr1{IFUDESIGN}}[$ifu]);
	    push(@ifuDesignSize,${$instr1{IFUDESIGNSIZE}}[$ifu]);
	    push(@manga_target1,${$instr1{MANGA_TARGET1}}[$ifu]);
	    push(@manga_target2,${$instr1{MANGA_TARGET2}}[$ifu]);
	  }
	$ifu++;
      }

    #now create empty mag array of arrays
    my @mag = dog(zeroes(5,$#ra+1));
    my @fib2mag;
    foreach my $i (0 .. $#mag)
      {
	my @m = list($mag[$i]);
	$fib2mag[$i] = \@m;
      }

    my %fibstr1;
    $fibstr1{RA} = \@ra;
    $fibstr1{DEC} = \@dec;
    $fibstr1{MANGAID} = \@mangaID;
    $fibstr1{IFUDESIGN} = \@ifuDesign;
    $fibstr1{IFUDESIGNSIZE} = \@ifuDesignSize;
    $fibstr1{FIBER2MAG} = \@fib2mag;
    $fibstr1{MANGA_TARGET1} = \@manga_target1;
    $fibstr1{MANGA_TARGET2} = \@manga_target2;

    #create hdr
    my %hdr; 
    $hdr{"TTYPE".1} = 'RA';
    $hdr{"TFORM".1} = '1D';
    $hdr{"TTYPE".2} = 'DEC';
    $hdr{"TFORM".2} = '1D';
    $hdr{"TTYPE".3} = 'MANGAID';
    $hdr{"TFORM".3} = '1D';
    $hdr{"TTYPE".4} = 'IFUDESIGN';
    $hdr{"TFORM".4} = '1D';
    $hdr{"TTYPE".5} = 'IFUDESIGNSIZE';
    $hdr{"TFORM".5} = '1D';
    $hdr{"TTYPE".6} = 'FIBER2MAG';
    $hdr{"TFORM".6} = '5D';
    $hdr{"TTYPE".7} = 'MANGA_TARGET1';
    $hdr{"TFORM".7} = '1D';
    $hdr{"TTYPE".8} = 'MANGA_TARGET2';
    $hdr{"TFORM".8} = '1D';

    #hash describing data types  
    my %fibstr1type = &structtype(\%fibstr1,\%hdr);
    #hash describing data format
    my %symbols1 = &symbols('STRUCT1',\%fibstr1,\%fibstr1type,\%hdr);

    #hash to output
    my %data = ('pointing',$par{pointing},'offset',$par{offset},'racen',$par{racen},'deccen',$par{deccen},'instrument',$par{instrument},'locationid',$par{locationid},'designid',$par{designid},'tileid',$par{tileid}, 'STRUCT1',\%fibstr1,'symbols',\%symbols1);

    my $outname = "plateInputFiber-".$par{designid}.".par";

    write_yanny($dir.$outname,\%data);
  }
