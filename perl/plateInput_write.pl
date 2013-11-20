#!/usr/bin/perl -w
#
#writes out plateInput yanny file
#written by David Wake 
#last update 11/14/13


use strict;
use PDL;
use efftickle;
require "/Users/daw/perl/sub/yanny_sub.pl";


my @ifus = (19,37,61,91,127);#ifu sizes


#define directorys
my $inputdir = "inputs/manga/2013.02.x.manga/";#inputs
 
#read in plate centers for current run
my ($rap,$decp,$designID,$plateID,$locationID,$fld) = rcols 'MANGA-plates.dat',1,2,3,4,5,{perlcol=>[0]};

#read in tiled galaxy catalog
#will replace this with tiled galaxy and std catalog
#this will include psf_mag for stds
my $catname = "/Users/daw/sdss_dr7/manga/sample_opt/MaNGA_targets_nsa_ran_tiles3_1.2_bun_2_4_4_2_5_Remaj_Ng_30000_np_1126.fits";
my $d = rfits($catname);
my %d = %$d;
my %hdr = %{$d{hdr}};
my $nhdr = (scalar keys(%d)) -2;#scalar keys(%d) counts the number of keys. Subtract 2 for tbl and hdr

#set up extra keywords
#currently these aren't in the catalog will add them later
if (defined $d{'MANGAID'}){}else
  {
    my @mangaid;
    foreach my $i (0 .. $d{CATIND}->nelem-1)
      {
	$mangaid[$i] = "1-".$d{CATIND}->at($i);
      }
    $d{'MANGAID'} = \@mangaid;
    $nhdr++;
    $hdr{"TTYPE".$nhdr} = 'MANGAID';
    $hdr{"TFORM".$nhdr} = '10A';
  }


if (defined $d{'MANGA_TARGET1'}){}else
  {
    $d{'MANGA_TARGET1'} = $d{PRIMARY} + 2*$d{SECONDARY};
    $nhdr++;
    $hdr{"TTYPE".$nhdr} = 'MANGA_TARGET1';
    $hdr{"TFORM".$nhdr} = '1D';
  }

if (defined $d{'MANGA_TARGET2'}){}else
  {
    $d{'MANGA_TARGET2'} = zeroes($d{RA}->nelem);;
    $nhdr++;
    $hdr{"TTYPE".$nhdr} = 'MANGA_TARGET2';
    $hdr{"TFORM".$nhdr} = '1D';
  }

my $priority =  1;
my $maxpriority =  2;
my $epoch =  2000;


#loop over each tile and write out plate input file
$d{IFUDESIGN} = zeroes($d{RA}->nelem-1);
$nhdr++;
$hdr{"TTYPE".$nhdr} = 'IFUDESIGN';
$hdr{"TFORM".$nhdr} = '1D';

 foreach my $i (0 .. $rap->nelem-1)
  {
    my $racen = $rap->at($i);
    my $deccen = $decp->at($i);

    my $locationid = $locationID->at($i);
    my $designid = $designID->at($i);#assigned by demitri

    my $name = "plateInput-".$designid.".par";
    my $in = which($d{LOCATIONID} == $locationid & $d{IFUDESIGNSIZE} > 0);

    #set up ifuDesign
    my $ifuid = zeroes($in->nelem);
    foreach my $size (@ifus)
      {
	if ( which($d{IFUDESIGNSIZE}->index($in) == $size)->nelem > 0)
	  {
	    $ifuid->where($d{IFUDESIGNSIZE}->index($in) == $size) .= sequence(which($d{IFUDESIGNSIZE}->index($in) == $size)->nelem)+1;
	  }
      }
    $d{IFUDESIGN}->index($in) .= $d{IFUDESIGNSIZE}->index($in)*100 + $ifuid;
   
    #define data structure
    my %struct1;

    #loop over hdr to keep in correct order
    foreach my $i (1 .. $nhdr)
      {
	my $key = $hdr{"TTYPE".$i};

	#check to see if $d{$key} is array
	my @d;
	if (ref $d{$key} eq 'ARRAY')
	  {
	    @d = @{$d{$key}}[list($in)];
	  }
	else
	  {
	    @d = list($d{$key}->index($in));
	  }
	$struct1{$key} = \@d;
      }    


    #hash describing data types  
    my %struct1type = &structtype(\%struct1,\%hdr);

    #hash describing data format
    my %symbols1 = &symbols('STRUCT1',\%struct1,\%struct1type,\%hdr);


    #hash to output
    my %data = ('pointing','1','offset',0,'racen',$racen,'deccen',$deccen,'instrument','MANGA','locationid',$locationid,'designid',$designid,'tileid','-1', 'STRUCT1',\%struct1,'symbols',\%symbols1);

    #check to see if input file already exists
    #if it does ask if we want to overwrite it
    if (-e $inputdir.$name)
      {
	print "$inputdir$name already exists. Do you wish to overwrite it (y/n)?\n";
	my $text = <STDIN>;
	chomp $text;
	if ($text =~ /y/){unlink $inputdir.$name;}else{next;}
      }

    print "writing $inputdir$name\n";
    write_yanny($inputdir.$name,\%data);
  }
