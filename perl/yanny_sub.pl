#!/usr/bin/perl -w
#
#creates the hashes required to write out yanny files using efftickle

use strict;
use Scalar::Util qw(looks_like_number);


sub structtype {
  #takes an input hash ref that contains the elements of a given structure
  #this hash has the format name and arrayref e.g. %struct1 = ('ra',\@ra,'dec',\@dec)
  #takes input hdr hash ref which gives order of hash elements
  #determins the type of variable in each of those elements
  #returns hash with the same keys and the variable type as the value
  #e.g. my %struct1type = &structtype(\%struct1,\%hdr);

  my %structtype;
  my ($struct,$hdr) = @_;
  my %struct = %$struct;
  my %hdr = %$hdr;
  my $nhdr = (scalar keys(%struct));#scalar keys(%struct) counts the number of keys.

  foreach my $i (1 .. $nhdr)
    {
      my $key = $hdr{"TTYPE".$i};

      #check to see if numeric
      my $num = looks_like_number(${$struct{$key}}[0]);
      if ($num != 0)
	{
	  #test to see if integer
	if (int(${$struct{$key}}[0]) ==  ${$struct{$key}}[0])
	  {
	    $structtype{$key} = 'int';
	  }
	else
	  {
	    $structtype{$key} = 'double';
	  }
      }
      #check to see if array
      elsif (ref ${$struct{$key}}[0] eq 'ARRAY')
	{
	  my $length = $#{${$struct{$key}}[0]} + 1;
	  $structtype{$key} = "double[".$length."]";
	}
      else{
	#find max size of character string
	my $maxlength = 0;
	foreach my $char (@{$struct{$key}})
	  {
	    my $len = length($char);
	    if ($len > $maxlength){$maxlength = $len;}
	  }
	$structtype{$key} = "char[".$maxlength."]";
      }
  }
%structtype;
}


sub symbols {
  #makes symbols hash which is used to define structure of yanny file
  #takes input of the name of the data structure and refs to the data structure and the structure type made by structype
  #takes input hdr hash ref which gives order of hash elements
  #e.g. my %symbols1 = &symbols('STRUCT1',\%struct1,\%struct1type,\%hdr);

  my ($structname,$struct,$structtype,$hdr) = @_;
  my %struct = %$struct;
  my %structtype = %$structtype;
  my %hdr = %$hdr;
  my $nhdr = (scalar keys(%struct));#scalar keys(%struct) counts the number of keys.

  #name of each element
  my @STRUCT;
  my $c = 0;
  foreach my $i (1 .. $nhdr)
    {
      my $key = $hdr{"TTYPE".$i};
      $STRUCT[$c] = $key;$c++;
    }

  #format of each element
  my @struct;
  $struct[0] = "typedef struct {\n"; 
  foreach my $i (1 .. $nhdr)
    {
      my $key = $hdr{"TTYPE".$i};
      $struct[0] = $struct[0]." $structtype{$key} $key;\n";
    }
  $struct[0] = $struct[0]."} ".$structname.";\n";

  my (@enum);
 
  my %symbols = ('struct',\@struct,$structname,\@STRUCT,'enum',\@enum);
  %symbols;
}

1
