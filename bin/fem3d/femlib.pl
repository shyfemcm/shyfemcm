#!/usr/bin/perl 
#
# replaces femlib and fembin with lib and bin
#
#---------------------------------------------------

$^I = '.bak';

while(<>) {

  if( /SHYFEMDIR/ ) {
    s/femlib/lib/g;
    s/fembin/bin/g;
  }

  print;
}

