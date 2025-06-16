#!/usr/bin/perl
#
# converts tabs to spaces
#
#---------------------------------------------------------

while(<>) {

  chomp;

  while( /\t/ ) {
    convert_tabs();
  }

  print "$_\n";
}

sub convert_tabs
{
  $i = index($_,"\t");

  if( $i < 0 ) {
    print STDERR "no tab found: $_\n";
  } else {
    $rel = $i % 8;
    $s = 8 - $rel;
    if( $s > 8 or $s < 0 ) { die "fatal error: $i $abs $rel $s $_\n"; }
    $spaces = " " x $s;
    s/\t/$spaces/;
  }
}

#---------------------------------------------------------

