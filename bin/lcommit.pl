#!/usr/bin/perl -sw
#
# computes local commit for version file
#
#---------------------------------------------

use strict;

$::branch = "" unless $::branch;
$::reset = "" unless $::reset;

while(<>) {

  if( /version\s*=(.*)/ ) {	# looks for line "version ="
    $::version = strip($1);
    print STDERR "version found: $::version\n";
  }

  if( /(.*)lcommit\s*=(.*)/ ) {	# looks for line "lcommit ="
    $::lcommit = strip($2);
    print STDERR "local commit found: $::lcommit\n";
    my $string = make_lcommit();
    print STDERR "local commit: $string\n";
    $_ = "$1" . "lcommit = '" . $string . "'\n";
  }

  print;
}

#---------------------------------------------

sub make_lcommit {

  my @s = split("-",$::lcommit);	# local commit is vers-branch-number

  my $vers = $s[0];
  my $branch = $s[1];
  my $number = $s[2];

  $number = 0 if $::reset;
  $number = 0 unless $number;
  $number = 0 if $vers ne $::version;	# if version has changed start from 0

  $number++;

  my $string = "$::version" . "-" . "$::branch" . "-" . "$number";

  return $string;
}

sub strip {				# strips white space and "'" from string

  my $string = shift;

  $string =~ s/^\s*\'\s*//;
  $string =~ s/\s*\'\s*$//;

  return $string;
}

#---------------------------------------------

