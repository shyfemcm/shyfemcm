#!/usr/bin/perl -ws

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# joins columns
#
#------------------------------------------------------------------------

use strict;

$::h = 0 unless $::h;
$::help = 0 unless $::help;
$::help = 1 if $::h;

$::notime = 0 unless $::notime;
$::newfile = 0 unless $::newfile;

if( $::help || !@ARGV ) {
  print "Usage: joincol.pl [-help|-h] files\n";
  exit 0;
}

if( $::notime ) {
  die "*** cannot yet join columns with option -notime\n";
}

#------------------------------------------------------------------------

my @tref = ();
my @dref = ();

while( my $file = shift ) {

  my ($tref,$dref) = read_file($file);

  push(@tref,$tref);
  push(@dref,$dref);
}

my $ref;
my $nfiles = @tref;
$ref = $tref[0];
my @tbasic = @$ref;
$ref = $dref[0];
my @dbasic = @$ref;

while( defined (my $t0 = shift(@tbasic)) ) {

  my $d = shift(@dbasic);
  my $line = "$t0 $d";

  for(my $i=1;$i<$nfiles;$i++) {
    my $ref = $tref[$i];
    my $t = shift(@$ref);
    if( $t ne $t0 ) {
	print STDERR "different times: $t0 $t\n";
	die "error\n";
    }
    $ref = $dref[$i];
    my $d = shift(@$ref);
    $line .= " $d";
  }

  print "$line\n";
}

#------------------------------------------------------------------------

sub read_file {

  my $file = shift;
  my ($t,$r);
  my @rest;
  my @t = ();
  my @r = ();

  open(FILE,"<$file");

  while(<FILE>) {

    s/^\s+//;

    ($t,@rest) = split;
    $r = join(" ",@rest);

    push(@t,$t);
    push(@r,$r);
  }

  close(FILE);

  return (\@t,\@r);
}

#------------------------------------------------------------------------

