#!/usr/bin/perl -s
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# splits xrange in lower and upper part
#
# year is needed, other parts can be obmitted
# format: YYYY-MM-DD::hh:mm:ss
#
#---------------------------------------------

use lib ("$ENV{SHYFEMDIR}/lib/perl","$ENV{HOME}/shyfem/lib/perl");

use date;

$::date = new date;

$::round = 0 unless $::round;

my $xrange = $ARGV[0];

#print STDERR "xrange (perl): $xrange\n";

my $xhigh = "";
my $xlow = "";

if( $xrange =~ /^\s*$/ ) {
  ;
} elsif( $xrange =~ /^:(.*)$/ ) {		# only high
  $xhigh = regular($1);
} elsif( $xrange =~ /^(.*):$/ ) {		# only low
  $xlow = regular($1);
} elsif( $xrange =~ /^(.*):(.*)::(.*)$/ ) {	# :: is from high
  $xlow = regular($1);
  $xhigh = regular("$2::$3");
} elsif( $xrange =~ /^(.*):(.*)-(.*)$/ ) {	# - is in date of high
  $xlow = regular($1);
  $xhigh = regular("$2-$3");
} else {
  die "cannot parse xrange: $xrange\n";
}

if( $round ) {
  print STDERR "rounding...\n";
  print STDERR "$xlow:$xhigh\n";
  ($xlow,$xhigh) = round_xrange($xlow,$xhigh);
  print "$xlow:$xhigh\n";
} else {
  print "$xlow=$xhigh\n";
}

#---------------------------------------------

sub round_xrange {

  my ($xl,$xh) = @_;

  my $sep = 0;

  my $ls = $::date->unformat_abs($xl);
  my $hs = $::date->unformat_abs($xh);

  my $secs = $hs - $ls;
  print STDERR "period $ls $hs - $secs\n";

  if( $secs < 60 ) {
    $sep = 0;
    return ($xl,$xh);
  } elsif( $secs < 600 ) {
    $sep = 0;
    ($xl,$xh) = round_secs($ls,$hs,60);
    return ($xl,$xh);
  } elsif( $secs < 3600 ) {
    $sep = 300;
    ($xl,$xh) = round_secs($ls,$hs,$sep);
    return ($xl,$xh);
  }

  my $hours = $secs / 3600;
  if( $hours < 24 ) {
    $sep = 3600;
    ($xl,$xh) = round_secs($ls,$hs,$sep);
    return ($xl,$xh);
  }

  my $days = $secs / 86400;
  if( $days < 4 ) {
    $sep = 43200;
    ($xl,$xh) = round_secs($ls,$hs,$sep);
    return ($xl,$xh);
  } elsif( $days < 10 ) {
    $sep = 86400;
    ($xl,$xh) = round_secs($ls,$hs,$sep);
    return ($xl,$xh);
  } elsif( $days < 30 ) {
    $sep = 2*86400;
    ($xl,$xh) = round_secs($ls,$hs,$sep);
    return ($xl,$xh);
  }

# month 6 12
# years 4 over

  return ($xl,$xh);
}

sub round_secs {

  my ($ls,$hs,$r) = @_;

  my $lsec = int($ls/$r)*$r;
  my $xl = $::date->format_abs($lsec);
  
  my $hsec = int($hs/$r)*$r;
  $hsec += $r if $hsec < $hs;
  my $xh = $::date->format_abs($hsec);
  
  my $sec = $hs - $ls;;
  my $rsec = $hsec - $lsec;;

  print STDERR "rounded secs: $sec - $rsec\n";

  return ($xl,$xh);
}

#---------------------------------------------

sub regular {

  my $string = shift;

  my ($dd,$tt);
  my ($y,$m,$d);
  my ($h,$M,$s);
  my ($out);

  if( $string =~ "::" ) {
    ($dd,$tt) = split(/::/,$string);
  } elsif( $string =~ /:/ ) {
    $tt = $string;
  } elsif( $string =~ /-/ ) {
    $dd = $string;
  } elsif( $string =~ /^(\d\d\d\d)$/ ) {
    $dd = $string;
  } else {
    die "cannot regularize date (something is missing): $string\n";
  }

  if( not $dd ) {
    die "missing date in string: $string\n";
  } elsif( $dd =~ /^(\d\d\d\d)$/ ) {
    $y = $1;
  } elsif( $dd =~ /^(\d\d\d\d)-(\d\d)$/ ) {
    $y = $1;
    $m = $2;
  } elsif( $dd =~ /^(\d\d\d\d)-(\d\d)-(\d\d)$/ ) {
    $y = $1;
    $m = $2;
    $d = $3;
  } else {
    die "cannot parse date: $dd ($string)\n";
  }

  $m = "01" unless $m;
  $d = "01" unless $d;

  if( not $tt ) {
    ;
  } elsif( $tt =~ /^(\d\d)$/ ) {
    $h = $1;
  } elsif( $tt =~ /^(\d\d):(\d\d)$/ ) {
    $h = $1;
    $M = $2;
  } elsif( $tt =~ /^(\d\d):(\d\d):(\d\d)$/ ) {
    $h = $1;
    $M = $2;
    $s = $3;
  } else {
    die "cannot parse time: $tt ($string)\n";
  }

  $h = "00" unless $h;
  $M = "00" unless $M;
  $s = "00" unless $s;

  #print STDERR "regular: $string $dd $tt\n";
  #print STDERR "regular: $y $m $d :: $h $M $s\n";

  return "$y-$m-$d\:\:$h:$M:$s";
}

#---------------------------------------------

