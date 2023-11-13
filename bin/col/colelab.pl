#!/usr/bin/perl -ws
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# elaborates file with time series
#
#------------------------------------------------------------------------

use strict;

$::h = 0 unless $::h;
$::help = 0 unless $::help;
$::help = 1 if $::h;

$::col = 0 unless $::col;
$::split = 0 unless $::split;
$::extract = 0 unless $::extract;
$::scale = 0 unless $::scale;
$::newfile = "col" unless $::newfile;
$::quiet = 0 unless $::quiet;
$::verbose = 0 unless $::verbose;

$::openfile = 0;
$::ncols = 0;

if( $::help || !@ARGV ) {
  print "Usage: elabcol.pl [-help|-h] [-options] file\n";
  exit 0;
}

#$::newfile = "col" unless $::newfile;

#------------------------------------------------------------------------

while(<>) {

  chomp;
  s/^\s+//;		#crop leading white space
  next if /^\#.*/;	#comment line -> ignore
  s/#.*//;		#comment -> delete
  #print STDERR "$_\n";

  my @data = split;

  store_cols(\@data);
}

if( $::split ) {
  split_cols();
} elsif( $::extract ) {
  extract_col($::extract);
} elsif( $::scale ) {
  my ($fact,$const) = make_scaling($::scale);
  scale_cols($::col,$fact,$const);
  write_cols();
} else {
  write_cols();
}

#close_files();

#------------------------------------------------------------------------

sub make_scaling {

  my $scale = shift;

  my ($fact,$const) = split(/,/,$scale);
  $fact = 1 unless defined $fact;
  $fact = eval( $fact );
  $const = 0 unless defined $const;
  $const = eval( $const );
  #print STDERR "scale: $::scale - $fact - $const\n";

  return ($fact,$const);
}

sub scale_cols {

  my ($col,$fact,$const) = @_;

  my $n = $::ncols;
  my $tcol = $::datacols[0];
  my $nt =  @$tcol;

  if( $col < 0 or $col >= $n ) {
    die "column number is not compatible: $col $n\n";
  }

  my $istart = 1;
  my $iend = $n;
  if( $col != 0 ) {
    $istart = $col;
    $iend = $col+1;
    unless( $::quiet ) {
      print STDERR "scaling column $col with fact=$fact and const=$const\n";
    }
  } else {
    unless( $::quiet ) {
      print STDERR "scaling all columns with fact=$fact and const=$const\n";
    }
  }

  for(my $i=$istart;$i<$iend;$i++) {
    my $rcol = $::datacols[$i];
    for(my $it=0;$it<$nt;$it++) {
      my $val = $rcol->[$it];
      $val = $val * $fact + $const;
      $rcol->[$it] = $val;
    }
  }
}

sub split_cols {

  my $n = $::ncols;
  my $tcol = $::datacols[0];
  my $nt =  @$tcol;

  for(my $i=1;$i<$n;$i++) {
    my $file = make_file_name($::newfile,$i,$n);
    open(FILE,">$file");
    write_col($nt,$tcol,$::datacols[$i]);
    close(FILE);
  }
}

sub extract_col {

  my $col = shift;

  my $n = $::ncols;
  my $tcol = $::datacols[0];
  my $nt =  @$tcol;

  my $file = make_file_name("extract",$col,$n);
  open(FILE,">$file");
  write_col($nt,$tcol,$::datacols[$col]);
  close(FILE);
}

sub write_col {

  my ($nt,$rt,$rv) = @_;

  for(my $it=0;$it<$nt;$it++) {
    my $t = $rt->[$it];
    my $v = $rv->[$it];
    my $line = "$t  $v\n";
    print FILE $line;
  }
}

sub make_file_name {

  my ($name,$number,$tot) = @_;

  my $np = 0;

  if( $tot >= 100 ) {
    $np = sprintf("%03d", $number);
  } elsif( $tot >= 10 ) {
    $np = sprintf("%02d", $number);
  } else {
    $np = sprintf("%d", $number);
  }

  my $file = "$name.$np.txt";

  #print STDERR "file name: $file\n";

  return $file;
}

#------------------------------------------------------------------------

sub store_cols {

  my $rdata = shift;

  my $n = @$rdata;
  $::ncols = $n unless $::ncols;

  return $n if $n == 0;

  unless( scalar @::datacols ) {
    @::datacols = () unless @::datacols;

    for(my $i=0;$i<$n;$i++) {
      my $rcol = [];
      $::datacols[$i] = $rcol;
    }
  }

  for(my $i=0;$i<$n;$i++) {
    my $rcol = $::datacols[$i];
    push(@$rcol,$rdata->[$i]);
  }
}

sub write_cols {

  my $n = @::datacols;
  my $tcol = $::datacols[0];
  my $nt =  @$tcol;

  my $n1 = $n - 1;
  print STDERR "writing file with $n1 data columns and $nt time records\n";

  for(my $it=0;$it<$nt;$it++) {
    my $line = "";
    for(my $i=0;$i<$n;$i++) {
      my $rcol = $::datacols[$i];
      my $val = $rcol->[$it];
      $line .= "$val  "
    }
    print "$line\n";
  }
}

sub open_files {

  my $rdata = shift;

  my $n = @$rdata;
  my $file;
  my $i;

  return $n if $n == 0;

  for(my $i=1;$i<$n;$i++) {
    my $file = "$::newfile.$i";
    open($i,">$file");
  }

  $::openfile = 1;

  return $n;
}
  
sub close_files {

  for(my $i=1;$i<$::ncols;$i++) {
    close($i);
  }
}

#------------------------------------------------------------------------

