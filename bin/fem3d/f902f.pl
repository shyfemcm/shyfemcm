#!/usr/bin/perl -s
#
# changes .f90 files back to .f format (only continuation lines)
#
#------------------------------------------

$::old = 0 unless $::old;
$old = $::old;

while(<>) {

  chomp;

  s/\s*\&\s*$//;

  if( $old ) {
    s/^\s*\&/     +/;
    next if check_bmpi_skip();
  }

  print "$_\n";
}

#------------------------------------------

sub check_bmpi_skip {

  if( /if\( bmpi_skip \) return/ ) {
    <>;
    return 1;
  }

  if( /if\( bmpi_skip \) then/ ) {
    while(<>) {
      last if /end\s*if/;
    }
    <>;
    return 1;
  }

  return 0;
}

