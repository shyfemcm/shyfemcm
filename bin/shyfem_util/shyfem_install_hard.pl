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
# installs shyfem through hard references
#
# called by shell script
#
#---------------------------------------------

$::verbose = 0;

$reset = 0 unless $reset;
$install = 0 unless $install;

$dir = shift;

if( $reset ) {
  ;
} elsif( $install ) {
  ;
} else {
  die "one of -reset or -install must be set... aborting\n";
}

$change = 0;

if( $reset ) {
  reset_fem();
} else {
  install_fem();
}

exit $change;

#-------------------------------------------------------------

sub install_fem {

 while(<>) {

  if( /SHYFEMDIR/ ) {
    if( /^FEMDIR=/ ) {
      $change = 1;
      print "SHYFEMDIR=$dir        # ggu_hard_install\n";
      print STDERR "$ARGV ($dir) (shell): $_" if $::verbose;
    } elsif( /^use lib/ and /lib\/perl/ ) {
      $change = 1;
      print "#ggu_hard_install $_";		# write old first
      print "use lib (\"$dir/lib/perl\");     # ggu_hard_install\n";
      print STDERR "$ARGV ($dir) (lib/perl): $_" if $::verbose;
      next;
    } elsif( /^use lib/ and /ENV/ and /bin/ ) {
      $change = 1;
      print "#ggu_hard_install $_";		# write old first
      print "use lib (\"$dir/bin\");     # ggu_hard_install\n";
      print STDERR "$ARGV ($dir) (perl): $_" if $::verbose;
      next;
    } else {
      print STDERR "$ARGV ($dir) (none): $_" if $::verbose;
    }
  }

  if( /^FEMDIR_INSTALL/ ) {
      $change = 1;
      print "SHYFEM_INSTALL=$dir        # ggu_hard_install\n";
      print STDERR "$ARGV ($dir) (shell_install): $_" if $::verbose;
  }

  print;
 }
}

sub reset_fem {

 while(<>) {

  if( /ggu_hard_install$/ ) {
    $change = 1;
    next;
  } elsif( s/^#ggu_hard_install // ) {
    $change = 1;
  }

  print;
 }
}

#-------------------------------------------------------------

