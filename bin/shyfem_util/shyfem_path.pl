#!/usr/bin/perl
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# cleans path from old bin entries
#
#-----------------------------------------------------

my $path = $ARGV[0];

my @path = split(/:/,$path);
my @new = ();

foreach my $dir (@path) {

  #next if $dir =~ /bin/;
  next if $dir =~ /shyfem/;

  push(@new,$dir);
}

my $path = join(":",@new);

print "$path\n";

