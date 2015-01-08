#!/usr/bin/perl

# AUTHOR: Joseph Fass
# LAST REVISED: June 2008
# 
# The Bioinformatics Core at UC Davis Genome Center
# http://bioinformatics.ucdavis.edu
# Copyright (c) 2008 The Regents of University of California, Davis Campus.
# All rights reserved.

use strict;

my $usage = "\nusage: fasta1line.pl <input file (fasta format)> <output file>\n\n".
            "For each sequence, puts nt/aa all on one line following header line.\n\n";

my $infile = shift or die $usage;
my $outfile = shift or die $usage;

open IN, "<$infile" or die $usage;
open OUT, ">$outfile" or die $usage;

my $header; my $sequence;
my $firstline = 1;

while (<IN>) {
  if (m/^>/) { # recognize a header line?
    if (!$firstline) { # output previous sequence and clear it, unless this is the first line of the file
      print OUT $sequence."\n";
      $sequence = "";
    } # if
    $firstline = 0; # we're not at the first line of the file anymore
    print OUT; # print out header line
  }
  else { # not a header line? - must be sequence
    chomp; # remove newline at end
    $sequence .= $_; # append additional sequence
  } # if-else
} # while
print OUT $sequence."\n"; # output final sequence

close IN; close OUT;
