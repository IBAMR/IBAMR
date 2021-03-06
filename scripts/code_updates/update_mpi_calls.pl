#! /usr/bin/perl
## ---------------------------------------------------------------------
##
## Copyright (c) 2019 - 2020 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

# This script replaces SAMRAI_MPI calls with the corresponding IBTK_MPI call
# To use this file, use the command
# perl update_mpi_calls.pl <directory>
# where <directory> is the base path to be changed.
use strict;

use File::Basename;
use File::Find;
use File::Compare;
use Cwd;

my $replaceDir = $ARGV[0];
my @fileExcludeList = ("IBTKInit.h","IBTKInit.cpp","IBTK_MPI.cpp","IBTK_MPI.h");

my $debug = 0;

my $end_of_line = $/;

my $filePattern = q/(.*\.[ChI]$)|(.*\.CPP$)|(.*\.cpp$)|(.*\.cxx$)|(.*\.CXX$)|(.*\.H$)|(.*\.hxx$)|(.*\.Hxx$)|(.*\.HXX$)|(.*\.txx$)|(.*\.ixx$)/;
my @filesToProcess = ();
sub selectAllSourceFiles {
    if ( $File::Find::name =~ m!/(build|\.git|CVS|contrib|scripts|config|configure|\{arch\})$!o ) {
        $File::Find::prune = 1;
    }
    elsif ( -f && m/$filePattern/ ) {
        push @filesToProcess, $File::Find::name;
        $filesToProcess[$#filesToProcess] =~ s|^\./||o;
   }
}

find( \&selectAllSourceFiles, $replaceDir );

print "@filesToProcess\n";

for my $file (@filesToProcess) {
    print "Replacing SAMRAI_MPI for source file $file\n";
    my $directory = dirname $file;
    my $filebasename = basename $file;
    my $tempFile = $filebasename . ".tmp";
    if ( grep( /$filebasename/, @fileExcludeList) )
    {
        print "Skipping file $file\n";
	next;
    }
    open FILE, "< $file" || die "Cannot open file $file";
    open TEMPFILE, "> $tempFile" || die "Cannot open temporary work file $tempFile";
    while ( my $str = <FILE> ) {
        # Replace SAMRAI_MPI header file
        $str =~ s/\#include\ [\"\<]tbox\/SAMRAI_MPI\.h[\"\>]/\#include\ "ibtk\/IBTK_MPI.h"/g;
        
        # Replace SAMRAI::tbox::SAMRAI_MPI::commWorld
        $str =~ s/(SAMRAI::)?(tbox::)?SAMRAI_MPI::commWorld/IBTK_MPI::getCommunicator()/g;

        # Replace SAMRAI_MPI with IBTK_MPI
        $str =~ s/(SAMRAI::)?(tbox::)?SAMRAI_MPI::/IBTK_MPI::/g;

        print TEMPFILE $str;
    }

    close FILE || die "Cannot close file $file\n";
    close TEMPFILE || die "Cannot close file $tempFile\n";
    if (compare($file,$tempFile) == 0) {
        unlink($tempFile);
    } else {
        unlink($file);
        rename($tempFile, $file);
    }

}

