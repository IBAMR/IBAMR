#! /usr/bin/perl

use strict;

use File::Basename;
use File::Find;
use File::Compare;
use Cwd;

my $replaceDir = $ARGV[0];

my $debug = 0;

my $end_of_line = $/;

my $filePattern = q/(.*\.[ChI]$)|(.*\.CPP$)|(.*\.cpp$)|(.*\.cxx$)|(.*\.CXX$)|(.*\.H$)|(.*\.hxx$)|(.*\.Hxx$)|(.*\.HXX$)|(.*\.txx$)|(.*\.ixx$)/;

my @filesToProcess = ();
sub selectAllSourceFiles {
    if ( $File::Find::name =~ m!/(build|\.git|contrib|config|configure|CVS|scripts|\{arch\}|CMakeLists)$!o ) {
        $File::Find::prune = 1;
    }
    elsif ( -f && m/$filePattern/ ) {
        push @filesToProcess, $File::Find::name;
        $filesToProcess[$#filesToProcess] =~ s|^\./||o;
   }
}

find( \&selectAllSourceFiles, $replaceDir );

my $newFile = $replaceDir . "/GriddingAlgorithm_calls.txt";
open NEWFILE, "> $newFile" || die "Cannot open output file $newFile\n";
for my $file (@filesToProcess) {
    print "Finding GriddingAlgorithm calls for source file $file\n";
    my $directory = dirname $file;
    my $filebasename = basename $file;
    open FILE, "< $file" || die "Cannot open file $file";
    my $lineNumber = 0;
    while ( my $str = <FILE> ) {
        $lineNumber = $lineNumber + 1;
        # If we find SAMRAI_MPI call, print file name, line number, line to newFile.
        if ( $str =~ /(makeCoarsestLevel\([\s\w\+\-\\.]+\))|(makeFinerLevel\([\s\w\+\-\\.]+\))|(regridAllFinerLevels\([\s\w\+\-\\.]+\))|(levelCanBeRefined\([\s\w\+\-\\.]+\))|(getMaxLevels\([\s\w\+\-\\.]*\))|(getRatioToCoarserLevel\([\s\w\+\-\\.]+\))|(getSmallestPatchSize\([\s\w\+\-\\.]+\))|(getLargestPatchSize\([\s\w\+\-\\.]+\))/ ) {
            print "Found possible GriddingAlgorithm call in $file\n";
            print NEWFILE "GriddingAlgorithm call needs replacement:\n";
            print NEWFILE "file: $file  on line: $lineNumber\n";
            print NEWFILE "line: $str\n";
        }
    }
    close FILE || die "Cannot close file $file\n";
}
close NEWFILE || die "Cannot close file $newFile\n";
