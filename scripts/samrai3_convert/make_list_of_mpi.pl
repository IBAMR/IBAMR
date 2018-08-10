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
    if ( $File::Find::name =~ m!/(build|\.git|contrib|config|configure|CVS|scripts|\{arch\})$!o ) {
        $File::Find::prune = 1;
    }
    elsif ( -f && m/$filePattern/ ) {
        push @filesToProcess, $File::Find::name;
        $filesToProcess[$#filesToProcess] =~ s|^\./||o;
   }
}

find( \&selectAllSourceFiles, $replaceDir );

my $newFile = $replaceDir . "/SAMRAI_MPI_calls.txt";
open NEWFILE, "> $newFile" || die "Cannot open output file $newFile\n";
for my $file (@filesToProcess) {
    print "Finding SAMRAI_MPI calls for source file $file\n";
    my $directory = dirname $file;
    my $filebasename = basename $file;
    open FILE, "< $file" || die "Cannot open file $file";
    my $lineNumber = 0;
    while ( my $str = <FILE> ) {
        $lineNumber = $lineNumber + 1;
        # If we find SAMRAI_MPI call, print file name, line number, line to newFile.
        if ( $str =~ /SAMRAI_MPI::/ ) {
            print "Found possible SAMRAI_MPI call in $file\n";
            print NEWFILE "SAMRAI_MPI call needs replacement:\n";
            print NEWFILE "file: $file  on line: $lineNumber\n";
            print NEWFILE "line: $str\n";
        }
    }
    close FILE || die "Cannot close file $file\n";
}
close NEWFILE || die "Cannot close file $newFile\n";
