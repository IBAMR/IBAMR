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

my $newFile = $replaceDir . "/BoxList_Tree_calls.txt";
open NEWFILE, "> $newFile" || die "Cannot open output file $newFile\n";
for my $file (@filesToProcess) {
    print "Finding BoxList|Tree calls for source file $file\n";
    my $directory = dirname $file;
    my $filebasename = basename $file;
    my $tempFile = $filebasename . ".tmp";
    open FILE, "< $file" || die "Cannot open file $file";
    open TEMPFILE, "> $tempFile" || die "Cannot open temporary work file $tempFile";
    my $lineNumber = 0;
    while ( my $str = <FILE> ) {
        $lineNumber = $lineNumber + 1;
        # If header file, replace with BoxContainer"
        $str =~ s/Box(Tree|List)\.h/BoxContainer\.h/g;
        # If we find BoxList or BoxTree call/constructor, print file name, line number, line to newFile.
        if ( $str =~ /BoxList/ ) {
            print "Found possible BoxList call in $file\n";
            print NEWFILE "BoxList call needs replacement:\n";
            print NEWFILE "file: $file  on line: $lineNumber\n";
            print NEWFILE "line: $str\n";
        }
        elsif ( $str =~ /BoxTree/ ) {
            print "Found possible BoxTree call in $file\n";
            print NEWFILE "BoxTree call needs replacement:\n";
            print NEWFILE "file: $file on line: $lineNumber\n";
            print NEWFILE "line: $str\n";
        }
        print TEMPFILE $str;
    }
    close FILE || die "Cannot close file $file\n";
    close TEMPFILE || die "Canot close file $tempFile\n";
    if (compare($file,$tempFile) == 0) {
        unlink($tempFile);
    } else {
        unlink($file);
        rename($tempFile, $file);
    }
}
close NEWFILE || die "Cannot close file $newFile\n";
