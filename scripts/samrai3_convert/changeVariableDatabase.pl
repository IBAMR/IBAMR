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


for my $file (@filesToProcess) {
    print "Adding zero/one IntVectors for default registers for source file $file\n";
    my $directory = dirname $file;
    my $filebasename = basename $file;
    my $tempFile = $filebasename . ".tmp";
    open FILE, "< $file" || die "Cannot open file $file";
    open TEMPFILE, "> $tempFile" || die "Cannot open temporary file $tempFile";
    while ( my $str = <FILE> ) {
        # Search for default ghost cell registers
        $str =~ s/registerVariableAndContext\(([.]*),([.]*)\)/registerVariableAndContext\($1,$2,IntVector::getZero\(NDIM\)\)/g;
        $str =~ s/registerVariableAndContext\((.*),(.*),.*0\);/registerVariableAndContext\($1,$2,IntVector::getZero\(NDIM\)\)/g;
        $str =~ s/registerVariableAndContext\((.*),(.*),.*1\);/registerVariableAndContext\($1,$2,IntVector::getOne\(NDIM\)\)/g;

        if($str =~ /registerVariableAndContext\((.*),(.*),.*0\);/) { print $str;}

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
