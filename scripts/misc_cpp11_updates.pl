#! /usr/bin/perl
# This script does a few miscellaneous C++11 updates such as:
#   -- replacing UniquePtr and AutoPtr with std::unique_ptr
#   -- replacing libMesh UniquePtr header with <memory> system header
#   -- replacing boost objects with corresponding std objects.
# To use this file, use the command
# perl misc_cpp11_updates.pl <directory>
# where <directory> is the base path to be changed.
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
    print "Replacing UniquePtr<> and AutoPtr<> for source file $file\n";
    my $directory = dirname $file;
    my $filebasename = basename $file;
    my $tempFile = $filebasename . ".tmp";
    open FILE, "< $file" || die "Cannot open file $file";
    open TEMPFILE, "> $tempFile" || die "Cannot open temporary work file $tempFile";
    while ( my $str = <FILE> ) {
        # Delete libMesh auto_ptr header include and replace with memory
        $str =~ s/\#include\ \"libmesh\/auto_ptr\.h\"/\#include\ <memory>/g;

        # Replace Pointer<STUFF> with unique_ptr<STUFF>
        $str =~ s/(libMesh::)?(UniquePtr|AutoPtr)\</std::unique_ptr\</g;

        # Replace boost::array<> with std::array<>
        $str =~ s/boost::array</std::array</g;

        # Replace boost::unordered_map<> with std::unordered_map<>
        $str =~ s/boost::unordered_map</std::unordered_map</g;

        # Replace boost::tuple<> with std::tuple<>
        $str =~ s/boost::tuple</std::tuple</g;

        # Replace boost::make_tuple
        $str =~ s/boost::make_tuple</std::make_tuple</g;

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

