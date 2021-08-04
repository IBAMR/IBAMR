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

print "@filesToProcess\n";

for my $file (@filesToProcess) {
    print "Replacing Pointer<> for source file $file\n";
    my $directory = dirname $file;
    my $filebasename = basename $file;
    my $tempFile = $filebasename . ".tmp";
    open FILE, "< $file" || die "Cannot open file $file";
    open TEMPFILE, "> $tempFile" || die "Cannot open temporary work file $tempFile";
    while ( my $str = <FILE> ) {
        # Delete SAMRAI::Pointer header include
        $str =~ s/#include\ "SAMRAI\/tbox\/Pointer\.h"//g;
        $str =~ s/#include\ [<"]tbox\/Pointer\.h[>"]//g;

        # Replace Pointer<STUFF> with shared_ptr<STUFF>
        $str =~ s/(SAMRAI::tbox::)?Pointer\</std::shared_ptr\</g;
        # Replace constructor with new STUFF to make_shared
        $str =~ s/([\=\s])new\ ([^\(]*)(\(.*)/$1std::make_shared<$2>$3/g;
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

