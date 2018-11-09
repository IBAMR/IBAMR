#! /usr/bin/perl

use strict;

use File::Basename;
use File::Find;
use File::Compare;
use Cwd;

my $replaceDir = $ARGV[0];

my $debug = 0;
my $end_of_line = $/;

my $filePattern = q/(.*\.f\.m4)/;
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

print "@filesToProcess\n";

for my $file (@filesToProcess) {
    print "Changing m4 includes for source file $file\n";
    my $directory = dirname $file;
    my $filebasename = basename $file;
    my $tempFile = $filebasename . ".tmp";
    open FILE, "< $file" || die "Cannot open file $file";
    open TEMPFILE, "> $tempFile" || die "Cannot open temporary work file $tempFile";
    while ( my $str = <FILE> ) {
        $str =~ s/SAMRAI_FORTDIR\/pdat_m4arrdim([23])d.i/SAMRAI_FORTDIR\/SAMRAI\/pdat\/fortran\/pdat_m4arrdim${1}d.i/g;
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

