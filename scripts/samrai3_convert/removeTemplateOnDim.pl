#! /usr/bin/perl
##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/tools/scripts/conversion2.0/renameXd.pl $
## Package:     SAMRAI scripts
## Copyright:   (c) 1997-2018 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1917 $
## Description: perl script to rename getNumber methods to be getNumberOf 
##

use strict;

use File::Basename;
use File::Find;
use File::Compare;
use Cwd;

my $replaceDir   = $ARGV[0];


# Flush I/O on write to avoid buffering
# $|=1;

my $debug=0;

my $end_of_line = $/;

#
# Remove duplicated values
#
sub unique {
    foreach my $test (@_){
	my $i = -1;
	my @indexes = map {$i++;$_ eq $test ? $i : ()} @_;
	shift @indexes;
	foreach my $index (@indexes){
	    splice(@_,$index,1);
	}
    }
    return @_;
}

my $pwd = cwd;

my @templatesOnDIMOnly = ();
my $filename="SAMRAI_classes_templated_on_DIM_only.txt";
open FILE, '<', $filename or die "Can't open $filename : $!";
while (<FILE>) {
    my $class=$_;
    chomp($class);
    push @templatesOnDIMOnly,$class;
}
close FILE;


my @templatesOnDIM = ();
my $filename="SAMRAI_classes_templated_on_DIM_and_other.txt";
open FILE, '<', $filename or die "Can't open $filename : $!";
while (<FILE>) {
    my $class=$_;
    chomp($class);
    push @templatesOnDIM,$class;
}
close FILE;

my @objectsNeedDIM = ();
my $filename="SAMRAI_classes_that_need_DIM.txt";
open FILE, '<', $filename or die "Can't open $filename : $!";
while (<FILE>) {
    my $class=$_;
    chomp($class);
    push @objectsNeedDIM,$class;
}
close FILE;

my @objectsNeedDIMOnly = ();
my $filename="SAMRAI_classes_that_need_DIM_only.txt";
open FILE, '<', $filename or die "Can't open $filename : $!";
while (<FILE>) {
    my $class=$_;
    chomp($class);
    push @objectsNeedDIMOnly,$class;
}
close FILE;

#=============================================================================
# Fixup for all source
#=============================================================================

#
# Excludes files that are in internal source code control directories.
#
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
print "$replaceDir\n";
find( \&selectAllSourceFiles, $replaceDir );

print "12 @filesToProcess\n" if ($debug);

for my $file (@filesToProcess) {
    print "Working on DIM fixes for source file $file\n";
    my $directory = dirname $file;

    my $filebasename = basename $file;

    my $tempFile = $filebasename . ".samrai.tmp";

    open FILE, "< $file" || die "Cannot open file $file";
    open TEMPFILE, "> $tempFile" || die "Cannot open temporary work file $tempFile";
    while ( my $str = <FILE> ) {

	# Replace variable definitions
	for my $classname (@templatesOnDIMOnly) {
#	    $str =~ s/$classname\<DIM\>/$classname/g;
	    $str =~ s/$classname\<NDIM\>/$classname/g;
#	    $str =~ s/$classname\<2\>/$classname/g;
#	    $str =~ s/$classname\<3\>/$classname/g;
	}

	for my $classname (@templatesOnDIM) {
#	    $str =~ s/$classname<DIM,\s*(.*)>/$classname<$1>/g;
#	    $str =~ s/$classname<DIM,(.*)>/$classname<$1>/g;
	    $str =~ s/$classname<NDIM,\s*(.*)>/$classname<$1>/g;
#	    $str =~ s/$classname<2,\s*(.*)>/$classname<$1>/g;
#	    $str =~ s/$classname<3,\s*(.*)>/$classname<$1>/g;
	}

    for my $classname (@objectsNeedDIM) {
        $str =~ s/(\s+)$classname(<.*>)?\((.*)\)/$1$classname$2\(Dimension\(NDIM\), $3\)/g;
    }


    for my $classname (@objectsNeedDIMOnly) {
        $str =~ s/(\s+)$classname(<.*>)?\(\)/$1$classname$2\(Dimension\(NDIM\)\)/g;
        $str =~ s/$classname(<.*>)?\ ([\w_]*);/$classname$1 $2(Dimension(NDIM));/g;
    }

    $str =~ s/template\ <int DIM>//g;
    $str =~ s/template\ <int DIM,\s*(.*)>/template <$1>/g;


    print TEMPFILE $str;
    }

    close FILE || die "Cannot close file $file";
    close TEMPFILE || die "Cannot close file $tempFile";

    # Only replace existing file if a replacement was done.
    if (compare($file,$tempFile) == 0) {
	unlink($tempFile);
    } else {
	unlink($file);
	rename( $tempFile, $file);
    }
}

