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


for my $file (@filesToProcess) {
    print "Converting SAMRAI_MPI calls to IBTK_MPI calls for source file $file\n";
    my $directory = dirname $file;
    my $filebasename = basename $file;
    my $tempFile = $filebasename . ".tmp";
    open FILE, "< $file" || die "Cannot open file $file";
    open TEMPFILE, "> $tempFile" || die "Cannot open temporary file $tempFile";
    my $found = 0;
    while ( my $str = <FILE> ) {

        # Search for SAMRAI_MPI calls
        if ($str =~ s/SAMRAI_MPI::/IBTK_MPI::/g) {
            $found = 1;
        }
        if ($str =~ /IBTK_MPI::setCommunicator/ || $str =~ /IBTK_MPI::setCallAbortInSerialInsteadOfExit/) {
            $str =~ s/IBTK_MPI::/SAMRAI_MPI::/g;
        }
        if ($str =~ s/IBTK_MPI::commWorld/IBTK_MPI::getSAMRAIWorld()/g) {
            $found = 1;
        }
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

    if ($found) {
        my $printed = 0;
        open FILE, "< $file" || die "Cannot open file $file";
        open TEMPFILE, "> $tempFile" || die "Cannot open temporary file $tempFile";
        while ( my $str = <FILE> ) {
            if ($printed == 0 && ($str =~ /#include\ [<"]ibtk/ || $str =~ /#include\ [<"]ibamr/)) {
                print TEMPFILE "#include \"ibtk/IBTK_MPI.h\"\n";
                $printed = 1;
            }
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

}
