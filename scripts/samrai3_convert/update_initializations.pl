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
    my $arg1 = "";
    my $arg2 = "";
    my $arg3 = "";
    my $arg4 = "";
    my $arg5 = "";
    my $libmesh_init = "\n";
    my $libmesh_found = 0;
    my $petsc_found = 0;
    my $samrai_found = 0;
    my $done = 0;
    while ( my $str = <FILE> ) {

        # Search for initialization calls.
        if ($str =~ /PetscInitialize\(&(.*),[\s]?&(.*),[\s]?(.*),[\s]?(.*)\);/)
        {
            $arg1 = $1;
            $arg2 = $2;
            $arg4 = $3;
            $arg5 = $4;
            $str = "\n";
            $petsc_found = 1;
            next;
        }
        if ($str =~ /SAMRAI_MPI::setCommunicator\((.*)\)/)
        {
            $arg3 = $1;
            $str = "\n";
            $samrai_found = 1;
            next;
        }
        if ($str =~ /LibMeshInit\ (.*)\((.*),\ (.*)\);/)
        {
            $libmesh_init = $1;
            $arg1 = $2;
            $arg2 = $3;
            $arg4 = "nullptr";
            $arg5 = "nullptr";
            $str = "\n";
            $libmesh_found = 1;
            next;
        }
        # Removed closing calls.
        if ($str =~ /SAMRAIManager::startup\(\)/ || $str =~ /SAMRAIManager::shutdown\(\)/ || $str =~ /PetscFinalize\(\)/)
        {
            $str = "\n";
            next;
        }
        if ($libmesh_found)
        {
            $str =~ s/$libmesh_init\./init.getLibMeshInit()->/g;
        }
        print TEMPFILE $str;
        if ($done == 0 && ($samrai_found && ($petsc_found || $libmesh_found)))
        {
            print TEMPFILE "IBTK_Init init($arg1, $arg2, $arg3, $arg4, $arg5);\n";
            $done = 1;
        }
    }
    close FILE || die "Cannot close file $file\n";
    close TEMPFILE || die "Cannot close file $tempFile\n";
    if (compare($file,$tempFile) == 0) {
        unlink($tempFile);
    } else {
        unlink($file);
        rename($tempFile, $file);
    }

    if ($done) {
        my $printed = 0;
        open FILE, "< $file" || die "Cannot open file $file";
        open TEMPFILE, "> $tempFile" || die "Cannot open temporary file $tempFile";
        while ( my $str = <FILE> ) {
            if ($printed == 0 && ($str =~ /#include\ [<"]IBAMR_config/))
            {
                my $str_new = "#include \"ibtk/IBTK_Init.h\"\n";
                $printed = 1;
                print TEMPFILE $str_new;
            }
            print TEMPFILE $str;
        }
        close FILE || die "Cannot close file $file\n";
        close TEMPFILE || die "Cannot close file $tempFile\n";
        if (compare($file,$tempFile) == 0) {
            unlink($tempFile);
        }
        else
        {
            unlink($file);
            rename($tempFile, $file);
        }
    }
}
