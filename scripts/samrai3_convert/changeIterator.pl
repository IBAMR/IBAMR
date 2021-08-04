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
    print "Updating Iterators for source file $file\n";
    my $directory = dirname $file;
    my $filebasename = basename $file;
    my $tempFile = $filebasename . ".tmp";
    open FILE, "< $file" || die "Cannot open file $file";
    open TEMPFILE, "> $tempFile" || die "Cannot open temporary work file $tempFile";
    while ( my $str = <FILE> ) {
        # Don't mess with includes...
        if ($str =~ /\#include/)
        {
            print TEMPFILE $str;
            next;
        }
        # Check what kind of iterator we are running into...
        # If it looks like PatchLevel::Iterator, we need to change to level->begin()
        # also change level->getLevel(p()) to p().
        if ($str =~ /PatchLevel::Iterator/)
        {
            $str =~ s/\(level\)/\ =\ level->begin\(\)/g;
            $str =~ s/([a-zA-Z_][\w]*);/$1\ !=\ level->end\(\);/g;
            
            #$str =~ /([\s]*)for\ \(PatchLevel::Iterator (.*)\((.*)\)/;
            #my $space = $1;
            #my $type = $3;
            #my $var = $2;
            #$str = join "", $space, "for (auto& ", $var, " : ", $type, ")\n";
        }
        if ($str =~ /getPatch\([a-zA-Z_][\w]*\(\)\)/)
        {
            $str =~ s/([a-zA-Z_][\w]*)->getPatch\(([a-zA-Z_][\w]*)\(\)\)/\*$2/g;
        }
        # If it looks like a Cell iterator, we need to add a line CellIterator end(box, false)
        # also change constructor and ending
        if ($str =~ /(Cell|Node)Iterator/)
        {
            my $type = $1;
            $str =~ /([ ]*)for/;
            my $spaces = $1;
            $str =~ s/([a-zA-Z_][\w]*)\(([a-zA-Z_][\w]*)\);/$1\(pdat::${type}Geometry::begin\($2\)\);/g;
            my $str_new = $spaces."${type}Iterator end(pdat::${type}Geometry::end($2));\n";
            print TEMPFILE $str_new;
            $str =~ s/([a-zA-Z_][\w]*);/$1\ \!=\ end;/g;
        }
        if ($str =~ /(Face|Cell|Side|Node)Index[\ &][\w]/)
        {
            print $str;
            my $idx_type = $1."Index";
            my $idx_1 = $idx_type;
            print $idx_type;
            print $idx_1;
            $str =~ s/(const\ )?$idx_1&?\ ?([\w][\w]*) = ([\w][\w]*)\(\)/const\ $idx_type\ $2\ =\ \*$3/g;
            print $str;
        }
        
        # If it looks like a Face iterator, we need to add a line FaceIterator end(box, axis, false)
        # also change constructor and ending
        if ($str =~ /(Side|Face)Iterator/)
        {
            my $type = $1;
            $str =~ /([ ]*)for/;
            my $spaces = $1;
            $str =~ s/([a-zA-Z_][\w]*)\(([a-zA-Z_][\w]*),([ ]?[a-zA-Z_][\w]*)\);/$1\(pdat::${type}Geometry::begin\($2,$3\)\);/g;
            my $str_new = $spaces."${type}Iterator end(pdat::${type}Geometry::end($2,$3));\n";
            print TEMPFILE $str_new;
            $str =~ s/([a-zA-Z_][\w]*);/$1\ \!=\ end;/g;
        }
        if ($str =~ /BoxIterator/)
        {
            $str =~ /([ ]*)for/;
            my $spaces = $1;
            $str =~ s/([a-zA-Z_][\w]*)\(([a-zA-Z_][\w]*)\);/$1\($2\.begin\(\)\);/g;
            my $box_name = $1;
            my $str_new = $spaces."BoxIterator end\($box_name\.end\(\)\);\n";
            print TEMPFILE $str_new;
            $str =~ s/([a-zA-Z_][\w]*);/$1\ \!=\ end;/g;
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

