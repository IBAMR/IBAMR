#! /usr/bin/perl
##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/scripts/source_manipulation/replaceIncludeGuards.pl $
## Package:     SAMRAI scripts
## Copyright:   (c) 1997-2018 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1917 $
## Description: perl script to insert the package prefix in #includes
##

use File::Basename;
use File::Find;
use File::Path;
use Cwd;
use Text::Wrap;

my $replaceDir   = $ARGV[1];

my $SAMRAISourceDir = $ARGV[0];

my $debug = 1;

print "$SAMRAISourceDir" if $debug > 0;

print "Replacing SAMRAI 2.x includes in directory \"$replaceDir\"\n\n";

my $pwd = cwd;

#
# File pattern to look for
#
my $filePattern;

my $excludePattern;

@allfiles = ();
sub selectFile {
    if ( $File::Find::dir =~ m!$excludePattern! ) {
	$File::Find::prune = true;
    }
    elsif ( -f && m/$filePattern/ ) {
	push @allfiles, $File::Find::name;
	$allfiles[$#allfiles] =~ s|^\./||;
    }
}

@packages=qw/tbox hier xfer pdat math mesh geom solv algs appu/;
#@packages= qw/tbox pdat/;

%packageDirectories=('tbox', 'toolbox',
		     'hier', 'hierarchy',
		     'xfer', 'transfer',
		     'pdat', 'patchdata',
		     'math', 'mathops',
		     'mesh', 'mesh',
		     'geom', 'geometry',
		     'solv', 'solvers',
		     'algs', 'algorithm',
		     'appu', 'apputils');

undef %fileToPackage;
# Create list of all source files
foreach $package (@packages) {
    print "Processing classs in package : $package\n" if ($debug); 
    $dir=$SAMRAISourceDir . "/source/SAMRAI/$package/";

    @allfiles = ();
    $filePattern = q|(.*\.[ChI]$)|;
    $excludePattern=q!/(.svn|CVS|templates)$!;
    find( \&selectFile, $dir );
    print "files=@allfiles" if ($debug > 1);
    foreach $file (@allfiles) {
	print "\tAdding file $file\n" if ($debug); 
	($pattern = basename $file) =~ s/(.*)\.[ChI]$/$1/;
	push @allPatterns, basename $pattern;
	$base=basename $file;
	$fileToPackage{$base}=$package;
    }
}



my $allSimplepatterns = join '|', @allPatterns;

print "$allSimplepatterns\n" if ($debug);

@allfiles = ();

$filePattern = q/(.*\.[ChI]$)|(.*\.CPP$)|(.*\.cpp$)|(.*\.cxx$)|(.*\.CXX$)|(.*\.H$)|(.*\.hxx$)|(.*\.Hxx$)|(.*\.HXX$)|(.*\.txx$)|(.*\.ixx$)/;
$excludePattern=q!/(build|\.git|contrib|config|configure|CVS|scripts|\{arch\})$!;
find( \&selectFile, $replaceDir );

print "@allfiles\n" if ($debug);

for $file (@allfiles) {
    print "\tWorking on $file\n";
    $tfile=$file . ".tmp";

    open OLDFILE, "< $file" || die "Cannot open file $file";
    open NEWFILE, "> $tfile" || die "Cannot open file $tfile";
    while ( $str = <OLDFILE> ) {

	# Fixup the non tbox includes
        $str =~ s/\#include(\s*)[<"]($allSimplepatterns)\.([hIC])[">]/#include$1\"SAMRAI\/$fileToPackage{"$2.$3"}\/$2.$3\"/go;

	# Fixup the tbox includes, these were already in namespace
        $str =~ s/\#include(\s*)["<]tbox\/(.*)[>"]/#include$1\"SAMRAI\/tbox\/$2\"/go;

	# Fixup configure header include
        $str =~ s/\#include(\s*)["<]SAMRAI_config.h[>"]/#include$1\"SAMRAI\/SAMRAI_config.h\"/go;
	
	print NEWFILE $str;
    }

    close OLDFILE || die "Cannot close file $file";
    close NEWFILE || die "Cannot close file $tfile";

    printf "unlink $file\n" if ( $debug > 1 );
    unlink($file);
    printf "rename $tfile $file\n" if ( $debug > 1 );
    rename( $tfile, $file);
}
