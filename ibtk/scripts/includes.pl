#!/usr/bin/perl

#
# Copy include information for all include files into the main include directory
#

# Make dirname the path to the scripts directory.
$_ = $0; s|[^/]+$||; s/^$/./; s|/$||; $dirname = $_;

#
# Check for "--link" flag to make links instead of copies.
#

if ( @ARGV && $ARGV[0] eq '--link' ) {
    shift;
    $link = 1;
}

$includedir="$dirname/../include/ibtk";

#
# Get relative path to and from source directory.
#

$pwd=`pwd`;
chop $pwd;

if ( $pwd =~ m'/src' ) {
    # We are in the source directory.
    ( $fromsource = $pwd ) =~ s|(.*/src)/?||g;

    $sourcepath = $1;
    ( $tosource = $fromsource ) =~ s/[^\/]+/../g;
    ( $filterdir = $fromsource ) =~ s/^$/./;
}
elsif ( -d 'src' ) {
    # We are directly above the source directory.
    $tosource = 'src';
    $fromsource = '..';
    $sourcepath = "$pwd/src";
    $filterdir = '.';
}
else {
    # We should not be running this script from here.
    die "This script must be run from somewhere in the IBTK working directory."
}
print "current dir: $pwd\n";
print "source path: $sourcepath\n";
print "from source: $fromsource\n";
print "to source: $tosource\n";
print "filter dir: $filterdir\n";

#
# Get all the files we will be looking for.
# Note that paths found are relative to the top level directory.
#

print "Scanning...\n";
open COM, "cd $sourcepath && find $filterdir -name '*.[hIC]' -print|";
@allfiles = <COM>;
print "Done scanning.\n";
for (@allfiles) { chop; s|^\./||; }
close COM;

print "Selecting files...\n";

#
# Gather common-case header files and inline files that should
# be in include.
#

@headers = grep /\.h$/, @allfiles;
@headers = grep !/(tests\/|examples\/)/, @headers;

@includes = grep /\.I$/, @allfiles;
@includes = grep !/(tests\/|examples\/)/, @includes;

#
#  Gather a few fortran include files that should be in included.
#

@fincludes = grep /\.i$/, @allfiles;
@fincludes = grep !/(tests\/|examples\/)/, @fincludes;

# All of the *.C files of type "template <class TYPE>"
# should be in include also!

@templates = grep /\.C$/, @allfiles;
@templates = grep {
    open IFILE, "< $sourcepath/$_" || die "Cannot open $_";
    $lines = join '', <IFILE>;
    close IFILE;
    $lines =~ m/(template\s?<\s?int DIM\s?>)|(template\s?<\s?int DIM\s?,\s?class\sTYPE\s?>)|(template\s<class\sTYPE>)|(template\s<>)/;
} @templates;
@templates = grep !/(tests\/|examples\/)/, @templates;

print "Done selecting files...\n";

#
# If a file should be in include and is not, link or copy it there.
# Warn if an existing link is to the wrong file.
#

print "Copying/linking files...\n";
for $path ( @headers, @includes, @fincludes, @templates ) {
    print "Checking $path\n";
    ( $file = $path ) =~ s:.*/::;  # File is the base name.

    if ( -l "$includedir/$file" ) {
	# The version in include is a link.  See if it is a correct link.
	$ltarget = readlink "$includedir/$file";
	print( "WARNING: File $includedir/$file appears to be a link to the wrong file.\n" )
	    if "$ltarget" ne "../../src/$path";
	unlink  "$includedir/$file"
	    || die "Cannot remove $includedir/$file before creating";
    }

    if ( ! -e "$includedir/$file" ) {
	# The version in include does not exist.  Create it.
	print "$includedir/$file -> ../$path\n";
	if ( $link ) {
	    symlink( "../../src/$path", "$includedir/$file" )
		|| die "Cannot create link $includedir/$file -> $path";
	}
	else {
	    &cp( "$tosource/$path", "$includedir/$file" )
		|| die "Cannot copy $path to $includedir/$file";
	}
    }
    else {
	# The version in include is a copy (exists and is not a link).
	# Update it if necessary.
	if ( &cmpfiles("./$tosource/$path", "$includedir/$file" ) ) {
	    unlink "$includedir/$file"
		|| die "Cannot remove $includedir/$file before creating";
	    if ( $link ) {
		print STDERR "WARNING: File $includedir/$file does not match $path\n";
	    }
	    else {
		&cp( "$tosource/$path", "$includedir/$file" )
		    || die "Cannot copy $path to $includedir/$file";
	    }
	}
    }
}
print "Done copying/linking files\n";

print "Creating convenience header\n";
open FILE, ">$includedir/ibtk.h" ||
    die "Cannot open file $includedir/ibtk.h\n";
for $path ( @headers, @includes ) {
    print "Checking $path\n";
    ( $file = $path ) =~ s:.*/::;  # File is the base name.
    print FILE "#include \"ibtk/$file\"\n";
}
close FILE ||
    die "Cannot close file $includedir/ibtk.h\n";

#
# Subroutine to check if two files are the same.
#
sub cmpfiles {
($ANAME,$BNAME) = @_;

open(AFILE, "$ANAME") || die "Cannot open input file $ANAME...";
open(BFILE, "$BNAME") || die "Cannot open input file $BNAME...";

while (!eof(AFILE) && !eof(BFILE)) {
    $ALINE = <AFILE>;
    $BLINE = <BFILE>;
    $_ = $ALINE;

    if (!/^(\/\/|c|C|#|##| \*)[ ]*(Release:[\t ]*\$Name|Revision:[\t ]*\$Revision|Modified:[\t ]*\$Date):[^\$]*\$/o) {
	    if ($ALINE ne $BLINE) {
		close AFILE;
		close BFILE;
		return 1;
	    }
	}
    }

    if (eof(AFILE) && eof(BFILE)) {
	$rvalue = 0;
    } else {
	$rvalue = 1;
    }

    close AFILE;
    close BFILE;
    return $rvalue;
}

#
# Subroutine to copy a file
#

sub cp {
    my ($fr,$to) = @_;
    open ( FR, "<$fr" ) || die "Cannot open $fr for reading";
    open ( TO, ">$to" ) || die "Cannot open $fr for writing";
    while ( $_ = <FR> ) { print TO $_; }
    close FR || die "Cannot close $fr";
    close TO || die "Cannot close $to";
    my($atime,$mtime) = (stat $fr)[8..9];
    utime $atime, $mtime, $to;
    return 1;
}
