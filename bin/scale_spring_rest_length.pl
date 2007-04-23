#!/usr/bin/perl -w
#
# filename: scale_spring_rest_length.pl
# author: Boyce Griffith
# usage: scale_spring_rest_length.pl <input filename> <scale factor> <output filename>
#
# A simple Perl script to scale the spring rest lengths in an IBAMR
# input file by a uniform, user-specified amount.

if ($#ARGV != 2) {
    die "incorrect number of command line arguments.\nusage:\n  scale_spring_rest_length.pl <input filename> <scale factor> <output filename>\n";
}

# parse the command line arguments
$input_filename = shift @ARGV;  chomp $input_filename;
$scale = shift @ARGV;  chomp $scale;
$output_filename = shift @ARGV;  chomp $output_filename;

print "input file: $input_filename\n";
print "scale factor: $scale\n";
print "output file: $output_filename\n";

# open the input and output files
if ($input_filename eq $output_filename) {
    die "error: input and output files must be different\n";
}
if (-e $output_filename) {
    print "warning: about to overwrite contents of $output_filename\n";
    print "press [Enter] to continue...";  <STDIN>;
}

open(IN, "$input_filename") || die "error: cannot open $input_filename for reading: $!";
open(OUT, ">$output_filename") || die "error: cannot open $output_filename for writing: $!";

# the first line in the input file has the format:
#
#   <number of edges> (comments)
#
# where items in ()'s are optional
$_ = <IN>;
chomp;
@line = split;
for ($i = 0; $i <= $#line; $i++) {
    $a = $line[$i];

    if (($i == 0) && ($a =~ /^[+-]?\d+$/)) {
	printf OUT "%6d", $a; # integer
    }
    else
    {
	print OUT $a;
    }

    if ($i < $#line) {
	print OUT " ";
    }
    else {
	print OUT "\n";
    }
}

# the remaining lines in the input file all have the format:
#
#   <first node> <second node> <stiffness> <rest length> (force fcn index) (comments)
#
# where items in ()'s are optional
while (<IN>) {
    chomp;
    @line = split;
    for ($i = 0; $i <= $#line; $i++) {
	$a = $line[$i];

	if ($i == 3) {
	    $a = $scale*$a; # rescale the rest length
	}

	if (($i == 0 || $i == 1 || $i == 4) && ($a =~ /^[+-]?\d+$/)) {
	    printf OUT "%6d", $a; # integer
	}
	elsif ($i == 2 || $i == 3) {
	    printf OUT "%1.16e", $a; # floating point value
	}
	else {
	    print OUT $a;
	}

	if ($i < $#line) {
	    print OUT " ";
	}
	else {
	    print OUT "\n";
	}
    }
}

# close the input and output files
close(IN) || die "error: cannot close $input_filename: $!";
close(OUT) || die "error: cannot close $output_filename: $!";
