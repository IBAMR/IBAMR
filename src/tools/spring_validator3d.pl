#!/usr/bin/perl -w
#
# filename: check_spring_rest_length.pl
# author: Boyce Griffith
# usage: check_spring_rest_length.pl <input filename> <stiffness> <output filename>
#
# A simple Perl script to determine whether the rest lengths of the
# springs in an IBAMR spring file are consistent with the initial
# lengths in the corresponding IBAMR vertex file.

use List::Util qw[min max];

if ($#ARGV != 0) {
    die "incorrect number of command line arguments.\nusage:\n  check_spring_rest_length.pl <input filename>\n";
}

# parse the command line arguments
$input_filename = shift @ARGV;  chomp $input_filename;

print "input file: $input_filename\n";

# open the input files
open(VERTEX_IN, "$input_filename.vertex") || die "error: cannot open $input_filename.vertex for reading: $!";
open(SPRING_IN, "$input_filename.spring") || die "error: cannot open $input_filename.spring for reading: $!";

# the first line in the vertex input file has the format:
#
#   <number of nodes> (comments)
#
# and the remaining lines in the input file all have the format:
#
#   <x> <y> <z> (comments)
#
# where items in ()'s are optional
<VERTEX_IN>;
$i = 0;
while (<VERTEX_IN>) {
    chomp;
    @line = split;
    $x[$i] = $line[0];
    $y[$i] = $line[1];
    $z[$i] = $line[2];
    $i += 1;
}
print "read $i vertices\n";

# the first line in the spring input file has the format:
#
#   <number of edges> (comments)
#
# and the remaining lines in the input file all have the format:
#
#   <first node> <second node> <stiffness> <rest length> (force fcn index) (comments)
#
# where items in ()'s are optional
<SPRING_IN>;
$i = 0;
while (<SPRING_IN>) {
    chomp;
    @line = split;
    $idx0 = $line[0];
    $idx1 = $line[1];
    $rst  = $line[3];

    $dx[0] = $x[$idx0] - $x[$idx1];
    $dx[1] = $y[$idx0] - $y[$idx1];
    $dx[2] = $z[$idx0] - $z[$idx1];

    $rst_actual = sqrt($dx[0]*$dx[0]+$dx[1]*$dx[1]+$dx[2]*$dx[2]);
    if (abs($rst_actual - $rst)/max($rst_actual,$rst,1.0) > 1.0e-6) {
        print "check spring on line: $i\n";
        print "rest length in spring file = $rst\n";
        print "initial displacement in vertex file = $rst_actual\n";
    }
    $i += 1;
}
print "read $i springs\n";
