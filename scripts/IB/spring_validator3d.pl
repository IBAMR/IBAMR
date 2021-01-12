#!/usr/bin/perl -w
## ---------------------------------------------------------------------
##
## Copyright (c) 2008 - 2020 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

#
# filename: spring_validator3d.pl
# author: Boyce Griffith
# usage: spring_validator3d.pl <input filename>
#
# A simple Perl script to determine whether the spring specifications
# in an IBAMR spring file are valid, and whether the specified spring
# resting lengths are consistent with the initial displacements in the
# corresponding IBAMR vertex file.
#
# NOTE: This utility will generate a large number of warning messages
# for sturctures where the initial lengths and resting lengths differ.

use List::Util qw[min max];

if ($#ARGV != 0) {
    die "incorrect number of command line arguments.\nusage:\n  spring_validator3d.pl <input filename>\n";
}

# boolean controling whether to check to see if the initial
# displacement equals the resting length.
$check_for_consistent_rest_lengths = 1;  # 1 = enabled; 0 = disabled
$mach_eps = 2.22044604925031e-16;
$rel_tol = sqrt($mach_eps);
$abs_tol = 20*$mach_eps;

print "\n";
if ($check_for_consistent_rest_lengths) {
    print "checks for consistent spring rest lengths are ENABLED!\n";
    print "warning messages WILL be printed if the initial displacements and resting lengths differ by more than a relative error tolerance of $rel_tol or an absolute error tolerance of $abs_tol.\n";
} else {
    print "checks for consistent spring rest lengths are DISABLED!\n";
    print "warning messages WILL NOT be printed if the initial displacements and resting lengths differ.\n";
}
print "\n";

# parse the command line arguments
$input_filename = shift @ARGV;  chomp $input_filename;

# the first line in the vertex input file has the format:
#
#   <number of nodes> (comments)
#
# and the remaining lines in the input file all have the format:
#
#   <x> <y> <z> (comments)
#
# where items in ()'s are optional
open(VERTEX_IN, "$input_filename.vertex") || die "error: cannot open $input_filename.vertex for reading: $!";
$_ = <VERTEX_IN>;
chomp;
@line = split;
$num_vertices_expected = $line[0];
$i = 0;
$num_vertices = 0;
$line_number = 1;
while (<VERTEX_IN>) {
    chomp;
    @line = split('#');
    $_ = $line[0];
    @line = split;

    $x[$i] = $line[0];
    $y[$i] = $line[1];
    $z[$i] = $line[2];

    $i += 1;
    $num_vertices += 1;
    $line_number += 1;
}
print "read $num_vertices vertices from vertex file: $input_filename.vertex\n";
if ($num_vertices != $num_vertices_expected) {
    print "WARNING: expected to read: $num_vertices_expected from vertex file $input_filename.vertex\n";
    print "         actually read:    $num_vertices\n";
}
print "\n";
close(VERTEX_IN);

# the first line in the spring input file has the format:
#
#   <number of edges> (comments)
#
# and the remaining lines in the input file all have the format:
#
#   <first node> <second node> <stiffness> <rest length> (force fcn index) (comments)
#
# where items in ()'s are optional
open(SPRING_IN, "$input_filename.spring") || die "error: cannot open $input_filename.spring for reading: $!";
$_ = <SPRING_IN>;
chomp;
@line = split;
$num_springs_expected = $line[0];
$i = 0;
$num_springs = 0;
$line_number = 1;
while (<SPRING_IN>) {
    chomp;
    @line = split('#');
    $_ = $line[0];
    @line = split;

    $idx0 = $line[0];
    $idx1 = $line[1];
    $stf  = $line[2];
    $rst  = $line[3];

    if ($idx0 < 0) {
	print "ERROR: index 0 on line $line_number in spring file $input_filename.spring is negative\n";
    }

    if ($idx1 < 0) {
	print "ERROR: index 1 on line $line_number in spring file $input_filename.spring is negative\n";
    }

    if ($idx0 >= $num_vertices) {
	print "ERROR: index 0 on line $line_number in spring file $input_filename.spring is greater than or equal to the number of vertices in vertex file $input_filename.vertex\n";
    }

    if ($idx1 >= $num_vertices) {
	print "ERROR: index 1 on line $line_number in spring file $input_filename.spring is greater than or equal to the number of vertices in vertex file $input_filename.vertex\n";
    }

    if ($stf < 0.0) {
	print "ERROR: spring constant on line $line_number in spring file $input_filename.spring is negative\n";
    }

    if ($rst < 0.0) {
	print "ERROR: rest length on line $line_number in spring file $input_filename.spring is negative\n";
    }

    if ($check_for_consistent_rest_lengths) {
	$dx[0] = $x[$idx0] - $x[$idx1];
	$dx[1] = $y[$idx0] - $y[$idx1];
	$dx[2] = $z[$idx0] - $z[$idx1];
	$rst_actual = sqrt($dx[0]*$dx[0]+$dx[1]*$dx[1]+$dx[2]*$dx[2]);
	$diff = $rst_actual - $rst;
	if (abs($diff)/max($rst_actual,$rst,$mach_eps) > $rel_tol && abs($diff) > $abs_tol) {
	    print "WARNING: rest length on line $line_number in spring file $input_filename.spring is different from initial displacement in vertex file $input_filename.vertex\n";
	    print "         rest length in spring file $input_filename.spring: $rst\n";
	    print "         initial displacement in vertex file $input_filename.vertex: $rst_actual\n";
	}
    }

    $i += 1;
    $num_springs += 1;
    $line_number += 1;
}
print "read $num_springs springs from spring file: $input_filename.spring\n";
if ($num_springs != $num_springs_expected) {
    print "WARNING: expected to read: $num_springs_expected springs from spring file $input_filename.spring\n";
    print "         actually read:    $num_springs springs\n";
}
print "\n";
close(SPRING_IN);
