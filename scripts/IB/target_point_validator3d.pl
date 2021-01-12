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
# filename: target_point_validator3d.pl
# author: Boyce Griffith
# usage: target_point_validator3d.pl <input filename>
#
# A simple Perl script to determine whether the target point
# specifications in an IBAMR target point file are valid.

use List::Util qw[min max];

if ($#ARGV != 0) {
    die "incorrect number of command line arguments.\nusage:\n  target_point_validator3d.pl <input filename>\n";
}

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

# the first line in the target point input file has the format:
#
#   <number of target point> (comments)
#
# and the remaining lines in the input file all have the format:
#
#   <node> <stiffness> (damping factor) (force fcn index) (comments)
#
# where items in ()'s are optional
open(TARGET_POINT_IN, "$input_filename.target") || die "error: cannot open $input_filename.target for reading: $!";
$_ = <TARGET_POINT_IN>;
chomp;
@line = split;
$num_target_points_expected = $line[0];
$i = 0;
$num_target_points = 0;
$line_number = 1;
while (<TARGET_POINT_IN>) {
    chomp;
    @line = split('#');
    $_ = $line[0];
    @line = split;

    $idx = $line[0];
    $stf = $line[1];

    if ($idx < 0) {
	print "ERROR: index on line $line_number in target point file $input_filename.target is negative\n";
    }

    if ($idx >= $num_vertices) {
	print "ERROR: index on line $line_number in target point file $input_filename.target is greater than or equal to the number of vertices in vertex file $input_filename.vertex\n";
    }

    if ($stf < 0.0) {
	print "ERROR: target point spring constant on line $line_number in target point file $input_filename.target is negative\n";
    }

    if (@line >= 3) {
	$rst = $line[2];
	if ($rst < 0.0) {
	    print "ERROR: target point damping constant line $line_number in target point file $input_filename.target is negative\n";
	}
    }

    $i += 1;
    $num_target_points += 1;
    $line_number += 1;
}
print "read $num_target_points target points from target point file: $input_filename.target\n";
if ($num_target_points != $num_target_points_expected) {
    print "WARNING: expected to read: $num_target_points_expected target points from target point file $input_filename.target\n";
    print "         actually read:    $num_target_points target points\n";
}
print "\n";
close(TARGET_POINT_IN);
