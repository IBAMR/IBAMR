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
# filename: beam_validator3d.pl
# author: Boyce Griffith
# usage: beam_validator3d.pl <input filename>
#
# A simple Perl script to determine whether the beam specifications in
# an IBAMR beam file are valid, and whether the specified beam
# curvatures are consistent with the initial curvatures in the
# corresponding IBAMR vertex file.
#
# NOTE: This utility will generate a large number of warning messages
# for sturctures where the initial curvatures and specified curvatures
# differ.

use List::Util qw[min max];

if ($#ARGV != 0) {
    die "incorrect number of command line arguments.\nusage:\n  beam_validator3d.pl <input filename>\n";
}

# boolean controling whether to check to see if the initial
# displacement equals the resting length.
$check_for_consistent_curvatures = 1;  # 1 = enabled; 0 = disabled
$mach_eps = 2.22044604925031e-16;
$rel_tol = sqrt($mach_eps);
$abs_tol = 20*$mach_eps;

print "\n";
if ($check_for_consistent_curvatures) {
    print "checks for consistent beam curvatures are ENABLED!\n";
    print "warning messages WILL be printed if the initial curvatures and specified curvatures differ by more than a relative error tolerance of $rel_tol or an absolute error tolerance of $abs_tol.\n";
} else {
    print "checks for consistent beam rest lengths are DISABLED!\n";
    print "warning messages WILL NOT be printed if the initial curvatures and specified curvatures differ.\n";
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

# the first line in the beam input file has the format:
#
#   <number of edges> (comments)
#
# and the remaining lines in the input file all have the format:
#
#   <first node> <second node> <second node> <stiffness> (<curvature x> <curvature y> <curvature z>) (comments)
#
# where items in ()'s are optional
open(BEAM_IN, "$input_filename.beam") || die "error: cannot open $input_filename.beam for reading: $!";
$_ = <BEAM_IN>;
chomp;
@line = split;
$num_beams_expected = $line[0];
$i = 0;
$num_beams = 0;
$line_number = 1;
while (<BEAM_IN>) {
    chomp;
    @line = split('#');
    $_ = $line[0];
    @line = split;

    $idx0 = $line[0];
    $idx1 = $line[1];
    $idx2 = $line[2];
    $bnd  = $line[3];

    if ($idx0 < 0) {
	print "ERROR: index 0 on line $line_number in beam file $input_filename.beam is negative\n";
    }

    if ($idx1 < 0) {
	print "ERROR: index 1 on line $line_number in beam file $input_filename.beam is negative\n";
    }

    if ($idx2 < 0) {
	print "ERROR: index 2 on line $line_number in beam file $input_filename.beam is negative\n";
    }

    if ($idx0 >= $num_vertices) {
	print "ERROR: index 0 on line $line_number in beam file $input_filename.beam is greater than or equal to the number of vertices in vertex file $input_filename.vertex\n";
    }

    if ($idx1 >= $num_vertices) {
	print "ERROR: index 1 on line $line_number in beam file $input_filename.beam is greater than or equal to the number of vertices in vertex file $input_filename.vertex\n";
    }

    if ($idx2 >= $num_vertices) {
	print "ERROR: index 2 on line $line_number in beam file $input_filename.beam is greater than or equal to the number of vertices in vertex file $input_filename.vertex\n";
    }

    if ($bnd < 0.0) {
	print "ERROR: beam constant on line $line_number in beam file $input_filename.beam is negative\n";
    }

    if ($check_for_consistent_curvatures && $#line >= 6) {
	$Lx = $line[4];
	$Ly = $line[5];
	$Lz = $line[6];

	$Lx_actual = $x[$idx0] - 2.0*$x[$idx1] + $x[$idx2];
	$Ly_actual = $y[$idx0] - 2.0*$y[$idx1] + $y[$idx2];
	$Lz_actual = $z[$idx0] - 2.0*$z[$idx1] + $z[$idx2];

	$diff = $Lx_actual - $Lx;
	if (abs($diff)/max(abs($Lx_actual),abs($Lx),$mach_eps) > $rel_tol && abs($diff) > $abs_tol) {
	    print "WARNING: x curvature on line $line_number in beam file $input_filename.beam is different from x curvature in vertex file $input_filename.vertex\n";
	    print "         x curvature in beam file $input_filename.beam: $Lx\n";
	    print "         initial x curvature in vertex file $input_filename.vertex: $Lx_actual\n";
	    print "         difference: $diff\n";
	    print "\n";
	}

	$diff = $Ly_actual - $Ly;
	if (abs($diff)/max(abs($Ly_actual),abs($Ly),$mach_eps) > $rel_tol && abs($diff) > $abs_tol) {
	    print "WARNING: y curvature on line $line_number in beam file $input_filename.beam is different from y curvature in vertex file $input_filename.vertex\n";
	    print "         y curvature in beam file $input_filename.beam: $Ly\n";
	    print "         initial y curvature in vertex file $input_filename.vertex: $Ly_actual\n";
	    print "         difference: $diff\n";
	    print "\n";
	}

	$diff = $Lz_actual - $Lz;
	if (abs($diff)/max(abs($Lz_actual),abs($Lz),$mach_eps) > $rel_tol && abs($diff) > $abs_tol) {
	    print "WARNING: z curvature on line $line_number in beam file $input_filename.beam is different from z curvature in vertex file $input_filename.vertex\n";
	    print "         z curvature in beam file $input_filename.beam: $Lz\n";
	    print "         initial z curvature in vertex file $input_filename.vertex: $Lz_actual\n";
	    print "         difference: $diff\n";
	    print "\n";
	}
    }

    $i += 1;
    $num_beams += 1;
    $line_number += 1;
}
print "read $num_beams beams from beam file: $input_filename.beam\n";
if ($num_beams != $num_beams_expected) {
    print "WARNING: expected to read: $num_beams_expected beams from beam file $input_filename.beam\n";
    print "         actually read:    $num_beams beams\n";
}
print "\n";
close(BEAM_IN);
