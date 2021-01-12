#!/usr/bin/perl -w
## ---------------------------------------------------------------------
##
## Copyright (c) 2007 - 2020 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

$start = 0;
$stop = 75775;
$skip = 2048;
$nprocs = 64;

$file_prefix = "amr0000";

$visit_file_name = "dumps.visit";
open(VISIT_OUT, ">$visit_file_name") || die "error: cannot open $visit_file_name for writing: $!";

for ($time_step_number = $start; $time_step_number <= $stop; $time_step_number += $skip) {
    printf VISIT_OUT "visit_dump.%05d/summary.samrai\n", $time_step_number;
}

close(VISIT_OUT) || die "error: cannot close $visit_file_name: $!";

$lag_data_file_name = "lag_data.visit";
open(LAG_DATA_OUT, ">$lag_data_file_name") || die "error: cannot open $lag_data_file_name for writing: $!";

for ($time_step_number = $start; $time_step_number <= $stop; $time_step_number += $skip) {
    printf LAG_DATA_OUT "lag_data.cycle_%06d/lag_data.cycle_%06d.summary.silo\n" , $time_step_number , $time_step_number;
}

close(LAG_DATA_OUT) || die "error: cannot close $lag_data_file_name: $!";

$meter_data_file_name = "meter_data.visit";
open(METER_DATA_OUT, ">$meter_data_file_name") || die "error: cannot open $meter_data_file_name for writing: $!";

for ($time_step_number = $start; $time_step_number <= $stop; $time_step_number += $skip) {
    printf METER_DATA_OUT "meter_data.cycle_%06d/meter_data.cycle_%06d.summary.silo\n" , $time_step_number , $time_step_number;
}

close(METER_DATA_OUT) || die "error: cannot close $meter_data_file_name: $!";

$list_file_name = "list.$file_prefix";
open(LIST_OUT, ">$list_file_name") || die "error: cannot open $list_file_name for writing: $!";

for ($time_step_number = $start; $time_step_number <= $stop; $time_step_number += $skip) {
    if ($time_step_number > 0) {
	print LIST_OUT "\n";  # we do not provide flow meter files
    }

    if ($time_step_number < 999999) {
	printf LIST_OUT "%s_%06d.mk\n", $file_prefix, $time_step_number;
	printf LIST_OUT "%s_%06d.xf\n", $file_prefix, $time_step_number;
    } else {
	printf LIST_OUT "%s%07d.mk\n", $file_prefix, $time_step_number;
	printf LIST_OUT "%s%07d.xf\n", $file_prefix, $time_step_number;
    }
}

close(LIST_OUT) || die "error: cannot close $list_file_name: $!";

$cat_file_name = "cat.$file_prefix.sh";
open(CAT_OUT, ">$cat_file_name") || die "error: cannot open $cat_file_name for writing: $!";

print CAT_OUT "#!/bin/sh\n";
print CAT_OUT "unset noclobber\n";
print CAT_OUT "echo \"About to concatenate files in directory \$PWD\"\n";
print CAT_OUT "echo\n";
print CAT_OUT "echo \"WARNING: this script *CAN* clobber existing files in the working directory\"\n";
print CAT_OUT "echo \"         but DOES NOT modify the constituent source data files\"\n";
print CAT_OUT "echo\n";
print CAT_OUT "read -p \"Press <Enter> to continue...\"\n";
print CAT_OUT "echo \"Concatenating files, please wait...\"\n";
print CAT_OUT "\n";

for ($time_step_number = $start; $time_step_number <= $stop; $time_step_number += $skip) {
    for ($rank = 0; $rank < $nprocs; $rank++) {
	if ($time_step_number < 999999) {
	    printf CAT_OUT "m3D_hdf5_marker_converter %s_%06d.mk.%04d.h5 %s_%06d.mk.%04d\n", $file_prefix, $time_step_number, $rank, $file_prefix, $time_step_number, $rank;
	} else {
	    printf CAT_OUT "m3D_hdf5_marker_converter %s%07d.mk.%04d.h5 %s%07d.mk.%04d\n", $file_prefix, $time_step_number, $rank, $file_prefix, $time_step_number, $rank;
	}
    }

    print CAT_OUT "\n";

    if ($time_step_number < 999999) {
	printf CAT_OUT "cat %s_%06d.mk.hdr", $file_prefix, $time_step_number;
    } else {
	printf CAT_OUT "cat %s%07d.mk.hdr", $file_prefix, $time_step_number;
    }

    for ($rank = 0; $rank < $nprocs; $rank++) {
	if ($time_step_number < 999999) {
	    printf CAT_OUT " %s_%06d.mk.%04d", $file_prefix, $time_step_number, $rank;
	} else {
	    printf CAT_OUT " %s%07d.mk.%04d", $file_prefix, $time_step_number, $rank;
	}
    }

    if ($time_step_number < 999999) {
	printf CAT_OUT " > %s_%06d.mk\n", $file_prefix, $time_step_number;
    } else {
	printf CAT_OUT " > %s%07d.mk\n", $file_prefix, $time_step_number;
    }

    print CAT_OUT "\n";

    for ($rank = 0; $rank < $nprocs; $rank++) {
	if ($time_step_number < 999999) {
	    printf CAT_OUT "rm -f %s_%06d.mk.%04d\n", $file_prefix, $time_step_number, $rank;
	} else {
	    printf CAT_OUT "rm -f %s%07d.mk.%04d\n", $file_prefix, $time_step_number, $rank;
	}
    }

    print CAT_OUT "\n";

    for ($rank = 0; $rank < $nprocs; $rank++) {
	if ($time_step_number < 999999) {
	    printf CAT_OUT "m3D_hdf5_fiber_converter %s_%06d.xf.%04d.h5 %s_%06d.xf.%04d\n", $file_prefix, $time_step_number, $rank, $file_prefix, $time_step_number, $rank;
	} else {
	    printf CAT_OUT "m3D_hdf5_fiber_converter %s%07d.xf.%04d.h5 %s%07d.xf.%04d\n", $file_prefix, $time_step_number, $rank, $file_prefix, $time_step_number, $rank;
	}
    }

    print CAT_OUT "\n";

    if ($time_step_number < 999999) {
	printf CAT_OUT "cat %s_%06d.xf.hdr", $file_prefix, $time_step_number;
    } else {
	printf CAT_OUT "cat %s%07d.xf.hdr", $file_prefix, $time_step_number;
    }

    for ($rank = 0; $rank < $nprocs; $rank++) {
	if ($time_step_number < 999999) {
	    printf CAT_OUT " %s_%06d.xf.%04d", $file_prefix, $time_step_number, $rank;
	} else {
	    printf CAT_OUT " %s%07d.xf.%04d", $file_prefix, $time_step_number, $rank;
	}
    }

    if ($time_step_number < 999999) {
	printf CAT_OUT " > %s_%06d.xf\n", $file_prefix, $time_step_number;
    } else {
	printf CAT_OUT " > %s%07d.xf\n", $file_prefix, $time_step_number;
    }

    print CAT_OUT "\n";

    for ($rank = 0; $rank < $nprocs; $rank++) {
	if ($time_step_number < 999999) {
	    printf CAT_OUT "rm -f %s_%06d.xf.%04d\n", $file_prefix, $time_step_number, $rank;
	} else {
	    printf CAT_OUT "rm -f %s%07d.xf.%04d\n", $file_prefix, $time_step_number, $rank;
	}
    }
}

close(CAT_OUT) || die "error: cannot close $cat_file_name: $!";
