#!/usr/bin/perl
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

use File::Copy;
use Env;

# import all possible environment variables.
Env::import();

# load the ConfigReader module.
use lib "$HOME/code/IBAMR/src/tools/modules";
use ConfigReader::Simple;

# announce that we are running the autostarter.
print "autostarter.pl:   \n";
print "autostarter.pl:   ************************************************************\n";
print "autostarter.pl:   \n";
print "autostarter.pl:   running the AUTO (re-)starter...\n";
print "autostarter.pl:   \n";
print "autostarter.pl:   ************************************************************\n";
print "autostarter.pl:   \n";

# determine the job settings by reading the input file.
$num_args = $#ARGV + 1;
if ($num_args != 1) {
    die "error: usage: autostarter.pl <input file>: $!";
}
$input_file = $ARGV[0];
print "autostarter.pl:   processing input file: $input_file\n";
print "autostarter.pl:   \n";

$config = ConfigReader::Simple->new($input_file);
die "error: could not read input file $input_file: $ConfigReader::Simple::ERROR" unless ref $config;

$working_directory = $config->get("working_dir");
$lock_file_name = $config->get("lock_file_name");
$executable = $config->get("executable");
$options = $config->get("options");
$viz_dir = $config->get("viz_dir");
$restart_dir = $config->get("restart_dir");

print "autostarter.pl:   using working directory: $working_directory\n";
print "autostarter.pl:   using lock file: $lock_file_name\n";
print "autostarter.pl:   using executable: $executable\n";
print "autostarter.pl:   using options: $options\n";
print "autostarter.pl:   using vizualization directory: $viz_dir\n";
print "autostarter.pl:   using restart directory: $restart_dir\n";
print "autostarter.pl:   \n";

chdir $working_directory || die "error: could not change working directory to $working_directory: $!";

# make sure that the lock file does not exist.
if (-e $lock_file_name) {
    die "error: lock file named $lock_file_name already exists: $!";
}

# set the restart log file name.
$restart_log_file = "$working_directory/restart.log";

# determine the timestep of the most recent previous restart.
$old_restart_num = -1;
if (-e $restart_log_file) {
    open(IN, "$restart_log_file") || die "error: cannot open $restart_log_file for reading: $!";
    while (<IN>) {
	chomp;
	$old_restart_num = $_;
    }
    close(IN) || die "error: cannot close $restart_log_file: $!";
}

if ($old_restart_num != -1) {
    print "autostarter.pl:   previous restart number: $old_restart_num\n";
}

# determine the most recent restart file.
$restart_num = -1;
if (-d $restart_dir) {
    while (defined($next = <$restart_dir/restore.*>)) {
	chomp($next);
	if (-d $next) {
	    @n = split(/$restart_dir\/restore\./,$next);
	    $num = $n[$n-1];
	    if ($num > $restart_num) {
		$restart_num = $num;
	    }
	}
    }
}

$from_restart = ($restart_num != -1);

if ($from_restart) {
    print "autostarter.pl:   current restart number: $restart_num\n";
} else {
    print "autostarter.pl:   initial program invokation (not from restart!)\n";
}

# make sure that the current restart number is larger than the
# previous restart number.
if ($from_restart) {
    ($restart_num > $old_restart_num) || die "error: restart numbers do not appear to be increasing: $!";
}

# if we are restarting, make backups of the VisIt master files.
if ($from_restart) {
    $visit_dumps_file = "$viz_dir/dumps.visit";
    $visit_dumps_backup_file = "$visit_dumps_file.bak";
    if (-e $visit_dumps_file) {
	move($visit_dumps_file, $visit_dumps_backup_file) || die "error: cannot move $visit_dumps_file to $visit_dumps_backup_file: $!";
    }

    $silo_dumps_file = "$viz_dir/lag_data.visit";
    $silo_dumps_backup_file = "$silo_dumps_file.bak";
    if (-e $silo_dumps_file) {
	move($silo_dumps_file, $silo_dumps_backup_file) || die "error: cannot move $silo_dumps_file to $silo_dumps_backup_file: $!";
    }

    $meter_dumps_file = "$viz_dir/meter_data.visit";
    $meter_dumps_backup_file = "$meter_dumps_file.bak";
    if (-e $meter_dumps_file) {
	move($meter_dumps_file, $meter_dumps_backup_file) || die "error: cannot move $meter_dumps_file to $meter_dumps_backup_file: $!";
    }
}

# if we are restarting, update the restart.log file.
if ($from_restart)
{
    open(LOGFILE, ">>$restart_log_file") || die "error: cannot open $restart_log_file for writing: $!";
    if ($old_restart_num != -1)
    {
	print LOGFILE "\n";
    }
    print LOGFILE $restart_num;
    close(LOGFILE) || die "error: cannot close $restart_log_file: $!";
}

# execute the command.
if ($from_restart) {
    $command = "$executable $restart_dir $restart_num $options";
} else {
    $command = "$executable $options";
}
$command =~ s/\s+/ /g; # remove any extra spaces

print "autostarter.pl:   about to execute: $command\n";
print "autostarter.pl:   \n";
print "autostarter.pl:   ************************************************************\n";
print "autostarter.pl:   \n";

system($command) == 0 || die "error: execution of $command failed";

if (-e $lock_file_name) {
    die "error: lock file named $lock_file_name exists following command execution: $!";
}

print "autostarter.pl:   \n";
print "autostarter.pl:   ************************************************************\n";
print "autostarter.pl:   \n";
print "autostarter.pl:   successfully executed: $command\n";

# if we started from a restart file, fix the VisIt master files.
if ($from_restart) {
    $visit_dumps_file = "$viz_dir/dumps.visit";
    $visit_dumps_backup_file = "$visit_dumps_file.bak";
    if (-e $visit_dumps_file) {
	$last = "";
	if (-e $visit_dumps_backup_file) {
	    open(IN, "<", $visit_dumps_backup_file) || die "error: cannot open $visit_dumps_backup_file for reading: $!";
	    while ($line = <IN>) {
		$last = $line;
	    }
	    close(IN) || die "error: cannot close $visit_dumps_backup_file: $!";
	}

	open(OUT, ">>", $visit_dumps_backup_file) || die "error: cannot open $visit_dumps_backup_file for writing: $!";
	open(IN, "<", $visit_dumps_file) || die "error: cannot open $visit_dumps_file for reading: $!";
	$first = <IN>;
	if (!($last eq $first)) {
	    print OUT $first;
	}
	while ($line = <IN>) {
	    print OUT $line;
	}
	close(OUT) || die "error: cannot close $visit_dumps_backup_file: $!";
	close(IN) || die "error: cannot close $visit_dumps_file: $!";
    }
    move($visit_dumps_backup_file,$visit_dumps_file);

    $silo_dumps_file = "$viz_dir/lag_data.visit";
    $silo_dumps_backup_file = "$silo_dumps_file.bak";
    if (-e $silo_dumps_file) {
	$last = "";
	if (-e $silo_dumps_backup_file) {
	    open(IN, "<", $silo_dumps_backup_file) || die "error: cannot open $silo_dumps_backup_file for reading: $!";
	    while ($line = <IN>) {
		$last = $line;
	    }
	    close(IN) || die "error: cannot close $silo_dumps_backup_file: $!";
	}

	open(OUT, ">>", $silo_dumps_backup_file) || die "error: cannot open $silo_dumps_backup_file for writing: $!";
	open(IN, "<", $silo_dumps_file) || die "error: cannot open $silo_dumps_file for reading: $!";
	$first = <IN>;
	if (!($last eq $first)) {
	    print OUT $first;
	}
	while ($line = <IN>) {
	    print OUT $line;
	}
	close(OUT) || die "error: cannot close $silo_dumps_backup_file: $!";
	close(IN) || die "error: cannot close $silo_dumps_file: $!";
    }
    move($silo_dumps_backup_file,$silo_dumps_file);

    $meter_dumps_file = "$viz_dir/meter_data.visit";
    $meter_dumps_backup_file = "$meter_dumps_file.bak";
    if (-e $meter_dumps_file) {
	$last = "";
	if (-e $meter_dumps_backup_file) {
	    open(IN, "<", $meter_dumps_backup_file) || die "error: cannot open $meter_dumps_backup_file for reading: $!";
	    while ($line = <IN>) {
		$last = $line;
	    }
	    close(IN) || die "error: cannot close $meter_dumps_backup_file: $!";
	}

	open(OUT, ">>", $meter_dumps_backup_file) || die "error: cannot open $meter_dumps_backup_file for writing: $!";
	open(IN, "<", $meter_dumps_file) || die "error: cannot open $meter_dumps_file for reading: $!";
	$first = <IN>;
	if (!($last eq $first)) {
	    print OUT $first;
	}
	while ($line = <IN>) {
	    print OUT $line;
	}
	close(OUT) || die "error: cannot close $meter_dumps_backup_file: $!";
	close(IN) || die "error: cannot close $meter_dumps_file: $!";
    }
    move($meter_dumps_backup_file,$meter_dumps_file);
}

# announce that we have successfully run the executable
print "autostarter.pl:   \n";
print "autostarter.pl:   ************************************************************\n";
print "autostarter.pl:   \n";
print "autostarter.pl:   ...AUTO (re-)starter complete!\n";
print "autostarter.pl:   \n";
print "autostarter.pl:   ************************************************************\n";
print "autostarter.pl:   \n";
