#!/usr/bin/perl
use Env;
use File::Copy;

# import all possible environment variables.
Env::import();
print "using mpiexec command: $MPIEXEC\n";
print "using executable: $EXEC\n";
print "using options: $OPTIONS\n";
print "using vizualization directory: $VIZ_DIR\n";
print "using restart directory: $RESTART_DIR\n";

# the restart log file name
$restart_log_file = "$PWD/restart.log";

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
    print "previous restart number: $old_restart_num\n";
}

# determine the most recent restart file.
$restart_num = -1;
if (-d $RESTART_DIR) {
    while (defined($next = <$RESTART_DIR/restore.*>)) {
	chomp($next);
	if (-d $next) {
	    @n = split(/$RESTART_DIR\/restore\./,$next);
	    $num = $n[$n-1];
	    if ($num > $restart_num) {
		$restart_num = $num;
	    }
	}
    }
}

$from_restart = ($restart_num != -1);

if ($from_restart) {
    print "current restart number: $restart_num\n";
}

# make sure that the current restart number is larger than the
# previous restart number.
if ($from_restart) {
    ($restart_num > $old_restart_num) || die "error: restart numbers do not appear to be increasing: $!";
}

# if we are restarting, make backups of the VisIt master files.
if ($from_restart) {
    $visit_dumps_file = "$VIZ_DIR/dumps.visit";
    $visit_dumps_backup_file = "$visit_dumps_file.bak";
    if (-e $visit_dumps_file) {
	move($visit_dumps_file, $visit_dumps_backup_file) || die "error: cannot move $visit_dumps_file to $visit_dumps_backup_file: $!";
    }

    $silo_dumps_file = "$VIZ_DIR/lag_data.visit";
    $silo_dumps_backup_file = "$silo_dumps_file.bak";
    if (-e $silo_dumps_file) {
	move($silo_dumps_file, $silo_dumps_backup_file) || die "error: cannot move $silo_dumps_file to $silo_dumps_backup_file: $!";
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
    $command = "$MPIEXEC $EXEC $RESTART_DIR $restart_num $OPTIONS";
} else {
    $command = "$MPIEXEC $EXEC $OPTIONS";
}
$command =~ s/\s+/ /g; # remove any extra spaces

print "about to execute: $command\n";
system($command) == 0 || die "error: $command failed: $!";

# if we started from a restart file, fix the VisIt master files.
if ($from_restart) {
    $visit_dumps_file = "$VIZ_DIR/dumps.visit";
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

    $silo_dumps_file = "$VIZ_DIR/lag_data.visit";
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
}
