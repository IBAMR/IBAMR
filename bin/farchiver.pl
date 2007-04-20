#!/usr/bin/perl
use Env;
use File::Path;

# import all possible environment variables.
Env::import();
$ARCHIVE_VIZ_DIR = "$ARCHIVE_DIR/$VIZ_DIR";
$ARCHIVE_RESTART_DIR = "$ARCHIVE_DIR/$RESTART_DIR";
print "using vizualization directory: $VIZ_DIR\n";
print "using restart directory: $RESTART_DIR\n";
print "using archive directory: $ARCHIVE_DIR\n";
print "using archive visualization directory: $ARCHIVE_VIZ_DIR\n";
print "using archive restart directory: $ARCHIVE_RESTART_DIR\n";

# setup the archiver commands.
$get_command = "far get";  # source-file target-file
$store_command = "far store";  # source-file target-file
$rget_command = "far rget";  # source-dir target-dir
$rstore_command = "far rstore";  # source-dir target-dir
$cp_command = "far cp";
$mkdir_command = "far mkdir";
$mv_command = "far mv";
$rm_command = "far rm";
$rmdir_command = "far rmdir";

# create the needed archive directories.
system("$mkdir_command -p $ARCHIVE_DIR") || die "error: cannot create archive directory $ARCHIVE_DIR: $!";
system("$mkdir_command -p $ARCHIVE_VIZ_DIR") || die "error: cannot create archive directory $ARCHIVE_VIZ_DIR: $!";
system("$mkdir_command -p $ARCHIVE_RESTART_DIR") || die "error: cannot create archive directory $ARCHIVE_RESTART_DIR: $!";

# archive all log files.
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime;
$year += 1900;
$mon += 1;
$datetime = sprintf "%04d%02d%02d%02d%02d%02d", $year, $mon, $mday, $hour, $min, $sec;
$logfile = "log_files.$datetime.tar.bz2";
@files = <*.log *.curve *.coords $VIZ_DIR/*.visit>;
system("tar cvfj $logfile @files") == 0 || die "error: cannot create tar file $logfile: $!";
system("$store_command $logfile $ARCHIVE_DIR/$logfile") || die "error: cannot store tar file $logfile to archive directory $ARCHIVE_DIR: $!";
system("$get_command $ARCHIVE_DIR/$logfile $logfile.tmp") || die "error: cannot retrieve tar file $logfile from archive directory $ARCHIVE_DIR: $!";
system("diff $logfile $logfile.tmp") == 0 || die "error: source tar file $logfile and archive tar file $ARCHIVE_DIR/$logfile appear to differ: $!";
unlink "$logfile" || die "error: cannot remove tar file $logfile: $!";
unlink "$logfile.tmp" || die "error: cannot remove tar file $logfile.tmp: $!";

# archive all viz files and remove all viz files.
if (-d $VIZ_DIR) {
    while (defined($next = <$VIZ_DIR/*>)) {
	chomp($next);
	if (-d $next) {
	    @split_name = split(/$VIZ_DIR\//,$next);
	    @split_name = split(/\//,$split_name[$split_name-1]);
	    $name = "$split_name[$split_name-1].tar.gz";

	    # archive the current viz file.
	    print "archiving: $next\n";
	    system("$rm_command -f $ARCHIVE_VIZ_DIR/$name");
	    system("tar cfz $name $next") == 0 || die "error: cannot create tar file $name: $!";
	    system("$store_command $name $ARCHIVE_VIZ_DIR/$name") || die "error: cannot store tar file $name to archive directory $ARCHIVE_VIZ_DIR: $!";
	    system("$get_command $ARCHIVE_VIZ_DIR/$name $name.tmp") || die "error: cannot retrieve tar file $name from archive directory $ARCHIVE_VIZ_DIR: $!";
	    system("diff $name $name.tmp") == 0 || die "error: source tar file $name and archive tar file $ARCHIVE_VIZ_DIR/$name appear to differ: $!";
	    unlink "$name" || die "error: cannot remove tar file $name: $!";
	    unlink "$name.tmp" || die "error: cannot remove tar file $name.tmp: $!";

	    # delete all viz files.
	    rmtree($next) || die "error: cannot remove directory $next: $!";
	}
    }
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

# archive all restart files and remove all but the most recent restart
# file.
if (-d $RESTART_DIR) {
    while (defined($next = <$RESTART_DIR/restore.*>)) {
	chomp($next);
	if (-d $next) {
	    @split_name = split(/$RESTART_DIR\//,$next);
	    @split_name = split(/\//,$split_name[$split_name-1]);
	    $name = "$split_name[$split_name-1].tar.gz";

	    # archive all restart files.
	    print "archiving: $next\n";
	    system("$rm_command -f $ARCHIVE_RESTART_DIR/$name");
	    system("tar cfz $name $next") == 0 || die "error: cannot create tar file $name: $!";
	    system("$store_command $name $ARCHIVE_RESTART_DIR/$name") || die "error: cannot store tar file $name to archive directory $ARCHIVE_RESTART_DIR: $!";
	    system("$get_command $ARCHIVE_RESTART_DIR/$name $name.tmp") || die "error: cannot retrieve tar file $name from archive directory $ARCHIVE_RESTART_DIR: $!";
	    system("diff $name $name.tmp") == 0 || die "error: source tar file $name and archive tar file $ARCHIVE_RESTART_DIR/$name appear to differ: $!";
	    unlink "$name" || die "error: cannot remove tar file $name: $!";
	    unlink "$name.tmp" || die "error: cannot remove tar file $name.tmp: $!";

	    # delete any older restart files.
	    @n = split(/$RESTART_DIR\/restore\./,$next);
	    $num = $n[$n-1];
	    if ($num != $restart_num) {
		rmtree($next) || die "error: cannot remove directory $next: $!";
	    }
	}
    }
}
