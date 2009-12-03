#!/usr/bin/perl
use Env;
use File::Path;

# announce that we are archiving files.
print "farchiver.pl:   \n";
print "farchiver.pl:   ************************************************************\n";
print "farchiver.pl:   \n";
print "farchiver.pl:   running the FARchiver...\n";
print "farchiver.pl:   \n";
print "farchiver.pl:   ************************************************************\n";
print "farchiver.pl:   \n";

# import all possible environment variables.
Env::import();
$ARCHIVE_VIZ_DIR = "$ARCHIVE_DIR/$VIZ_DIR";
$ARCHIVE_RESTART_DIR = "$ARCHIVE_DIR/$RESTART_DIR";
print "farchiver.pl:   using lock file: $LOCK_FILE_NAME\n";
print "farchiver.pl:   using vizualization directory: $VIZ_DIR\n";
print "farchiver.pl:   using restart directory: $RESTART_DIR\n";
print "farchiver.pl:   using archive directory: $ARCHIVE_DIR\n";
print "farchiver.pl:   using archive visualization directory: $ARCHIVE_VIZ_DIR\n";
print "farchiver.pl:   using archive restart directory: $ARCHIVE_RESTART_DIR\n";
print "farchiver.pl:   \n";

# make sure that the lock file does not exist.
if (-e $LOCK_FILE_NAME) {
    die "error: lock file named $LOCK_FILE_NAME already exists: $!";
}

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

# recursively create the needed archive directories.
print "farchiver.pl:   creating archive directories: $ARCHIVE_DIR\n";
print "farchiver.pl:                                 $ARCHIVE_VIZ_DIR\n";
print "farchiver.pl:                                 $ARCHIVE_RESTART_DIR\n";
system("$mkdir_command -p $ARCHIVE_DIR") == 0 || die "error: cannot create archive directory $ARCHIVE_DIR: $!";
system("$mkdir_command -p $ARCHIVE_VIZ_DIR") == 0 || die "error: cannot create archive directory $ARCHIVE_VIZ_DIR: $!";
system("$mkdir_command -p $ARCHIVE_RESTART_DIR") == 0 || die "error: cannot create archive directory $ARCHIVE_RESTART_DIR: $!";
print "farchiver.pl:   \n";

# archive all log files.
print "farchiver.pl:   archiving: log files\n";

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime;
$year += 1900;
$mon += 1;
$datetime = sprintf "%04d%02d%02d%02d%02d%02d", $year, $mon, $mday, $hour, $min, $sec;

print "farchiver.pl:   creating gzipped tar file...\n";
$logfile = "log_files.$datetime.tar.gz";
@files = <*.log *.curve *.coords $VIZ_DIR/*.visit>;
system("tar cfz $logfile @files") == 0 || die "error: cannot create tar file $logfile: $!";

print "farchiver.pl:   archiving gzipped tar file...\n";
system("$store_command $logfile $ARCHIVE_DIR/$logfile") == 0 || die "error: cannot store tar file $logfile to archive directory $ARCHIVE_DIR: $!";

print "farchiver.pl:   ensuring file was successfully archived...\n";
system("$get_command $ARCHIVE_DIR/$logfile $logfile.tmp") == 0 || die "error: cannot retrieve tar file $logfile from archive directory $ARCHIVE_DIR: $!";
system("diff $logfile $logfile.tmp") == 0 || die "error: source tar file $logfile and archive tar file $ARCHIVE_DIR/$logfile appear to differ: $!";

print "farchiver.pl:   cleaning up temporary files...\n";
unlink "$logfile" || die "error: cannot remove tar file $logfile: $!";
unlink "$logfile.tmp" || die "error: cannot remove tar file $logfile.tmp: $!";

print "farchiver.pl:   \n";

# archive all viz files and remove all viz files.
if (-d $VIZ_DIR) {
    while (defined($next = <$VIZ_DIR/*>)) {
	chomp($next);
	if (-d $next) {
	    @split_name = split(/$VIZ_DIR\//,$next);
	    @split_name = split(/\//,$split_name[$split_name-1]);
	    $name = "$split_name[$split_name-1].tar.gz";

	    print "farchiver.pl:   archiving: $next\n";
	    system("$rm_command -f $ARCHIVE_VIZ_DIR/$name");

	    print "farchiver.pl:   creating gzipped tar file...\n";
	    system("tar cfz $name $next") == 0 || die "error: cannot create tar file $name: $!";

	    print "farchiver.pl:   archiving gzipped tar file...\n";
	    system("$store_command $name $ARCHIVE_VIZ_DIR/$name") == 0 || die "error: cannot store tar file $name to archive directory $ARCHIVE_VIZ_DIR: $!";

	    print "farchiver.pl:   ensuring file was successfully archived...\n";
	    system("$get_command $ARCHIVE_VIZ_DIR/$name $name.tmp") == 0 || die "error: cannot retrieve tar file $name from archive directory $ARCHIVE_VIZ_DIR: $!";
	    system("diff $name $name.tmp") == 0 || die "error: source tar file $name and archive tar file $ARCHIVE_VIZ_DIR/$name appear to differ: $!";

	    print "farchiver.pl:   cleaning up temporary files...\n";
	    unlink "$name" || die "error: cannot remove tar file $name: $!";
	    unlink "$name.tmp" || die "error: cannot remove tar file $name.tmp: $!";

	    print "farchiver.pl:   deleting source files...\n";
	    rmtree($next) || die "error: cannot remove directory $next: $!";

	    print "farchiver.pl:   \n";
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

	    print "farchiver.pl:   archiving: $next\n";
	    system("$rm_command -f $ARCHIVE_RESTART_DIR/$name");

	    print "farchiver.pl:   creating gzipped tar file...\n";
	    system("tar cfz $name $next") == 0 || die "error: cannot create tar file $name: $!";

	    print "farchiver.pl:   archiving gzipped tar file...\n";
	    system("$store_command $name $ARCHIVE_RESTART_DIR/$name") == 0 || die "error: cannot store tar file $name to archive directory $ARCHIVE_RESTART_DIR: $!";

	    print "farchiver.pl:   ensuring file was successfully archived...\n";
	    system("$get_command $ARCHIVE_RESTART_DIR/$name $name.tmp") == 0 || die "error: cannot retrieve tar file $name from archive directory $ARCHIVE_RESTART_DIR: $!";
	    system("diff $name $name.tmp") == 0 || die "error: source tar file $name and archive tar file $ARCHIVE_RESTART_DIR/$name appear to differ: $!";

	    print "farchiver.pl:   cleaning up temporary files...\n";
	    unlink "$name" || die "error: cannot remove tar file $name: $!";
	    unlink "$name.tmp" || die "error: cannot remove tar file $name.tmp: $!";

	    @n = split(/$RESTART_DIR\/restore\./,$next);
	    $num = $n[$n-1];
	    if ($num != $restart_num) {
		print "farchiver.pl:   deleting source files...\n";
		rmtree($next) || die "error: cannot remove directory $next: $!";
	    }

	    print "farchiver.pl:   \n";
	}
    }
}

# announce that we have successfully archived all files.
print "farchiver.pl:   \n";
print "farchiver.pl:   ************************************************************\n";
print "farchiver.pl:   \n";
print "farchiver.pl:   ...FARchiver complete!\n";
print "farchiver.pl:   \n";
print "farchiver.pl:   ************************************************************\n";
print "farchiver.pl:   \n";
