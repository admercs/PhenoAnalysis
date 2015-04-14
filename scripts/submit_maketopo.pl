#!/usr/bin/perl

# 2013/03/27 Reto Stockli (Blue Marble Research)
# This perl script submits a batch job to generate topography maps

# Arguments
# projdir : provides the project top level diretory where output boundary condition and input maps are located
# sitecase : provides a site case of a site list where sites are defined
# sitestring : provide a site to process out of the site case list

use warnings;
use strict;

use read_list;
use Time::Local;

# parameters not defined in initialization file
my $email = "reto.stoeckli\@meteoswiss.ch";

# for single point sites, choose a default grid size
my $deltax_default = 0.05;
my $deltay_default = 0.05;

# distribute IDL-based boundary processing on nxprocs x nyprocs jobs
my $interactive = 0;
my $nxprocs = 1;
my $nyprocs = 1;
my $nprocess = 1;
my $corespernode = 16;
my $account = "msclim";
my $queue = "normal";
my $memory = "2GB";
my $walltime = "12:00:00";

my $histogram = 0;

my $topocase = "test";

# define model executable name
my $exename="maketopo";
my $configname="phenoanalysis-config";

# --------------------
# START OF MAIN SCRIPT 
# --------------------

# read arguments
my $nargs = $#ARGV + 1;

# define runcase, site/domain case and sitestring (provide with command line argument)
my $projdir;
my $sitecase;
my $sitestring;
if ($nargs == 3) {
    $projdir = $ARGV[0];
    $sitecase = $ARGV[1];
    $sitestring = $ARGV[2];
} else {
    die "Please provide three arguments: (1) project directory (2) site/domain case (3) site/domain number";
}

# find out current and base directory (script directory)
my $basedir;
my $install;
my $curdir = `pwd`;
chomp $curdir;
if (substr($curdir,-13) eq "share\/scripts") {
    print "Running from installation location\n";
    $basedir = `cd ..\/..\/; pwd `;
    $install = 1;
} else {
    print "Running from source tree\n";
    $basedir = `cd ..\/; pwd `;
    $install = 0;
}
chomp $basedir;

my $sitedir;
my $scriptdir;
my $bindir;
my $configdir;
my $idldir;
if ($install == 1) {
    # directory for site / domain files 
    $sitedir = $basedir."\/share\/sites";
    
    # directory for scripts
    $scriptdir = $basedir."\/share\/scripts";

    # directory for compiled model code
    $bindir = $basedir."\/bin";

    # directory for configuration code
    $configdir = $basedir."\/bin";

    # directory for compiled idl code
    $idldir = $basedir."\/bin";
} else {
    # directory for site / domain files 
    $sitedir = $basedir."\/data\/sites";
    
    # directory for scripts
    $scriptdir = $basedir."\/scripts";

    # directory for compiled model code
    $bindir = $basedir."\/src";

    # directory for configuration code
    $configdir = $basedir."\/config";

    # directory for source idl code
    $idldir = $basedir."\/analysis\/idl";
}

# define system type
my $os = `uname -s`;
chomp $os;
my $host = `uname -n`;
chomp $host;

# read system configuration

# check if configuration program exists
if (! -e $configdir."\/".$configname) {
    if ($install == 1) {
        die "Please compile and install (make install) before executing this script! \n";
    } else {
        die "Please compile (make) before executing this script! \n";
    }
}

# check for availability of job submission system
my $pbs = `$configdir\/$configname --has-pbs`;
my $poe = `$configdir\/$configname --has-poe`;
my $slurm = `$configdir\/$configname --has-slurm`;

chomp $pbs;
chomp $poe;
chomp $slurm;

# check for availability of MPI
my $have_mpi = `$configdir\/$configname --has-mpi`;

chomp $have_mpi;

# read code version
my $version = `$configdir\/$configname --version`;

chomp $version;

# define site file
my $sitefile = $sitedir."\/".$sitecase.".dat";
if (! -e $sitefile) {
    die "Site file not found: $sitefile \n";
}

# sitestring is an array of sites

# sitestring is an array of sites
my @sites = split(",",$sitestring);
if ($#sites == 0) {
    # sitestring is a sequence of sites
    eval '@sites = '.$sitestring;
}
my $nsites = $#sites+1;

# read site data
my @ref = read_list::read_sitelist($sitefile,@sites);

my @siteshort=@{$ref[0]};
my @sitename=@{$ref[1]};
my @sitexmin=@{$ref[2]};
my @sitexmax=@{$ref[3]};
my @siteymin=@{$ref[4]};
my @siteymax=@{$ref[5]};
my @sitedeltax=@{$ref[6]};
my @sitedeltay=@{$ref[7]};
my @sitetype=@{$ref[8]};
my @siteprojection=@{$ref[9]};
my @siteprojparam=@{$ref[10]};
my @siteelevation=@{$ref[11]};

# set and create processing directories
if (! -d $projdir) {die "Top level project directory does not exist: $projdir \n"};

my $topodir = $projdir."\/topo";
if ($topocase ne "") {
    $topodir .= "-".$topocase;
}

if (! -d $topodir) {
    if (`mkdir $topodir` ne "") {
        die "Cannot create $topodir \n";
    }
    print "Created directory: $topodir \n";
}

my $batchdir = $topodir."\/batch";
if (! -d $batchdir) {
    if (`mkdir $batchdir` ne "") {
        die "Cannot create $batchdir";
    }
    print "Created directory: $batchdir \n";
}

my $logdir = $topodir."\/log";
if (! -d $logdir) {
    if (`mkdir $logdir` ne "") {
        die "Cannot create $logdir \n";
    }
    print "Created directory: $logdir \n";
}

# edit job submission parameters
my $nnodes = int($nprocess/$corespernode + 0.999); # number of nodes needed : ceiling of fractional number
my $nprocess2 = $corespernode * $nnodes; # update number of processes to full load of all nodes

my $suffix;
my $logfile;
my $errfile;
my $jobname;
my $batchfile;
my $jobfile;

for (my $s=0;$s<$nsites;$s++) {

    if (($sitedeltax[$s] eq "NA") || ($sitedeltax[$s] == 0.0))  {
	$sitedeltax[$s] = $deltax_default;
	$sitexmin[$s] -= 0.5*$deltax_default;
	$sitexmax[$s] += 0.5*$deltax_default;
    }
    if (($sitedeltay[$s] eq "NA") || ($sitedeltay[$s] == 0.0)) {
	$sitedeltay[$s] = $deltay_default;
	$siteymin[$s] -= 0.5*$deltay_default;
	$siteymax[$s] += 0.5*$deltay_default;
    }

    # loop through tiles and submit one processing job per tile
    for (my $yproc=1;$yproc<=$nyprocs;$yproc++) {
	for (my $xproc=1;$xproc<=$nxprocs;$xproc++) {
#    for (my $yproc=15;$yproc<=15;$yproc++) {
#	for (my $xproc=18;$xproc<=18;$xproc++) {

	    $suffix = sprintf("%02dx%02d",$xproc,$yproc);

	    if (($nxprocs > 1) || ($nyprocs > 1)) {
		$logfile = $logdir."\/".$exename.".".$suffix.".log";
		$errfile = $logdir."\/".$exename.".".$suffix.".err";
		$jobname = $exename.".".$suffix;
		$batchfile = $batchdir."\/".$exename.".".$suffix.".idl";
		$jobfile = $batchdir."\/".$exename.".".$suffix.".sh";
	    } else {
		$logfile = $logdir."\/".$exename.".log";
		$errfile = $logdir."\/".$exename.".err";
		$jobname = $exename;
		$batchfile = $batchdir."\/".$exename.".idl";
		$jobfile = $batchdir."\/".$exename.".sh";
	    }

	    if ($install == 0) {
		# write idl batch script
		print "Writing $batchfile \n";

		# create idl batch file
		open OUT,">$batchfile" or die "cannot open $batchfile for writing :$!";
		print OUT <<EOF;
.r $exename
$exename,basedir="$projdir",topocase="$topocase", \$
  lonmin0=$sitexmin[$s]d,lonmax0=$sitexmax[$s]d, \$
  latmin0=$siteymin[$s]d,latmax0=$siteymax[$s]d, \$
  dlon0=$sitedeltax[$s]d, dlat0=$sitedeltay[$s]d, region="$sitename[$s]", \$
  xproc=$xproc,yproc=$yproc,nxprocs=$nxprocs,nyprocs=$nyprocs,histo=$histogram
exit
EOF
		close OUT;
	    }

	    print "Writing $jobfile \n";
	    my $runcmd="";
	    open OUT,">$jobfile" or die "cannot open $jobfile for writing :$!";
	    print OUT <<EOF;
#!/bin/sh
#
EOF

	    if ($slurm == 1) {
# SLURM job submission parameters
		print OUT <<EOF;
#SBATCH --job-name=$jobname 
#SBATCH --nodes=$nnodes 
#SBATCH --ntasks=$nprocess 
##SBATCH --cpus-per-task=1 
##SBATCH --ntasks-per-node=$corespernode 
##SBATCH --mem=$memory
#SBATCH --mem-per-cpu=$memory
#SBATCH --time=$walltime 
#SBATCH --partition=$queue 
#SBATCH --account=$account 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=$email 
#SBATCH --output=$logfile
#SBATCH --error=$errfile

EOF
#		$runcmd="srun";
	    }

	    if ($pbs == 1) {
## PBS job submission parameters
## For instance useful for the NASA/NCCS discover Linux cluster
## or on the CSCS buin&dole cray machine or rigi linux cluster

## CSCS Dole
		print OUT <<EOF;
#PBS -M $email
#PBS -N $jobname
#PBS -l walltime=$walltime
#PBS -l select=1:ncpus=$nprocess:mpiprocs=$nprocess
#PBS -l mppwidth=$nprocess
#PBS -e $errfile
#PBS -o $logfile
#PBS -W group_list=$account

## NASA NCCS Discover
##PBS -V
##PBS -q $queue
##PBS -l mem=$memory
##PBS -l select=$nnodes:ncpus=$corespernode
##PBS -W group_list=$account

EOF
	    }

	    if ($poe == 1) {
## IBM POE LSF job submission parameters
## for instance useful for the NCAR/UCAR bluefire IBM Power6
		print OUT <<EOF;
#BSUB -J $jobname
#BSUB -n $nprocess
#BSUB -e $errfile
#BSUB -o $logfile
#BSUB -q $queue
#BSUB -W $walltime
#BSUB -P $account
#BSUB -R "span[ptile=$corespernode]"

EOF
	    }

	    if ($install == 1) {
		print OUT <<EOF;
cd $idldir; $runcmd idl -queue -rt="$exename.sav" -args basedir="$projdir" topocase="$topocase" lonmin0=$sitexmin[$s]d lonmax0=$sitexmax[$s]d latmin0=$siteymin[$s]d latmax0=$siteymax[$s]d dlon0=$sitedeltax[$s]d dlat0=$sitedeltay[$s]d region="$sitename[$s]" xproc=$xproc yproc=$yproc nxprocs=$nxprocs nyprocs=$nyprocs histo=$histogram
EOF
	    } else {
		print OUT <<EOF;
cd $idldir; $runcmd idl -queue $batchfile
EOF
	    }
	    close OUT;

	    # make job file executable
	    if (`chmod +x $jobfile` ne "") {die "Cannot make $jobfile executable \n"};

	    my $submissioncommand;
	    if ((($pbs == 0) && ($slurm == 0) && ($poe == 0)) || ($interactive == 1)) {
		# regular execution
		$submissioncommand="$jobfile";
	    } else {
		if ($pbs == 1) {
		    # submit job to PBS queue
		    $submissioncommand="qsub $jobfile";
		}
		if ($poe == 1) {
		    # submit job to POE queue
		    $submissioncommand="bsub < $jobfile";
		}
		if ($slurm == 1) {
		    # submit job to SLURM queue
		    $submissioncommand="sbatch $jobfile";
		}
	    }

	    # submit job
	    print "Executing $jobfile \n";
	    system("$submissioncommand \n");

	} # x tile loop

    } # y tile loop

} # site loop

exit;
