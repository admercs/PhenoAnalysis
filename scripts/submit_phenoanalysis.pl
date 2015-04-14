#!/usr/bin/perl

# 2013/03/27 Reto Stockli (Blue Marble Research)
# This perl script submits a Pheno Analysis batch job

# Arguments
# runcase : provide the Pheno Analysis case name as a double quoted string
# sitecase: provides the Site list case name as a double quoted string
# sitestring : provide a list of sites to process inside a double quoted string

use warnings;
use strict;

use read_list;
use Time::Local;

# parameters not defined in initialization file
my $email = "reto.stoeckli\@meteoswiss.ch";

# define model executable name
my $exename="phenoanalysis";
my $configname="phenoanalysis-config";

# --------------------
# START OF MAIN SCRIPT 
# --------------------

# read arguments
my $nargs = $#ARGV + 1;

# define runcase, site/domain case and sitestring (provide with command line argument)
my $runcase;
my $sitecase;
my $sitestring;
if ($nargs == 3) {
    $runcase = $ARGV[0];
    $sitecase = $ARGV[1];
    $sitestring = $ARGV[2];
} else {
    die "Please provide three arguments: (1) model case (2) site/domain case (3) site/domain number(s)";
}

# find out current and base directory (script directory)
my $basedir;
my $install;
my $interactive;
my $curdir = `pwd`;
chomp $curdir;
if (substr($curdir,-13) eq "share\/scripts") {
    print "Running from installation location\n";
    $basedir = `cd ..\/..\/; pwd `;
    $install = 1;
    $interactive = 0; # run job interactively (no submission, even on systems with job submission)
} else {
    print "Running from source tree\n";
    $basedir = `cd ..\/; pwd `;
    $install = 0;
    $interactive = 1; # run job interactively (no submission, even on systems with job submission)
}
chomp $basedir;

my $inidir;
my $sitedir;
my $scriptdir;
my $exedir;
my $configdir;
if ($install == 1) {
    # directory for initialization files
    $inidir = $basedir."\/share\/ini\/".$runcase;
    
    # directory for site / domain files 
    $sitedir = $basedir."\/share\/sites";
    
    # directory for scripts
    $scriptdir = $basedir."\/share\/scripts";

    # directory for compiled model code
    $exedir = $basedir."\/bin";

    # directory for configuration code
    $configdir = $basedir."\/bin";
} else {
    # directory for initialization files
    $inidir = $basedir."\/data\/ini\/".$runcase;
    
    # directory for site / domain files 
    $sitedir = $basedir."\/data\/sites";
    
    # directory for scripts
    $scriptdir = $basedir."\/scripts";

    # directory for compiled model code
    $exedir = $basedir."\/src";

    # directory for configuration code
    $configdir = $basedir."\/config";
}


# read model namelist file
my $namelistfile = "pheno.nml";
if (! -e $inidir."\/".$namelistfile) {
    die "Model namelist file not found: $inidir\/$namelistfile \n";
}

my $regional = read_list::read_namelist($inidir."\/".$namelistfile,"regional");
my $nrens = read_list::read_namelist($inidir."\/".$namelistfile,"nrens");
my $resubmit = read_list::read_namelist($inidir."\/".$namelistfile,"resubmit");
my $restart = read_list::read_namelist($inidir."\/".$namelistfile,"restart");

my $analysis =  read_list::read_namelist($inidir."\/".$namelistfile,"analysis");
my $saveobs =  read_list::read_namelist($inidir."\/".$namelistfile,"saveobs");
my $simulator = read_list::read_namelist($inidir."\/".$namelistfile,"simulator");

my $nprocess = read_list::read_namelist($inidir."\/".$namelistfile,"nprocess");

my $date_start = read_list::read_namelist($inidir."\/".$namelistfile,"date_start");
my $date_end = read_list::read_namelist($inidir."\/".$namelistfile,"date_end");
my $loop = read_list::read_namelist($inidir."\/".$namelistfile,"loop");

my $walltime = read_list::read_namelist($inidir."\/".$namelistfile,"walltime");
my $account = read_list::read_namelist($inidir."\/".$namelistfile,"account");
my $queue = read_list::read_namelist($inidir."\/".$namelistfile,"queue");
my $memory = read_list::read_namelist($inidir."\/".$namelistfile,"memory");
my $cpupernode =  read_list::read_namelist($inidir."\/".$namelistfile,"cpupernode");

my $siteloop = read_list::read_namelist($inidir."\/".$namelistfile,"siteloop");

my $meteorologytype = read_list::read_namelist($inidir."\/".$namelistfile,"meteorologytype");
my $pfttype = read_list::read_namelist($inidir."\/".$namelistfile,"pfttype");
my $landcovertype = read_list::read_namelist($inidir."\/".$namelistfile,"landcovertype");
my $topotype = read_list::read_namelist($inidir."\/".$namelistfile,"topotype");
my $satellitetype = read_list::read_namelist($inidir."\/".$namelistfile,"satellitetype");

my $projectdir = read_list::read_namelist($inidir."\/".$namelistfile,"projectdir");
my $datadir = read_list::read_namelist($inidir."\/".$namelistfile,"datadir");
my $meteorologydir = read_list::read_namelist($inidir."\/".$namelistfile,"meteorologydir");
my $pftdir = read_list::read_namelist($inidir."\/".$namelistfile,"pftdir");
my $topodir = read_list::read_namelist($inidir."\/".$namelistfile,"topodir");
my $landcoverdir = read_list::read_namelist($inidir."\/".$namelistfile,"landcoverdir");
my $satellitedir = read_list::read_namelist($inidir."\/".$namelistfile,"satellitedir");

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
my $mpi = `$configdir\/$configname --has-mpi`;

chomp $mpi;

# read code version
my $version = `$configdir\/$configname --version`;

chomp $version;

# make sure not distributed run is selected without mpi support
if (($mpi == 0) && ($nprocess gt 1)) {
    die "MPI support is not compiled, please therefore turn subregion processing off.\n";
}

# only use the observation simulator for analysis runs
if (($analysis eq ".false.") && ($simulator eq ".true.")) {
    die "Please set Analysis to .true. when simulator is set to .true. \n";
}

# make sure we choose enough ensemble members for the analysis
if (($analysis eq ".true.") && ($nrens < 5)) {
    die "Please select more than 5 ensemble members when the analysis is set to .true. \n";
}

# determine number of parallel processes
if (($mpi == 0) && ($nprocess gt 1)) {
    die "Please set nprocess to 1 when MPI support is not compiled \n";
}

# make sure restart is set when resubmitting jobs
if ((($resubmit == 1) || ($loop > 1)) && ($restart eq ".false.")) {
    die "Set restart to .true. if resubmission or loop > 1 is chosen"
}

# set and create processing directories
if (! -d $projectdir) {die "Top level project directory does not exist: $projectdir \n"};

my $rundir = $projectdir."\/".$runcase;

if (! -d $rundir) {
    if (`mkdir $rundir` ne "") {
        die "Cannot create $rundir \n";
    }
    print "Created directory: $rundir \n";
}

my $batchdir = $rundir."\/batch";
if (! -d $batchdir) {
    if (`mkdir $batchdir` ne "") {
        die "Cannot create $batchdir";
    }
    print "Created directory: $batchdir \n";
}

my $logdir = $rundir."\/log";
if (! -d $logdir) {
    if (`mkdir $logdir` ne "") {
        die "Cannot create $logdir \n";
    }
    print "Created directory: $logdir \n";
}

my $outdir = $rundir."\/output";
if (! -d $outdir) {
    if (`mkdir $outdir` ne "") {
        die "Cannot create $outdir \n";
    }
    print "Created directory: $outdir \n";
}

my $restartdir = $rundir."\/restart";
if (! -d $restartdir) {
    if (`mkdir $restartdir` ne "") {
        die "Cannot create $restartdir \n";
    }
    print "Created directory: $restartdir \n";
}

# set input data directories
if (! -d $datadir) {die "Top level data directory does not exist: $datadir \n"};
if (substr($datadir,-1,1) ne "\/") {
    $datadir = $datadir."\/";
}

# set directory containing meteorological data
if (substr($meteorologydir,0,1) ne "\/") {
    $meteorologydir = $datadir.$meteorologydir;
}
if (substr($meteorologydir,-1,1) ne "\/") {
    $meteorologydir = $meteorologydir."\/";
}
if (! -d $meteorologydir) {die "Meteorological data directory does not exist: $meteorologydir \n"};

# set directory containing plant functional type data
if (substr($pftdir,0,1) ne "\/") {
    $pftdir = $datadir.$pftdir;
}
if (substr($pftdir,-1,1) ne "\/") {
    $pftdir = $pftdir."\/";
}
if (! -d $pftdir) {die "PFT data directory does not exist: $pftdir \n"};

# set directory containing MODIS landcover data
if (substr($landcoverdir,0,1) ne "\/") {
    $landcoverdir = $datadir.$landcoverdir;
}
if (substr($landcoverdir,-1,1) ne "\/") {
    $landcoverdir = $landcoverdir."\/";
}
if (! -d $landcoverdir) {die "Landcover data directory does not exist: $landcoverdir \n"};

# set directory containing elevation data
if (substr($topodir,0,1) ne "\/") {
    $topodir = $datadir.$topodir;
}
if (substr($topodir,-1,1) ne "\/") {
    $topodir = $topodir."\/";
}
if (! -d $topodir) {die "Topo data directory does not exist: $topodir \n"};

# set directory containing MODIS FPAR&LAI assimilation data
if (substr($satellitedir,0,1) ne "\/") {
    $satellitedir = $datadir.$satellitedir;
}
if (substr($satellitedir,-1,1) ne "\/") {
    $satellitedir = $satellitedir."\/";
}
if (! -d $satellitedir) {die "Input satellite data directory does not exist: $satellitedir \n"};

# check data files
my $forcingsfile = $inidir."\/forcings.dat";
if (! -e $forcingsfile) {
    die "Forcings file not found : $forcingsfile \n";
}

my $observationsfile = $inidir."\/observations.dat";
if (! -e $observationsfile) {
    die "Observations file not found : $observationsfile \n";
}

my $parametersfile = $inidir."\/parameters.dat";
if (! -e $parametersfile) {
    die "Parameters file not found : $parametersfile \n";
}

my $statesfile = $inidir."\/states.dat";
if (! -e $statesfile) {
    die "States file not found : $statesfile \n";
}

my $predictionparametersfile;
if ($analysis eq ".false.") {
# copy prediction parameter files to runtime directory
    if ($regional eq ".true.") {
	$predictionparametersfile = $inidir."\/prediction_parameters_regional.dat";
    } else {
	$predictionparametersfile = $inidir."\/prediction_parameters_local.dat";
    }

    if (! -e $predictionparametersfile) {
	die "Prediction Parameter file not found : $predictionparametersfile \n";
    }
}

my $simulatorparametersfile;
if ($simulator eq ".true.") {
# copy  simulator parameter files to runtime directory
    if ($regional eq ".true.") {
	$simulatorparametersfile = $inidir."\/simulator_parameters_regional.dat";
    } else {
	$simulatorparametersfile = $inidir."\/simulator_parameters_local.dat";
    }

    if (! -e $simulatorparametersfile) {
	die "Simulator Parameter file not found : $simulatorparametersfile \n";
    }
}


# edit job submission parameters
my $nodes = int($nprocess/$cpupernode + 0.999); # number of nodes needed : ceiling of fractional number

my $logfile = $logdir."\/".$exename.".log";
my $errfile = $logdir."\/".$exename.".err";
my $jobname = $exename."-".$runcase;


# ---------------------------------------------
# loop over sites or simulate all sites at once
# ---------------------------------------------
my $sitefile = $sitedir."\/".$sitecase.".dat";
if (! -e $sitefile) {
    die "Site file not found: $sitefile \n";
}

# sitestring is an array of sites
my @sites = split(",",$sitestring);
if ($#sites == 0) {
    # sitestring is a sequence of sites
    eval '@sites = '.$sitestring;
}
my $nsites = $#sites+1;

my ($site0, $site1);
if ($siteloop == 1) {
    $site0 = 0;
    $site1 = $nsites-1;
} else {
    $site0 = 0;
    $site1 = 0;
}
for (my $site=$site0;$site<=$site1;$site++) {

    # if we loop over time series, always resubmit the job for each loop
    if ($loop > 1) { $resubmit = 1 }

    # write batch script for job submission
    my $jobfile = $batchdir."\/".$exename.".sh";

    # generate submission command
    my $submissioncommand;
    if ($interactive == 1) {
	$submissioncommand="$jobfile";
    } else {
	if ($pbs == 1) {
    	    # submit job to PBS queue
    	    $submissioncommand="qsub $jobfile";
    	} elsif ($poe == 1) {
    	    # submit job to POE queue
    	    $submissioncommand="bsub < $jobfile";
        } elsif ($slurm == 1) {
    	    # submit job to SLURM queue
    	    $submissioncommand="sbatch $jobfile";
        } else {
    	    # regular execution
	    $submissioncommand="nohup $jobfile > $logfile 2> $errfile < /dev/null &";
	}
    }

    print "Writing $jobfile \n";
    open OUT,">$jobfile" or die "cannot open $jobfile for writing :$!";
    print OUT <<EOF;
#!/bin/sh
#
EOF

    my $runcmd;
    if ($pbs == 1) {
	# add PBS job parameters
	print OUT <<EOF;
#PBS -M $email
#PBS -N $jobname
#PBS -m a
#PBS -l cput=$walltime
#PBS -l select=1:ncpus=1:mem=$memory
#PBS -q $queue
#PBS -r n 
#PBS -o $logfile
#PBS -e $errfile
#PBS -V

## NASA NCCS Discover
##PBS -V
##PBS -q $queue
##PBS -l mem=$memory
##PBS -l select=$nodes:ncpus=$cpupernode
##PBS -W group_list=$account

EOF
	if (($mpi == 1) && ($nprocess > 1)) {
	    $runcmd="mpirun -np $nprocess";
	} else {
	    $runcmd="";
	}
    } elsif ($slurm == 1) {
	# add SLURM job parameters
	if (substr($host,0,4) eq "lema") {
		print OUT <<EOF;
#SBATCH --mail-user=$email
#SBATCH --job-name=$jobname
#SBATCH --mail-type=FAIL
#SBATCH --time=$walltime
#SBATCH --nodes=$nodes
#SBATCH --ntasks=$nprocess
#SBATCH --ntasks-per-node=$cpupernode
#SBATCH --mem=$memory
#SBATCH --partition=$queue
#SBATCH --account=$account
#SBATCH --output=$logfile
#SBATCH --error=$errfile

EOF
	} else {
		print OUT <<EOF;
#SBATCH --mail-user=$email
#SBATCH --job-name=$jobname
#SBATCH --mail-type=FAIL
#SBATCH --time=$walltime
#SBATCH --nodes=1
#SBATCH --ntasks=$nprocess
#SBATCH --mem=$memory
#SBATCH --partition=$queue
#SBATCH --account=$account
#SBATCH --output=$logfile
#SBATCH --error=$errfile

EOF
	}

	if (substr($host,0,4) eq "lema") {
	    $runcmd = "aprun -B";
	} else {
	    if (($mpi == 1) && ($nprocess > 1)) {
		    $runcmd="srun -n $nprocess";
		} else {
		    $runcmd="";
		}
	}
    } elsif ($poe == 1) {
	# add POE job parameters
	print OUT <<EOF;
#BSUB -a poe
#BSUB -J $jobname
#BSUB -n $nprocess
#BSUB -o $logfile
#BSUB -e $errfile
#BSUB -q $queue
#BSUB -W $walltime
#BSUB -P $account
#BSUB -R "span[ptile=$cpupernode]"

EOF
	if (($mpi == 1) && ($nprocess > 1)) {
	    $runcmd = "mpirun.lsf";
	} else {
	    $runcmd = "mpirun.lsf";
	}
    } else {
	if (($mpi == 1) && ($nprocess > 0)) {
	    $runcmd = "mpirun -np $nprocess";
	} else {
	    $runcmd = "";
	}
    }

    # reset loop, start date and end date at start of submission 
    if ($resubmit == 1) {
	print OUT <<EOF;

resubmit=$resubmit
date_start=$date_start
date_curr=$date_start
date_end=$date_end
l=1

EOF
    }

    # remove successful flag file
    print OUT <<EOF;
cd $batchdir

if [ -e "model_successful" ]
then
    rm "model_successful"
fi

EOF

    my $execmd = "$runcmd $exedir\/$exename $inidir\/$namelistfile $sitefile";
    if ($siteloop == 1) {
	$execmd .= " $sites[$site]";
    } else {
	$execmd .= " $sitestring";
    }
    
    if ($resubmit == 1) {
	$execmd .= " \$date_curr";
    }

	print OUT <<EOF;
$execmd

EOF

# abort batch script if success flag is not set after execution
# why so complicated? Unlike in C, there is no standard exit status for Fortran
    print OUT <<EOF;
if [ ! -e "model_successful" ]
then
    echo "Model terminated. Please check Log File."
    exit 1
fi

EOF

# for resubmission, let the submission shell script modify new resubmission dates in namelist
    if ($resubmit == 1) {

	print OUT <<EOF;
case \${#date_curr} in
4)
    date_curr=\$((\$date_curr + 1))
    ;;
6)
    year_curr=\$((\$date_curr / 100))
    month_curr=\$((\$date_curr - ( \$year_curr * 100 )))

    month_curr=\$((\$month_curr + 1))
    if [ \$month_curr -gt 12 ]
    then
       month_curr=1
       year_curr=\$((\$year_curr + 1))
    fi
    date_curr=\$((\$year_curr * 100 + \$month_curr))
    ;;
8)
    year_curr=\$(( \$date_curr / 10000 ))
    month_curr=\$(( (\$date_curr - (\$year_curr * 10000)) / 100 ))
    day_curr=\$(( \$date_curr - (\$year_curr * 10000) - (\$month_curr * 100) ))

    ndayofmonth=31
    if ([ \$month_curr == 4 ] || [ \$month_curr == 6 ] || [ \$month_curr == 9 ] || [ \$month_curr == 11 ])
    then
       ndayofmonth=30
    fi
    if [ \$month_curr == 2 ]
    then
        if (([ \$(( (\$year_curr % 4) )) -eq 0 ] && [ \$(( (\$year_curr % 100) )) -ne 0 ]) || [ \$(( (\$year_curr % 400) )) -eq 0 ])
        then
           ndayofmonth=29
        else
           ndayofmonth=28
        fi
    fi

    (( day_curr++ ))
    if [ \$day_curr -gt \$ndayofmonth ]
    then
       day_curr=1
       month_curr=\$(( \$month_curr + 1 ))
       if [ \$month_curr -gt 12 ]
       then
          month_curr=1
          year_curr=\$(( \$year_curr + 1 ))
       fi
    fi
    date_curr=\$(( \$year_curr * 10000 + \$month_curr * 100 + \$day_curr ))
    ;;
*)
    echo "sub-daily resubmission not implemented"
    exit 1
    ;;
esac

if [ \$date_curr -gt \$date_end ]
then
    if [ \$l -lt $loop ]
    then
       l=\$(( \$l + 1 ))
       date_curr=\$date_start
    else
       resubmit=0
    fi
fi

if [ \$resubmit -eq 1 ]
then
    echo "new resubmission date: \$date_curr loop: \$l"
    sed -e '1,/date_curr=.*/s/date_curr=.*/date_curr='\$date_curr'/' -e '1,/l=.*/s/l=.*/l='\$l'/' $jobfile > $jobfile.tmp
    mv -f $jobfile.tmp $jobfile
    chmod +x $jobfile

    $submissioncommand
fi

EOF
    }

    close OUT;

    # make job file executable
    system("chmod +x $jobfile \n");

    # submit job
    system("$submissioncommand \n");

} # site loop

exit;
