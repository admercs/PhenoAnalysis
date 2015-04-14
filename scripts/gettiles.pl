#!/usr/bin/perl
#
#SBATCH --job-name=gettiles2013
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=24:00:00 
#SBATCH --partition=normal 
#SBATCH --account=msclim 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=reto.stoeckli@meteoswiss.ch 
#SBATCH --output=gettiles2013.log
#SBATCH --error=gettiles2013.err
#
#
# 2013/07/15 Reto Stockli (Blue Marble Research)
#
# Perl script fetches MODIS LAND (TERRA satellite) tiles the edc ftp data pool 
# and put them in date-specific directories if they don't exist yet locally 
# also check for dupliates and remove them if needed
use strict;
use warnings;
use diagnostics;
use LWP::Simple;

# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

my $forceget = 0;      # get file even if it is stored locally
my $datapool = 2;      # get from http data pool (2), ftp data pool OBSOLETE (1) or order system WIST/ECHO (0)
my $sensor   = "MOLT"; # type of sensor archive (MOLT, MOLA, MOTA, SRTM, WELD, ASTT)

#my $basedir = "/Users/stockli/modis/";
my $basedir = "/store/msclim/stockli/modis/";
#my $basedir = "/glade/p/cgd/tss/people/stockli/modis/";
my $year_start = 2013;
my $year_end = 2013;
my $prefix = "MOD15A2";
#my $prefix = "MOD44B";
#my $prefix = "MOD44W";
#my $prefix = "MCD12Q1";
my $version = "005";
my $tile = ""; # look for a specific tile? (set to "" if you want global search!)
#my $tile = "h18v04";
#my $tile = "h12v04";
#my $tile = "h20v07";
#my $tile = "h22v01";
#my $tile = "h11v05";
#my $tile = "h12v09";
#my $tile = "h12v03";
#my $tile = "h08v05";

my $remoteget;
my $remotels;
my $server;
my $serverdir;
if ($datapool == 0) {
    $remoteget = "ncftpget";
    $remotels = "ncftpls";
    $server = "ftp://e4ftl01u.ecs.nasa.gov";
    $serverdir = "/PullDir/0301780854hsdlL";
} elsif ($datapool == 1) {
    $remoteget = "ncftpget";
    $remotels = "ncftpls -x -F";
    $server = "ftp://e4ftl01.cr.usgs.gov";
    $serverdir = "/".$sensor."/".$prefix.".".$version;
} else {
    $server = "http://e4ftl01.cr.usgs.gov";
    $serverdir = "/".$sensor."/".$prefix.".".$version;
}

my $dt = 8; # number of days: MODIS product time step (1, 8 or 16 days)

my @allserverfiles;
my @allserverdates;
if ($datapool == 0) {
    # fetch EDC FTP Pull directory content
    @allserverfiles = grep (!/xml/, grep (/.hdf/, `$remotels $server$serverdir/`));    
} elsif ($datapool == 1) {
} else {
    # fetch http directory content
    my $dirs = &GetDirectories($server,$serverdir);
    @allserverdates = @$dirs;
}

my $proddir = $basedir.$prefix.".".$version;
if (-d $proddir) {} else { `mkdir $proddir`; }

foreach (@allserverdates) {

    (my $year, my $month, my $day) = split(/\./);

    if (($year >= $year_start) & ($year <= $year_end)) {

	my @monthsum;
	my $leapyear;
	if (((($year % 4) == 0) & (($year % 100) != 0)) || (($year % 400) == 0)) {
	    $leapyear = 1;
	    @monthsum = (0,31,60,91,121,152,182,213,244,274,305,335,366);
	} else {
	    $leapyear = 0;
	    @monthsum = (0,31,59,90,120,151,181,212,243,273,304,334,365);
	}

	# Generate YYYYDDD
	my $doy = $monthsum[$month-1] + $day;

        my $date1 = "$year".sprintf("%3d",$doy);
	$date1=~s/\s/0/g;

	# Generate YYYY.MM.DD
        my $date2 = "$year.".sprintf("%2d",$month).".".sprintf("%2d",$day);
	$date2=~s/\s/0/g;
		
	my $outdir = $proddir."\/".$date1."\/";
	
	if (-d $outdir) {} else { `mkdir $outdir`; }
	
	# fetch list of tiles available for that date from server
	my @serverfiles;
	if ($datapool == 0) {
    	    if ($tile ne "") {
		@serverfiles = grep (/A$date1/, grep (/$tile/, @allserverfiles));
	    } else {	    
		@serverfiles = grep (/A$date1/, @allserverfiles);
	    }
	} elsif ($datapool == 1) {
    	    if ($tile ne "") {
		@serverfiles = grep (/A$date1/, grep (/$tile/, grep (!/xml/, grep (/.hdf/, `$remotels $server$serverdir/$date2/`))));
	    } else {	    
		@serverfiles = grep (/A$date1/, grep (!/xml/, grep (/.hdf/, `$remotels $server$serverdir/$date2/`)));
	    }
	} else {
	    my $files = &GetFiles($server, $serverdir."/".$date2);

    	    if ($tile ne "") {
		
		@serverfiles = grep (/A$date1/, grep (/$tile/, grep (!/xml/, grep (/.hdf/,@$files))));
	    } else {	    
		@serverfiles = grep (/A$date1/, grep (!/xml/, grep (/.hdf/,@$files)));
	    }
	}

	my $nserverfiles = $#serverfiles + 1;
	
	for (my $f=0; $f<$nserverfiles; $f++) {
	    my $serverfile = trim($serverfiles[$f]);
	    
	    my $servertile = substr($serverfile,length($prefix)+10,6);

	    # check if we already have this tile
	    opendir(DIR,$outdir);
	    my @alllocalfiles = readdir(DIR);
	    closedir DIR;
	    my @localfiles = grep(/$prefix.A$date1.$servertile/,@alllocalfiles);
	    my $nlocalfiles = $#localfiles + 1;
	    		
	    # if we have multiple copies stored locally, remove the copies
	    for (my $l=1;$l<$nlocalfiles;$l++) {
		my $localfile = trim($localfiles[$l]);
		chomp $localfile;
		
		print "Removing local copy No. $l+1: $outdir/$localfile \n"; 

		`\\rm -f $outdir/$localfile `;
	    }

	    if (($nlocalfiles == 0) || ($forceget == 1)) {
		print "downloading: $serverfile \n";
	    
		if ($datapool == 0) {
		    `cd $outdir; $remoteget $server$serverdir/$serverfile `;
		} elsif ($datapool == 1) {
		    `cd $outdir; $remoteget $server$serverdir/$date2/$serverfile `;
		} else {
		    `cd $outdir; curl $server$serverdir/$date2/$serverfile -OL -s `;
		}
	    } else {
		print "$serverfile already locally available \n";
	    }
	} # next file on server
	
    } # date in requested year range
} # next date

###############################################################################
#
#  GetDirectories( <host>, <directory> )
#  
# Gets a list of sub-directories in a given directory. The sub directories
# are returned as relative paths, rather than full paths.
#
###############################################################################
sub GetDirectories
{
  my $host = shift;
  my $dir = shift;

  # Check that the passed in directory begins with a '/'
    
  $dir = '/' if !$dir;
  $dir = '/' . $dir if ( $dir !~ m/^\//o );


  # Build the URL for the resource

#  my $url = 'http://' . $host . $dir;
  my $url = $host . $dir;
  $url = $url . '/' if ( $dir !~ m/\/$/o );

  
  # Retrieve the content
   
  my $content = get($url);
  my @subdirs = ();

  if (defined $content) {
      while ($content =~ m/<a href=\"([\w\-\/\.]*?)\">/go)
	{
	    # If this is a directory, chop the ending '/' off and push it on the list.
	    # Note that we exclude anything begining with '/' to eliminate the 'Parent Directory'
	    # entry 
	    
	    push @subdirs, $1 if( $1 =~ m/(^[^\/].*)\/$/o );
	}
  }
  
  
  # Return the list
  
  return \@subdirs;
}

###############################################################################
#
#  GetFiles( <host>, <directory> )
#  
# Gets a list of files in a given directory. The files
# are returned as relative paths, rather than full paths.
#
###############################################################################
sub GetFiles
{
  my $host = shift;
  my $dir = shift;

  # Check that the passed in directory beings with a '/'
    
  $dir = '/' if !$dir;
  $dir = '/' . $dir if ( $dir !~ m/^\//o );


  # Build the URL for the resource

#  my $url = 'http://' . $host . $dir;
  my $url = $host . $dir;
  $url = $url . '/' if ( $dir !~ m/\/$/o );

  
  # Retrieve the content
   
  my $content = get($url);
  my @files = ();

  if (defined $content) {
      while ($content =~ m/<a href=\"([\w\-\/\.]*?)\">/go)
	{
	    # Files should not begin or end with a '/'. Anything that does is considered
	    # a directory
	    
	    push @files, $1 if( $1 =~ m/(^[^\/].*[^\/])$/o );
	}
  }
  
  # Return the list
  
  return \@files;
}

