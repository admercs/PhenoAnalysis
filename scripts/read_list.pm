package read_list;

# 2013/03/28 Reto Stockli (Blue Marble Research)
# This perl script reads values from autoconf config.h, from ascii namelist and from a geographical site list

use warnings;
use strict;

# read argument value from config.h file
sub read_config {

    my ($configfile, $argument) = @_;
    my $value = "";
    my (@array, @array2, @array3);

    open(IN,trim($configfile));
    my $line=<IN>;
    while (!eof) {
	$line=<IN>;
	chomp $line;
	
	if (trim($line) ne "") {
	    @array = split("#define",trim($line));

	    if ($#array > 0) {
		@array2 = split(" ",trim($array[1]));
		
		if (trim($array2[0]) eq trim($argument)) {
		    # found argument, check if it has a value
		    if ($#array2 > 0) {
			# remove double quotes of value if present
			@array3 = split(/[",']/,trim($array2[1]));
			if (trim($array3[0]) eq "") {
			    $value = trim($array3[1]);
			} else {
			    $value = trim($array3[0]);
			}
		    } else {
			# or assume default value of #defined statement
			$value = 1
		    }
		}
	    }
	}
    }
    close IN;
    
    return $value;
}    

# read argument value from namelist file
sub read_namelist {

    my ($namelistfile, $argument) = @_;
    my $value = "";
    my (@array, @array2, @array3);

    open(IN,$namelistfile);
    my $line=<IN>;
    while (!eof) {
	$line=<IN>;
	chomp $line;
	
	if (trim($line) ne "") {    
	    @array = split("=",$line);
	    
	    if (trim($array[0]) eq $argument) {
		# found value to argument, now remove double quotes and comments if present
		@array2 = split(/[#,!]/,trim($array[1]));

		if (trim($array2[0]) ne "") {
		    @array3 = split(/[",']/,trim($array2[0]));
		    if ($#array3 == 1) {
			if (trim($array3[1]) ne "") {
			    $value = trim($array3[1]);
			}
		    } elsif ($#array3 == 0) {
			if (trim($array3[0]) ne "") {
			    $value = trim($array3[0]);
			}
		    } else {
			$value = "";
		    }
		} else {
		    $value = "";
		}
	    } 	
	}
    }
    close IN;
    
    return $value;
}    

# get site data
sub read_sitelist {

    my ($sitefile, @sites) = @_;

    # clear arrays
    my @short = ();
    my @name = ();
    my @xmin = ();
    my @xmax = ();
    my @ymin = ();
    my @ymax = ();
    my @deltax = ();
    my @deltay = ();
    my @type = ();
    my @projection = ();
    my @projparam = ();
    my @elevation = ();


    my $nsites = $#sites+1;
    my $i=0;
    my $count=0;
    my @array;
    my $id;
    my $line;

    while (($i<$nsites) && ($count < 25)) {

	open(IN,$sitefile ) || die "file $sitefile not found \n";
	$line=<IN>;
	$line=<IN>;
	while (!eof) {
	    $line=<IN>;

	    if ($i<$nsites) {
		chomp $line;
		@array=split(",",$line);
		
		$id=trim($array[0]);
		
		if ($id == $sites[$i]) {
		    push(@short,trim($array[1]));
		    push(@name,trim($array[2]));
		    push(@xmin,trim($array[3]));
		    push(@xmax,trim($array[4]));
		    push(@ymin,trim($array[5]));
		    push(@ymax,trim($array[6]));
		    push(@deltax,trim($array[7]));
		    push(@deltay,trim($array[8]));
		    push(@type,trim($array[9]));
		    push(@projection,trim($array[10]));
		    push(@projparam,trim($array[11]));
		    push(@elevation,trim($array[12]));
		    $i++;
		}
	    }
	}
	close(IN);
	$count++;
    }

    if ($i<$nsites) {
	print "Not all sites were found in $sitefile . \n";
	die;
    }
    
    return(\@short,\@name,\@xmin,\@xmax,\@ymin,\@ymax,\@deltax,\@deltay,\@type,\@projection,\@projparam,\@elevation);
    
}

# trim white spaces from string
sub trim($_) {
    my $str = $_[0];
    $str=~s/^\s+|\s+$//g;
    return $str;
}

1;
