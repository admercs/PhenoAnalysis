#!/usr/bin/perl

$logfile = "/Users/stockli/pheno_analysis/data/parameters/prediction_parameters_regional_Global-Prediction.dat";

@lines = `more $logfile | grep "Total obs used:"`;

$nlines = $#lines;

$nobs = 0;

for ($l=0;$l<$nlines;$l++) {

    @temp=trim(split(" ",$lines[$l]));
    print "$lines[$l]\n";
    $nobs += $temp[9];
    

}

print "$nobs\n";

sub trim {
    my @out = @_;
    for (@out) {
	s/^\s+//;
	s/\s+$//;
    }
    return wantarray ? @out : $out[0];
}
				

