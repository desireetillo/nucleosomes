#!/usr/bin/perl -w
use strict;

# Desiree Tillo
# 4 May 2010
# email: desiree.tillo[at]utoronto[dot]ca
#
# calcInvitroLassoModelScore.pl
# Script to calculate intrinsic nucleosome occupancy based on the linear model given in: 
# Tillo D and  Hughes TR.  G+C content dominates intrinsic nucleosome occupancy.  BMC Bioinformatics. 2009.  10:442.
#
# Usage: ./calcInvitroLassoModelScore.pl  <input sequence file in fasta format, one sequence per file> <output file name>
# Outputs zero-based midpoint coordinate of 150 base window<tab>Intrinsic nucleosome occupancy prediction
# Note:  This script is not optimized for speed, and time to run is dependent on sequence length and step size

my $stepSize = 1;  # currently hardcoded as 1 (predicts every base in a sequence), but feel free to change this as you see fit.
my $yIntercept = 1.788022;

if(@ARGV == 2){
    my ($seqFile, $outFileName) = @ARGV;
    my $seq = getSeq($seqFile);
    my $lengthSeq = length($seq);
    print "$lengthSeq\n";
    if($lengthSeq < 150){
	die "Sorry, this script only predicts nucleosome occupancy for sequences that are at least 150bp in length\n";
    }
    else{
	open(OUT, ">$outFileName") || die;
	for (my $i = 0; $i <= ($lengthSeq - 150); $i += $stepSize){ 
	    my $window150 = substr($seq, $i, 150);
	    my $window75 = substr($window150, 37, 75);
	    my $midpoint = $i + 74;
	    if($window150 !~ /N/){ # will not compute predictions for sequences with ambiguous bases (i.e. Ns).
		my $gc = calcWeightedGC($window75); 
		my $slide = calcWeightedSlide($window75);
		my $propeller = calcWeightedPropeller($window75);
		my $motifs = calcMotifScore($window150);
		my $modelScore = $gc + $slide + $propeller + $motifs + $yIntercept;
		print OUT "$midpoint\t",  sprintf("%.5f", $modelScore),"\n";
	    }
	}
	close OUT;
    }
}
else{
    print "calcInvitroLassoModelScore.pl script to calculate intrinsic nucleosome occupancy based on the Tillo and Hughes linear model (PMID: 20028554)\n";
    print "usage: ./calcInvitroLassoModelScore.pl  <input sequence file in fasta format> <output file name>\n";
}



# Subroutines #

sub getSeq{
    my ($input) = @_;
    my ($seq);
    open(INPUT, "$input") || die "Could not open $input:  $!\n";
    while (my $line = <INPUT>){
        chomp ($line);
        if ($line =~ /^\>/){
            next;
        }
        else {
            $seq .= uc($line);
        }
    }
    close INPUT;
    return $seq;
}


sub calcWeightedGC{
    # inputs: sequence window
    # output:  weighted proportion G+C
    my ($window) = @_;
    my ($gcBase) = 0;
    my ($gcCont);
    my ($lSeq) = length($window);
    while($window =~ /(?=[CG])/g){
        $gcBase++;
    }
    $gcCont = 1.67175 * ($gcBase / $lSeq);
    return $gcCont;
}

sub calcWeightedSlide{
    # inputs:  sequence window;
    # output:  slide value weighted by in vitro lasso model beta.  
    my ($window) = @_; 
    my ($tmp) = 0;
    my ($weightedSlide);
    my %slideParam = qw(AA      -0.03
			AT      -0.37
			AG      0.47
			AC      -0.13
			TA      0.74
			TT      -0.03
			TG      1.46
			TC      -0.07
			GA      -0.07
			GT      -0.13
			GG      0.6
			GC      0.29
			CA      1.46
			CT      0.47
			CG      0.63
			CC      0.6
			);
    for my $i(0..(length($window) - 2)){
	my $dinuc = substr($window, $i, 2);

	$tmp += $slideParam{$dinuc};
    }

    $weightedSlide = 1.31928 * ($tmp / 75);
    return $weightedSlide;
}

sub calcWeightedPropeller{
    #inputs: sequence window
    #output: propeller value weighted by in vitro lasso model weight
    my ($window) = @_;
    my ($tmp) = 0;
    my ($weightedPropeller);
    my (%propellerParam) = qw(AA      -17.3
			      AT      -16.9
			      AG      -14.3
			      AC      -6.7
			      TA      -11.1
			      TT      -17.3
			      TG      -8.6
			      TC      -15.1
			      GA      -15.1
			      GT      -6.7
			      GG      -12.8
			      GC      -11.7
			      CA      -8.6
			      CT      -14.3
			      CG      -11.2
			      CC      -12.8
			      );
    for my $i(0..(length($window) - 2)){
        my $dinuc = substr($window, $i, 2);
        $tmp += $propellerParam{$dinuc};
    }
    $weightedPropeller = 0.145742 * ($tmp / 75);
    return $weightedPropeller;

}

# 4-mer score 
sub calcMotifScore{
    #inputs:  sequence window
    #output:  weighted score for all 4-mers in model
    my ($window) = @_;
    my ($rc) = scalar reverse $window;
    $rc =~ tr/ATGC/TACG/;
    my ($totalMotifScore) = 0;
    my %motifWeights = qw(AAAA -0.10549
			  AAAT -0.07628
			  AAGT -0.03006
			  AATA -0.05055
			  AATT -0.02564
			  AGAA -0.02154
			  ATAA -0.03949
			  ATAT -0.02354
			  ATTA -0.03214
			  GAAA -0.03314
			  TATA -0.0334
			  );
    
    
    foreach my $motif(keys %motifWeights){
	chomp $motif;
	my $motifCount = 0;
	# scan both forward and reverse complement for motif
	while($window =~ /(?=($motif))/g){
	    $motifCount++;
	}
	while($rc =~ /(?=($motif))/g){
	    $motifCount++;
	}
	$totalMotifScore += ($motifCount * $motifWeights{$motif});
	
    }
    return $totalMotifScore;
}

