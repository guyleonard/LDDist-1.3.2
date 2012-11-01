#!/usr/bin/perl -w

use Bio::AlignIO;
use Getopt::Std;
use LDDist;

use vars qw ($VERSION $Revision);
$VERSION = '1.3';
$Revision = qq /Protein LogDet PERL program and wrapper, version $VERSION from 2004-07-13/;

#
#	Program constants and defaults
#

$format='clustalw'; #Use ClustalW format for input
$BootTmpFileName = "pLDtmp".(int (rand 10000)).".trees";
$InvText="";


#
#	Program switches
#
#	-c n	number of classes when using rate heterogeneity. Default is 1, i.e. no heterogeneity
#	-u t	filename of user defined classes
#	-i n	percentage of invariant sites to exclude
#	-r d	estimate i using Sidow's capture-recapture
#	-r t	estimate i using Steel's capture-recapture
#	-b n	do n bootstrap replicates
#	-s t	source alignment format (default clustalw)
#	-d		override type to DNA
#	-v		print the version of LDDist.pm
#	-h		print this list
#

getopts ('hvc:i:b:r:s:du:', \%ProgSwitches);

if (exists $ProgSwitches{"v"}) {
	print "\n", LDDist::GetVersion(),"\n";
	print $Revision, "\n\n";
	exit;
}
if (exists $ProgSwitches{"h"}) {
	print "\nUsage: ProtLogDet.pm <InputFile >OutputFile\nAdditionally, the following program switches can be used\n	
	-s t	input alignment file format. Default clustalw; supported formats are those in AlignIO, currently
	
           fasta       FASTA format
           selex       selex (hmmer) format
           stockholm   stockholm format
           prodom      prodom (protein domain) format
           clustalw    clustalw (.aln) format
           msf         msf (GCG) format
           mase        mase (seaview) format
           bl2seq      Bl2seq Blast output
           nexus       Swofford et al NEXUS format
           pfam        Pfam sequence alignment format
           phylip      Felsenstein's PHYLIP format
           emboss      EMBOSS water and needle format
           mega        MEGA format
           meme        MEME format
           psi         PSI-BLAST format

	-c n	use rate heterogenity with n classes, calculated with Shannon-Wiener. Default is 1, i.e. no heterogeneity
	-u t	filename with vector of user defined character rate classes
	-i n	exclude n percent of invariant sites
	-r d	calculate (and apply) invariant fraction using the capture-recapture method of Sidow et al. Overrides -i
	-r t	calculate (and apply) invariant fraction using the capture-recapture method of Steele et al. Overrides -i, -r d
	-b n	do n bootstrap replicates in addition to the ordinary analysis
	-d  	override datatype to be DNA rather than Protein
	-v	print the version number and then exit
	-h	print this list and then exit\n\n\n";
	exit;
}

unless (exists $ProgSwitches{"c"}) {
	$ProgSwitches{"c"} = 1;
} else {
	if ($ProgSwitches{"c"}<1) {
		die "Can not have less than 1 rate class!";
	}
}


if (exists $ProgSwitches{"i"} && $ProgSwitches{"i"}>100) {die "Can not exclude more than 100% of invariant sites!"}
if (exists $ProgSwitches{"i"} && exists $ProgSwitches{"r"}) {die "Make up your mind: explicit fraction to exclude OR calculate by capture-recapture!"}
if (exists $ProgSwitches{"s"}) {
	$format=$ProgSwitches{"s"};
}

#if (exists $ProgSwitches{"i"} && $ProgSwitches{"i"}<1) {die "Percentage of invariant sites should be between 1 and 100, but be carefulÉ"}

#
#	Inits
#


#
#	Slurp in the alignment and do the stuffÉ
#

$InData  = Bio::AlignIO->new('-fh' => \*STDIN, '-format' => $format);

while ($Alignment = $InData->next_aln() ) {
	$Alignment->uppercase();
	$Alignment->unmatch();
	if (! $Alignment->is_flush()) {
	   $Alignment->throw("These sequences do not seem to be aligned properly - not the same length");
	}
	
	print STDERR "The file contains ",$Alignment->no_sequences()," sequences and ",$Alignment->length()," sites\n";
	print STDERR "The alphabet appears to be ",$Alignment->get_seq_by_pos(1)->alphabet,"\n";

	if ($Alignment->get_seq_by_pos(1)->alphabet eq "protein") {
		$Alignment->missing_char('X');
		$Alignment->gap_char('-');
		$DNA=0;
	} else {
#		die ("Sorry, can not do DNA/RNA sequences yet\n")
		$Alignment->missing_char('X');
		$Alignment->gap_char('-');
		$DNA=1;
	}
	if (exists $ProgSwitches{"d"}) {
		$DNA=1;
		print STDERR "Datatype overridden to DNA\n";
	}
	
	
#
#Create a new Taxon-Character matrix
#

	$theMatrix= LDDist::TCM->new();
	
#Now, fill the TCM with the sequences

	for (my $Taxon=1;$Taxon<=$Alignment->no_sequences();$Taxon++) {
		push @SeqArr, $Alignment->get_seq_by_pos($Taxon)->seq;
	}
	$theMatrix->setTCM(\@SeqArr, $DNA);
	
#
#	Using invariant sites?
#
	
	if (exists $ProgSwitches{"i"}) {
		$theMatrix->SetInvariants($ProgSwitches{"i"}/100)
	}

#
#	Estimate invariant sites from capture-recapture?
#
	
	if (exists $ProgSwitches{"r"}) {
		if (lc ($ProgSwitches{"r"}) eq "d") {
			$InvariantFraction = $theMatrix->SidowCRCInvariants();
			$EstType ="Sidow et al.";
		} else {
			$InvariantFraction = $theMatrix->SteelCRCInvariants();
			$EstType ="Steele et al.";
		}
		$theMatrix->SetInvariants($InvariantFraction);
		$InvText = "[!pinvar= $InvariantFraction, estimated according to $EstType];";
		print STDERR "Fraction of constant sites that are estimated to be invariant is ", $InvariantFraction, " (using $EstType)\n",
	}
		
#
#	Partition into rate classes if more than one, or use user defined classes (precedence)
#

	if (exists $ProgSwitches{"u"}) {
		open (RATEFILE, $ProgSwitches{"u"});
		chomp (my $theRateVec = <RATEFILE>);
		close RATEFILE;
		my @RateVec = split //,$theRateVec;
		if ($Alignment->length() == scalar @RateVec) {
			foreach my $theCode (@RateVec) { $ClassCodes{$theCode}++ }
			my $theClass=0;
			foreach my $theCode (keys %ClassCodes) {
				$Classes{$theCode}=$theClass++;
			}
			$theMatrix->NewRateClasses (scalar keys %ClassCodes);
			for (my $theSite=0;$theSite<@RateVec;$theSite++) { 
				$theMatrix->SetRateClass($theSite,$Classes{$RateVec[$theSite]}) 
			}
			
		} else {
			die "Length of rate class vector does not match the number of aligned sites!\n";
		}
	} elsif (exists $ProgSwitches{"c"} && $ProgSwitches{"c"}>1) {
		$theMatrix->SWgrouping($ProgSwitches{"c"});
		my $theRateVec="";
		for (my $theSite=0;$theSite<$Alignment->length();$theSite++) { 
				$theRateVec .= $theMatrix->GetRateClass($theSite); 
			}
		print STDERR "RateVector:\n", $theRateVec,"\n";
	}

#
#	Create a pairwise distance matrix for all pairs of taxa
#

	
	$NumberOfPairs = $theMatrix->LogDetDistances(0);
	@Distances = $theMatrix->getDistance(0);
	for (my $thePair=1;$thePair<$NumberOfPairs;$thePair++) {
		push @Distances, $theMatrix->getDistance($thePair);
	}

#
#	Create and print the NEXUS file header
#		

$PAUPCOMMANDBLOCK = "
begin paup;
	$InvText
	log file=slask.log replace=yes;
	set increase=auto;
	nj showtree=yes;
end;\n\n";
$BOOTCOMMANDBLOCK = "
begin paup;
	log file=boot.log replace=yes;
	nj treefile=$BootTmpFileName append=yes showtree=no;
end;\n\n";
$POSTBOOTCOMMANDBLOCK = "
begin paup;
	gettree file=$BootTmpFileName allblocks=yes;
	contree /majrule=yes strict=no;
	quit;
end;\n\n";
$LASTBLOCKCOMMAND="
begin paup;
	quit warntsave=no;
end;\n\n";


	my $TaxList ="";
	for (my $Taxon=1;$Taxon<=$Alignment->no_sequences();$Taxon++) {
		$TaxList.="\t\t'".$Alignment->get_seq_by_pos($Taxon)->display_id()."'\n";
	}
	print "#NEXUS\n\nbegin taxa;\n\tdimensions ntax=",$Alignment->no_sequences(),";\n\ttaxlabels\n$TaxList;\nend;\n\n";

#
#	Build the NEXUS distance matrix
#


	$DistanceMxText="begin distances;\nformat triangle = lower nodiagonal;\nmatrix\n";

	$DistanceMxText.="'".$Alignment->get_seq_by_pos(1)->display_id()."'\n";
	for ($FirstIndex=2;$FirstIndex<=$Alignment->no_sequences();$FirstIndex++) {
		$DistanceMxText.="'".$Alignment->get_seq_by_pos($FirstIndex)->display_id()."'\t";
		for ($SecondIndex=1;$SecondIndex<$FirstIndex;$SecondIndex++) {
			my $Distance = shift @Distances;	 
			$DistanceMxText.=$Distance." ";
		}
		$DistanceMxText.="\n"
	}   	 

	$DistanceMxText.=";\nend;\n\n";
	
	print $DistanceMxText;
	print $PAUPCOMMANDBLOCK;

#
#	If -b switch, do -b n bootstrap replicates
#

	if (exists $ProgSwitches{"b"}) {
	#Do bootstrap
		for (my $replicate=0;$replicate<$ProgSwitches{"b"};$replicate++) {
				$NumberOfPairs = $theMatrix->LogDetDistances(1);
				@Distances = $theMatrix->getDistance(0);
				for (my $thePair=1;$thePair<$NumberOfPairs;$thePair++) {
					push @Distances, $theMatrix->getDistance($thePair);
				}

				$DistanceMxText="begin distances;\nformat triangle = lower nodiagonal;\nmatrix\n";
				$DistanceMxText.="'".$Alignment->get_seq_by_pos(1)->display_id()."'\n";
				for ($FirstIndex=2;$FirstIndex<=$Alignment->no_sequences();$FirstIndex++) {
					$DistanceMxText.="'".$Alignment->get_seq_by_pos($FirstIndex)->display_id()."'\t";
					for ($SecondIndex=1;$SecondIndex<$FirstIndex;$SecondIndex++) {
						my $Distance = shift @Distances;	 
						$DistanceMxText.=$Distance." ";
					}
					$DistanceMxText.="\n"
				}   	 

				$DistanceMxText.=";\nend;\n\n";
				print $DistanceMxText;
				print $BOOTCOMMANDBLOCK;

		}
		print $POSTBOOTCOMMANDBLOCK;
	}
	print $LASTBLOCKCOMMAND;
$theMatrix->DESTROY();
}
