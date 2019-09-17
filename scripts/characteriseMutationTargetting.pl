#!/usr/bin/perl

use warnings;
use strict;

use IO::File;
use IO::Uncompress::Bunzip2 qw(bunzip2 $Bunzip2Error);

if (@ARGV != 2) {
	print "usage: $0 <input> <output>\n";
	print "where:\n\t<input> igblast, AIRR fmt, bzip2 compressed\n";
	print "\t<output> tab-delim, summary of mutation targetting in each sequence\n";
	exit;
}

chomp (my $inputFile = $ARGV[0]);
chomp (my $outputFile = $ARGV[1]);

#open and uncompress the input file
my $input = new IO::File $inputFile, or die "There was a problem opening the input file: $inputFile\n$!\n";
my $in = new IO::Uncompress::Bunzip2 $input, or die "There was a problem uncompressing the input file: $inputFile\n$Bunzip2Error\n";

#get the header from the input
my $hdr = <$in>;
chomp $hdr;
my %fields = ();
my @fieldLst = split(/\t/, $hdr, -1);
for (my $i = 0; $i < scalar @fieldLst; $i++) {
	my $field = $fieldLst[$i];

	$fields{$field} = $i;
}

#prepare the output file
open (OUT, ">$outputFile"), or die "There was a problem opening the output file: $outputFile\n$!\n";

#make a header for the output file t$vdjNt\t$glAA\t$vdjAA\t$aaChange
my $outHdr = "id\tclone\tsample\talignedPos\tregion\tfrmNt\ttoNt\tfrmAA\ttoAA\tSorR\tntLabel\ttnOrTv\thotspot\t5mer\n";
print OUT $outHdr;

#counter to give a progress update
my $cnt = 0;

#work through the igblast file
while (my $l = <$in>) {
	chomp $l;

	#print status update every 5k reads
	$cnt ++;
	if (($cnt % 5000) == 0) {
		print "Processed mutation details from $inputFile: $cnt\n";
	}

	my @data = split(/\t/, $l, -1);
	#get the ID
	my $id = $data[($fields{sequence_id})];
	my $clone = $data[($fields{cloneLabel})];
	my $sample = $data[($fields{sample})];
	#print OUT "$id\t$clone\t$sample\t";

	#get the germline and rearranged sequences to compare and parse mutations from
	my $ntVDJSeq = $data[($fields{sequence_alignment})];
	my $ntGLSeq = $data[($fields{germline_alignment})];
	my $aaVDJSeq = $data[($fields{sequence_alignment_aa})];
	my $aaGLSeq = $data[($fields{germline_alignment_aa})];

	#get offsets that delim the different sequence regions to keep track of which we are in (these are relative to the aligned sequence, not the input VDJ sequence)
	my $vStart = $data[($fields{v_alignment_start})];
	my $vEnd = $data[($fields{v_alignment_end})];
	my $dStart = $data[($fields{d_alignment_start})];
	my $dEnd = $data[($fields{d_alignment_end})];
	my $jStart = $data[($fields{j_alignment_start})];
	my $jEnd = $data[($fields{j_alignment_end})];

	my $fr1Start = $data[($fields{fwr1_start})];
	my $fr1End = $data[($fields{fwr1_end})];
	my $fr2Start = $data[($fields{fwr2_start})];
	my $fr2End = $data[($fields{fwr2_end})];
	my $fr3Start = $data[($fields{fwr3_start})];
	my $fr3End = $data[($fields{fwr3_end})];
	my $cdr1Start = $data[($fields{cdr1_start})];
	my $cdr1End = $data[($fields{cdr1_end})];
	my $cdr2Start = $data[($fields{cdr2_start})];
	my $cdr2End = $data[($fields{cdr2_end})];
	my $cdr3Start = $data[($fields{cdr3_start})];
	my $cdr3End = $data[($fields{cdr3_end})];

	#some additional offsets need for correcting the start/ends to aligned sequence versus to the input sequence
	my $alignedOffset = $data[($fields{v_sequence_start})];
	my $glStartOffset = $data[($fields{v_germline_start})];
	my $startOffset = 0;
	if ($glStartOffset < 3) {
		if ($glStartOffset == 1) {
			$startOffset = 2;
		} elsif ($glStartOffset == 2) {
			$startOffset = 1;
		}
	} else {
		if (($glStartOffset % 3) == 1) {
			$startOffset = 2;
		} elsif (($glStartOffset % 3) == 2) {
			$startOffset = 1;
		}
	}
	
	#build a hash with the region postions to make it easier to track which region a position is part of
	my %regions = ();
	for (my $i = 0; $i < $jEnd; $i++) {
		my $region = "";
		if ($i <= ($fr1End - $alignedOffset)) {
			$region = "fr1";
		} elsif ($i <= ($cdr1End - $alignedOffset)) {
			$region = "cdr1";
		} elsif ($i <= ($fr2End - $alignedOffset)) {
			$region = "fr2";
		} elsif ($i <= ($cdr2End - $alignedOffset)) {
			$region = "cdr2";
		} elsif ($i <= ($fr3End - $alignedOffset)) {
			$region = "fr3";
		} elsif ($i <= ($cdr3End - $alignedOffset)) {
			$region = "cdr3";
		} else {
			$region = "fr4";
		}

		$regions{$i} = $region;
	}

	#convert the strings to character arrays
	my @glNts = split(/\s*/, $ntGLSeq, -1);
	my @vdjNts = split(/\s*/, $ntVDJSeq, -1);
	my @glAAs = split(/\s*/, $aaGLSeq, -1);
	my @vdjAAs = split(/\s*/, $aaVDJSeq, -1);

	#hash for tracking the nature of the mutations within the current sequence
	my %muts = ();

	#confirm that the germline and vdj arrays contain the same number of chars
	#skipping any sequences that include indels
	if (($ntGLSeq !~ m/\-/) && ($ntVDJSeq !~ m/\-/ )) {
		#need to step through the nt and aa sequence simulateously
		my $j = 0;
		my $end = length $ntGLSeq;
		if ($end > ((scalar @glAAs) * 3) ) {
			$end = (scalar @glAAs) * 3;
		}
		for (my $i = $startOffset; $i < $end; $i++) {
			#get the nt at current position in the germline and VDJ
			my $glNt = $glNts[$i];
			my $vdjNt = $vdjNts[$i];
			#get the AA for germline and VDJ for the codon which the current nt is part of
			my $glAA = $glAAs[$j];
			my $vdjAA = $vdjAAs[$j];

			#determine if the change at the AA level was replacement or silent
			my $aaChange = "";
			if ($glAA eq $vdjAA) {
				$aaChange = "silent";
			} elsif ($glAA !~ m/X/) {
				$aaChange = "replacement";
			} else {
				$aaChange = "unknown";
			}

			my $region = "";
			if (exists $regions{$i}) {
				$region = $regions{$i};
			} 

			#compare the nts to determine if there was a mutation at this position
			if ($glNt eq $vdjNt) {
				#no mismatch at this position

			} elsif ($glNt =~ m/[A|G|C|T]/) {
				my $outStr = "$id\t$clone\t$sample\t$i\t$region\t$glNt\t$vdjNt\t$glAA\t$vdjAA\t$aaChange";

				#mismatch at this position, simply to variable naming
				my $frm = $glNt;
				my $to = $vdjNt;

				#keep a tally of mutation targeting to each nt regardless of the outcome (ie transition or transversion)
				if (exists $muts{$region}{changes}{$frm}) {
					$muts{$region}{changes}{$frm}++;
				} else {
					$muts{$region}{changes}{$frm} = 1;
				}

				#track whether or not this was a silent or replacement change at the AA level (note that codons with more than one nt change will get counted multiple times)
				if (exists $muts{$region}{aaChange}{$aaChange}) {
					$muts{$region}{aaChange}{$aaChange}++;
				} else {
					$muts{$region}{aaChange}{$aaChange} = 1;
				}

				#keep tally of each distinct mutation change
				my $change = $frm . ">" . $to;
				$outStr .= "\t$change";
				if (exists $muts{$region}{changes}{$change}) {
					$muts{$region}{changes}{$change}++;
				} else {
					$muts{$region}{changes}{$change} = 1;
				}

				#determine if transition (A<>G or C<>T) or transversion (A<>T or G<>C or C<>A or G<>T)
				if ((($frm =~ m/A|G/) && ($to =~ m/A|G/)) || (($frm =~ m/C|T/) && ($to =~ m/C|T/))) {
					#transition
					if (exists $muts{$region}{changes}{transition}) {
						$muts{$region}{changes}{transition}++;
					} else {
						$muts{$region}{changes}{transition} = 1;
					}
					$outStr .= "\ttransition";	
				} elsif (($frm =~ m/A|C|G|T/) && ($to =~ m/A|G|C|T/)) {
					#transverion, excluding indels by requiring A/C/G/T, also not counting any degenerate nt symbols like N/S/W etc.
					if (exists $muts{$region}{changes}{transversion}) {
						$muts{$region}{changes}{transversion}++;
					} else {
						$muts{$region}{changes}{transversion} = 1;
					}
					$outStr .= "\ttransversion";	
				} else {
					$outStr .= "\t";
				}

				#determine if at hotspot
				##can only do this if have two upstream and two downstream nts available (well can do with less but this is max requirement)
				if (($i > 2) && ($i < ($jEnd - 2))) {
					#build the 5mer with the current position at the central position, eg nninn
					#this is from the germline sequence, ignoring the potential alternation of the local sequence context by 
					#mutation as don't know the order in which muts were added to the current sequence
					my $fivemer = $glNts[($i-2)] . $glNts[($i-1)] . $glNt . $glNts[($i+1)] . $glNts[($i+2)];

					#test the 5mer for whether or not it includes any of the reported hotspots
					my $hotspot = "none";
					if ($fivemer =~ m/[A|T][A|G]C[C|T][A|C|G|T]/) {
						#[at][ag][c][ct][acgt]
						$hotspot = "WRCY";
					} elsif ($fivemer =~ m/[A|C|G|T][A|G]G[A|T][C|T]/) {
						#[acgt][ag]g[at][ct]
						$hotspot = "RGYW";
					} elsif ($fivemer =~ m/[A|C|G|T][A|T]A[A|C|G|T][A|C|G|T]/) {
						#[acgt][at]a[acgt][acgt]
						$hotspot = "WAN";
					} elsif ($fivemer =~ m/[A|C|G|T][A|C|G|T]T[A|T][A|T|C|G]/) {
						#[acgt][acgt]t[at][acgt]
						$hotspot = "NTW";
					}
					$outStr .= "\t$hotspot";

					#save both the 5mer and the hotspot
					if (exists $muts{$region}{fivemers}{$fivemer}) {
						$muts{$region}{fivemers}{$fivemer}++;
					} else {
						$muts{$region}{fivemers}{$fivemer} = 1;
					}
					if (exists $muts{$region}{hotspots}{$hotspot}) {
						$muts{$region}{hotspots}{$hotspot}++;
					} else {
						$muts{$region}{hotspots}{$hotspot} = 1;
					}
					$outStr .= "\t$fivemer\n";

				} else {
					#too close the the sequence start/end
					$outStr .= "\t\t\n";
				}

				#print the string to the output file
				print OUT $outStr;
			}
		
			#check if the index for the AA sequences needs to be incremented
			if (( ($i - $startOffset + 1) % 3) == 0) {
				$j++;
			}
		}	
	} else {
		#print "The germline and VDJ sequences differ in length.\n";
	}

	#check if any muts have been tracked for the current sequence and build an output string with the summary of the counts of different mutation types
	#foreach my $region (keys %muts) {
	#	foreach my $change (keys %{$muts{$region}{changes}}) {
	#		print "$id\t$clone\t$sample\t$region\t$change\t" . $muts{$region}{changes}{$change} . "\n";
	#	}
	#	foreach my $hotspot (keys %{$muts{$region}{hotspots}}) {
	#		print "$id\t$clone\t$sample\t$region\t$hotspot\t" . $muts{$region}{hotspots}{$hotspot} . "\n";
	#	}
	#	foreach my $mer (keys %{$muts{$region}{fivemers}}) {
	#		print "$id\t$clone\t$sample\t$region\t$mer\t" . $muts{$region}{fivemers}{$mer} . "\n";
	#	}
	#	foreach my $aaChange (keys %{$muts{$region}{aaChange}}) {
	#		print "$id\t$clone\t$sample\t$region\t$aaChange\t" . $muts{$region}{aaChange}{$aaChange} . "\n";
	#	}
	#}
	
}
#close the input and output file
close $in;
close OUT;


#done
exit;
