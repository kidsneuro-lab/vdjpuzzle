#/usr/bin/perl 

use strict;
use warnings;

my $file1 = $ARGV[0];
my $file2 = $ARGV[1];
open(FH, "<$file1") or "die cannot open file";
open(FH1, "<$file2") or "die cannot open file";
my @data1 = <FH>; my @data2=<FH1>;  
my @da1=(); my @da2=(); my @Ids=(); my @fr1=(); my @fr2=(); my @fr3=(); my @cdr1=(); my @cdr2=(); my @cdr3=();
print "mutations.nt.FR1\tmutations.aa.FR1\tmutations.nt.FR2\tmutations.aa.FR2\tmutations.nt.FR3\tmutations.aa.FR3\tmutations.nt.CDR1\tmutations.aa.CDR1\tmutations.nt.CDR2\tmutations.aa.CDR2tmutations.nt.CDR3\tmutations.aa.CDR3\n";
my @mutations=(); 
for(my $i=1; $i<(@data1); $i++)
{
	chomp($data1[$i]);
	@da1 = split('\t',$data1[$i]);
	push(@Ids, $da1[0]);
}
for(my $j=0; $j<(@Ids); $j++)
{
	chomp($Ids[$j]);
	for(my $k=1;$k<(@data2);$k++)
	{
                      chomp($data2[$k]);
                      @da2 = split('\t', $data2[$k]);
		      if($Ids[$j] =~ /$da2[0]/ &&  $da2[4] =~ /fr1/)
		      {
				print "$da2[3],$da2[10],$da2[7],$da2[8]";
		      }
		else
		{
			print "\t";
		}
	@mutations=();
	}
print "\n"
}
close(FH); close(FH1);
