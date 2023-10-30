#!/usr/bin/perl -w

use lib "$ENV{\"HOME\"}/Perl_Codes/Jun_lib/";
use NA_SEQ;
#Add other packages here with syntax: use package;

$Usage="==========================================
|| Jun Lu, Apr 2014
|| this code removes short reads from a fa file, with a cutoff set
||
|| Usage: perl $0 -f <input_fa> -o <output_fa> -L <length_cutoff>
==========================================
||-f <input_fa>:	Input fa file
||-o <output_fa>:	Output fa file
||-L <length_cutoff>:	length of cutoff; optional; default 15 (meaning retaining everything >=15nt
==========================================
";

#default values for parameters
$input_fa=0;
$output_fa=0;
$length_cutoff=15;

for($i=0;$i<@ARGV;$i+=2)
{
	if($ARGV[$i] eq "-f")
	{
		$input_fa=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-o")
	{
		$output_fa=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-L")
	{
		$length_cutoff=$ARGV[$i+1];
	}
	else
	{
		
	}
}
if(!$input_fa || !$output_fa) {die "$Usage\n";}
#============= start code here ============================
$below_cutoff_count=0;

$seqs=NA_SEQ::read_FASTA_file($input_fa);
open($outfile,">$output_fa") || die "couldn't open $output_fa outfile\n";

for($i=0;$i<@{$seqs};$i++)
{
    if($seqs->[$i]->length()>=$length_cutoff)
    {
	$seqs->[$i]->write_FASTA_seq($outfile);
    }
    else
    {
	$below_cutoff_count++;
    }
}

close($outfile);
print "Seqs below $length_cutoff Nt: $below_cutoff_count.\n";


