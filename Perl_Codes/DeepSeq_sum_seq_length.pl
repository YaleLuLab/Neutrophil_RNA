#!/usr/bin/perl -w

use lib "/home/jun/Perl_Codes/Jun_lib/";
use NA_SEQ;

$Usage="==========================================
|| Jun Lu, Feb 2012
|| Usage: perl $0 -f <fasta_file> -A <adaptor> -s <second_adaptor> -u <use_2nd_adaptor> -o <outfilename>
||
|| This program uses collapsed fasta file to calculate the length distribution of
||	clipped reads. The use of adaptor sequences are optional to calculate primer dimers
==========================================
||-f <fasta_file>:	Name of fasta file here
||-A <adaptor>:	Illumina adaptor seq; default: TGGAATTCTCGGGTGCCAAGG
||-s <second_adaptor>:	Minor form of adaptor; default: TGGACTTCTCGGGTGCCAAGG
||-u <use_2nd_adaptor>:	whether to use 2nd adaptor seq; default 1 (use)
||-o <outfilename>:	output fa file name
==========================================
";

#default values for parameters
$fasta_file=0;
$adaptor="TGGAATTCTCGGGTGCCAAGG";
$second_adaptor="TGGACTTCTCGGGTGCCAAGG";
$use_2nd_adaptor=1;
$outfilename=0;

for($i=0;$i<@ARGV;$i+=2)
{
	if($ARGV[$i] eq "-f")
	{
		$fasta_file=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-A")
	{
		$adaptor=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-s")
	{
		$second_adaptor=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-u")
	{
		$use_2nd_adaptor=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-o")
	{
		$outfilename=$ARGV[$i+1];
	}
	else
	{
		
	}
}
if(!$fasta_file || !$outfilename) {die "$Usage\n";}

$seqs=NA_SEQ::read_FASTA_file($fasta_file);

%size_hash=();

($total_reads, $reads_with_N, $reads_start_N, $adaptor_dimers)=(0,0,0,0);

for($i=0; $i<@{$seqs};$i++)
{
	$id=$seqs->[$i]->name;
	$seq=$seqs->[$i]->seq;
	$seq_len=length($seq);

	if($id=~/^(\S+)_\S+_x(\d+)$/)
	{
		$temp_count=$2;
		$total_reads+=$temp_count;
	}
	else
	{
		$temp_count=0;
		print STDERR "WARNing: $id is of incorrect format\n";
		exit;
	}

	# now calculate if adaptor dimer
	if($seq=~/^$adaptor/)
	{
		$adaptor_dimers+=$temp_count;
	}
	elsif($use_2nd_adaptor && $seq=~/^$second_adaptor/)
	{
		$adaptor_dimers+=$temp_count;
	}
	else
	{
		# now calculate if containing N
		if($seq=~/N/)
		{	
			$reads_with_N+=$temp_count;
			if($seq=~/^N/)
			{
				$reads_start_N+=$temp_count;
			}
		}
	
		# now parse into hash with length
		if(exists($size_hash{$seq_len}))
		{
			$size_hash{$seq_len}+=$temp_count;
		}
		else
		{
			$size_hash{$seq_len}=$temp_count;
		}
	}
}
	
#output summary
open(OUT,">$outfilename") || die "couldn't open output file $outfilename\n";
print OUT "Data File\t$fasta_file\n";
print OUT "Total Reads\t$total_reads\n";
print OUT "PrimerDimer\t$adaptor_dimers\n";
print OUT "Reads with N\t$reads_with_N\n";
print OUT "Reads Start with N\t$reads_start_N\n";	

print OUT "Length\tReads\n";



foreach $templen (sort {$a<=>$b} keys %size_hash)
{
	print OUT "$templen\t$size_hash{$templen}\n";
}
close(OUT);

exit;






	




