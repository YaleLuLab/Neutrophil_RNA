#!/usr/bin/perl -w

use lib "/home/jun/Perl_Codes/Jun_lib/";
#Add other packages here with syntax: use package;

$Usage="==========================================
|| Jun Lu, Feb 2012; modified May2017; increased mode2 adaptor match length, and uses non-greedy match to get rid of adaptor concatemers
|| Usage: perl $0 -f <fasta_file> -A <adaptor> -M <mode> -B <bases> -L <minLength> -u <use_2nd_adaptor> -s <second_adaptor> -o <outfile>
|| 
|| This code was modified from clip_adaptors.pl from miRDeep2 package
|| This code aims to get rid of 3'adaptor sequences from small RNA seq data
|| The original miRDeep code uses 6 bases of the input adatpor and restricts length of small 
||                 RNA to be at least 18 bases, or cutoff any base match to the adpator at the very end
||This program provides the flexibiity of parameter setting using mode of 1
|| In addition, mode 2 is especially suited for data from 50bp reads or more from HiSeq
||                 Specifically, in mode 2, the first pass match uses a defined length of match (default 15 bases of adaptor), and only requiring seq to have at least 1 base in front of adatpor, then in second pass, try to match at least
||                 6 bases after a defined minimal length for the small RNA (default to 40);
||Another modification is the addition of choice to use a second adaptor (could be the major unclipped form of sequence variation of the original adaptor)
==========================================
|| -f <fasta_file>: name of fa file, with deep seq data in fa format
||-A <adaptor>: Adaptor seq; default is Illumina TGGAATTCTCGGGTGCCAAGG
||-M <mode>:   1 or 2; default 2; mode 1 uses miRDeep2 algorithm, with specific cutoffs with certain number of bases for the
||                      adaptor used in match (-B), and restrict the length of seq to have a minimum (-L) in first pass match
||                      mode 2 uses HiSeq friendly match, see above
||-B <bases>:	number of bases of adaptor used, default 10
||-L <minLength>:	number of bases of small RNA length, default 18, see -M description
||-u <use_2nd_adaptor>: default 0; set to 1 if want to use
||-s <second_adaptor>: sequence of second adatpor; default: TGGACTTCTCGGGTGCCAAGG
||-o <outfile>: Output file name
==========================================
";

#default values for parameters
$adaptor="TGGAATTCTCGGGTGCCAAGG";
$mode=2;
$bases=15;
$minLength=18;
$second_adaptor="TGGACTTCTCGGGTGCCAAGG";
$use_2nd_adaptor=0;

for($i=0;$i<@ARGV;$i+=2)
{
	if($ARGV[$i] eq "-A")
	{
		$adaptor=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-f")
	{
		$fasta_file=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-M")
	{
		$mode=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-B")
	{
		$bases=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-L")
	{
		$minLength=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-u")
	{
		$use_2nd_adaptor=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-s")
	{
		$second_adaptor=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-o")
	{
	    $outfile=$ARGV[$i+1];
	}
	else
	{
		
	}
}

if(!$fasta_file||!$outfile) {die "$Usage\n";}



my $prefix=substr($adaptor,0,$bases);
$prefix=~tr/[acgtun\.]/[ACGTTNN]/; 

my $prefix2=substr($second_adaptor,0,$bases);
$prefix2=~tr/[acgtun\.]/[ACGTTNN]/; 

remove_adapters($fasta_file,$prefix,$minLength, $mode,$outfile,$use_2nd_adaptor,$prefix2);


exit;


sub remove_adapters{

    my ($file_fasta,$prefix,$minLength,$mode,$outfile,$use_2nd_adaptor,$prefix2) = @_;
    my ($id,$seq) = ();

    open($OUTF,">$outfile") or die "cannot open $outfile\n";


    open (FASTA, "<$fasta_file") or die "can not open $fasta_file\n";
    while (<FASTA>)
    {
        chomp;
        if (/^>(\S+)/)
        {
            $id  = $1;
	    $seq = "";
            while (<FASTA>){
                chomp;
                if (/^>(\S+)/){
		    remove_adapter($id,$seq,$prefix,$minLength,$mode,$OUTF,$use_2nd_adaptor,$prefix2);
                    $id    = $1;
                    $seq   = "";
                    next;
                }
                $seq .= $_;
            }
        }
    }
    remove_adapter($id,$seq,$prefix,$minLength,$mode,$OUTF,$use_2nd_adaptor,$prefix2);
    close FASTA;
    close $OUTF;
    return;
}





sub remove_adapter{

    my($id,$seq,$prefix,$minLength,$mode,$OUTF,$use_2nd_adaptor,$prefix2)=@_;
    $seq=~tr/[acgtun\.]/[ACGTTNN]/;   
    my $seq_clipped;

    if($mode==1)
    {
	
	if($seq=~/(\w{$minLength,})$prefix/){
	    
	    $seq_clipped=$1;
	    
	}else{
	    
	    my $finish=0;
	    
	    while(not $finish and (length($prefix)>0)){
		
		chop $prefix;
		
		if($seq=~/(\w+)$prefix$/){
		    
		    $seq_clipped=$1;
		    $finish=1;
		}
	    }
	    
	}
	
	if(not $seq_clipped){
	    
	    if($use_2nd_adaptor && $seq=~/(\w{$minLength,})$prefix2/){
		
		$seq_clipped=$1;
		
	    }else{
		
		$seq_clipped=$seq;
	    }
	}
	
	print $OUTF ">$id\n$seq_clipped\n";
    }
    elsif($mode==2)
    {
	#note: uses non-greedy match below
	if($seq=~/(\w+?)$prefix/){
	    
	    $seq_clipped=$1;
	    
	}
	else{
	    
	    my $finish=0;
	    
	    while(not $finish and (length($prefix)>=6)){
		
		chop $prefix;
		#note, again uses greedy match, after the first $minLength bases
		if($seq=~/(\w{$minLength,})$prefix/){
		    
		    $seq_clipped=$1;
		    $finish=1;
		}
	    }
	    
	}
	
	if(not $seq_clipped){
	    if($use_2nd_adaptor && $seq=~/(\w+?)$prefix2/){
		
		$seq_clipped=$1;
		
	    }else{
		
		$seq_clipped=$seq;
	    }
	}
	print $OUTF ">$id\n$seq_clipped\n";    
    }
}

