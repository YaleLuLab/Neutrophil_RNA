#!/usr/bin/perl -w

use lib "$ENV{\"HOME\"}/Perl_Codes/Jun_lib/";
use READ_FILE;
#Add other packages here with syntax: use package;

$Usage="==========================================
|| Jun Lu, July 2022, tested
|| v13 is modified from v11, to be compatible with major changes in the v13 of map_to_known program.
|| v13 also adds support for outputting a length summary for mapped reads
|| This code is revised from v1, to be compatible with v6 of map to known. Also added a parameter to input folder surfix
||
|| The v1 code summarizes mapping results to known RNA species
|| It takes directories of mapping for multiple samples and outputs two summary files
|| First summary file summarizes mapping read counts to different classes
|| Second summarizes isoform information
||
|| Usage: perl $0 -d <InDir> -s <SummaryFile> -o <IsoformFile> -f <ThreshReads> -b <BowtieV> -p <Species> -i <sum_isoform> -fs <FolderSurfix> -m <sumLen> -a <Assembly>
==========================================
||-d <InDir>:	Name of directory; this option can be used mulitple times to input multiple samples
||-s <SummaryFile>:	Readcount summary file name
||-o <IsoformFile>:	Isoform summary file name
||-f <ThreshReads>:	Threshold level; sequences with fewer than this for total reads across all samples will be elimited in isoform reporting; default is 10
||-b <BowtieV>:         Defines bowtie version to work on in case there are both existing; 1 for bowtie and 2 for bowtie2; default 2
||-p <Species>:         Optional, use hsa for human and mmu for mouse; default hsa
||-i <sum_isoform>:     Optional, set to 1 to perform isoform summary; set to 0 to avoid; default 1
||-fs <FolderSurfix>:   Optional, defines the string for FolderSurfix; default is _Map_to_Known_v13
||-m <sumLen>:          Optional, set to 0 to disable; set to file name to define the output for length summary of mapped reads; default 0; this task needs sum_isoform to be set to 1 to run. 
||-a <Assembly>:        Optional, goes with sumLen--if summing length, Assembly refers to the genome that was aligned against.
==========================================
";

#default values for parameters
@InDir=();
$SummaryFile=0;
$IsoformFile=0;
$ThreshReads=10;
$BowtieV=2;
$Species="hsa";
$sum_isoform=1;  #make sure to change back to 1 
$FolderSurfix="_Map_to_Known_v13";
$sumLen=0;
$Assembly=0;

for($i=0;$i<@ARGV;$i+=2)
{
	if($ARGV[$i] eq "-d")
	{
		push @InDir, $ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-s")
	{
		$SummaryFile=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-o")
	{
		$IsoformFile=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-f")
	{
		$ThreshReads=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-b")
	{
		$BowtieV=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-p")
	{
		$Species=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-i")
	{
		$sum_isoform=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-m")
	{
		$sumLen=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-fs")
	{
		$FolderSurfix=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-a")
	{
		$Assembly=$ARGV[$i+1];
	}
	else
	{
		
	}
}
if(!@InDir || !$SummaryFile || !$IsoformFile) {die "$Usage\n";}
#============= start code here ============================

#---Setup parameters

$MapSummarySurfix="_Mapping_Summary_bt".$BowtieV.".txt";

#@MappingTypes=("AdaptorSpikeIn","ncRNAs","Genome");
@TypeSurfix=(qr/_ncRbt\d\.sam/, qr/_GENOME.+?bt\d\.sam/);
$genomeIdx=1; #defines that the index for genome in @TypeSurfix

@LenSumRange=(1..151);
$InitialReads_Surfix="_IaSbt".$BowtieV.".fa";
$UnMappedReads_Surfix="_IaSbt".$BowtieV."_ncRbt".$BowtieV."_GENOME".$Assembly."bt".$BowtieV.".fa";


#---first verify directory names and get the sample names
@Samples=();
@ParentFolders=();

for($i=0;$i<@InDir;$i++)
{
    if(! -d $InDir[$i])
    {
	print STDERR "Known Mapping Summary: Error: $InDir[$i] does not exist\n";
    }

    $tempFolder=$InDir[$i];
    $tempFolder=~s/\/$//; # remove last splash if existing
    $InDir[$i]=$tempFolder.'/';
    if($tempFolder=~/^(.+)\/(.+)$FolderSurfix/)
    {
	$ParentFolders[$i]=$1.'/';
	$Samples[$i]=$2;
    }
    else
    {
	$tempFolder=~/^(.+)$FolderSurfix/;
	$ParentFolders[$i]='./';
	$Samples[$i]=$1;
    }
}

#------Next summary mapping read counts for different classes

%ReadCat=();

for($i=0;$i<@Samples;$i++)
{
    $curFile=$InDir[$i].$Samples[$i].$MapSummarySurfix;
    open($infile,$curFile) || die "Couldn't open $curFile\n";

    @wholefile=<$infile>;
    close($infile);
    chomp(@wholefile);
    shift(@wholefile);
    for($j=0;$j<@wholefile;$j++)
    {
	@temparray=split("\t",$wholefile[$j]);
	if(!exists($ReadCat{$temparray[0]})) {
	    my @cat_counts=(0)x@Samples; #creates a zero count array for each sample
	    $ReadCat{$temparray[0]}=\@cat_counts;
	    $ReadCat{$temparray[0]}->[$i]=$temparray[3]
	}
	else {
	    $ReadCat{$temparray[0]}->[$i]=$temparray[3];
	}
    }
}

open($outfile,">$SummaryFile") || die "couldn't open $SummaryFile\n";
print $outfile "Step";
for($i=0;$i<@Samples;$i++)
{
    print $outfile "\t",$Samples[$i];
}
print $outfile "\n";

foreach $j (keys %ReadCat)
{
    for($i=0;$i<@Samples;$i++)
    {
	if($i==0)
	{
	    print $outfile $j;
	}
	print $outfile "\t",$ReadCat{$j}->[$i];
    }
    print $outfile "\n";
}

close($outfile);

print STDERR "\tFinished summarizing read counts\n";

########Summarize Mapped Read Length

if($sumLen)
{
    print STDERR "\tStart summarizing length of mapped reads\n";

    #set up 0 array for tallying
    @LenSumData=();
    for($i=0;$i<@LenSumRange;$i++)
    {
	my @tempsams=(0)x@Samples;
	$LenSumData[$i]=\@tempsams;
    }	    

    #now start to read each samples' fa files to determine which ones are mapped

    for($i=0;$i<@Samples;$i++)
    {
	print STDERR "\t\tStart summarizing sample $i, $Samples[$i]\n";
	my %UnMappedReads=();
	my $UnMappedSeqs=NA_SEQ::read_FASTA_file($InDir[$i].$Samples[$i].$UnMappedReads_Surfix);
	my $InitialSeqs=NA_SEQ::read_FASTA_file($InDir[$i].$Samples[$i].$InitialReads_Surfix);
	
	for($j=0;$j<@{$UnMappedSeqs};$j++)
	{
	    if(!exists($UnMappedReads{$UnMappedSeqs->[$j]->name()})) {
		$UnMappedReads{$UnMappedSeqs->[$j]->name()}=1;
	    }
	}

	for($j=0;$j<@{$InitialSeqs};$j++)
	{
	    if(!exists($UnMappedReads{$InitialSeqs->[$j]->name()})) {
		my $tempSeqName=$InitialSeqs->[$j]->name();
		if($tempSeqName=~/.+_x(\d+)$/) {
		    my $tempSeqLength=$InitialSeqs->[$j]->length();
		    $LenSumData[$tempSeqLength-1][$i]+=$1;
		}
	    }
	}
    }

    #output the results
    print STDERR "\tStart producing output of length summary\n";

    open(OUT,">$sumLen") || die "couldn't open length summary file $sumLen\n";
    print OUT "Length";
    for($i=0;$i<@Samples;$i++)
    {
	print OUT "\t",$Samples[$i];
    }
    print OUT "\n";
        
    for($j=0;$j<@LenSumRange;$j++)
    {
	print OUT $LenSumRange[$j];
	print OUT "\t",join("\t",@{$LenSumData[$j]}),"\n";
    }
    close OUT;
}






##########################------Next summarize isoforms
if(!$sum_isoform) {exit;}

%SeqTarget=(); #This hash stores combined Seqeunce and Target as a unique identifier
%UniqueSeqs=(); #This has stores unique sequences to check whether or not there are multiple matches

@curAlignFile=(); #this stores all file names to be used

#find all files with alignment
for($i=0;$i<@Samples;$i++)
{
    opendir($dh, $InDir[$i]) || die "couldn't open dir $InDir[$i]\n";
    @curdirfiles=readdir($dh);
    closedir $dh;

    for($j=0;$j<@TypeSurfix;$j++)
    {
	@curRNAfiles=grep {/$TypeSurfix[$j]/} @curdirfiles;
	if(!@curRNAfiles)
	{
	    die "cannot find file for $Samples[$i]: $RNAtypes[$j]\n";
	}
	elsif(@curRNAfiles>1)
	{
	    die "found more than one file for $Samples[$i]: $RNAtypes[$j]:",join(';',@curRNAfiles),"\n";
	}
	else
	{
	    $curAlignFile[$i][$j]=$curRNAfiles[0];
	}
    }
}

#start reading alignment files
@OutputData=();
@OutTable=();
@TotalCount=();
@OutputSeqs=();
$counter=0;

for($j=0;$j<@TypeSurfix;$j++)
{
    for($i=0;$i<@Samples;$i++)
    {
	print STDERR "\tprocessing sample $i, file $curAlignFile[$i][$j]\n";

	my $curSamData=READ_FILE::read_bowtie2_sam_alignments($InDir[$i].$curAlignFile[$i][$j], 1);

	for($k=0;$k<@{$curSamData};$k++)
	{
	    if(!$curSamData->[$k]->{mapped}) {next;}
	    
	    if($k % 20000==0) {print STDERR "\t\tprocessing entry line $k\n";}

	    if($curSamData->[$k]->{reverse})
	    {
		$curSamData->[$k]->{seq}=reverse($curSamData->[$k]->{seq});
		$curSamData->[$k]->{seq}=~tr/ATCGN/TAGCN/;
	    }
	    
	    $tempSeqTar=$curSamData->[$k]->{seq}."|".$curSamData->[$k]->{rname}."|".$curSamData->[$k]->{pos};
	    $curSamData->[$k]->{qname}=~/.+_x(\d+)$/;
	    $tempCount=$1;
	    if($j==$genomeIdx) {
		$tempRNAtype="Genome";
	    } elsif ($curSamData->[$k]->{rname}=~/^(.+?)__/) {
		$tempRNAtype=$1;
	    } else {
		$tempRNAtype="NA";
	    }
		
	    if(!exists($SeqTarget{$tempSeqTar}))
	    {
		my @tempCountArray=(0)x@Samples;
		$tempCountArray[$i]=$tempCount;

		$OutputData[$counter]=$curSamData->[$k]->{rname};
		$OutputData[$counter].="\t".$tempRNAtype;
		$OutputData[$counter].="\t".$curSamData->[$k]->{seq};
		$OutputData[$counter].="\t".$curSamData->[$k]->{pos};
		$OutputData[$counter].="\t".($curSamData->[$k]->{pos}+length($curSamData->[$k]->{seq})-1);
		$OutputData[$counter].="\t".$curSamData->[$k]->{reverse};
		$OutputData[$counter].="\t".length($curSamData->[$k]->{seq});

		$OutTable[$counter]=\@tempCountArray;
		$OutputSeqs[$counter]=$curSamData->[$k]->{seq};
		
		#debug
		#if($counter<5)
		#{    #print STDERR join("\t",@tempCountArray),"\n";
		#    print STDERR $OutputData[$counter],"\n";
		#    print STDERR join("\t",@{$OutTable[$counter]}),"\n";
		#}
		#end debug

		$TotalCount[$counter]=$tempCount;
		
		$SeqTarget{$tempSeqTar}=$counter;

		$counter++;

		if(!exists($UniqueSeqs{$curSamData->[$k]->{seq}}))
		{
		    $UniqueSeqs{$curSamData->[$k]->{seq}}=1;
		}		
		else
		{
		    $UniqueSeqs{$curSamData->[$k]->{seq}}+=1;
		}
	    
		#if($sumLen) #this way is count length summary is not right--the seq from sam file may only reflect the matched seq, rather than raw read seq.
		#{
		#    if(!exists($SamUniqueSeqs[$i]->{$curSamData->[$k]->{seq}})) {
		#	$SamUniqueSeqs[$i]->{$curSamData->[$k]->{seq}}=1;
		#	$LenSumData[length($curSamData->[$k]->{seq})-1][$i]+=$tempCount;
		#    }
		#}
	    }
	    else
	    {
		if((!($curSamData->[$k]->{flag} & 256) || $OutTable[$SeqTarget{$tempSeqTar}]->[$i]==0) && !($curSamData->[$k]->{flag} & 512) && !($curSamData->[$k]->{flag} & 1024) ) # avoid counting secondary alignments etc.
		#if( !($curSamData->[$k]->{flag} & 512) && !($curSamData->[$k]->{flag} & 1024) )
		{
		    $OutTable[$SeqTarget{$tempSeqTar}]->[$i]+=$tempCount;
		    $TotalCount[$SeqTarget{$tempSeqTar}]+=$tempCount;
		}
	    }
	}
    }
}

#---output
print STDERR "\tStart producing output of isoform summary\n";

open(OUT,">$IsoformFile") || die "couldn't open $IsoformFile\n";
print OUT "RNAname\tRNAtype\tReadSeq\tStart\tEnd\tReverseStrand\tReadLength";
for($i=0;$i<@Samples;$i++)
{
    print OUT "\t",$Samples[$i];
}
print OUT "\tTotalCount\tNumOfMatch\n";

for($j=0;$j<@OutputData;$j++)
{
    if($TotalCount[$j]>=$ThreshReads)
    {
	print OUT $OutputData[$j];
	
	print OUT "\t",join("\t",@{$OutTable[$j]});

	print OUT "\t",$TotalCount[$j];
	print OUT "\t",$UniqueSeqs{$OutputSeqs[$j]},"\n";
    }
}

close OUT;



