#!/usr/bin/perl -w

use lib "$ENV{\"HOME\"}/Perl_Codes/Jun_lib/";
use READ_FILE;
use NA_SEQ;

$Usage="==========================================
|| Jun Lu, Feb 2017; mod July 2017
|| Usage: perl $0 -p <precursorFA> -m <matureArf> [-f <samArfs>] -o <outfilename>
||
|| This program summarizes the mapping result from miRDeep2 in terms of isoforms
|| In the output file, the field of BestMisMatch was previously defined as: if a sequence is a mismatch, this is the best mismatch among all mismatches, even there is perfect match for this sequence
|| To avoid confusion and increase practical usability, the field BestMisMatch is now replaced with two fields: IsBestMatch and NumBestMatch
||    IsBestMatch tells if a match is the best match for this sequence. If the sequence is mapped equally well to two or more different reference sequences, all will be recorded as 1
||    NumBestMatch tells the total number of best matches for a given sequence 
==========================================
||-p <precursorFA>:     miR precursor fa file
||-m <matureArf>:	mature miRNA mapping result in Arf format
||-f <samArfs>:	repeatable; sample mapping result in Arf format
||-o <outfilename>: output file name
==========================================
";

#default values for parameters
$precursorFA=0;
$matureArf=0;
@samArfs=();
$outfilename=0;

for($i=0;$i<@ARGV;$i+=2)
{
	if($ARGV[$i] eq "-m")
	{
		$matureArf=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-p")
	{
		$precursorFA=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-f")
	{
		push @samArfs,$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-o")
	{
	    $outfilename=$ARGV[$i+1];
	}
	else
	{
		
	}
}
if(!$precursorFA || !$matureArf || !@samArfs ||!$outfilename) {die "$Usage\n";}
#============= start code here ============================
%Precursor=();
$NoMatchLoopLen=-70;
$wiggleroom=5;


#--first read mature miRNA mapping to Pre-miRNA and annotate pre-miRNA

#before doing anything else, read in precursor FA file to get the length of premiR

$premiRs=NA_SEQ::read_FASTA_file($precursorFA);
%premiRLen=();
for($i=0;$i<@{$premiRs};$i++)
{
    if(!exists($premiRLen{$premiRs->[$i]->name()}))
    {
	$premiRLen{$premiRs->[$i]->name()}=$premiRs->[$i]->length();
    }
}

#then read mature mapping file

$miRMapping=READ_FILE::read_arf_file($matureArf);

for($i=0;$i<@{$miRMapping};$i++)
{
    $matureMapping=$miRMapping->[$i];
    if($matureMapping->{MisMatch}>0) {next;}
    
    if(!exists($Precursor{$matureMapping->{RefName}}))
    {
	if($matureMapping->{RefStart}<($premiRLen{$matureMapping->{RefName}}/2-5))
	{
	    $Precursor{$matureMapping->{RefName}}={
		Name5p=>$matureMapping->{ReadName},
		Len5p=>$matureMapping->{ReadLen},
		Start5p=>$matureMapping->{RefStart},
		End5p=>$matureMapping->{RefEnd},
	    }
	}
        else
	{
	    $Precursor{$matureMapping->{RefName}}={
		Name3p=>$matureMapping->{ReadName},
		Len3p=>$matureMapping->{ReadLen},
		Start3p=>$matureMapping->{RefStart},
		End3p=>$matureMapping->{RefEnd},
	    }
	}
    }
    else
    {
	if($matureMapping->{RefStart}<($premiRLen{$matureMapping->{RefName}}/2-5))
	{
	    if(!exists($Precursor{$matureMapping->{RefName}}->{Name5p}))
	    {
		$Precursor{$matureMapping->{RefName}}->{Name5p}=$matureMapping->{ReadName};
		$Precursor{$matureMapping->{RefName}}->{Len5p}=$matureMapping->{ReadLen};
		$Precursor{$matureMapping->{RefName}}->{Start5p}=$matureMapping->{RefStart};
		$Precursor{$matureMapping->{RefName}}->{End5p}=$matureMapping->{RefEnd};
	    }
	    else
	    {
		print STDERR "ERROR: Found another instance of 5p miRNA, ",$matureMapping->{RefName},", ",$matureMapping->{ReadName},", ",$matureMapping->{RefStart},"\n";
	    }
	}
	else
	{
	    if(!exists($Precursor{$matureMapping->{RefName}}->{Name3p}))
	    {
		$Precursor{$matureMapping->{RefName}}->{Name3p}=$matureMapping->{ReadName};
		$Precursor{$matureMapping->{RefName}}->{Len3p}=$matureMapping->{ReadLen};
		$Precursor{$matureMapping->{RefName}}->{Start3p}=$matureMapping->{RefStart};
		$Precursor{$matureMapping->{RefName}}->{End3p}=$matureMapping->{RefEnd};
	    }
	    else
	    {
		print STDERR "ERROR: Found another instance of 3p miRNA, ",$matureMapping->{RefName},", ",$matureMapping->{ReadName},", ",$matureMapping->{RefStart},", ", $Precursor{$matureMapping->{RefName}}->{Name3p},"\n";
	    }
	}
    }
}

@AllPrecursors=keys %Precursor;

for($i=0;$i<@AllPrecursors;$i++)
{
    if(exists($Precursor{$AllPrecursors[$i]}->{Name5p}) && exists($Precursor{$AllPrecursors[$i]}->{Name3p}))
    {
	$Precursor{$AllPrecursors[$i]}->{LoopLen}=$Precursor{$AllPrecursors[$i]}->{Start3p}-$Precursor{$AllPrecursors[$i]}->{End5p}-1;
    }
    elsif(exists($Precursor{$AllPrecursors[$i]}->{Name5p}))
    {
	$Precursor{$AllPrecursors[$i]}->{LoopLen}=$NoMatchLoopLen;
	$Precursor{$AllPrecursors[$i]}->{Name3p}=0;
	$Precursor{$AllPrecursors[$i]}->{Len3p}=0;
	$Precursor{$AllPrecursors[$i]}->{Start3p}=0;
	$Precursor{$AllPrecursors[$i]}->{End3p}=0;
    }
    elsif(exists($Precursor{$AllPrecursors[$i]}->{Name3p}))
    {
	$Precursor{$AllPrecursors[$i]}->{LoopLen}=$NoMatchLoopLen;
	$Precursor{$AllPrecursors[$i]}->{Name5p}=0;
	$Precursor{$AllPrecursors[$i]}->{Len5p}=0;
	$Precursor{$AllPrecursors[$i]}->{Start5p}=0;
	$Precursor{$AllPrecursors[$i]}->{End5p}=0;
    }
}

print "Finished mapping mature miRNAs onto Precursor\n";

#----Now go through reads mapping files

%Reads=(); #this hash stores all possible precursors that the reads map to; whereas %Precursor stores other info for reads that map to a specific precursor

#-first process the sample names
@SamNames=();
for($i=0;$i<@samArfs;$i++)
{
    $samArfs[$i]=~/([^\/]+)\.arf/;
    $SamNames[$i]=$1;
}

#-second read in the arf files
for($i=0;$i<@samArfs;$i++)
{
    $curArf=READ_FILE::read_arf_file($samArfs[$i]);
    for($j=0;$j<@{$curArf};$j++)
    {
	$curArf->[$j]->{ReadName}=~/_x(\d+)$/;
	$curReadNum=$1;
	$curSeq=$curArf->[$j]->{ReadSeq};
	$curPrecur=$curArf->[$j]->{RefName};

	if(!exists($Reads{$curSeq}))
	{
	    $Reads{$curSeq}={$curPrecur=>$curArf->[$j]->{MisMatch},
	    };
	}
	elsif(!exists($Reads{$curSeq}->{$curPrecur}))
	{
	    $Reads{$curSeq}->{$curPrecur}=$curArf->[$j]->{MisMatch};
	}

	if(!exists($Precursor{$curPrecur}->{Seqs}))
	{
	    my @curCounts=();
	    for($k=0;$k<@samArfs;$k++)
	    {
		$curCounts[$k]=0;
	    }
	    $curCounts[$i]=$curReadNum;


	    $Precursor{$curPrecur}->{Seqs}={$curSeq=>{
		MappedSeq=>$curArf->[$j]->{MappedSeq},
		ReadLen=>$curArf->[$j]->{ReadLen},
		ReadStart=>$curArf->[$j]->{ReadStart},
		ReadEnd=>$curArf->[$j]->{ReadEnd},
		RefStart=>$curArf->[$j]->{RefStart},
		RefEnd=>$curArf->[$j]->{RefEnd},
		MisMatch=>$curArf->[$j]->{MisMatch},
		MatchString=>$curArf->[$j]->{MatchString},
		Counts=>\@curCounts,
					    },
	    };
	    
	    if($curArf->[$j]->{RefEnd}<=($Precursor{$curPrecur}->{Start5p}+$wiggleroom))
	    {
		$Precursor{$curPrecur}->{Seqs}->{$curSeq}->{ReadRegion}="5p_flank";
	    }
	    elsif($curArf->[$j]->{RefStart}>=($Precursor{$curPrecur}->{Start5p}-$wiggleroom) &&
		  $curArf->[$j]->{RefEnd}<=($Precursor{$curPrecur}->{End5p}+$wiggleroom))
	    {
		$Precursor{$curPrecur}->{Seqs}->{$curSeq}->{ReadRegion}="5p";
	    }
	    elsif($curArf->[$j]->{RefStart}>=($Precursor{$curPrecur}->{Start3p}-$wiggleroom) &&
		  $curArf->[$j]->{RefEnd}<=($Precursor{$curPrecur}->{End3p}+$wiggleroom))
	    {
		$Precursor{$curPrecur}->{Seqs}->{$curSeq}->{ReadRegion}="3p";
	    }
	    elsif($curArf->[$j]->{RefStart}>=($Precursor{$curPrecur}->{End5p}-$wiggleroom) &&
		  $curArf->[$j]->{RefEnd}<=($Precursor{$curPrecur}->{Start3p}+$wiggleroom) &&
		  $Precursor{$curPrecur}->{End5p}>0)
	    {
		$Precursor{$curPrecur}->{Seqs}->{$curSeq}->{ReadRegion}="loop";
	    }
	    elsif($curArf->[$j]->{RefStart}>=($Precursor{$curPrecur}->{End3p}-$wiggleroom) &&
		  $Precursor{$curPrecur}->{End3p}>0)
	    {
		$Precursor{$curPrecur}->{Seqs}->{$curSeq}->{ReadRegion}="3p_flank";
	    }
	    else
	    {
		$Precursor{$curPrecur}->{Seqs}->{$curSeq}->{ReadRegion}="Other";
	    }
	}
	elsif(!exists($Precursor{$curPrecur}->{Seqs}->{$curSeq}))
	{
	    my @curCounts=();
	    for($k=0;$k<@samArfs;$k++)
	    {
		$curCounts[$k]=0;
	    }
	    $curCounts[$i]=$curReadNum;

	    $Precursor{$curPrecur}->{Seqs}->{$curSeq}={
		MappedSeq=>$curArf->[$j]->{MappedSeq},
		ReadLen=>$curArf->[$j]->{ReadLen},
		ReadStart=>$curArf->[$j]->{ReadStart},
		ReadEnd=>$curArf->[$j]->{ReadEnd},
		RefStart=>$curArf->[$j]->{RefStart},
		RefEnd=>$curArf->[$j]->{RefEnd},
		MisMatch=>$curArf->[$j]->{MisMatch},
		MatchString=>$curArf->[$j]->{MatchString},
		Counts=>\@curCounts,
	    };

	    if($curArf->[$j]->{RefEnd}<=($Precursor{$curPrecur}->{Start5p}+$wiggleroom))
	    {
		$Precursor{$curPrecur}->{Seqs}->{$curSeq}->{ReadRegion}="5p_flank";
	    }
	    elsif($curArf->[$j]->{RefStart}>=($Precursor{$curPrecur}->{Start5p}-$wiggleroom) &&
		  $curArf->[$j]->{RefEnd}<=($Precursor{$curPrecur}->{End5p}+$wiggleroom))
	    {
		$Precursor{$curPrecur}->{Seqs}->{$curSeq}->{ReadRegion}="5p";
	    }
	    elsif($curArf->[$j]->{RefStart}>=($Precursor{$curPrecur}->{Start3p}-$wiggleroom) &&
		  $curArf->[$j]->{RefEnd}<=($Precursor{$curPrecur}->{End3p}+$wiggleroom))
	    {
		$Precursor{$curPrecur}->{Seqs}->{$curSeq}->{ReadRegion}="3p";
	    }
	    elsif($curArf->[$j]->{RefStart}>=($Precursor{$curPrecur}->{End5p}-$wiggleroom) &&
		  $curArf->[$j]->{RefEnd}<=($Precursor{$curPrecur}->{Start3p}+$wiggleroom) &&
		  $Precursor{$curPrecur}->{End5p}>0)
	    {
		$Precursor{$curPrecur}->{Seqs}->{$curSeq}->{ReadRegion}="loop";
	    }
	    elsif($curArf->[$j]->{RefStart}>=($Precursor{$curPrecur}->{End3p}-$wiggleroom) &&
		  $Precursor{$curPrecur}->{End3p}>0)
	    {
		$Precursor{$curPrecur}->{Seqs}->{$curSeq}->{ReadRegion}="3p_flank";
	    }
	    else
	    {
		$Precursor{$curPrecur}->{Seqs}->{$curSeq}->{ReadRegion}="Other";
	    }
	}
	else
	{
	    if($Precursor{$curPrecur}->{Seqs}->{$curSeq}->{RefStart}!=$curArf->[$j]->{RefStart})
	    {
		#print "WARNING: multiple mapping on the same precursor detected: $curPrecur, $curSeq\n";
	    }
	    elsif($Precursor{$curPrecur}->{Seqs}->{$curSeq}->{Counts}->[$i]>0)
	    {
		print "WARNING: same seq mapped twice to the same precursor detected: $curPrecur, $curSeq\n";
	    }
	    else
	    {
		$Precursor{$curPrecur}->{Seqs}->{$curSeq}->{Counts}->[$i]+=$curReadNum;
	    }
	}
    }
}

print "Ready to produce output\n";


#======print out output===================

open($outfile,">$outfilename") || die "couldn't open output file $outfilename\n";

print $outfile "Hairpin\tmiR5p\t5pLen\t5pStart\t5pEnd\tmiR3p\t3pLen\t3pStart\t3pEnd\tLoopLen\tReadSeq\tReadLen\tReadStart\tReadEnd\tMappedSeq\tPrecurStart\tPrecurEnd\tReadRegion\tMisMatches\tMatchCode\tMultipleMatch\tIsBestMatch\tNumBestMatch";
for($i=0;$i<@samArfs;$i++)
{
    print $outfile "\t",$SamNames[$i];
}
print $outfile "\n";

for($i=0;$i<@AllPrecursors;$i++)
{
    $curPrecur=$AllPrecursors[$i];
    #debug
    #print $curPrecur,"\n";
    #end debug

    @curMapped=keys %{$Precursor{$curPrecur}->{Seqs}};
    #debug
    #print join("\n",@curMapped),"\n";
    #end debug

    if(!@curMapped) {next;}
    
    for($j=0;$j<@curMapped;$j++)
    {
	print $outfile $curPrecur,"\t",$Precursor{$curPrecur}->{Name5p},"\t",$Precursor{$curPrecur}->{Len5p},"\t",$Precursor{$curPrecur}->{Start5p},"\t",$Precursor{$curPrecur}->{End5p},"\t",$Precursor{$curPrecur}->{Name3p},"\t",$Precursor{$curPrecur}->{Len3p},"\t",$Precursor{$curPrecur}->{Start3p},"\t",$Precursor{$curPrecur}->{End3p},"\t",$Precursor{$curPrecur}->{LoopLen};

	print $outfile "\t",$curMapped[$j],"\t",$Precursor{$curPrecur}->{Seqs}->{$curMapped[$j]}->{ReadLen},"\t",$Precursor{$curPrecur}->{Seqs}->{$curMapped[$j]}->{ReadStart},"\t",$Precursor{$curPrecur}->{Seqs}->{$curMapped[$j]}->{ReadEnd},"\t",$Precursor{$curPrecur}->{Seqs}->{$curMapped[$j]}->{MappedSeq},"\t",$Precursor{$curPrecur}->{Seqs}->{$curMapped[$j]}->{RefStart},"\t",$Precursor{$curPrecur}->{Seqs}->{$curMapped[$j]}->{RefEnd},"\t",$Precursor{$curPrecur}->{Seqs}->{$curMapped[$j]}->{ReadRegion},"\t",$Precursor{$curPrecur}->{Seqs}->{$curMapped[$j]}->{MisMatch},"\t",$Precursor{$curPrecur}->{Seqs}->{$curMapped[$j]}->{MatchString};

	@tempMappedPre=keys %{$Reads{$curMapped[$j]}};
	$tempMultiMatch=@tempMappedPre;
	print $outfile "\t",$tempMultiMatch;

	@tempBestIdx=(0);
	$tempMinMisMatch=200;
	for($k=0;$k<@tempMappedPre;$k++)
	{
	    if($Reads{$curMapped[$j]}->{$tempMappedPre[$k]}<$tempMinMisMatch)
	    {
		$tempMinMisMatch=$Reads{$curMapped[$j]}->{$tempMappedPre[$k]};
		@tempBestIdx=($k);
	    }
	    elsif($Reads{$curMapped[$j]}->{$tempMappedPre[$k]}==$tempMinMisMatch && $k>0)
	    {
		push @tempBestIdx,$k;
	    }
	}
	
	#if($tempMultiMatch==1)
	#{
	#   print $outfile "\t0"; 
	#}
	$tempBestIsMe=0;
	for ($k=0;$k<@tempBestIdx;$k++)
	{
	    if($tempMappedPre[$tempBestIdx[$k]] eq $curPrecur)
	    {
		$tempBestIsMe=1;
	    }
	}
	
	if($tempBestIsMe)
	{
	    print $outfile "\t1";
	}
	else
	{
	    print $outfile "\t0";
	}
	    
	$tempNumBestMatches=@tempBestIdx;
	print $outfile "\t$tempNumBestMatches";


	for($k=0;$k<@samArfs;$k++)
	{
	    print $outfile "\t",$Precursor{$curPrecur}->{Seqs}->{$curMapped[$j]}->{Counts}->[$k];
	}
	print $outfile "\n";
    }
}
