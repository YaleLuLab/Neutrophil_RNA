#!/usr/bin/perl -w

use lib "$ENV{\"HOME\"}/Perl_Codes/Jun_lib/";
#Add other packages here with syntax: use package;

$Usage="==========================================
|| Jun Lu, July 2022
|| V13 from v12--modified to use the map-to-known of v13. 
|| This v12 is modified from v11, to add the ability to define two or more data files per sample. 
||
|| This v11 is modified from v10, to allow compatilbity on HPC
|| It establishes a standard folder structure in the project folder, and also creates run log for tracking the procedures
|| The following folders will be created in the project folder, if not already there: RawSeq (for storing raw fastq and fasta files), ProcessedSeq (intermediate files), Quantified (stores quantification results)
||
|| Sample Definition File (v8): txt file with a single header line
||      Column 1: ProjectName; this will be used for certain outputs
||      Column 2: ProjectFolder; please add path, ideally from root folder
||      Column 3: DataFileName (without surfix) or SRR ID, more than one datafile can be provided, separated by comma (no space after comma)
||      Column 4: ConvertedFileName; these names will be used after fa conversion instead of the original DataFileName
||      Column 5: Adaptor sequence to remove; if standard (TGGAATTCTCGGGTGCCAAGG), please specify \'standard\'; otherwise provide sequence; Note that CLIP data often use another 3'adaptor (AGATCGGAAGAGCGGTTC... check irCLIP paper)
||      Column 6: Species; please specify human (or hsa) or mouse (or mmu); case insensitive
||      Column 7: miRBase version; specify v18,v18_FPCM,v20,v21; if custom files, please specify path and prefix of files _hairpin.fa and _mature.fa will be appended to the prefix to find the files.
||      Column 8: Length Cutoff; number used for filtering reads; normally we use 15 or 16 as a cutoff
||      Column 9: Mismatch; tolerated mismatch during mapping to pre-miRNAs databases. (The info in this column is not used in v10 pipeline).
||      Column 10: Genome Assembly, for mapping to genome; choose between hg18, hg19, hg38, mm9, mm10 and mm39 (mm10 and mm39 only use bowtie 2)
||      Column 11: Mapping algorithm; 1 for bowtie, 2 for bowtie 2; use 2 recommended
||
|| STEP 1: download and/or convert fastq to fasta
|| STEP 2: clip off adatpor sequence
|| STEP 3: collapse sequences to unique ones
|| STEP 4: summarize seq length distribution
|| STEP 5: filter out short reads
|| STEP 6: quantify the reads using miRDeep2
|| STEP 7: summarize miRNA quantification results
|| STEP 8: quantify miRNA isoform information
|| STEP 9: map to known RNAs or the genome
|| STEP 10: quantify other RNA isoforms
|| STEP 11: map directly to mitochondria without map to known
|| STEP 12: quantify mitochondria reads
||
|| Usage: perl $0 -switch <parameter>
==========================================
||-d <def_file>:	Folder and file definition file; this is a text file with a single header line, and four tab-delimited fields: 
||                      ProjectName, ProjectFolder(including path), DataFileName or SRR ID (without surfix), ConvertedFileName, and Custom Ligation Adaptor Seq (use \'standard\' if the standard TGGAATTCTCGGGTGCCAAGG is used 
||-t <start_step>:	Optional: The step of analysis to begin with, 1-based numbering
||-e <end_step>:	Optional: The step of analysis to end; 1-based; use 0 for running to the last step
||-sra <is_sra>         Optional: if 1, specifies to treat FileName as SRR ID and download from SRR; if 0, specifies not to download from SRR; default 0
||-fa <is_fa>           Optional: if 1, skips the fastq to fasta conversion, and if download is set to 1, specifies that the data to be downloaded should be in fa format; default is 0
||-irt <IsoformThresh>  Optional: defines the threshold for the number of reads of a given isoform (summed acrossed all samples) to report; default 0; only applies to map to known classes
||-rl <runlog>          Optional: name of runlog file to save all commands; default will use ProjectName
||-g <MisMatchN>        Optioanl: number of mismatches for mapping to miRNA database. default 0.
=========================================
";

#default values for parameters
$def_file=0;
$start_step=1;
$end_step=0; #end_step of 0 means end at the very end
$is_sra=0;
$is_fa=0;
$runlog=0;
$IsoformThresh=0;
$MisMatchN=0;

for($i=0;$i<@ARGV;$i+=2)
{
	if($ARGV[$i] eq "-d")
	{
		$def_file=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-t")
	{
		$start_step=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-e")
	{
		$end_step=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-sra")
	{
	    $is_sra=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-fa")
	{
	    $is_fa=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-irt")
	{
	    $IsoformThresh=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-rl")
	{
	    $runlog=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-g")
	{
	    $MisMatchN=$ARGV[$i+1];
	}
	else
	{
		
	}
}
if(!$def_file) {die "$Usage\n";}
#============= start code here ============================

$index_dir_base="~/LuLab/DeepSeq/databases/";

@StepSurfix=('NA','.fa','_clipped_v2.fa','_collapsed_v2.fa','_read_len_summary_v2.txt','_collapsed_v2_L','','','');

$RawFolder="RawSeq/";
$ProcessedFolder="ProcessedSeq/";
$QuantifiedFolder="Quantified/";

$fastqdumpFolder='';

#----read in the folder/file definition file-----
@samples=();
open(IN,$def_file) || die "couldn't read definition file $def_file\n";
@wholefile=<IN>;
close(IN);
chomp(@wholefile);
shift(@wholefile);
for($i=0;$i<@wholefile;$i++)
{
    $wholefile[$i]=~s/\r//g;
    @temparray=split("\t",$wholefile[$i]);

    $tempfolder=$temparray[1];
    
    if($tempfolder!~/\/$/)
    {
	$tempfolder=$tempfolder.'/';
    }

    my @tempfiles=split("\,",$temparray[2]);

    #print join("\n",@tempfiles),"\n";

    $samples[$i]={Project=>$temparray[0],
		  Folder=>$tempfolder,
		  File=>\@tempfiles,
		  NewName=>$temparray[3],
		  Adaptor=>$temparray[4],
		  Species=>$temparray[5],
		  miRBase=>$temparray[6],
		  LenCutoff=>$temparray[7],
		  misMatch=>$temparray[8],
		  Assembly=>$temparray[9],
		  BowtieV=>$temparray[10],
    };
    #print $samples[$i]->{File}->[0],"\n";
}

#---setup analysis parameters
if(lc($samples[0]->{Species}) eq 'mouse' || lc($samples[0]->{Species}) eq 'mmu')
{
    $curSpecies="mmu";
}
else
{
    $curSpecies="hsa";
}

$miRBaseVer='';
if($samples[0]->{miRBase} eq 'v18')
{
    $PreFile="miRBase/v18_hairpin.fa";$MatureFile="miRBase/v18_mature.fa";$miRBaseVer="v18";}
elsif($samples[0]->{miRBase} eq 'v18_FPCM')
{
    $PreFile="miRBase/v18_FPCM_spiked_hairpin.fa";$MatureFile="miRBase/v18_FPCM_spiked_mature.fa";$miRBaseVer="v18FPCM";}
elsif($samples[0]->{miRBase} eq 'v20')
{
    $PreFile="miRBase/v20_hairpin.fa";$MatureFile="miRBase/v20_mature.fa";$miRBaseVer="v20";}
elsif($samples[0]->{miRBase} eq 'v21')
{
    $PreFile="miRBase/v21_hairpin.fa";$MatureFile="miRBase/v21_mature.fa";$miRBaseVer="v21";}
else 
{
    $PreFile=$samples[$i]->{miRBase}."_hairpin.fa";$MatureFile=$samples[$i]->{miRBase}."_mature.fa";$miRBaseVer="custom";
}

#--
$GenomeAssembly=$samples[0]->{Assembly};
$BowtieVersion=$samples[0]->{BowtieV};

$KnownFolder="Map_to_Known_v13_".$GenomeAssembly."_bt".$BowtieVersion.'/';
$MapToKnownSummaryFileSurfix="_mapping_summary_file_bt".$BowtieVersion."_".$curSpecies."_".$GenomeAssembly."_v13.txt";
$MapToKnownIsoformFileSurfix="_nonmiR_isoform_summary_bt".$BowtieVersion."_".$curSpecies."_".$GenomeAssembly."_f".$IsoformThresh."_v13.txt";
$MapToKnownLenFileSurfix="_Length_Summary_Mapped_Reads_bt".$BowtieVersion."_".$curSpecies."_".$GenomeAssembly."_v13.txt";
$MapToKnownFolderSurfix="_Map_to_Known_v13/";

$MitoFolder="Map_to_Mito_".$GenomeAssembly."_bt".$BowtieVersion.'/';
$MapToMitoSummaryFileSurfix="_mapping_to_Mito_summary_file_bt".$BowtieVersion."_".$curSpecies."_".$GenomeAssembly.".txt";
$MapToMitoIsoformFileSurfix="_Mito_isoform_summary_bt".$BowtieVersion."_".$curSpecies."_".$GenomeAssembly."_f".$IsoformThresh.".txt";
$MapToMitoFolderSurfix="_Map_to_Mito_v6/";

%MitoSize=("hg38",16569,"mm10",16299);

#---setup runlog
if(!$runlog)
{
    $runlog=$samples[0]->{Folder}.$samples[0]->{Project}.'_runlog.txt';
}
open($outfile,">>$runlog") || die "couldn't open runlog file $runlog\n";

print $outfile "\n------------------------------------------\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
$year+=1900;
$mon+=1;
$wday=$yday=$isdst=0;
print $outfile "\t\tTime: $year-$mon-$mday $hour:$min:$sec\n";
print $outfile $0." ".(join " ", @ARGV), "\n\n";


#---setup Folders
if(! -d $RawFolder) {mkdir($RawFolder);}
if(! -d $ProcessedFolder) {mkdir($ProcessedFolder);}
if(! -d $QuantifiedFolder) {mkdir($QuantifiedFolder);}


#------------
$cur_step=1;
if($start_step<=$cur_step)
{
    for($i=0;$i<@samples;$i++)
    {
	for($j=0;$j<@{$samples[$i]->{File}};$j++)
	{
	    $temp_cat_file_command='cat ';
	    if ($is_sra)
	    {
		if($is_fa)
		{
		    $tempcommand=$fastqdumpFolder."fastq-dump --fasta -O ".$samples[$i]->{Folder}.$RawFolder." ".$samples[$i]->{File}->[$j];
		    $temp_cat_file_command.=$samples[$i]->{Folder}.$RawFolder.$samples[$i]->{File}->[$j].".fasta ";
		}
		else
		{
		    $tempcommand=$fastqdumpFolder."fastq-dump -O ".$samples[$i]->{Folder}.$RawFolder." ".$samples[$i]->{File}->[$j];
		    $temp_cat_file_command.=$samples[$i]->{Folder}.$RawFolder.$samples[$i]->{File}->[$j].".fastq ";
		}
	    
		print "Start download for sample $i: $samples[$i]->{NewName}\n";
		print $outfile $tempcommand,"\n";
		`$tempcommand`;
	    }
	    else
	    {
		if($is_fa)
		{
		    $temp_cat_file_command.=$samples[$i]->{Folder}.$RawFolder.$samples[$i]->{File}->[$j].".fa";
		}
		else
		{
		    $temp_cat_file_command.=$samples[$i]->{Folder}.$RawFolder.$samples[$i]->{File}->[$j].".fastq";
		}
	    }
	}
	print "\tFinished downloading, sample $i\n";	

	#concatenate files

	if($is_fa)
	{
	    $tempcommand=$temp_cat_file_command.'> '.$samples[$i]->{Folder}.$RawFolder.$samples[$i]->{NewName}.$StepSurfix[$cur_step];
	    print $outfile $tempcommand,"\n";
	    `$tempcommand`;
	}
	else
	{
	    $tempcommand=$temp_cat_file_command.'> '.$samples[$i]->{Folder}.$RawFolder.$samples[$i]->{NewName}.".fastq";
	    print $outfile $tempcommand,"\n";
	    `$tempcommand`;
	}
    }

    for($i=0;$i<@samples;$i++)
    {
	if (!$is_fa)
	{
	    $tempinfile=$samples[$i]->{Folder}.$RawFolder.$samples[$i]->{NewName}.'.fastq';
	    $tempoutfile=$samples[$i]->{Folder}.$RawFolder.$samples[$i]->{NewName}.$StepSurfix[$cur_step];
	
	    $tempcommand="DeepSeq_fastq2fasta.pl ".$tempinfile.' > '.$tempoutfile;
	    print "Start STEP 1 for file $i: $samples[$i]->{NewName}\n";
	    print $outfile $tempcommand,"\n";
	    `$tempcommand`;
	    print "\tFinished. STEP1, file $i\n";
	}	
    }
}

#---STEP 2: clip off adatpor sequence---------
$cur_step=2;
if($start_step<=$cur_step && ($end_step>=$cur_step || $end_step==0))
{
    for($i=0;$i<@samples;$i++)
    {
	$tempinfile=$samples[$i]->{Folder}.$RawFolder.$samples[$i]->{NewName}.$StepSurfix[$cur_step-1];
	$tempoutfile=$samples[$i]->{Folder}.$ProcessedFolder.$samples[$i]->{NewName}.$StepSurfix[$cur_step];
	
	if(lc($samples[$i]->{Adaptor}) eq "standard") # this comparison is case insensitive
	{
	    $tempcommand="DeepSeq_clip_adaptors.pl -f ".$tempinfile.' -o '.$tempoutfile;
	}
	else
	{
	    $tempcommand="DeepSeq_clip_adaptors.pl -f ".$tempinfile.' -o '.$tempoutfile.' -A '.$samples[$i]->{Adaptor};
	}
	print "Start STEP 2 for file $i: $samples[$i]->{NewName}\n";
	print $outfile $tempcommand,"\n";
	`$tempcommand`;
	print "\tFinished. STEP2, file $i\n";


    }
}

#---STEP 3: collapse sequences to unique ones---------
$cur_step=3;
if($start_step<=$cur_step && ($end_step>=$cur_step || $end_step==0))
{
    for($i=0;$i<@samples;$i++)
    {
	$tempinfile=$samples[$i]->{Folder}.$ProcessedFolder.$samples[$i]->{NewName}.$StepSurfix[$cur_step-1];
	$tempoutfile=$samples[$i]->{Folder}.$ProcessedFolder.$samples[$i]->{NewName}.$StepSurfix[$cur_step];

	$tempcommand="DeepSeq_collapse_reads.pl ".$tempinfile.' '.$curSpecies.' > '.$tempoutfile;
	
	print "Start STEP 3 for file $i: $samples[$i]->{NewName}\n";
	print $outfile $tempcommand,"\n";
	`$tempcommand`;
	print "\tFinished. STEP3, file $i\n";

    }
}

#---STEP 4: summarize seq length distribution---------
$cur_step=4;
if($start_step<=$cur_step && ($end_step>=$cur_step || $end_step==0))
{
    for($i=0;$i<@samples;$i++)
    {
	$tempinfile=$samples[$i]->{Folder}.$ProcessedFolder.$samples[$i]->{NewName}.$StepSurfix[$cur_step-1];
	$tempoutfile=$samples[$i]->{Folder}.$QuantifiedFolder.$samples[$i]->{NewName}.$StepSurfix[$cur_step];
	$tempcommand="DeepSeq_sum_seq_length_v11.pl -f ".$tempinfile.' -o '.$tempoutfile;

	print "Start STEP 4 for file $i: $samples[$i]->{NewName}\n";
	print $outfile $tempcommand,"\n";
	`$tempcommand`;
	print "\tFinished. STEP4, file $i\n";
    }
}

#---STEP 5: filter out short reads---------
$cur_step=5;
if($start_step<=$cur_step && ($end_step>=$cur_step || $end_step==0))
{
    for($i=0;$i<@samples;$i++)
    {
	$tempinfile=$samples[$i]->{Folder}.$ProcessedFolder.$samples[$i]->{NewName}.$StepSurfix[$cur_step-2];
	$tempoutfile=$samples[$i]->{Folder}.$ProcessedFolder.$samples[$i]->{NewName}.$StepSurfix[$cur_step].$samples[$i]->{LenCutoff}.'.fa';
	$tempcommand="DeepSeq_filter_fa_with_length_v11.pl -L ".$samples[$i]->{LenCutoff}." -f ".$tempinfile.' -o '.$tempoutfile;

	print "Start STEP 5 for file $i: $samples[$i]->{NewName}\n";
	print $outfile $tempcommand,"\n";
	`$tempcommand`;
	print "\tFinished. STEP5, file $i\n";
    }
}

#---STEP 6: quantify the reads using miRDeep2---------
$cur_step=6;
if($start_step<=$cur_step && ($end_step>=$cur_step || $end_step==0))
{
    for($i=0;$i<@samples;$i++)
    {
	$tempinfile=$samples[$i]->{Folder}.$ProcessedFolder.$samples[$i]->{NewName}.$StepSurfix[$cur_step-1].$samples[$i]->{LenCutoff}.'.fa';
	
	$tempoutfile=$samples[$i]->{NewName}.'_v2_L'.$samples[$i]->{LenCutoff}.'g'.$MisMatchN.'_'.$curSpecies.'_'.$miRBaseVer;
	$tempcommand="DeepSeq_quantify_miR_reads.pl -p ".$index_dir_base.$PreFile." -m ".$index_dir_base.$MatureFile." -r ".$tempinfile." -t ".$curSpecies." -y ".$tempoutfile." -g ".$MisMatchN;

	chdir($samples[$i]->{Folder}) || die "error executing cd command\n";

	
	print "Start STEP 6 for file $i: $samples[$i]->{NewName}\n";
	print $outfile $tempcommand,"\n";
	`$tempcommand`;
	print "\tFinished. STEP6, file $i\n";

    }
}

#---STEP 7: summarize quantification data---------
$cur_step=7;
if($start_step<=$cur_step && ($end_step>=$cur_step || $end_step==0))
{
    print "Start STEP $cur_step\n";

    $tempinprefix="expression_analyses/expression_analyses_";

    $tempoutfile=$samples[0]->{Folder}.$QuantifiedFolder.$samples[0]->{Project}."_miRNA_Quantification_Summary_".$curSpecies."_".$miRBaseVer."_L".$samples[0]->{LenCutoff}."g".$MisMatchN.".txt";
  
    @allfiles=();


    #read in data
    for($i=0;$i<@samples;$i++)
    {
	$tempinfile=$samples[$i]->{Folder}.$tempinprefix.$samples[$i]->{NewName}.'_v2_L'.$samples[$i]->{LenCutoff}.'g'.$MisMatchN.'_'.$curSpecies.'_'.$miRBaseVer."/miRNA_expressed.csv";

	open($infile,$tempinfile) || die "couldn't open $tempinfile\n";
	my @tempwholefile=<$infile>;
	close $infile;

	shift @tempwholefile;
	chomp @tempwholefile;

	my @tempExpression=();
	for($j=0;$j<@tempwholefile;$j++)
	{
	    @temparray=split("\t",$tempwholefile[$j]);
	    $tempExpression[$j]={Mature=>$temparray[0],
				 Reads=>$temparray[1],
				 Precursor=>$temparray[2]
	    };
	}
	push @allfiles,\@tempExpression;

    }

    #summarize data
    $titleline='';
    @sum_data=();
    for($i=0;$i<@samples;$i++)
    {
	if($i==0)
	{
	    $titleline="Mature\tPrecursor\t".$samples[$i]->{NewName};
	    for($j=0;$j<@{$allfiles[$i]};$j++)
	    {
		$sum_data[$j]=$allfiles[$i]->[$j]->{Mature}."\t".$allfiles[$i]->[$j]->{Precursor}."\t".$allfiles[$i]->[$j]->{Reads};
	    }
	}
	else
	{
	    $titleline=$titleline."\t".$samples[$i]->{NewName};

	    for($j=0;$j<@{$allfiles[$i]};$j++)
	    {
		if($allfiles[$i]->[$j]->{Mature} eq $allfiles[0]->[$j]->{Mature})
		{
		    $sum_data[$j]=$sum_data[$j]."\t".$allfiles[$i]->[$j]->{Reads};
		}
		else
		{
		    print "Error encoutered $i\n";
		}
	    }
	}
    }

    #output data

    open($outfile2,">$tempoutfile") || die "couldn't open output $tempoutfile\n";

    print $outfile2 $titleline,"\n";
    print $outfile2 join("\n",@sum_data),"\n";
    
    close($outfile2);

    print "\tFinished. STEP $cur_step, file $i\n";

}

#---STEP 8: quantify miRNA isoforms---------
$cur_step=8;
if($start_step<=$cur_step && ($end_step>=$cur_step || $end_step==0))
{
    $tempinprefix="expression_analyses/expression_analyses_";
    $MatureFile=~/([^\/]+)\.fa/;
    $tempMatureMapped_prefix=$1;

    $tempSummaryFile=$samples[0]->{Folder}.$QuantifiedFolder.$samples[0]->{Project}."_miRNA_isoform_summary_".$curSpecies."_".$miRBaseVer."_L".$samples[0]->{LenCutoff}."g".$MisMatchN.".txt";

    $tempcommand="DeepSeq_summarize_miRNA_isoforms_v11.pl -o ".$tempSummaryFile." -p ".$index_dir_base.$PreFile." -m ".$samples[0]->{Folder}.$tempinprefix.$samples[0]->{NewName}.'_v2_L'.$samples[0]->{LenCutoff}.'g'.$MisMatchN.'_'.$curSpecies.'_'.$miRBaseVer."/".$tempMatureMapped_prefix."_mapped.arf";

    for($i=0;$i<@samples;$i++)
    {
	$tempcommand=$tempcommand." -f ".$samples[$i]->{Folder}.$tempinprefix.$samples[$i]->{NewName}.'_v2_L'.$samples[$i]->{LenCutoff}.'g'.$MisMatchN.'_'.$curSpecies.'_'.$miRBaseVer."/".$samples[$i]->{NewName}.$StepSurfix[$cur_step-3].$samples[$i]->{LenCutoff}."_mapped.arf";
    }

    print "Start STEP $cur_step\n";
    print $outfile $tempcommand,"\n";
    `$tempcommand`;
    print "\tFinished. STEP$cur_step\n";

}


#---STEP 9: map to known RNAs or the genome---------
$cur_step=9;
if($start_step<=$cur_step && ($end_step>=$cur_step || $end_step==0))
{
    if(! -d $KnownFolder) {mkdir($KnownFolder);}

    for($i=0;$i<@samples;$i++)
    {
	$tempinfile=$samples[$i]->{Folder}.$ProcessedFolder.$samples[$i]->{NewName}.$StepSurfix[5].$samples[$i]->{LenCutoff}.'.fa';
	$tempoutfile=$samples[$i]->{NewName};

	chdir($samples[$i]->{Folder}) || die "error executing cd command\n";


	$tempcommand="DeepSeq_map_small_RNAs_to_Known_v13.pl -f ".$tempinfile.' -d '.$samples[$i]->{Folder}.$KnownFolder.' -b '.$tempoutfile." -v ".$BowtieVersion." -p 6 -L ".$samples[$i]->{LenCutoff}.' -s '.$curSpecies.' -a '.$GenomeAssembly;


	print "Start STEP 9 for file $i: $samples[$i]->{NewName}\n";
	print $outfile $tempcommand,"\n";
	`$tempcommand`;
	print "\tFinished. STEP9, file $i\n";

    }
}

#---STEP 10: quantify  known RNA or genome isoforms---------
$cur_step=10;
if($start_step<=$cur_step && ($end_step>=$cur_step || $end_step==0))
{
    $tempSummaryFile=$samples[0]->{Folder}.$QuantifiedFolder.$samples[0]->{Project}.$MapToKnownSummaryFileSurfix;

    $tempIsoformFile=$samples[0]->{Folder}.$QuantifiedFolder.$samples[0]->{Project}.$MapToKnownIsoformFileSurfix;
    $tempLenFile=$samples[0]->{Folder}.$QuantifiedFolder.$samples[0]->{Project}.$MapToKnownLenFileSurfix;

    $tempcommand="DeepSeq_summarize_isoforms_for_known_classes_v13.pl -s ".$tempSummaryFile." -o ".$tempIsoformFile." -b ".$BowtieVersion." -f ".$IsoformThresh." -p ".$curSpecies." -m ".$tempLenFile." -a ".$GenomeAssembly;

    for($i=0;$i<@samples;$i++)
    {
	$tempcommand=$tempcommand." -d ".$samples[$i]->{Folder}.$KnownFolder.$samples[$i]->{NewName}.$MapToKnownFolderSurfix;
    }

    print "Start STEP 10\n";
    print $outfile $tempcommand,"\n";
    `$tempcommand`;
    print "\tFinished. STEP10\n";

}

#---STEP 11: map directly to mitochondria genome---------
$cur_step=11;
if($start_step<=$cur_step && ($end_step>=$cur_step || $end_step==0))
{
    if(! -d $MitoFolder) {mkdir($MitoFolder);}

    for($i=0;$i<@samples;$i++)
    {
	$tempinfile=$samples[$i]->{Folder}.$ProcessedFolder.$samples[$i]->{NewName}.$StepSurfix[5].$samples[$i]->{LenCutoff}.'.fa';
	$tempoutfile=$samples[$i]->{NewName};

	chdir($samples[$i]->{Folder}) || die "error executing cd command\n";


	$tempcommand="DeepSeq_map_small_RNAs_to_Mito_v11.pl -f ".$tempinfile.' -d '.$samples[$i]->{Folder}.$MitoFolder." -v 2"." -p 6 ".' -s '.$curSpecies.' -a '.$GenomeAssembly.' -b '.$tempoutfile;


	print "Start STEP 11 for file $i: $samples[$i]->{NewName}\n";
	print $outfile $tempcommand,"\n";
	`$tempcommand`;
	print "\tFinished. STEP11, file $i\n";

    }
}

#---STEP 12: quantify  mitochondria mapping isoforms---------
$cur_step=12;
if($start_step<=$cur_step && ($end_step>=$cur_step || $end_step==0))
{
    $tempSummaryFile=$samples[0]->{Folder}.$QuantifiedFolder.$samples[0]->{Project}.$MapToMitoSummaryFileSurfix;

    $tempIsoformFile=$samples[0]->{Folder}.$QuantifiedFolder.$samples[0]->{Project}.$MapToMitoIsoformFileSurfix;

    $tempcommand="DeepSeq_summarize_mito_isoforms_v11.pl -s ".$tempSummaryFile." -o ".$tempIsoformFile." -b ".$BowtieVersion." -f ".$IsoformThresh." -p ".$curSpecies." -bed ".$MitoSize{$GenomeAssembly}." -bf "."_bt".$BowtieVersion.'.bg';

    for($i=0;$i<@samples;$i++)
    {
	$tempcommand=$tempcommand." -d ".$samples[$i]->{Folder}.$MitoFolder.$samples[$i]->{NewName}.$MapToMitoFolderSurfix;
    }

    print "Start STEP $cur_step\n";
    print $outfile $tempcommand,"\n";
    `$tempcommand`;
    print "\tFinished. STEP $cur_step\n";

}


#----------------------
close($outfile);
