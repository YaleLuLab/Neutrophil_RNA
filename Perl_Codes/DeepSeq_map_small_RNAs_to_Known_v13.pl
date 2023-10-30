#!/usr/bin/perl -w

use lib "$ENV{\"HOME\"}/Perl_Codes/Jun_lib/";
use READ_FILE;
use NA_SEQ;

$Usage="==========================================
|| Jun Lu, June 2022, tested
|| v13 is modified from v11--use new database, based on RNAcentral Release 20, for mapping
|| v13 also made the program more efficient by putting multiple categories of RNA into the same mapping step rather than separately.
|| v13 added support for mm39
|| v13 outputs up to 10 mappings for the same input sequence.
|| v6 is modified on v1
|| v6 maps to Illumina adaptors first (uses bowtie 1 with 1 mismatch, regardless of naming) to remove artifacts of eventual mapping to genome
|| Currently allowed genomes: hg18, hg19, hg38, mm9, mm10(mm10 only with bowtie2), mm39
||
|| Usage: perl $0 -f <fa_file> -d <dir> -b <base> -L <len_cutoff> -v <bowtie_version> -p <processes> -r <re_run> -a <Assembly> -s <run_summary>
==========================================
||-f <fa_file>:	DeepSeq reads clipped and collapsed and after length cutoff
||-d <dir>:	directory to hold outputs
||-b <base>:	base name to add to the output files
||-L <len_cutoff>:	not used
||-v <bowtie_version>: specify 1 for bowtie, specify 2 for bowtie2; default 2
||-p <processes>: specify the number of threads to run using bowtie; default 6
||                         Note: no more than 7 on our server
||-r <re_run>:      set to 1 will run through whole script again, regardless of whether some intermediate files already exist or not
||                         default is 0, which means do not re-run to generate files that were already generated
||-s <Species>:     Define human or hsa, or mouse or mmu; default hsa
||-a <Assembly>: the assembly of genome that the last step of mapping should map to; e.g. hg19; default hg19
||-rs <run_summary>: default is 1; which runs summary; change to 0 if not; currently only running with bowtie 2 results
==========================================
";

#default values for parameters
$fa_file=0;
$dir=0;
$base=0;
$len_cutoff=16;
$bowtie_version=2;
$processes=6;
$re_run=0; # this specifies if an output file exists, whether to rerun the relevant code; default not to rerun
$run_summary=1;
$Species="hsa";

for($i=0;$i<@ARGV;$i+=2)
{
	if($ARGV[$i] eq "-f")
	{
		$fa_file=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-d")
	{
		$dir=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-b")
	{
		$base=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-L")
	{
		$len_cutoff=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-v")
	{
		$bowtie_version=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-p")
	{
		$processes=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-r")
	{
		$re_run=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-a")
	{
		$Assembly=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-rs")
	{
		$run_summary=$ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-s")
	{
		$Species=$ARGV[$i+1];
	}
	else
	{
		
	}
}
if(!$fa_file || !$dir || !$base) {die "$Usage\n";}
#============= start code here ============================

if(lc($Species) eq "human" || lc($Species) eq "hsa")
{
    $Species="hsa";
}
else
{
    $Species="mmu";
}

#------------------------------------
$run_Illumina=1;
$run_ncRNA=1;
$run_genome=1;

#------------------------------------
$output_dir_base="Map_to_Known_v13";
$output_dir="";

$outfile_base="";
$infilename="";
#$codedir="/gpfs/ysm/project/jun_lu/jl894/LuLab/HPC_codes/";
#$index_dir_base="/gpfs/ysm/project/jun_lu/jl894/LuLab/DeepSeq/databases/";
$index_dir_base="$ENV{\"HOME\"}/LuLab/DeepSeq/databases/"; #change to cope with McCleary migration

@len_summary_range=(1,100);

if($bowtie_version==1)
{
    $bowtieV="bt1";
    $bowtieSwitches="-f -p $processes -v 1 -k 3 --best ";
    $bowtieProgram="bowtie ";
    $bowtieGenomeSwitches="-f -p $processes -v 1 -k 10 ";
}
else
{
    $bowtieV="bt2";
    $bowtieSwitches="-f -p $processes -D 20 -R 3 -N 0 -L 18 -i S,1,0.50 -k 3 --score-min L,0,-0.6 ";
    $bowtieProgram="bowtie2 ";
    $bowtieGenomeSwitches="-f -p $processes -D 20 -R 3 -N 0 -L 18 -i S,1,0.50 -k 10 --score-min L,0,-0.6 ";
}

if($Assembly eq 'hg19')
{
   $index_base_genome=$index_dir_base."genomes/hg19/hg19";
   #$TSS_file_name=$index_dir_base."genomes/hg19_annotations/RefSeq_NM_TSS_hg19.txt";
}
elsif($Assembly eq 'hg18')
{
    $index_base_genome=$index_dir_base."genomes/hg18/hg18";
}
elsif($Assembly eq 'hg38')
{
    if($bowtie_version==1) {$index_base_genome=$index_dir_base."genomes/hg38/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set";}
    else {$index_base_genome=$index_dir_base."genomes/hg38/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index";}
}
elsif($Assembly eq 'mm9')
{
    $index_base_genome=$index_dir_base."genomes/mm9/mm9";
}
elsif($Assembly eq 'mm10')
{
    $index_base_genome=$index_dir_base."genomes/mm10/mm10/mm10";
}
elsif($Assembly eq 'mm39')
{
    $index_base_genome=$index_dir_base."genomes/mm39/GRCm39/GRCm39";
}
#------first set output parameters----------------------

$dir=~s/\/$//;  #remove end slash

if (! -d $dir) {mkdir $dir;}

$output_dir=$dir.'/'.$base."_".$output_dir_base;

if(! -d $output_dir) {mkdir $output_dir};

$outfile_base=$output_dir.'/'.$base;
$infilename=$fa_file;

$summary_file=$outfile_base."_Mapping_Summary_".$bowtieV.'.txt';

print STDERR "\n";

if($run_summary)
{
    print STDERR "Summarizing original sequence file...\n";
    open(SUMMARY,">$summary_file") || die "couldn't open summary file $summary_file\n";
   
    print SUMMARY "Step\tFile\tTotal Unique Seq\tTotal Reads";
    
    foreach $i ($len_summary_range[0]..$len_summary_range[1])
    {
	print SUMMARY "\t$i";
    }
    print SUMMARY "\n";
 
    $tempsummary=summarize_fa($infilename,\@len_summary_range);
    print SUMMARY "Input File\t$tempsummary\n";

    print STDERR "\t\tSummarizing step 1 completed\n";
}

#-------------run mapping to Illumina Adaptors------------------
if($run_Illumina)
{
    print STDERR "Bowtie mapping to Illumina Adaptors and Spike-In ...\n";

    $outfile_base=$outfile_base."_IaS".$bowtieV;

    my $outfilename=$outfile_base.".fa";

    my $samfile=$outfile_base.".sam";
    my $bwtfile=$outfile_base.".bwt";

    my $index_base=$index_dir_base."RNAcentral/Genome_Coordinates/Release_20/subselected_gff3/bowtie_index/IlluminaAdaptor_SpikeIn";

    my $command;
    
    #only run bowtie1; modified to fit the lastest bowtie syntax.
    
    $command="bowtie -f -p $processes -v 1 -k 3 --best "."--norc "."--un $outfilename "."-x $index_base "."$infilename "."$bwtfile";
    
    #print "$command\n";
    
    if($re_run || ! -e $outfilename) {
	`$command`;
	print STDERR "\t\tFinished. $outfilename generated.\n";
    }
    else
    {
	print STDERR "\t\tSkipped. $outfilename exists.\n";
    }

    #.................start summarize mapped...............
    if($run_summary) {
	print STDERR "\t\tSummarizing alignments...\n";
	my $summary;
	
	$summary=summarizeAlignments(1,$bwtfile,\@len_summary_range,0);
		
	print SUMMARY "$summary";
	print STDERR "\t\tSummarizing completed\n";
    }

    #.............end summarize.................
    $infilename=$outfilename;
}

#-------------run mapping to ncRNA------------------
if($run_ncRNA)
{
    print STDERR "Bowtie mapping to ncRNAs...\n";

    $outfile_base=$outfile_base."_ncR".$bowtieV;

    my $outfilename=$outfile_base.".fa";

    my $samfile=$outfile_base.".sam";
    my $bwtfile=$outfile_base.".bwt";

    my $index_base;

    if($Species eq 'hsa')
    {
	$index_base=$index_dir_base."RNAcentral/Genome_Coordinates/Release_20/subselected_gff3/bowtie_index/RNAcenV20_hsa_noLPM_50flank_with45S";
    }
    elsif($Species eq 'mmu')
    {
	$index_base=$index_dir_base."RNAcentral/Genome_Coordinates/Release_20/subselected_gff3/bowtie_index/RNAcenV20_mmu_noLPM_50flank_with45S";
    }

    my $command;
    if($bowtie_version==2)
    {
	$command=$bowtieProgram.$bowtieSwitches."--norc "."--un $outfilename "."-x $index_base "."-U $infilename "."-S $samfile";
    }
    elsif($bowtie_version==1)
    {
	$command=$bowtieProgram.$bowtieSwitches."--norc "."--un $outfilename "."$index_base "."$infilename "."$bwtfile";
    }

    #print "$command\n";
    
    if($re_run || ! -e $outfilename) {
	`$command`;
	print STDERR "\t\tFinished. $outfilename generated.\n";
    }
    else
    {
	print STDERR "\t\tSkipped. $outfilename exists.\n";
    }

    #................start summarize mapped...............
    if($run_summary) {
	print STDERR "\t\tSummarizing alignments...\n";
	my $summary;
	
	if($bowtie_version==1)
	{
	    $summary=summarizeAlignments($bowtie_version,$bwtfile,\@len_summary_range,0);
	}
	elsif($bowtie_version==2)
	{
	    $summary=summarizeAlignments($bowtie_version,$samfile,\@len_summary_range,0);
	}
	
	print SUMMARY "$summary";
	print STDERR "\t\tSummarizing complete\n";
    }

    #..............end summarize............
    $infilename=$outfilename;
}


#-------------run mapping to the genome with the remaining sequences------------------

if($run_genome)
{
    print STDERR "Bowtie mapping to the genome...\n";

    $outfile_base=$outfile_base."_GENOME$Assembly".$bowtieV;

    my $outfilename=$outfile_base.".fa";

    my $samfile=$outfile_base.".sam";
    my $bwtfile=$outfile_base.".bwt";

    my $command;
    if($bowtie_version==2)
    {
	$command=$bowtieProgram.$bowtieGenomeSwitches."--un $outfilename "."-x $index_base_genome "."-U $infilename "."-S $samfile";
    }
    elsif($bowtie_version==1)
    {
	$command=$bowtieProgram.$bowtieGenomeSwitches."--un $outfilename "."$index_base_genome "."$infilename "."$bwtfile";
    }

    #print "$command\n";
    
    if($re_run || ! -e $outfilename) {
	`$command`;
	print STDERR "\t\tFinished. $outfilename generated.\n";
    }
    else
    {
	print STDERR "\t\tSkipped. $outfilename exists.\n";
    }

    #................start summarize mapped...............
    if($run_summary) {
	print "\t\tSummarizing alignments...\n";
	my $summary;
	
	if($bowtie_version==1)
	{
	    $summary=summarizeAlignments($bowtie_version,$bwtfile,\@len_summary_range,1);
	}
	elsif($bowtie_version==2)
	{
	    $summary=summarizeAlignments($bowtie_version,$samfile,\@len_summary_range,1);
	}
	
	print SUMMARY "$summary";
	print STDERR "\t\tSummarizing complete\n";
    }

    #..............end summarize............
    $infilename=$outfilename;
    #$output_sam=$samfile;
}


#-----------------------------------------------------------------------------
if($run_summary)
{
    print "Summarizing unmapped fa file...\n";
    my $summary=summarize_fa($infilename,\@len_summary_range);
    print SUMMARY "UnMapped\t$summary\n";
    print "\t\tSummarizing complete\n";    
    close SUMMARY;
}

exit;
#----------------------------------------------------------------------------
sub summarize_fa{

    my ($infile,$len_summary_range)=@_;
    my $summary;

    my $seqs=NA_SEQ::read_FASTA_file($infile);

    my %size_hash=();
    my $total_seqs=@{$seqs};
    my $total_reads=0;
    my ($i,$id,$seq,$seq_len,$tempcount);

    for($i=0; $i<$total_seqs;$i++)
    {
        $id=$seqs->[$i]->name();
        $seq=$seqs->[$i]->seq;
        $seq_len=length($seq);

        if($id=~/.+_x(\d+)$/)
        {
                $temp_count=$1;
                $total_reads+=$temp_count;

		if(!exists($size_hash{$seq_len}))
		{
		    $size_hash{$seq_len}=$temp_count;
		}
		else
		{
		    $size_hash{$seq_len}+=$temp_count;
		}
	}
        else
        {
                $temp_count=0;
		print STDERR "WARNing: $id is of incorrect format\n";
                exit;
        }
    }
	
    $summary="$infile\t$total_seqs\t$total_reads";

    foreach $i ($len_summary_range->[0]..$len_summary_range->[1])
    {
	if(exists($size_hash{$i})) {
	    $summary.="\t$size_hash{$i}";
	}
	else {
	    $summary.="\t0";
	}
    }       

    return $summary;
}


sub summarizeAlignments{
    my ($bowtie_version,$infile,$len_summary_range,$isgenome)=@_;
    my $summary;

    my %RNAtypes=();

    #my %size_hash=();
    my %unique_seqs;
    my $i=0;
    #my @wholefile;
    #my @temparray;
    my $occurance;
    my $seq_len;

    my $align;

    if($bowtie_version==1)
    {
	$align=READ_FILE::read_bowtie1_bwt_alignments($infile);
    }
    elsif($bowtie_version==2)
    {
	$align=READ_FILE::read_bowtie2_sam_alignments($infile,1);
    }

    my $outname=$infile;
    if($bowtie_version==2) {
	$outname=~s/\.sam/\.mapped_reads\.txt/;
    } elsif ($bowtie_version==1) {
	$outname=~s/\.bwt/\.mapped_reads\.txt/;
    }
    
    my $tempout;
    open($tempout,">$outname") || die "couldn't open output file $outname\n";
    print $tempout "SeqName\tSequence\tLength\n";
    
    my $tempRNAtype;
    
    for($i=0;$i<@$align;$i++)
    {
	if(!$isgenome) {
	    if($align->[$i]->{rname}=~/(.+?)__/) {
		$tempRNAtype=$1;
	    } else {
		$tempRNAtype="NA";
	    }
	    
	} else {
	    $tempRNAtype="Genome";
	}
	
	if(!exists($RNAtypes{$tempRNAtype})){
	    my %size_hash=();
	    $RNAtypes{$tempRNAtype}=\%size_hash;
	}
	
	if(!exists($unique_seqs{$align->[$i]->{qname}}) && $align->[$i]->{mapped})
	{
	    $unique_seqs{$align->[$i]->{qname}}=1;
	    
	    print $tempout "$align->[$i]->{qname}\t";
	    if(!$align->[$i]->{reverse}) {
		print $tempout $align->[$i]->{seq},"\t",length($align->[$i]->{seq}),"\n";
	    } else {
		my $temprc=$align->[$i]->{seq};
		$temprc=reverse($temprc);
		$temprc=~tr/ATCG/TAGC/;
		print $tempout $temprc,"\t",length($temprc),"\n";
	    }
	    
	    if($align->[$i]->{qname}=~/.+_x(\d+)/)
	    {
		$occurance=$1;
		
		if(!exists($RNAtypes{$tempRNAtype}->{TotalSeqs})) {
		    $RNAtypes{$tempRNAtype}->{TotalSeqs}=1;
		} else {
		    $RNAtypes{$tempRNAtype}->{TotalSeqs}+=1;
		}
		
		if(!exists($RNAtypes{$tempRNAtype}->{TotalReads})) {
		    $RNAtypes{$tempRNAtype}->{TotalReads}=$occurance;
		} else {
		    $RNAtypes{$tempRNAtype}->{TotalReads}+=$occurance;
		}
		
		$seq_len=length($align->[$i]->{seq});
		
		if(!exists($RNAtypes{$tempRNAtype}->{$seq_len}))
		{
		    $RNAtypes{$tempRNAtype}->{$seq_len}=$occurance;
		}
		else
		{
		    $RNAtypes{$tempRNAtype}->{$seq_len}+=$occurance;
		}
	    }
	    else
	    {
		print STDERR "$align->[$i]->{qname} is of wrong format\n";
		exit;
	    }
	}
    }
  
    #output to $summary
    $summary="";
    foreach my $tempRNAtype2 (sort keys(%RNAtypes)) {
	if(!exists($RNAtypes{$tempRNAtype2}->{TotalSeqs})) {
	    $RNAtypes{$tempRNAtype2}->{TotalSeqs}=0;
	}
	if(!exists($RNAtypes{$tempRNAtype2}->{TotalReads})) {
	    $RNAtypes{$tempRNAtype2}->{TotalReads}=0;
	}

	$summary.="$tempRNAtype2\t$infile\t".$RNAtypes{$tempRNAtype2}->{TotalSeqs}."\t".$RNAtypes{$tempRNAtype2}->{TotalReads};
	foreach $i ($len_summary_range->[0]..$len_summary_range->[1])
	{
	    if(exists($RNAtypes{$tempRNAtype2}->{$i})) {
		$summary.="\t".$RNAtypes{$tempRNAtype2}->{$i}; 
	    }
	    else {
		$summary.="\t0";
	    }
	}
	$summary.="\n";
    }
    return $summary;
}

