% add code to enter the direcotry with mapped data

% Define data specific parameters
PlotFolder='Plots/';

if(~exist(PlotFolder,'dir'))
    mkdir(PlotFolder);
end

Samples={'HeLa_Input','HeLa_Biotin','H9_Input','H9_Biotin'};
SampleFields={'HeLa_Input','HeLa_Biotin','H9_Input','H9_Biotin'};
Samples_for_p=[1 2;3 4];

Mapped_ReadLength_File='GlycoRNA_Length_Summary_Mapped_Reads_bt2_hsa_hg19_v13.txt';

PlotLines={'-k','-b','-g','--r'};

IsoformFileName='GlycoRNA_nonmiR_isoform_summary_bt2_hsa_hg19_f10_v13.txt';
IsoformOutFileName='GlycoRNA_nonmiR_isoform_Unique_Seq_summary_f100_v13.txt';
NormOutFileName='GlycoRNA_nonmiR_isoform_Normalized_Unique_Seq_summary_f100.txt';

save_workspace='GlycoRNA_newMap_data.mat';

ReadSummaryFile='GlycoRNA_mapping_summary_file_bt2_hsa_hg19_v13.txt';

Min_floor_level=1; %set RPM of 1 as floor



%================Process Alignment Summary File to get rid of redundancy=====================

DataTable=Jun_read_tdf_MultiCPU(IsoformFileName);

%Find unique sequences and assemble counts for the unqiue sequences

[UniqueSeq tempIdx1 tempIdx2]=unique(DataTable.ReadSeq);

%UniqueData stores only the unique sequences and their counts etc
UniqueData=DataTable;
TableFields=fieldnames(UniqueData);
for i=1:length(TableFields)
    UniqueData.(TableFields{i})=UniqueData.(TableFields{i})(1:length(UniqueSeq));    
end

for i=1:length(UniqueSeq)
    cur_idx=find(tempIdx2==i);
    
    UniqueData.RNAname{i}=DataTable.RNAname{cur_idx(1)};
    UniqueData.RNAtype{i}=DataTable.RNAtype{cur_idx(1)};
    UniqueData.ReadSeq{i}=DataTable.ReadSeq{cur_idx(1)};
    for j=4:7
        UniqueData.(TableFields{j})(i)=DataTable.(TableFields{j})(cur_idx(1));
    end
    
    tmpTotalCount=0;
    for j=8:(length(TableFields)-2) %for all sample count fields
        UniqueData.(TableFields{j})(i)=max(DataTable.(TableFields{j})(cur_idx));
        tmpTotalCount=tmpTotalCount+UniqueData.(TableFields{j})(i);
    end
    
    UniqueData.TotalCount(i)=tmpTotalCount;
    UniqueData.NumOfMatch(i)=DataTable.NumOfMatch(cur_idx(1));
end


%Filter UniqueData
UniqueDataOrig=UniqueData;

ThresholdCount=100;
ThresholdLen=16;

SequencesForAnalysis=(UniqueDataOrig.TotalCount>ThresholdCount) &...
    (UniqueDataOrig.ReadLength>=ThresholdLen);

UniqueData=reorder_tdf_structure(UniqueDataOrig,SequencesForAnalysis,1); 

save(save_workspace,'DataTable','UniqueDataOrig','UniqueData','UniqueSeq');


%Output data
TableFields=fieldnames(UniqueData);

fid=fopen(IsoformOutFileName,'w');
fprintf(fid,[strjoin(TableFields,'\t') '\n']);
for i=1:length(UniqueData.RNAname)
    for j=1:3
        fprintf(fid,'%s\t',UniqueData.(TableFields{j}){i});
    end
    for j=4:length(TableFields)
        fprintf(fid,'%f\t',UniqueData.(TableFields{j})(i));
    end
    fprintf(fid,'\n');
end
fclose(fid);

%======Calculate and Plot Enrichment for genome mapping=====================
% read and normalize data to RPM
load(save_workspace);
MapSummary=Jun_read_tdf(ReadSummaryFile);
[tempData MapSumIdx]=setdiff(MapSummary.Step,{'Input File','IlluminaAdaptor','UnMapped'});

for i=1:length(SampleFields)
   TotalMappedReads(i)=sum(MapSummary.(SampleFields{i})(MapSumIdx));
end


NormUnique=UniqueData;
for i=1:length(SampleFields)
   NormUnique.(SampleFields{i})=UniqueData.(SampleFields{i})/sum(MapSummary.(SampleFields{i})(MapSumIdx))*1000000;
end


%calculate P_value using Fisher Exact test
for i=1:size(Samples_for_p,1)
   for j=1:length(UniqueData.(SampleFields{i}))
       if mod(j,2000)==0
           disp(num2str(j));
       end
       [htemp, ptemp]=fishertest([UniqueData.(SampleFields{Samples_for_p(i,2)})(j),...
           TotalMappedReads(Samples_for_p(i,2))-UniqueData.(SampleFields{Samples_for_p(i,2)})(j);...
           UniqueData.(SampleFields{Samples_for_p(i,1)})(j),...
           TotalMappedReads(Samples_for_p(i,1))-UniqueData.(SampleFields{Samples_for_p(i,1)})(j)],'Tail','right');
       NormUnique.([SampleFields{Samples_for_p(i,2)} '_P'])(j,1)=ptemp;
   end
   NormUnique.([SampleFields{Samples_for_p(i,2)} '_FDR'])=calc_fdr_value(NormUnique.([SampleFields{Samples_for_p(i,2)} '_P']));
end


save(save_workspace,'DataTable','UniqueDataOrig','UniqueData','UniqueSeq','NormUnique','Samples','SampleFields');

%Output data
TableFields=fieldnames(NormUnique);

fid=fopen(NormOutFileName,'w');
fprintf(fid,[strjoin(TableFields,'\t') '\n']);
for i=1:length(NormUnique.RNAname)
    for j=1:3
        fprintf(fid,'%s\t',NormUnique.(TableFields{j}){i});
    end
    for j=4:length(TableFields)
        fprintf(fid,'%f\t',NormUnique.(TableFields{j})(i));
    end
    fprintf(fid,'\n');
end
fclose(fid);


%------plot Input vs Biotin
curSamples=[1,2];

plot(log2(NormUnique.(SampleFields{curSamples(1)})+Min_floor_level), log2(NormUnique.(SampleFields{curSamples(2)})+Min_floor_level),'.b','MarkerSize',8);
xlim([0, 9]);
ylim([0 9]);
set(gca,'Visible','off')
fig_prop=get(gca);
axes('Position',fig_prop.Position,...
       'XAxisLocation','bottom',...
       'YAxisLocation','left',...
       'Color','none',...
       'XTickLabel',fig_prop.XTickLabel,...
       'YTickLabel',fig_prop.YTickLabel,...
       'XColor','k','YColor','k',...
       'LineWidth',1,...
       'XTick',fig_prop.XTick,...
       'YTick',fig_prop.YTick,...
       'XLim',fig_prop.XLim,...
       'YLim',fig_prop.YLim,...
       'TickDir','out');
set(gca,'FontSize',16,'FontWeight','bold','FontName','Arial');
saveas(gca, [PlotFolder 'Hela_Biotin_vs_Input_log2_scatter_plot.png'],'png');

close
curSamples=[3,4];
plot(log2(NormUnique.(SampleFields{curSamples(1)})+Min_floor_level), log2(NormUnique.(SampleFields{curSamples(2)})+Min_floor_level),'.b','MarkerSize',8);
xlim([0, 8]);
ylim([0 8]);
set(gca,'Visible','off')
fig_prop=get(gca);
axes('Position',fig_prop.Position,...
       'XAxisLocation','bottom',...
       'YAxisLocation','left',...
       'Color','none',...
       'XTickLabel',fig_prop.XTickLabel,...
       'YTickLabel',fig_prop.YTickLabel,...
       'XColor','k','YColor','k',...
       'LineWidth',1,...
       'XTick',fig_prop.XTick,...
       'YTick',fig_prop.YTick,...
       'XLim',fig_prop.XLim,...
       'YLim',fig_prop.YLim,...
       'TickDir','out');
set(gca,'FontSize',16,'FontWeight','bold','FontName','Arial');
saveas(gca, [PlotFolder 'H9_Biotin_vs_Input_log2_scatter_plot.png'],'png');


