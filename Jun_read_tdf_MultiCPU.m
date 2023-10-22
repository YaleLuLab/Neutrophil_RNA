function S=Jun_read_tdf_MultiCPU(fname,PoolSize,dlm,nlines)
% function S=Jun_read_tdf_MultiCPU(fname,PoolSize,dlm,nlines)
%    This code is modified based on Jun_read_tdf, with the following changes
%          1. use parallel computing (uses parpool and may not be
%          compatible with old matlab versions
%          2. if name of fields have space to start, eliminate such spaces
% 
%    reads a deliminated text file and partitions it to a structure S.
%    
%    functions similar to tdfread function in Matlab, with first line as column header for becoming names of fields for S
%      Compared to tdfread, it has the following differences
%      1. for character strings fields, use cell array rather than
%      character array as storage
%      2. reads the first three data line to determine the nature of the
%      columns; later entried in column is assumed to be the same as the
%      first three data lines
%      3. hopefully is faster on larger files
%
%    Input: 
%      fname: specifies the file name; or specifies an already open file id
%      dlm: deliminator; the default is \t (9)
%      nlines: when specified, read this many data lines from the file in total,
%                nlines need to be >3.
%      PoolSize, number of cpu cores to use; default CPU number or 2, whichever is larger;
%
%
%Jun Lu, Feb 2015
%

if nargin==1
  dlm=9; %means \t
  PoolSize=max(2,feature('numCores'));
elseif nargin==2
    dlm=9;
end

if ischar(fname)
  fid=fopen(fname,'r');
else
  fid=fname;
end

ln=1;
if ~exist('nlines','var')
  nlines=Inf;
  do_close=1;
else
  do_close=0;
end

TestLineCount=3;
DummyFirstLetter='F';

%=======Read Header Line and Prepares for Use as Field Name
headerLine=fgetl(fid);
Fields=dlmsep(headerLine,dlm); % separate line into a cell array of strings
% doublecheck whether headers are good fo field names; if not, modify
longLineCounter=0;
for i=1:length(Fields)
    Fields{i}=regexprep(Fields{i},'^\s+',''); %eliminate beginning spaces
    
    if ~isstrprop(Fields{i}(1),'alpha') % if first character is a letter
        Fields{i}=[DummyFirstLetter Fields{i}];
    end
    
    Fields{i}=regexprep(Fields{i},'\s','_'); %change all white space character into '_'
    Fields{i}=regexprep(Fields{i},'\W','_'); %change any non alphabetic, numeric or underscore into '_'
    
    if length(Fields{i})>namelengthmax
        longLineCounter=longLineCounter+1;
        tempNumStr=num2str(longLineCounter);
        Fields{i}=[Fields{i}(1:(namelengthmax-length(tempNumStr))) tempNumStr];
    end
end

%=======Now read in all lines from file and store in AllData
%read first line
disp('Start Reading File...');

AllData=[];

ln=1;
while(ln<=nlines)
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    
    AllData{ln}=dlmsep(tline,dlm);
    
    ln=ln+1;
    if mod(ln,1000)==0
        fprintf(1,['...' num2str(ln)]);
        if mod(ln, 10000)==0
            fprintf(1,'\n');
        end
    end
end
if do_close
  fclose(fid);
  fid=-1;
end
disp('Finish Reading File.');

%=======Read the rest and parse==================
disp('Start Parsing Lines...');

parpool(PoolSize);

for i=1:length(Fields)
    disp(['parsing field ' num2str(i)]);
    LineIsNum=1;
    for j=1:TestLineCount %determine if a specific field is num
        %if isempty(str2num(AllData{j}{i}))
        if ~all(ismember(AllData{j}{i}, ' 0123456789+-.eEdD'))
            LineIsNum=0;
        end
    end
    
    if ~LineIsNum
        parfor j=1:length(AllData)
            tempArray{j}=AllData{j}{i};
        end
        %S.(Fields{i})=(cellfun(@(x) x{i},AllData,'UniformOutput',false))';
    else
        parfor j=1:length(AllData)
            tempArray(j)=str2num(AllData{j}{i});
        end
        %S.(Fields{i})=(cellfun(@(x) str2num(x{i}),AllData))';
    end
    S.(Fields{i})=tempArray';
    clear tempArray;
end

disp('Completed.');

poolobj = gcp('nocreate');
delete(poolobj);
