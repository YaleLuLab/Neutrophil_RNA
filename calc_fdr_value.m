function fdr=calc_fdr_value(p)
%function fdr=calc_fdr_value(p)
[sp,ord]=sort(p);

fdr=sp*length(p)./((1:length(p))'); % make sure identical values
                                    % are taken care of 

fdr=min([ fdr ones(size(fdr,1),1) ],[],2);

fdr=[fdr; 1];
for i=length(p):-1:1
  fdr(i)=min(fdr(i+1),fdr(i));
end
fdr=fdr(1:(end-1));
fdr(ord)=fdr;
