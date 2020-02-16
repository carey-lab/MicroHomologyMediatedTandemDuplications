%% ScerHapVSDip__ModelToPredictMHD
% load dup_sites_found____genes_with_essentiality.txt
% & the genome ( to get the sequences )
%
% predict the # of MA lines in which an event was observed, 
%   given MHlen, DistBetween, %GC, expression of gene, etc
%
%% load data
T = readtable('~/Downloads/dup_sites_found____genes_with_essentiality.txt');
T.Properties.VariableNames = {'chr' 'GeneStart' 'GeneEnd' 'type' 'GENE' 'ORF' 'chrx' 's1' 'e2' 'SRA' 'NDupreads' 'e1'} ; 
T.chrx = [] ; 

A = readtable('~/Downloads/Scer_MHRs.txt' , 'ReadVariableNames',false);
A.Properties.VariableNames= {'chr' 's1' 'e1' 's2' 'e2'} ;
A.MHlen = A.e1 - A.s1 ; 

A = outerjoin( A , T(:,{'s1' 'e1' 'e2' 'NDupreads' 'ORF' 'SRA' } ), 'Key' ,{'s1' 'e1' 'e2'} ) ;
A.NDupreads(isnan(A.NDupreads))=0; 
A.HasDup = A.NDupreads>0 ; 

A.MHlen(A.MHlen>=10)=10;
%%
G = grpstats( A , {'MHlen' }, 'mean' , 'DataVars' , {'HasDup' 'NDupreads'} );
%%

plot(  G.MHlen , 100*G.mean_HasDup  )
