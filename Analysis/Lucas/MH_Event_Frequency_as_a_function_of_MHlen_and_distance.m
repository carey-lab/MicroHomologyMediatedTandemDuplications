function MH_Event_Frequency_as_a_function_of_MHlen_and_distance( FILENAME_sign_count_tsv )
%% MH_Event_Frequency_as_a_function_of_MHlen_and_distance( FILENAME_sign_count_tsv )
%   calculate and plot duplication & collapse frequencies as a function of
%    the length of the MH repeats, and as a function of the distance between repeats
%
%    Input : the output of catch_signatures.awk
%       I	4048	4051	4060	4063	0	5050
%       I	4054	4057	4066	4069	3	1
%       I	4048	4051	4066	4069	0	2
%      chr   st      end     st      en   iCount  dCount
%    iCount = Insertions      dCount = Deletions / Collapse
%
%    output is a pair of tables with the event frequencies calculated,
%    grouped by MHlen or by inter-MH distance
% 
% October 2019, LBC
%% load data
if ~exist('FILENAME_sign_count_tsv' , 'var')
    FILENAME_sign_count_tsv = '~/Downloads/alltsvs.txt' ; 
end
XLS_FILE_NAME = '~/Downloads/MH_Experiment4.xlsx'; delete(XLS_FILE_NAME);


T = readtable( FILENAME_sign_count_tsv , 'FileType','text','Format','%s%s%d%d%d%d%d%d%f%d%f' , 'TreatAsEmpty','-' );

T.Properties.VariableNames = {'lib' 'chr' 's1' 'e1' 's2' 'e2' 'ReadDepth' 'DupCounts' 'DupFreq' 'CollapseCounts' 'CollapseFreq' };
T.HasDup = T.DupCounts > 0 ; 
T.HasColl = T.CollapseCounts > 0 ; 
%T.lib = regexprep( T.lib , '.sign.norm.tsv' ,'');
T.lib = regexprep( T.lib , '^lib_' ,'');
T.lib = str2double(T.lib) ; 

T.MHlen = T.e1-T.s1+1 ;
T.TotalLen = T.e2 - T.s1 ; 
T.NTBetweenRepeats = T.s2 - T.e1 ; 
%T.chr = categorical(T.chr) ; 
T.s1 = uint32(T.s1) ; T.e1 = uint32(T.e1) ; T.s2 = uint32(T.s2) ; T.e2 = uint32(T.e2) ; 
% get rid of NaNs
T.DupCounts( isnan(T.DupCounts) ) = 0 ; 
T.DupFreq( isnan(T.DupFreq) ) = 0 ; 
T.CollapseCounts( isnan(T.CollapseCounts) ) = 0 ; 
T.CollapseFreq( isnan(T.CollapseFreq) ) = 0 ; 

T.MHlen(T.MHlen>10) = 11 ; 
T.e1 = [] ; T2.e2 = [] ; % after we calculate MHlen we don't need these columns
T = sortrows(T ,{'lib' 'chr' 's1' 's2'}, 'ascend');
writetable( T , XLS_FILE_NAME , 'Sheet' ,'raw' );

Q = sortrows( grpstats(T,{'chr' 'lib'},'max','datavars','s2') , 'max_s2' ) ;

G = grpstats( T ,{ 'lib' 'chr' 'MHlen'} , {'mean' 'median' 'max'}  , 'DataVars' , {'ReadDepth' 'DupFreq' 'CollapseFreq' 'HasDup' 'HasColl'  's2' } );

% remove summary stats that aren't interesting;
G.median_HasDup = [] ; G.median_HasColl=[];
G.max_HasColl=[] ; G.max_HasDup=[]; 
G.max_DupFreq=[] ; G.max_CollapseFreq=[];
G.mean_s2=[] ; G.median_s2=[];
G.max_ReadDepth=[] ; G.mean_ReadDepth=[];
G.ChemSyn  = ismember(G.chr , Q.chr(Q.max_s2<600)) ;
% remove if < 10
%G = G( G.GroupCount >= 10 , :);

% normalize by MHL
G.mean_DupFreq_norm_by_MHL = NaN(height(G),1);
G.mean_CollapseFreq_norm_by_MHL = NaN(height(G),1);
for I = 1:height(G)
    G.mean_DupFreq_norm_by_MHL(I) = log2( G.mean_DupFreq(I) ./ median(G.mean_DupFreq(G.MHlen==G.MHlen(I))));
    G.mean_CollapseFreq_norm_by_MHL(I) = log2( G.mean_CollapseFreq(I) ./ median(G.mean_CollapseFreq(G.MHlen==G.MHlen(I))));
end
G.mean_DupFreq_norm_by_MHL(G.mean_DupFreq_norm_by_MHL < -5) = -5 ; % deals with -Inf (log2(zero))
G.mean_CollapseFreq_norm_by_MHL(G.mean_CollapseFreq_norm_by_MHL < -5) = -5 ; % deals with -Inf (log2(zero))

G = sortrows(G ,{'lib' 'chr' 'MHlen'} , 'ascend');

writetable( G , XLS_FILE_NAME , 'Sheet' ,'grouped' );


%% now we can plot individual experiments


%G.mean_CollapseFreq( G.mean_CollapseFreq>50 ) = 50 ; 
idx = regexpcmp(G.chr,'ssp1') & G.max_s2>200 & ~regexpcmp(G.chr,'PCR') & G.median_ReadDepth>1000 ;
lbl = strcat(string(G.lib)  , "--" , string(G.chr) ) ; 
G = sortrows(G,'mean_CollapseFreq_norm_by_MHL','ascend');
clrgrp1 = cellfun(@(X)~isempty(regexp(X,'mrc1')),lbl(idx))  ; 
clrgrp2 = cellfun(@(X)~isempty(regexp(X,'rad52')),lbl(idx))  ; 
clrgrp3 = cellfun(@(X)~isempty(regexp(X,'rad50')),lbl(idx))  ; 
clrgrp = zeros(sum(idx),1); 
clrgrp(clrgrp1)=1;clrgrp(clrgrp2)=2;clrgrp(clrgrp3)=3;

Y =  G.mean_HasDup(idx ) ; 
%Y(Y>3)=3;
figure;
bh = boxplot( Y , { lbl(idx ) }  , 'Colors' , 'rgbm' , 'ColorGroup' , clrgrp );
xtickangle(45)
ylabel('Duplication Freq (normalized by MHlen)');
line( xlim , [0 0 ] ,'LineStyle','--')

%%
fh = figure('units','centimeters','position',[5 5 5 5]) ;
Y =  log2( double(1+T.CollapseCounts( strcmp(T.chr,'mrc1.ssp1.dup.3'))) ./ double( 1+T.CollapseCounts(strcmp(T.chr , 'ssp1.dup.3')))) ; 
Y(Y<-4)=-4; Y(Y>4)=4;
ecdf( Y ); grid on 
xlabel('WT more    mrc1 more')
set(gca,'ytick',0:.1:1);
grid on ; 

fh = figure('units','centimeters','position',[5 5 5 5]) ;
Y =  log2( double(1+T.CollapseCounts( strcmp(T.chr,'mrc1.ssp1.dup.2'))) ./ double( 1+T.CollapseCounts(strcmp(T.chr , 'ssp1.dup.2')))) ; 
Y(Y<-4)=-4; Y(Y>4)=4;
ecdf( Y ); grid on 
xlabel('WT more    mrc1 more')
set(gca,'ytick',0:.1:1);
grid on ; 

fh = figure('units','centimeters','position',[5 5 5 5]) ;
Y =  log2( double(1+T.CollapseCounts( strcmp(T.chr,'mrc1.ssp1.dup.4'))) ./ double( 1+T.CollapseCounts(strcmp(T.chr , 'ssp1.dup.1')))) ; 
Y(Y<-4)=-4; Y(Y>4)=4;
ecdf( Y ); grid on 
xlabel('WT more    mrc1 more')
set(gca,'ytick',0:.1:1);
grid on ; 

%% group by MH-length ,and by distance-between-repeats, for calculating statistics
%G.mean_DupFreq(G.mean_DupFreq>200) = 200 ; 
%G.mean_CollapseFreq(G.mean_CollapseFreq>50) = 50 ; 

%%
idx = regexpcmp( G.chr , 'fragment')  ;
fh = figure('units','centimeters','position',[5 5 40 30]) ;
    Y = G.mean_DupFreq_norm_by_MHL; Y(Y>35)=35;
gh = gscatter( G.MHlen(idx) , Y(idx) , {G.lib(idx) , G.chr(idx)} , ['bbrrmmm'] )
ylabel('Collapse freq')
xlabel('MH length')
set(gca,'yscale','lin')
%%
ul = unique(G.lib);
fh = figure('units','centimeters','position',[5 5 40 30]) ;
for I = 1:numel(ul)
    Q = G( strcmp(G.lib, ul{I}) , :);
    subplot(5,5,I)
    gscatter( Q.MHlen , Q.mean_CollapseFreq+0.01 , Q.chr)
    legend('location','NorthOutside')
    set(gca,'yscale','log')
end
