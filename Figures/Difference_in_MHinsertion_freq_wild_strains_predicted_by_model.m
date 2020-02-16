% Fig_NaturalIsolates_TandemDuplications_WithWithout_MH
%   what % of >10bp insertions in wild S. pombe isolates have MH? 
FIGDIR  = '~/Nutstore Files/Microhomology shared folder/Figures/Fig4/' ;
FN = [ FIGDIR '57 nature isolates_sequencing summary.xlsx' ] ;

T = readtable( FN , 'Sheet' , 2);
T = T( : , 1:13) ; 
T.Properties.VariableNames{6} = 'REF_N' ; 
T.Properties.VariableNames{7} = 'ALT_N' ; 
T = T(2:end,:);
T.REF_N = str2double(T.REF_N) ; 
T.ALT_N = str2double(T.ALT_N) ; 
T.HasMH = ~cellfun(@isempty,T.MH);

T.REF_len = cellfun(@length,T.REF);
T.ALT_len = cellfun(@length,T.ALT);
T.IsInsertion = T.ALT_len > T.REF_len ; 
T.Location = categorical(T.Location);

N = T ; 
clear 'T' ; 


allins = readtable(  FN , 'Sheet' , 1);

genes_have_an_insert = unique(allins.Gene) ; 
% per gene
DATADIR = '~/SynologyDrive/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/Sarah/MH_project/Manuscript-todo/processeddata/' ;
FIGDIR  = '~/Nutstore Files/Microhomology shared folder/Figures/Fig2 - cis-determinants of MTD through ultra-deep sequencing/' ;
D = readtable([DATADIR 'MHRinGene.txt' ],'TreatAsEmpty','.');
D.Properties.VariableNames = {'chr' 'start1' 'end1' 'type' 'systematic_name' 'N_MHPs'} ;
D.systematic_name = regexprep( regexprep( D.systematic_name , 'ID=','') , ';Name.*' , '') ;

T = readtable([DATADIR 'MHRSumPreinGene.txt' ],'TreatAsEmpty','.');
T.Properties.VariableNames = { 'chr'	'start1'	'end2' 'dunno'	'gene'	'systematic_name'	'sumObs'	'sumPre'	'sumFloat' 'dunno2'}; 
T.GeneLength = abs(T.start1 - T.end2) ./ 1000 ; 
T.SumFloat_N = T.sumFloat ./ T.GeneLength ; 
T.sumObs_N = T.sumObs ./ T.GeneLength ; 
T.sumPre_N = T.sumObs ./ T.GeneLength ; 
T = sortrows(T,'sumFloat','ascend');
T.sumFloat_Rank = (1:height(T))' ;

T = sortrows(T,'SumFloat_N','ascend');
T.sumFloat_N_Rank = (1:height(T))' ;

T = innerjoin(D , T  , 'key' , 'systematic_name'); 

G = T ; 
clear 'T' 'D' ; 


%%
idxN1 =  ~cellfun(@isempty,N.MH) & strcmp(N.REF,'.') ; 
idxN0 =  cellfun(@isempty,N.MH)  ; 

Y  = G.sumFloat ; 
Y(Y>250)=250;
G.FoundInNature = (ismember(G.systematic_name,N.Gene(idxN1)) | ismember(G.gene,N.Gene(idxN1)) ) ;
G.NoMHFoundInNature = (ismember(G.systematic_name,N.Gene(idxN0)) | ismember(G.gene,N.Gene(idxN0)) ) ;
G.AnyInsertFoundInNature = (ismember(G.systematic_name,genes_have_an_insert) | ismember(G.gene,genes_have_an_insert) ) ;

fh = figure('units','centimeters','position',[5 5 8 5]);
hold on ; 
histogram( Y(G.FoundInNature) , 15, 'Normalization','Probability')
histogram( Y(G.AnyInsertFoundInNature) , 15 , 'Normalization','Probability')
lh = legend({'Found in nature' 'not found'}); 
set(lh,'box','off');
xlabel('Predicted MTD score')
ylabel('Fraction of genes')
[~,p] = ttest2( Y(G.FoundInNature)  , Y(G.AnyInsertFoundInNature) ) ; 
[~,p,~] = kstest2( G.sumFloat(G.FoundInNature)  , G.sumFloat(G.AnyInsertFoundInNature) ) ; 
fc = log2(  nanmean(Y(G.FoundInNature))  / nanmean(Y(G.AnyInsertFoundInNature)) ) ; 
title( sprintf('p = %0.05f fold-change = %0.04f' , p , fc) )
print('-dpng','~/Downloads/MH_Mediated_inserts_found_in_nature_vs_PredictorScore', '-r600');

% %%
% 
% N0 = N( ~cellfun(@isempty,N.MH) & N.Location=='reading_frame' , :) ; 
% 
% N1 = innerjoin( N0 , G , 'LeftKey', 'Gene' , 'RightKey', 'systematic_name') ; 
% N2 = innerjoin( N0 , G , 'LeftKey', 'Gene' , 'RightKey', 'gene') ; 
% 
% N1 = sortrows(N1 ,  'sumFloat_Rank' , 'ascend') ; 
% N2 = sortrows(N2 ,  'sumFloat_Rank' , 'ascend') ; 
% 
% 
