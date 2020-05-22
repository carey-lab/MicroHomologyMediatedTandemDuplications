C = readtable('~/SynologyDrive/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/DataFromCluster/wt_polh_subset.MTD.counts.txt');
C.HasMTD = C.Var6 > 0 ; C.Var6 = [] ;
C.HasMTC = C.Var7 > 0 ; C.Var7 = [] ;
C.CHR = regexprep(C.Var1 , 'CHROMOSOME_'  , ''); C.Var1 = [];

R = readtable('/Volumes/MicroscopyAndSequencing/CElegans_MutationalSignatures/mapped_reads_per_experiment.txt' , 'FileType','text');
R.Properties.VariableNames = { 'ReadsMapped'  'FileName'};
R.ID = regexprep(  regexprep( R.FileName , '^.*/' ,'') , '\.[bc][ra].*m' , '');
R = grpstats(R ,'ID' ,'max' , 'DataVars','ReadsMapped') ;
R.MillionsReadsMapped = R.max_ReadsMapped ./ 1e6 ; 

A = readtable('~/CareyLab/ExternalData/Volkova2020/Volkova2020.txt');
A.Var20 = []; A.Var21 = []; A.Var22 = [];A.Var23 = [];
%
A.ID = regexprep( regexprep( regexprep(A.Var10,'.*\/','') ,'\.bam','') ,'\.cram' ,'') ;
for I = 1:height(A)
    if isempty(A.ID{I})
        A.ID{I} = regexprep(  regexprep( regexprep(A.Var11{I},'.*\/','') ,'\.bam',''),'\.cram' ,'') ;
    end
end
A.Var10 = []; A.Var11 = []; A.Var12 = [];A.Var13 = [];A.Var14 = [];

%%
T = outerjoin( C , A , 'LeftKey','Var8' ,'RightKey','ID');
T = T( ~isnan(T.Var2_C),:);


G = grpstats( T , {'Var3_A' 'ID'} , 'sum' , 'DataVars' , {'HasMTD' 'HasMTC'});
G = G( G.sum_HasMTD>0,:);

G = innerjoin( G , R(:,{'ID' 'MillionsReadsMapped'}) , 'Key', 'ID') ;

G.MTDs_per_MillionReads = G.sum_HasMTD ./ G.MillionsReadsMapped ;
G.MTCs_per_MillionReads = G.sum_HasMTC ./ G.MillionsReadsMapped ;
%G.MTDs_per_MillionReads = G.sum_HasMTD ./ 1 ;
%G.MTCs_per_MillionReads = G.sum_HasMTC ./ 1 ;

%Gc = grpstats( grpstats( T , {'Var3_A' 'Var8_C'} , 'sum' , 'DataVars' , {'HasMTD' 'HasMTC'}) , 'Var3_A' ,'mean', 'DataVars','GroupCount');
%Gc.mean_GroupCount = []; 

figure; 
tiledlayout(2,2)

nexttile;
boxplot(G.MTDs_per_MillionReads , G.Var3_A);
ylabel('MT DUPs (per million reads)')
grid on ;
title('MH Duplications')

nexttile;
boxplot(G.MTCs_per_MillionReads , G.Var3_A);
ylabel('MT COLLs (per million reads)')
grid on ;
title('MH Collapses')


nexttile;
boxplot(log2( G.sum_HasMTD ./ G.sum_HasMTC ) , G.Var3_A);
ylabel('log_2( Dup / Collapse )')
grid on ;
title('ratio')


nexttile;
gscatter(G.MTDs_per_MillionReads , G.MTCs_per_MillionReads , G.Var3_A );
ylabel('MT COLLs (per million reads)')
xlabel('MT DUPs(per million reads)')
lh = line(xlim,xlim);set(lh,'HandleVisibility','off')