%% human MTD probability vs MH length
% Figure 2 - show that human behaves somewhat similarly to yeast
%  MHLen vs % MHPs with an MTD
%



%% load S. pombe data
DATADIR = '~/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/Sarah/MH_project/ProcessedData/';
T= readtable( [DATADIR '10k.sign.count.tsv.sorted_full_MHR__interect100pct__intersectANY_intersectEssential_CDS.txt'],'FileType','text','Format','%s%d%d%d%d%d%d%d%d');
T.Properties.VariableNames = {'chr' 's1' 'e2' 'MHLen' 'DupCounts' 'IntersectCDS_100pct' 'IntersectCDS_any' 'IntersectBothEssential100pct'  'IntersectAnyEssential100pct' };
T.HasDup = T.DupCounts > 0 ; 
T.MHLen(T.MHLen>12)=12 ; 
G = grpstats( T , 'MHLen' , {'mean'} , 'DataVars' , 'HasDup');
% Gstd = NaN(height(G),1);
% parfor I = 1:height(G)
%     v = double(T.HasDup( T.MHLen==G.MHLen(I))) ;
%     Gstd(I) = std(bootstrp( 10 , @mean , v )) ;
% end
FIGDIR = '~/Nutstore Files/Microhomology shared folder/Figures/Fig2 - cis-determinants of MTD through ultra-deep sequencing/' ; 
FIGNAME = [FIGDIR 'Human_MHLen_vs_MTDfreq' ] ;

%% load human data

DATADIR = '~/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/DataFromCluster/counts_per_chr/' ;
BEDFILES = dir( [ DATADIR 'chr*bed' ]) ;
G  = cell( numel(BEDFILES) , 1);
parfor chrI = 1:numel(BEDFILES)
    tic; 
    fn = [ BEDFILES(chrI).folder filesep BEDFILES(chrI).name ]  ; 
    T = readtable( fn , 'FileType','text','Delimiter','\t','Format','%s%d%d%d%d');
    T.Properties.VariableNames = {'chr' 's' 'e' 'MHlen' 'NDupReads'} ; 
    T.HasDup = T.NDupReads>0 ; 
    fprintf('%d/%d\tRead %s\t' , chrI , numel(BEDFILES) , BEDFILES(chrI).name); 
    G{chrI} =  grpstats( T , {'chr' 'MHlen' } , 'mean' , 'DataVars' , 'HasDup')   ; 
    fprintf('read & proc took %0.02f min\n' , toc / 60 );
end
clear 'T' ; 


%%
R = table();
for I = 1:numel(BEDFILES)
    R = vertcat(R,G{I});
end
R.MHlen(R.MHlen>12)=12 ; 
save('~/Downloads/counts_per_chr.mat' , 'R' , 'G');
Q = grpstats( R(R.GroupCount>1000,:) , 'MHlen' , {'mean' 'sem' 'std'} , 'DataVars' , 'mean_HasDup'); 

%%

fh = figure('units','centimeters','position',[5 5  4 7]) ;
hold on ;
clrs = get(gca,'ColorOrder');
errorbar( Q.MHlen , 100*Q.mean_mean_HasDup , Q.std_mean_HasDup , 'o-' , 'LineWidth',2 ,'Color',clrs(2,:))
plot( G.MHLen , 100*G.mean_HasDup ,'o-' , 'LineWidth',2 ,'Color',clrs(1,:))

xlim([3.5 12.5])
ylabel('% of MHPs with an MTD')
set(gca,'yscale','log')
set(gca,'xtick',4:2:20)
xlabel('MH length (nt)')
print('-depsc2' , FIGNAME );
close all ; 