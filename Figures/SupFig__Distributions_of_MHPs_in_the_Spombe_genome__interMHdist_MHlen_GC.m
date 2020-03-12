% SupFig__Distributions_of_MHPs_in_the_Spombe_genome__interMHdist_MHlen_GC
% LLD : It may be useful to have a histogram of predicted MTDs used for sequencing read mapping. 
% LBC Feb 2020

% load data
FIGDIR = '/Users/lcarey/Nutstore Files/Microhomology shared folder/Figures/Fig2 - cis-determinants of MTD through ultra-deep sequencing/' ;
FIGNAME = [FIGDIR 'SupFig__Distributions_of_MHPs_in_the_Spombe_genome__interMHdist_MHlen_GC'];

%DATADIR = '~/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/DataFromCluster/' ;

DATADIR = '/Users/lcarey/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/Sarah/MH_project/ProcessedData/';
T= readtable( [DATADIR '10k.sign.count.tsv.sorted_full_MHR__interect100pct__intersectANY_intersectEssential_CDS.txt'],'FileType','text','Format','%s%d%d%d%d%d%d%d%d');

T.Properties.VariableNames = {'chr' 's1' 'e2' 'MHLen' 'DupCounts' 'IntersectCDS_100pct' 'IntersectCDS_any' 'IntersectBothEssential100pct'  'IntersectAnyEssential100pct' };


T.InterMHDist = (T.e2 - T.s1 + 1);  
T.DuplicatatedDistance = T.InterMHDist - T.MHLen ;  
 

T.HasDup = T.DupCounts > 0 ; 

T.HasIntersect100pct = T.IntersectCDS_100pct > 0 ;
T.HasIntersectAny = T.IntersectCDS_any > 0 ;

T.MHLenR2 = T.MHLen ; 
T.MHLenR2( T.MHLenR2 >17 ) = 17  ;

T.IntersectNonessentialCDS = T.IntersectAnyEssential100pct==0 &  T.IntersectCDS_100pct>0; 
%% calc set categories for the legend

T.MHLenR = T.MHLen ; 
T.MHLenR( T.MHLenR >7 ) = 7  ;
T.MHLenR = categorical(T.MHLenR);
T.MHLenR(T.MHLenR=='4')='MH len = 4' ; 
T.MHLenR(T.MHLenR=='5')='MH len = 5' ; 
T.MHLenR(T.MHLenR=='6')='MH len = 6' ; 
T.MHLenR(T.MHLenR=='7')='MH len >= 7' ; 
%T.MHLenR(T.MHLenR=='7')='MH len = 7' ; 
%T.MHLenR(T.MHLenR=='8')='MH len = 8' ; 
%T.MHLenR(T.MHLenR=='9')='MH len >= 9' ; 
T.MHLenR=categorical(string(T.MHLenR));
%%

%
%G = grpstats(T,{'DuplicatatedDistance' 'MHLenR'},{'sum'},'DataVars','HasDup');

G = grpstats(T,{'DuplicatatedDistance' 'MHLenR' 'HasIntersect100pct' },{'mean' 'sum'},'DataVars','HasDup');

% G = grpstats(T,{'MHLenR2' 'HasIntersect100pct' },{'mean' 'sum'},'DataVars','HasDup');


%% depletion of MHPs likely to cause an MTD


R = table();
R.MHL = unique(T.MHLenR2) ; 

for I = 1:height(R)
    MHLR = R.MHL(I) ; 
    [a,b]=count_unique(mod(T.DuplicatatedDistance(T.MHLenR2 == MHLR & T.IntersectBothEssential100pct),3));
    v = 100* b(a==0) ./ sum(b) ; 
    R.InEssentialGene(I) = v; 

    [a,b]=count_unique(mod(T.DuplicatatedDistance(T.MHLenR2 == MHLR & T.IntersectNonessentialCDS),3));
    v = 100* b(a==0) ./ sum(b) ; 
    R.InNonEssentialGene(I) = v;     
    
    [a,b]=count_unique(mod(T.DuplicatatedDistance(T.MHLenR2 == MHLR & ~T.HasIntersect100pct),3));
    v = 100* b(a==0) ./ sum(b) ; 
    R.out_of_gene(I) = v; 
end
%%
fh = figure('units','centimeters','position',[5 5 25 7]) ; 
clrs = get(gca,'ColorOrder');
clrs(3,:) = [0.3 0.3 0.3] ; 
%clrs(3,:) = [0.5 0 0] ; 

xl = [142.5 155.5] ;
yl = [3e2 4e4] ;
t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');


nexttile
idx = G.DuplicatatedDistance>130 & G.DuplicatatedDistance<160    &  ~G.HasIntersect100pct  ;
gh = gscatter( G.DuplicatatedDistance(idx), G.GroupCount(idx) ,  G.MHLenR(idx) , clrs);
arrayfun(@(I)set(gh(I),'LineStyle','-') , 1:numel(gh))
set(gca,'yscale','log')
%axis tight; 
ylim(yl) ;
set(gca,'xtick',3:3:1e4) ;
xlim(xl)
ylabel('# of MH Pairs')
set(gca,'XGrid','on')
xlabel('BP that would be duplicated by an MTD');
lh = legend('location','nw');
title('intergenic MHPs')

nexttile
idx = G.DuplicatatedDistance>130 & G.DuplicatatedDistance<160    &  G.HasIntersect100pct  ;
gh = gscatter( G.DuplicatatedDistance(idx), G.GroupCount(idx) ,  G.MHLenR(idx) , clrs);
arrayfun(@(I)set(gh(I),'LineStyle','-') , 1:numel(gh))
set(gca,'yscale','log')
%axis tight; 
ylim(yl) ;
set(gca,'xtick',3:3:1e4) ;
xlim(xl)
ylabel('# of MH Pairs')
set(gca,'XGrid','on')
xlabel('BP that would be duplicated by an MTD')
legend('off')
title('within gene MHPs')
% enrichment with increasing MHLen ? 


%
nexttile
hold on ; 
clrs = cbrewer( 'qual' , 'Set1' , 6) ;
plot( R.MHL , R.InEssentialGene    , '.-k' , 'MarkerFaceColor','k'  , 'DisplayName' , 'in essential genes' ,'LineWidth',2,'MarkerSize',15)
plot( R.MHL , R.InNonEssentialGene , '.-' , 'MarkerFaceColor',clrs(1,:) ,'Color',clrs(1,:) , 'DisplayName' , 'in non-essential'  ,'LineWidth',2,'MarkerSize',15)
plot( R.MHL , R.out_of_gene , '.-' ,  'MarkerFaceColor',clrs(3,:) ,'Color','b' , 'DisplayName' , 'integenic MHPs' ,'LineWidth',2,'MarkerSize',15)

ylabel('% of MHPs that are in-frame')
xlabel('MH length (nt)')
xlim([4,11.1])
ylim([30 80])
lh = line( xlim , [33 33] , 'LineStyle','--','Color', clrs(5,:) ,'DisplayName','random expectation','LineWidth',2);
title('MHPs: intergenic  vs  within gene')
legend('location','nw')


print('-dpng',FIGNAME,'-r300');
%%
legend('location','EastOutside')
print('-dpng',[FIGNAME '_legend'],'-r300');

