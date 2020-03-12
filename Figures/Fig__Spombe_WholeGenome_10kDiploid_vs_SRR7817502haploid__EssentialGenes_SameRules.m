% Compare E. coli vs S.pombe haploid vs S.pombe diploid
% we have two high-coverage whole-genome S. pombe datasets: 10k, & SRR7817502
% same MTDfrequency rules? 
% 
% SRR7817502 is liklely haploid. Fewer MTDs in essential genes? 
%
%

%% load data
D = '10k_rm.sign.count.tsv' ;
H = 'SRR7817502_rm.sign.count.tsv' ;
E = 'Ecoli.sign.count.tsv' ;
DATADIR = '~/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/DataFromCluster/' ;

FIGDIR = '~/Nutstore Files/Microhomology shared folder/Figures/Supplementary Figures/' ; 
FIGNAME = [ FIGDIR 'SupFig__Haploid_vs_Diploid_vs_Ecoli_MTDfreq__InterMH_Distance'  ] ; 

vn = {'chr' 's1' 'e1' 's2' 'e2' 'NDup' 'NCol'}; 
vt = {'string' 'uint32' 'uint32' 'uint32' 'uint32' 'uint32' 'uint32'}; 

opts = delimitedTextImportOptions( 'Delimiter','\t', 'VariableNames', vn , 'VariableTypes', vt ) ;
D = readtable( [DATADIR D]  , opts );
H = readtable( [DATADIR H]  , opts );
E = readtable( [DATADIR E]  , opts );

D.strain = repmat( {'diploid'} , height(D),1);
H.strain = repmat( {'haploid'} , height(H),1);
E.strain = repmat( {'E. coli'} , height(E),1);


D.HasDup = D.NDup>0;
D.MHLen = D.e1 - D.s1 + 1 ; 
D.InterMHDist = D.s2 - D.e1 + 1 ; 


H.HasDup = H.NDup>0;
H.MHLen = H.e1 - H.s1 + 1 ; 
H.InterMHDist = H.s2 - H.e1 + 1 ; 

E.HasDup = E.NDup>0;
E.MHLen = E.e1 - E.s1 + 1 ; 
E.InterMHDist = E.s2 - E.e1 + 1 ; 


%%
n=5; 
D.InterMHDistR = round(double(D.InterMHDist)./n)*n ;
H.InterMHDistR = round(double(H.InterMHDist)./n)*n ;
E.InterMHDistR = round(double(E.InterMHDist)./n)*n ;


T = vertcat(D,H);
T = vertcat(T,E);

T.strain = categorical(T.strain); 

T.InterMHDistR(T.InterMHDistR>300)=300; 
T.MHLen(T.MHLen>6)=6; 

G = grpstats( T , {'MHLen' 'InterMHDistR' 'strain'} , {'mean' 'sum'} , 'DataVars' , 'HasDup');
Go=G;

%%
G=Go;
G = G(G.GroupCount>1e3,:);

e = 1e-2 ;
yl = [0.01  1.6] ; 
fh = figure('units','centimeters','position',[5 5  20 7]) ;
t = tiledlayout(1,3 , 'Padding' , 'compact' , 'TileSpacing' , 'compact' );

nexttile;
gh = gscatter(  G.InterMHDistR(G.MHLen==4) , 100*G.mean_HasDup(G.MHLen==4)+e , G.strain(G.MHLen==4) );
arrayfun(@(I)set(gh(I),'LineStyle','-') , 1:numel(gh)) ;
title('MH length = 4')
set(gca,'yscale','log')
ylim(yl)
set(gca,'yticklabel',{'0' '0.1' '1'});

nexttile;
gh =gscatter(  G.InterMHDistR(G.MHLen==5) , 100*G.mean_HasDup(G.MHLen==5)+e , G.strain(G.MHLen==5) ) ;
arrayfun(@(I)set(gh(I),'LineStyle','-') , 1:numel(gh));
title('MH length = 5')
set(gca,'yscale','log')
ylim(yl)
set(gca,'yticklabel',{'0' '0.1' '1'});

nexttile;
gh = gscatter(  G.InterMHDistR(G.MHLen==6) , 100*G.mean_HasDup(G.MHLen==6)+e , G.strain(G.MHLen==6) );
arrayfun(@(I)set(gh(I),'LineStyle','-') , 1:numel(gh));
title('MH length >= 6')
set(gca,'yscale','log')
ylim(yl)
set(gca,'yticklabel',{'0' '0.1' '1'});

xlabel( t , 'Inter MH distance' , 'FontSize' , 15 )
ylabel( t, '% of MHPs with an MTD' , 'FontSize' , 15 )

print( '-dpng' , FIGNAME , '-r200');
close ; 

%% look at <50nt

n=2; 
D.InterMHDistR = round(double(D.InterMHDist)./n)*n ;
H.InterMHDistR = round(double(H.InterMHDist)./n)*n ;
E.InterMHDistR = round(double(E.InterMHDist)./n)*n ;


T = vertcat(D,H);
T = vertcat(T,E);

T.strain = categorical(T.strain); 

T.InterMHDistR(T.InterMHDistR>300)=300; 
T.MHLen(T.MHLen>6)=6; 
T = T(T.InterMHDist<50,:);

G = grpstats( T , {'MHLen' 'InterMHDistR' 'strain'} , {'mean' 'sum'} , 'DataVars' , 'HasDup');
%%

e = 1e-2 ;
yl = [0.01  1.6] ; 
fh = figure('units','centimeters','position',[5 5  20 7]) ;
t = tiledlayout(1,3 , 'Padding' , 'compact' , 'TileSpacing' , 'compact' );

nexttile;
gh = gscatter(  G.InterMHDistR(G.MHLen==4) , 100*G.mean_HasDup(G.MHLen==4)+e , G.strain(G.MHLen==4) );
arrayfun(@(I)set(gh(I),'LineStyle','-') , 1:numel(gh)) ;
title('MH length = 4')
set(gca,'yscale','log')
ylim(yl)
set(gca,'yticklabel',{'0' '0.1' '1'});
xlim([3 50])

nexttile;
gh =gscatter(  G.InterMHDistR(G.MHLen==5) , 100*G.mean_HasDup(G.MHLen==5)+e , G.strain(G.MHLen==5) ) ;
arrayfun(@(I)set(gh(I),'LineStyle','-') , 1:numel(gh));
title('MH length = 5')
set(gca,'yscale','log')
ylim(yl)
set(gca,'yticklabel',{'0' '0.1' '1'});
xlim([3 50])

nexttile;
gh = gscatter(  G.InterMHDistR(G.MHLen==6) , 100*G.mean_HasDup(G.MHLen==6)+e , G.strain(G.MHLen==6) );
arrayfun(@(I)set(gh(I),'LineStyle','-') , 1:numel(gh));
title('MH length >= 6')
set(gca,'yscale','log')
ylim(yl)
set(gca,'yticklabel',{'0' '0.1' '1'});
xlim([3 50])
legend('off')
xlabel( t , 'Inter MH distance' , 'FontSize' , 15 )
ylabel( t, '% of MHPs with an MTD' , 'FontSize' , 15 )

print( '-dpng' , [FIGNAME '_zoom'], '-r200');
close ; 

