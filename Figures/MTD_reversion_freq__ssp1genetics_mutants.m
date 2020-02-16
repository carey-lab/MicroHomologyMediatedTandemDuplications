% load data
SOURCEDATA = '/Users/lcarey/Nutstore Files/Microhomology shared folder/Figures/supplemental data/reversion frequency.xlsx' ;
T = readtable( SOURCEDATA , 'Sheet' , 2);

V = stack(T , 1:size(T,2));
V.Properties.VariableNames = {'Genotype' , 'MTDReversionRate'};
V.Genotype = string(V.Genotype) ;
V.Genotype = regexprep( V.Genotype , '_d' , '\\Delta') ; 
V.Genotype = regexprep( V.Genotype , '_' , ',') ; 
V = V( ~isnan(V.MTDReversionRate),:);
G = grpstats( V , 'Genotype' , {'mean' 'sem'  'std' 'median'}) ;


%%
G = sortrows(G,'median_MTDReversionRate');
fh = figure('units','centimeters','position',[5 5 70 8]);
hold on ;
bh = bar(G.mean_MTDReversionRate , 'FaceColor' , [.85 .85 .85] );

yWT = V.MTDReversionRate(strcmp(V.Genotype,'wt'));
idxWT = find(strcmp(G.Genotype,'wt'));
bh = bar(idxWT , G.mean_MTDReversionRate(idxWT) , 'FaceColor' , [.5 .5 .5] );
for I = 1:height(G)
    y = V.MTDReversionRate(strcmp(V.Genotype,G.Genotype{I}));
    [~,p] = ttest2( y , yWT);
    if p>0.01 & p<0.05
        text( I , 0.0075 , '*')
    elseif p>0.001 & p<0.05
        text( I , 0.0075 , '**')
    elseif p<0.05
        text( I , 0.0075 , '***')
    end
end
        
 errorbar( 1:height(G) , G.mean_MTDReversionRate , G.sem_MTDReversionRate , '.k','LineWidth',3);

set(gca,'xtick',1:height(G));
set(gca,'xticklabels',G.Genotype)
ylabel('MTD Reversion Rate')
set(gca,'yscale','log')

for I = 1:height(G)
    y = V.MTDReversionRate( strcmp(V.Genotype , G.Genotype{I})) ; 
    plot( linspace(I-0.2 , I+0.2 , numel(y))  ,  y  , 'ok');
end

for I = 1:height(G)
    y = V.MTDReversionRate( strcmp(V.Genotype , G.Genotype{I})) ; 
    plot( linspace(I-0.2 , I+0.2 , numel(y))  ,  y  , 'ok');
end

yl = G.mean_MTDReversionRate(strcmp(G.Genotype,'wt')) - 2*G.std_MTDReversionRate(strcmp(G.Genotype,'wt')) ; 
yh = G.mean_MTDReversionRate(strcmp(G.Genotype,'wt')) + 2*G.std_MTDReversionRate(strcmp(G.Genotype,'wt')) ; 
%line( xlim , [yl yl],'LineStyle','--')
%line( xlim , [yh yh],'LineStyle','--')

print('-dpng','~/Downloads/bars_with_points.png','-r500')
close ; 

%% for cytoscape
