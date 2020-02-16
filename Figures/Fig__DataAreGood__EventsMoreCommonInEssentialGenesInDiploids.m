% 2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/ScerHaploidVSDiploid_PRJNA449976
DATADIR = '~/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/DataFromCluster/' ; 
FIGUREDIR = '~/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/FIGURES/' ; 

T = readtable([ DATADIR 'Scer_HapVSDip__dup_sites_found____genes_with_essentiality.txt' ], 'ReadVariableNames' ,false,'Delimiter','\t');
%T = readtable([ DATADIR 'Scer_HapVSDip__dup_sites_found____genes_with_essentiality_overlap=1.0.txt' ], 'ReadVariableNames' ,false,'Delimiter','\t');
A = readtable([ DATADIR 'Scer_HapVSDip__MetaData_Ploidy_and_MAlines.csv' ]);

haploid_SRAs = A.SRA( strcmp(A.Ploidy , 'hap')) ; 
diploid_SRAs = A.SRA( strcmp(A.Ploidy , 'dip')) ; 

hap_idx = ismember(T.Var10 , haploid_SRAs); 
dip_idx = ismember(T.Var10 , diploid_SRAs); 
T.haploid = hap_idx ; 
T.intergenic = strcmp(T.Var4,'intergenic'); 
T.SizeDivByThree = mod(T.Var9-T.Var8,3)==0 ; 
T.Desc = T.Var4; 
T.Var4 = [] ;
T.Essential = regexpcmp(T.Desc , '^Essential gene');
T.NonEssentialGene = ~T.Essential & ~T.intergenic ; 

% BS = readtable('bamstats.tab','FileType','text');
% BS.Properties.VariableNames = {'SRA' 'Nreadspaired'} ;
% BS = innerjoin(A(:,{'Ploidy' 'SRA'}),BS,'Key','SRA');
% BSG = grpstats( BS , 'Ploidy' , 'median' , 'DataVars' , 'Nreadspaired')
% grpstats( grpstats(T,'Var10','mean','datavars','haploid') , 'mean_haploid' ,'mean' ,'datavars','GroupCount')
%
% statistics for number of times each found MHR Dup was found in each set
% of genes + experiments
G = grpstats( T , {'haploid'  'intergenic' 'Essential'} , {'mean' 'sum'} , 'DataVars' , {'Var11'});
G.NonEssentialGene = ~G.Essential & ~G.intergenic ; 

fh = figure('units','centimeters','position',[5 5 15 15]) ;
clrs = get(gca,'ColorOrder');
subplot(2,2,1)
txt = arrayfun(@(X)sprintf('%d',X),G.GroupCount(G.haploid),'UniformOutput',false) ; 
ph = pie( G.GroupCount(G.haploid) , txt );
set(ph(1),'FaceColor',clrs(1,:)) ; 
set(ph(3),'FaceColor',clrs(2,:)) ; 
set(ph(5),'FaceColor',clrs(3,:)) ;
title('haploid')
subplot(2,2,2)
txt = arrayfun(@(X)sprintf('%d',X),G.GroupCount(~G.haploid),'UniformOutput',false) ; 
ph = pie( G.GroupCount(~G.haploid) , txt );
set(ph(1),'FaceColor',clrs(1,:)) ; 
set(ph(3),'FaceColor',clrs(2,:)) ; 
set(ph(5),'FaceColor',clrs(3,:)) ;
title('diploid')

legend({'intergenic' 'essential'  'non-essential' } , 'location','East')


%subplot(3,2,3)
%x = [ G.GroupCount(G.haploid & G.intergenic) G.GroupCount(~G.haploid & G.intergenic) ; ...
%      G.GroupCount(G.haploid & G.Essential) G.GroupCount(~G.haploid & G.Essential) ; ...
%      G.GroupCount(G.haploid & G.NonEssentialGene) G.GroupCount(~G.haploid & G.NonEssentialGene)  ] ; 
%bar(  (x./sum(x))'  , 'stacked')
%set(gca,'xticklabel',{'Haploid' 'Diploid'})
%ylabel('Fraction of Duplication loci')
%xlim([0.4 2.6])

x = [ G.GroupCount(G.haploid & G.Essential) sum(G.GroupCount(G.haploid & ~G.Essential)) ; ...
    G.GroupCount(~G.haploid & G.Essential) sum(G.GroupCount(~G.haploid & ~G.Essential))]
[~,p,s]=fishertest( flipud( x ) )

x = [ G.GroupCount(G.haploid & G.Essential) sum(G.GroupCount(G.haploid & G.NonEssentialGene)) ; ...
    G.GroupCount(~G.haploid & G.Essential) sum(G.GroupCount(~G.haploid & G.NonEssentialGene))]
[~,p,s]=fishertest( flipud( x ) )



subplot(2,2,3); hold on; 
bh = bar( 0.5 , x(1,1) , 'FaceColor',clrs(4,:));
bh = bar( 1.5 , x(2,1) , 'FaceColor',clrs(5,:));
bh = bar( 3 , x(1,2) , 'FaceColor',clrs(4,:));
bh = bar( 4 , x(2,2) , 'FaceColor',clrs(5,:));
set(gca,'yscale','log')
set(gca,'ytick',[1 2 5 10 25 50 100 250 500 1000])
ylim([ 0.5 max(x(:))])
set(gca,'xtick',[1 3.5])
set(gca,'xticklabel', {'Essential' 'Non-Essential'});
legend( {'haploid' 'diploid'},'location','nw')
ylabel('# of MHRs w/a Duplication')


% bootstrap loci locations
subplot(2,2,4); hold on; 
hap_fraction_in_essential = bootstrp( 10000 , @mean , T.Essential(T.haploid )) ; 
dip_fraction_in_essential = bootstrp( 10000 , @mean , T.Essential(~T.haploid)) ; 
Y = 100 * [ (hap_fraction_in_essential) (dip_fraction_in_essential)] ; 
bh = bar( mean(Y) , 'FaceColor', [.7 .7 .7] ) ;
errorbar( bh.XData , mean(Y) , std(Y) , 'ok' , 'LineWidth',2); 
ylabel('% of Duplications that are in essential genes')
xlim([0.4 2.6])
set(gca,'xtick',1:2)
set(gca,'ytick',0:10)
set(gca,'xticklabel', {'Haploid' 'Diploid'});


print('-dpng', [ FIGUREDIR 'Scer_HapVSDip__dup_sites_found____genes_with_essentiality' ] , '-r300');





%G = sortrows(G,{'haploid' 'Essential' 'GroupCount'}  );
%G.Properties.RowNames = strcat(G.Var5,'-',num2str(G.haploid) , '-' , num2str(G.intergenic)); 


idx_E = G.Essential ; % & ~G.SizeDivByThree ;
idx_notE = ~G.Essential ; % & ~G.SizeDivByThree ;

x = [sum(G.GroupCount(idx_E & ~G.haploid)) sum(G.GroupCount(idx_notE & ~G.haploid)) ;...
    sum(G.GroupCount(idx_E & G.haploid)) sum(G.GroupCount(idx_notE & G.haploid)) ]

[~,p,s]=fishertest( x) 

pct_in_essential = x(:,1) ./ sum(x,2) * 100 


% %%
% load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features_new.mat');
% idx = find(strcmp(DS.TYPE , 'ORF')); D = DS(idx , :);
% D = D(: , {'ORF' , 'dist_to_the_end_kb' , 'percent_underreplicated_cdc20' , ...
%     'percent_underreplicated_dbf2' , 'Phenotype_Observable' , 'Essentiality'});
% T = dataset('file' , '~/Develop/Mendoza__ReplicationEvolution/Data/ExternalData/Baryshnikova10.tab');
% T = T(: , {'ORF' , 'fitness'}); T = unique(T);
% D = join(D , T , 'Type' , 'left' , 'Keys' , 'ORF' , 'MergeKeys' , true);
% D = dataset2table(D);
% 
% Q=innerjoin(G,D(:,{'ORF' 'Essentiality' 'fitness'}) , 'RightKey','ORF','LeftKey','Var6') ;
% 
% %%
% 
% idx_bad = Q.Essentiality & ~Q.SizeDivByThree ;
% x = [sum(idx_bad & ~Q.haploid) sum(~idx_bad & ~Q.haploid) ;...
%     sum(idx_bad & Q.haploid) sum(~idx_bad & Q.haploid) ]
% 
% [~,p,s]=fishertest( x) 
% 
% pct_in_essential = x(:,1) ./ sum(x,2) * 100 