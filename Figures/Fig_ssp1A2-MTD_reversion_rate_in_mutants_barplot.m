% Fig_ssp1A2-MTD_reversion_rate_in_mutants_barplot.m
% barplot of MTD reversion rates for the ssp1-A2-MTD

%% load data and sort by MTD reversion rate
XLS = '~/Nutstore Files/Microhomology shared folder/Figures/supplemental data/reversion frequency.xlsx';
T = readtable( XLS , 'Sheet' , 'MTD reversion raw data');
datamean = nanmean(table2array(T));
[~,o] = sort( datamean , 'ascend');
T = T(: , o);

%%
data = (table2array(T));
vn = T.Properties.VariableNames ; 
vn = regexprep( vn , '_d','\\Delta');
vn = regexprep( vn , '_',' ');
wtIDX = find(strcmp(vn,'wt'));
wtDATA = data( ~isnan(data( : , wtIDX )) , wtIDX );

sem = arrayfun(@(I)nansem(data(:,I)),1:size(data,2)) ; 

fh = figure('units','centimeters','position',[5 5  35 6]) ;
hold on ;
bar( nanmean(data) ,'FaceColor',[.8 .8 .8] );
bar( wtIDX , nanmean(data(:,wtIDX)) ,'FaceColor',[.4 .4 .4] );
errorbar( 1:size(data,2) , nanmean(data) , sem ,'.k','LineWidth',2)

% single data & p-value ***s
for I = 1:size(data,2)
    d = data( ~isnan(data(:,I)) , I ); 
    xl = linspace(I-0.2 , I+0.2 , numel(d) );
    plot( xl , d ,'ok' ,'MarkerSize',4)
    [~,p] = ttest2( d , wtDATA );
    if p<0.001
        text( I , 0.008 , '***');
    elseif p<0.01
        text( I , 0.008 , '**');
    elseif p<0.05
        text( I , 0.008 , '**');
    end
end

set(gca,'yscale','log')
ylabel('MTD revertant frequency')
xlim([0.4 size(data,2)+0.6])
set(gca,'xtick',1:size(data,2),'xticklabels',vn);
fix_xticklabels(gca)

%print('-dpng','~/Downloads/Fig_ssp1A2-MTD_reversion_rate_in_mutants_barplot','-r600');
%close ; 



%%
data = (table2array(T));
wtDATA = data( ~isnan(data( : , wtIDX )) , wtIDX );

data =  log2(data ./ mean(wtDATA) ) ; 
sem = arrayfun(@(I)nansem(data(:,I)),1:size(data,2)) ; 

fh = figure('units','centimeters','position',[5 5  30 7]) ;
hold on ;
bar( nanmean(data) ,'FaceColor',[.8 .8 .8] );
bar( wtIDX , nanmean(data(:,wtIDX)) ,'FaceColor',[.4 .4 .4] );
errorbar( 1:size(data,2) , nanmean(data) , sem ,'.k','LineWidth',2)

% single data & p-value ***s
for I = 1:size(data,2)
    d = data( ~isnan(data(:,I)) , I ); 
    xl = linspace(I-0.2 , I+0.2 , numel(d) );
    plot( xl , d ,'ok' ,'MarkerSize',4)
    [~,p] = ttest2( d , wtDATA );
    if p<0.001
        text( I , 5 , '***');
    elseif p<0.01
        text( I , 5 , '**');
    elseif p<0.05
        text( I , 5 , '**');
    end
end

ylabel('MTD reversion rate (log_2(mutant/WT))')
xlim([0.4 size(data,2)+0.6])
set(gca,'xtick',1:size(data,2),'xticklabels',vn);
fix_xticklabels(gca)

print('-dpng','~/Downloads/Fig_ssp1A2-MTD_reversion_rate_in_mutants_barplot_log2','-r600');
close ; 

%% epistasis

mat = table2array(T);

R = table();
R.strain = T.Properties.VariableNames';
R.rev_frq = nanmedian( table2array(T))' ; 

for I = 1:numel(vn)
    v = log(mat(:,1) ./ mat(:,strcmp(R.strain,'wt'))' );
    v = v(~isnan(v)) ; 
    R.sem(I) = sem(v);
    R.std(I) = std(v);
end



R.rel_to_wt_log2R = log2( R.rev_frq ./ R.rev_frq(strcmp(R.strain,'wt')));
R.rel_to_mrc1_log2R = log2( R.rev_frq ./ R.rev_frq(strcmp(R.strain,'mrc1_d')));
R.rel_to_wt_diff = ( R.rev_frq - R.rev_frq(strcmp(R.strain,'wt')));
R.rel_to_mrc1_diff = ( R.rev_frq - R.rev_frq(strcmp(R.strain,'mrc1_d')));
R.rel_to_wt_R = ( R.rev_frq ./ R.rev_frq(strcmp(R.strain,'wt')));
R.rel_to_mrc1_R = ( R.rev_frq ./ R.rev_frq(strcmp(R.strain,'mrc1_d')));

genes = { 'rad50_d' 'pds5_d'  'rik1_d' 'rad2_d' 'cds1_d' 'ssb3_d'} ;

fh = figure('units','centimeters','position',[5 5 8 8]);
hold on ; 
for I = 1:numel(genes)
    idx = strcmp(R.strain,genes{I}) ; 
    s = R.rel_to_wt_log2R(idx);
    m = R.rel_to_wt_log2R(strcmp(R.strain,'mrc1_d'));
    d = R.rel_to_wt_log2R(strcmp(R.strain,[genes{I} '_mrc1_d']));
    e = s+m ; 
    if d>3
        lw=1;
    else
        lw=3;
    end
    ph = errorbar(e,d,R.std(idx),'o' ,'LineWidth',lw,'DisplayName',  regexprep(genes{I},'_d','')  ) ; 
end
lh = line(xlim,xlim,'Color','k','LineStyle','--');
lh.HandleVisibility = 'off';
xlabel('Expected fold change (single + mrc1\Delta)')
ylabel('Measured fold change (double mutant)')

y = rev_frq(strcmp(vn,'mrc1_d')) ; 
%lh = plot( m ,m , 'ok' ,'MarkerFaceColor',[.7 .7 .7]) ; lh.HandleVisibility = 'off';
%lh = line(xlim,[m m],'Color',[.7 .7 .7]) ; lh.HandleVisibility = 'off';
rectangle( 'Position' , [ min(xlim) m-0.5 range(xlim) 1] ,'EdgeColor',[.7 .7 .7])
%lh = line([m m],ylim,'Color',[.7 .7 .7]) ; lh.HandleVisibility = 'off';
legend('location','best','box','off')
ylim([-0.9 5.8])

print('-dpng' ,'~/Downloads/Fig_ssp1A2-MTD_reversion_rate__epistasis','-r600');

