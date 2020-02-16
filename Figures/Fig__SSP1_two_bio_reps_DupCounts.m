%% Figure to plot the frequency of duplication events across ssp1

% To identify the complete spectrum of MTD events in a single gene we performed 106x coverage sequencing of ssp1,
% identifying many MTDs that are supported dozens of overlapping sequencing reads (Fig. @B).
% We recovered all three MTD events found in the genetic screen, as well as @ additional events 
% spread across the gene (Fig.@C,D). 

%% load data
DATADIR = '~/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/DataFromCluster/' ; 
FIGUREDIR = '~/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/FIGURES/' ; 

T = readtable( [ DATADIR  'PombeAmpliconSeq_E2_alltsvs.txt'] ,'TreatAsEmpty','-'); 
T.Properties.VariableNames = {'lib' 'chr' 's1' 'e1' 's2' 'e2' 'ReadDepth' 'DupCounts' 'DupFreq' 'CollapseCounts' 'CollapseFreq' };
T.DupCounts(isnan(T.DupCounts))=0 ;
T.DupFreq(isnan(T.DupFreq))=0 ;
T.HasDup = T.DupCounts>0 ;
T.MHlen = T.e1-T.s1+1 ;

idx_to_keep = strcmp( T.chr , 'ssp1') ;
T = T( idx_to_keep , :) ; 
T.midpt = mean([T.s1 T.e2],2);
T = sortrows(T ,'midpt');

U = unstack( T , {'DupCounts' 'DupFreq'} , 'lib' ,  'GroupingVariables'  , {'s1' 's2' 'e1' 'e2' 'midpt'} ) ; 

ssp1 = fastaread( [ DATADIR 'ssp1_wt_sequenced_amplicon.fasta' ] ) ;



%%
clrs = parula(5);
idx1 = strcmp(T.lib,'lib_1');
fh = figure('units','centimeters','position',[5 5 20 7]) ;
hold on ; 
plot( T.midpt(idx1) , T.DupCounts(idx1) , '-','Color',clrs(1,:))
plot( 2+T.midpt(~idx1) , T.DupCounts(~idx1) , '-','Color',clrs(4,:))
%plot( T.midpt(~idx1) , -1*T.DupCounts(~idx1) , '-')

axis tight ;
legend({'rep 1' 'rep 2'})
%ylim([-20 20])
ylim([0 20])
ylabel('# of reads supporting each MTD')
set(gca,'ytick',-20:10:20)
set(gca,'yticklabel',[20 10 0 10 20])
xlabel('Position along ssp1 (nt)')

A1_GTCGTCCG_loc = 582.5; %midpt 512 start in ORF , 533 in seq , diff=21 ;  533    540    625    632
A2_AGGCA_loc = 778.5  ; %s1 730 in ORF , 749
A3_GCTTT_loc = 1069; % 1013 in ORF 

text(A1_GTCGTCCG_loc,10,'A1')
text(A2_AGGCA_loc,10,'A2')
text(A3_GCTTT_loc,10,'A3')
print('-dpng','~/Downloads/Figure2__ssp1_MTDs_across_gene_barplot.png','-r600')
close ;
%% correlation between the two
X = U.DupCounts_lib_1 ; Y = U.DupCounts_lib_2 ; 
%X = U.DupFreq_lib_1 ;  Y = U.DupFreq_lib_2 ; 

X(X>50)=50;
Y(Y>50)=50;
%X(X==0)=random('uniform',-0.05,0.05,sum(X==0),1);
%Y(Y==0)=random('uniform',-0.05,0.05,sum(Y==0),1);

X(X==0)=0.5;
Y(Y==0)=0.5;
fh = figure('units','centimeters','position',[5 5 7 7]) ;
hold on ; 
sh = scatter( X , Y , 'ok' , 'MarkerFaceColor',[.7 .7 .7]);
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'xtick',[0.5 1 2 5 10 20 50])
set(gca,'ytick',get(gca,'xtick'))
set(gca,'xticklabel',[0 1 2 5 10 20 50])
set(gca,'yticklabel',get(gca,'xticklabel'))
axis tight ;
xlabel('# of reads in rep 1')
ylabel('# of reads in rep 2')
[rS,pS]=corr(U.DupCounts_lib_1,U.DupCounts_lib_2,'type','Spearman')
[rP,pP]=corr(U.DupCounts_lib_1,U.DupCounts_lib_2,'type','Pearson')
text(0.6 ,40 , sprintf( 'r = %0.02f',rP));
text(0.6 ,25 , sprintf( '\\rho = %0.02f',rS));
print('-dpng','~/Downloads/Figure2__ssp1_MTDs_across_gene_correlation_scatterplot.png','-r600')
close ;