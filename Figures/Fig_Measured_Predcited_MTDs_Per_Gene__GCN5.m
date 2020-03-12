% We identified three MTDs in the SAGA complex histone acetyltransferase catalytic subunit gcn5, far higher than expected by 
% chance based on the gene length (p=@) or the number of MHPs (p=@). The model predicts gcn5 to be in the top @% of genes for MTD frequency, 
%  suggesting that MTDs in gcn5 should be frequent. Indeed, whole-genome sequencing of 16 suppressors of htb1G52D identified MTDs in gcn5, as well as in 
%  ubp8, where we also observed a single MTD in our high-coverage sequencing data. 
%  These results suggest that MTDs likely exist at high frequency in most genes, 
%  and are frequently the raw material on which natural selection acts. 
%
% LBC Febuary 2020

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
T = sortrows(T,'sumFloat','descend');
T.sumFloat_Rank = (1:height(T))' ;

T = sortrows(T,'SumFloat_N','ascend');
T.sumFloat_N_Rank = (1:height(T))' ;


T = innerjoin(D , T  , 'key' , 'systematic_name'); 
 
%% save for supplementary table
T.N_observed_MTDs = T.sumObs ; 
T.N_predicted_MTDs = T.sumPre ; 
T.sum_predicted_MHP_score = T.sumFloat ; 
T.sum_predicted_MHP_score_Rank = T.sumFloat_Rank ; 
T.sum_predicted_MHP_score_Rank( isnan(T.sum_predicted_MHP_score) ) = NaN ; 

writetable(  T(:,{'gene'  'systematic_name' 'N_observed_MTDs' 'N_predicted_MTDs' 'sum_predicted_MHP_score' 'GeneLength' 'sum_predicted_MHP_score_Rank' 'N_MHPs'} )  ,  '~/Downloads/MTDs_per_gene.xlsx' );
 
%% histogram distribution of MTDs per gene
fh = figure('units','centimeters','position',[5 5 5 6]);
hold on ; 

Y=T.sumObs;
Y(Y>=10)=10.5;

Yn=T.sumObs_N;
Yn(Yn>=10)=10.5;

histogram(Y,-0.5:11);
%histogram(Yn,-0.5:11);

set(gca,'yscale','log')
ylabel('# of genes')
xlabel('Obs. MTDs per gene')
set(gca,'xtick',[0 5 9.5])
set(gca,'xticklabel',{'0' '5' '>=10'})
set(gca,'ytick',[1 10 1e2 1e3 1e4 1e5])
axis tight; 
yl = ylim ;

print('-dpng', [ FIGDIR 'Distribution_of_MTDs_per_gene' ] ,'-r300');
close

% normalized by gene length
fh = figure('units','centimeters','position',[5 5 5 6]);
hold on ; 

Yn=T.sumObs_N;
Yn(Yn>=10)=10.5;
histogram(Yn,-0.5:11);

set(gca,'yscale','log')
ylabel('# of genes')
xlabel('Obs. MTDs/gene/kb')
set(gca,'xtick',[0 5 9.5])
set(gca,'xticklabel',{'0' '5' '>=10'})
set(gca,'ytick',[1 10 1e2 1e3 1e4 1e5])
axis tight; 
ylim(yl)

print('-dpng', [ FIGDIR 'Distribution_of_MTDs_per_gene_per_kb' ] ,'-r300');
close

%% histogram  of MHPs per gene, and MTDs normalized by MHP
fh = figure('units','centimeters','position',[5 5 6 6]);
hold on ; 

G = T.sumObs ; 
G(G>=6)=6;
boxplot( (T.N_MHPs) , G ,'symbol','','notch','on' , 'Positions' , [ 0:5 6.2] ) ;
ylim([0 10500])
xlabel('# of MTDs per gene')
ylabel('# of MHPs per gene')
set(gca,'xtick',0:6.2)
set(gca,'xticklabel', {'0' '1' '2' '3' '4' '5' '>=6'})

%set(gca,'ytick',[0 1e3 2e3 5e3 1e4])
print('-dpng', [ FIGDIR 'Distribution_of_MTDs_per_gene_boxplot' ] ,'-r300');
close

%%
Yobs = T.sumObs; 
Ypred = T.sumPre ; 
idx = find(strcmp(T.gene,'gcn5'));

pctile_pred = mean( T.sumPre>T.sumPre(idx))*100
pctile_obs = mean( T.sumObs>T.sumObs(idx))*100

pctile_pred_N = mean( T.sumPre_N>T.sumPre_N(idx))*100
pctile_obs_N = mean( T.sumObs_N>T.sumObs_N(idx))*100
%% ecdf showing $ of MTDs, with gcn5


gcn5_lg_txt = sprintf('\\it{gcn5} , top %0.0f%% obs, , %0.0f%% pred' , pctile_obs , pctile_pred) ; 

fh = figure('units','centimeters','position',[5 5 8 8]);
hold on ;
[f,x]=ecdf( T.sumObs );
x(x==0)=0.5;
plot(x,f*100,'-k','LineWidth',2,'DisplayName','Observed MTDs')

[f,x] = ecdf( T.sumPre ) ; 
x(x==0)=0.5;
plot(x,f*100,'-r','LineWidth',2,'DisplayName','Predicted MTDs')

line( [ T.sumObs(idx) T.sumObs(idx)] , ylim ,'Color','b','DisplayName', gcn5_lg_txt ,'LineWidth',3)
xlabel('# of MTDs per gene')

ax = gca ; 
grid on ;
set(gca,'ytick',0:10:100)
xlim([0 10])
legend('location','se');
ylabel('Cummulative % of genes')
set(gca,'xscale','log')
set(gca,'xtick',[0.5 1:10] )
set(gca,'xticklabel',[0 1:10] )
ax.XMinorTick =  'off'  ; 
ax.XMinorGrid = 'off' ;

print('-dpng', [ FIGDIR 'Distribution_of_MTDs_per_gene_ecdf_gcn5' ] ,'-r300');
close ; 