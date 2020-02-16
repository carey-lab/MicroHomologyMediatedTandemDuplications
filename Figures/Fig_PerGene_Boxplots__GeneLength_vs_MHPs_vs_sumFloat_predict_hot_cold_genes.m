% show that the model can predict the # of MTDs in each gene, and that
% different genes have very different predicted #s of MTDs


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
clear 'D' ; 

%T = T( ~isnan(T.N_MHPs) , :);
%T = T( ~isnan(T.sumObs) , :);
%T = T( ~isnan(T.sumFloat) , :);


%% anything else with higher correlation ? 
idx = T.sumObs<=20 ; 
for I = 1:30
    try
    r(I) = corr( table2array(T(idx,I)) , T.sumObs(idx) ,'rows','complete');
    catch
    end
end
R = table();
R.r = r';
R.vn = T.Properties.VariableNames' ; 
R = sortrows(R,'r','descend')
%%
clr_len = [ 0.85 0.85 0.85 ] ; 
clr_mdl = [77	175	74	]./255 ; 

T.sumObs_THRESH = T.sumObs ; 
    T.sumObs_THRESH(T.sumObs_THRESH>10) = 12 ; 

fh = figure('units','centimeters','position',[5 5 30 7]) ;
subplot(1,4,1) 
boxplot( T.GeneLength , T.sumObs_THRESH ,'symbol','' , 'Positions' , unique(T.sumObs_THRESH) ,'PlotStyle','compact' ,'Color' , clr_len )
ylabel('Gene Length (kb)')
%xlabel('# of observed MTDs per gene')
ylim([0.5 prctile( T.GeneLength,99.9)])
[r1,~] = bootstrp(1000,@corr,T.sumObs,T.GeneLength);
%title( sprintf( 'r=%0.03f' , rP1 ) )
set(gca,'xtick',[0 2 4 6 8 10 12]) ; set(gca,'xticklabel',{'0' '2' '4' '6' '8' '10' '>10'}) 


subplot(1,4,2) 
boxplot( T.N_MHPs , T.sumObs_THRESH ,'symbol','', 'Positions' , unique(T.sumObs_THRESH) ,'PlotStyle','compact' ,'Color','r')
ylabel('# of MHPs')
xlabel( '# of observed MTDs per gene')
ylim([1000 prctile( T.N_MHPs,99.9)])
[r2,~] = bootstrp(1000,@corr,T.sumObs,T.N_MHPs);
set(gca,'xtick',[0 2 4 6 8 10 12]) ; set(gca,'xticklabel',{'0' '2' '4' '6' '8' '10' '>10'}) 


subplot(1,4,3) 
boxplot( T.sumFloat , T.sumObs_THRESH ,'symbol','' , 'Positions' , unique(T.sumObs_THRESH) ,'PlotStyle','compact' ,'Color',clr_mdl)
ylabel('Predicted MTD score')
%xlabel('# of observed MTDs per gene')
ylim([5 prctile( T.sumFloat,99.9)])
[r3,~] = bootstrp(1000,@corr,T.sumObs,T.sumFloat);
set(gca,'xtick',[0 2 4 6 8 10 12]) ; set(gca,'xticklabel',{'0' '2' '4' '6' '8' '10' '>10'}) 

subplot(1,4,4)  ; hold on ;
r1 = corr(T.sumObs,T.GeneLength,'rows','complete');
r2 = corr(T.sumObs,T.N_MHPs,'rows','complete');
r3 = corr(T.sumObs,T.sumFloat,'rows','complete');
bar( 1,  mean(r1) ,'FaceColor',clr_len) ; %errorbar(1,mean(r1),std(r1),'.k');
bar( 2,  mean(r2) ,'FaceColor','r') ; %errorbar(2,mean(r2),std(r2),'.k');
bar( 3,  mean(r3) ,'FaceColor',clr_mdl) ; %errorbar(3,mean(r3),std(r3),'.k');
xlim([0.4 3.6])
set(gca,'xticklabels',{'Length' '# MHPs' 'model'})
ylabel('Correlation')

print('-dpng',[ FIGDIR '_PerGene_Boxplots__GeneLength_vs_MHPs_vs_sumFloat_predict_hot_cold_genes'] ,'-r300');
close ; 