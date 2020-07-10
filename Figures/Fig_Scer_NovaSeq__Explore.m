FN = '/Volumes/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/Sarah/new_MH_project_S.cer/ProcessedData/';
FN = [FN 'MHfeatures.Scer_union.sign.count.tsv.txt' ];
opts = detectImportOptions( FN )
opts = setvartype(opts,{'MHP_start_position' 'ntclosestMHR' 'MHP_end_position' 'dis_to_telomere' 'duplication'},'uint32');
opts = setvartype(opts,{'MHlen' 'interMH'},'uint16');
%opts = setvartype(opts,{'duplication' 'dup' 'bintclosestMHR'},'logical');

T = readtable( FN , opts);
T.dup = logical(T.dup); 
% G = grpstats(T ,'MHlen' ,'sum','DataVars',{'dup'})
%% more data
%  cut.pl -n chr,MHP_start_position,MHP_end_position,dup < MHfeatures.Scer_union.sign.count.tsv.txt > ~/Downloads/Scer.bed
%   convert2bed -i gff < David_2006_total_RNA_transcripts_V64.gff3 | bedtools map -o mean -a ~/Downloads/Scer.bed -b - | perl -pne 's/\.$/0/'  | cut -f 5 > ~/Downloads/Scer.david06.txt
%   convert2bed -i wig < Steinmetz_2006_PolII_occupancy_WT_V64.wig | bedtools map -o mean -a ~/Downloads/Scer.bed -b - | perl -pne 's/\.$/NaN/' | cut -f 5 > ~/Downloads/Steinmetz_2006_PolII.txt

T.DAVID = dlmread('~/Downloads/Scer.david06.txt');
T.DAVID(T.DAVID<0)=0;

T.POLII = dlmread('~/Downloads/Steinmetz_2006_PolII.txt');

%%
KD = '~/ExternalData/Kaplan09/';
T.Kpred = dlmread([KD 'Kaplan_2009_predicted_average_nucleosome_occupancy_V64.txt']);
%T.KpredNorm = dlmread([KD 'Kaplan_2009_predicted_average_nucleosome_occupancy_normalized_V64.txt']);
T.KpredScore = dlmread([KD 'Kaplan_2009_predicted_nucleosome_positioning_model_score_V64.txt']);
T.KoccLog2 = dlmread([KD 'Kaplan_2009_YPD_nucleosome_occupancy_map_dMean_log2_sMOL_V64.txt']);
T.KoccInVitro = dlmread([KD 'Kaplan_2009_InVitro_nucleosome_occupancy_map_dMean_log2_sMOL_V64.txt']);

%%
% grep 'CDS' ~/ExternalData/SGD/saccharomyces_cerevisiae.gff | bedtools intersect  -c   -a ~/Downloads/Scer.bed -b - | cut -f 5  > ~/Downloads/Scer_anyf.txt
% grep 'CDS' ~/ExternalData/SGD/saccharomyces_cerevisiae.gff | bedtools intersect -f 1 -c   -a ~/Downloads/Scer.bed -b - | cut -f 5  > ~/Downloads/Scer_f1.txt

T.InGeneANY = dlmread( '~/Downloads/Scer_anyf.txt' );
T.InGeneALL = dlmread( '~/Downloads/Scer_f1.txt' );
T.InGeneANY = T.InGeneANY>0 ; 
T.InGeneALL = T.InGeneALL>0 ; 

%% classifier
N_with_dup = sum(T.dup)
Q = T( vertcat( find(T.dup==1) , randsample(find(T.dup==0),N_with_dup) ) , :);


%%
T.subtel = T.dis_to_telomere<50 ; 

T.d2t = round(T.dis_to_telomere./5)*5;
T.d2t(T.d2t>150)=150;
G1 = grpstats( T , 'd2t' ,'mean' ,'DataVars',{'MHlen' 'dup'});

T.rept = round(T.replication_time./4)*4;
T.rept(T.rept>50)=50;
G2 = grpstats( T , 'rept' ,'mean' ,'DataVars',{'MHlen' 'dup'});
G3 = grpstats( T , {'rept' 'subtel'} ,'mean' ,'DataVars',{'MHlen' 'dup'});



T.txn = round(T.DAVID*5)./5 ;
T.txn(T.txn>3)=3;
T.txnC = NaN(height(T),1);
T.txnC(T.DAVID==0) = 0 ; 
T.txnC(T.DAVID>0 & T.DAVID<3) = 1 ; 
T.txnC(T.DAVID>=2.5) = 2 ; 


T.txn2C = ones(height(T),1);
T.txn2C(T.POLII < prctile(T.POLII,25)) = 0 ; 
T.txn2C(T.POLII > prctile(T.POLII,75)) = 2 ; 


G4 = grpstats( T , {'txnC' 'd2t'},'mean' ,'DataVars',{'MHlen' 'dup'}) ; 
G5 = grpstats( T , {'txn2C' 'd2t'},'mean' ,'DataVars',{'MHlen' 'dup'}) ; 

%G4 = grpstats( T , {'rept' 'subtel'} ,'mean' ,'DataVars',{'MHlen' 'dup'});

%%
lgtxt = 'MTD frequency' ; 

figure; tiledlayout(3,2);

nexttile; hold on ;
plot( G1.d2t,100*G1.mean_dup,'-ok','LineWidth',2)
set(gca,'yscale','log')
ylabel('% of MHPs w/MTD')
xlabel('kb to closest telomere')
legend(lgtxt , 'location','ne');

nexttile; hold on ;
plot( G2.rept,100*G2.mean_dup,'-ok','LineWidth',2);
set(gca, 'XDir','reverse')
set(gca,'yscale','log')
ylabel('% of MHPs w/MTD')
xlabel('replication timing (minutes after release from G1)')
legend(lgtxt , 'location','ne');

nexttile; hold on ;
plot( G3.rept(G3.subtel==1),100*G3.mean_dup(G3.subtel==1),'-ob','LineWidth',2);
plot( G3.rept(G3.subtel==0),100*G3.mean_dup(G3.subtel==0),'-or','LineWidth',2);
set(gca, 'XDir','reverse')
set(gca,'yscale','log')
ylabel('% of MHPs w/MTD')
xlabel('replication timing (minutes after release from G1)')
legend({'subtelomeric' ,'>50kb from tel' } , 'location','ne');


nexttile; hold on ;
plot( G4.d2t(G4.txnC==0),100*G4.mean_dup(G4.txnC==0),'-ok','LineWidth',2,'DisplayName','not expressed');
plot( G4.d2t(G4.txnC==1),100*G4.mean_dup(G4.txnC==1),'-ob','LineWidth',2,'DisplayName','mid-expressed');
plot( G4.d2t(G4.txnC==2),100*G4.mean_dup(G4.txnC==2),'-or','LineWidth',2,'DisplayName','high expressed');
set(gca,'yscale','log')
ylabel('% of MHPs w/MTD')
xlabel('kb to nearest telomere')
legend( 'location','ne');


nexttile; hold on ;
plot( G5.d2t(G5.txn2C==0),100*G5.mean_dup(G5.txn2C==0),'-ok','LineWidth',2,'DisplayName','low PolII occ');
plot( G5.d2t(G5.txn2C==1),100*G5.mean_dup(G5.txn2C==1),'-ob','LineWidth',2,'DisplayName','mid PolII occ');
plot( G5.d2t(G5.txn2C==2),100*G5.mean_dup(G5.txn2C==2),'-or','LineWidth',2,'DisplayName','high PolII occ');
set(gca,'yscale','log')
ylabel('% of MHPs w/MTD')
xlabel('kb to nearest telomere')
legend( 'location','ne');
%%
[~,tbl,stats] = anovan( Q.dup , [Q.MHlen Q.dis_to_telomere_log10],'Display','off');
multcompare(stats,'Dimension',2);