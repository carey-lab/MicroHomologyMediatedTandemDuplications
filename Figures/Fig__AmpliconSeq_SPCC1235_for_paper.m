DATADIR = '/Users/lcarey/SynologyDrive/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/DataFromCluster/';
T = readtable( [ DATADIR 'PombeAmpliconSeq_E3_alltsvs.txt' ],'TreatAsEmpty','-');T.Properties.VariableNames = {'chr' 's1' 'e1' 's2' 'e2' 'ReadDepth' 'Ndup' 'Ncol' 'DupFrq' 'ColFrq'};

T.DupFrq(isnan(T.DupFrq))=0;
T.Ndup(isnan(T.Ndup))=0;
T.Ncol = []; T.ColFrq = [] ;

T.MHlen = T.e1 - T.s1 + 1 ;
T.InterMHDistance = T.s2 - T.e1 + 1 ;
T.InterMHDistanceR = round(T.InterMHDistance./10)*10 ; 

%T=T(T.MHlen>=6,:);
for I = 1:height(T),  if strcmp(T.chr{I},'SPC1235.01.rad27d.1'), T.chr{I}='SPCC1235.01.rad27d.1'; end; end 
T = sortrows(T,{'chr'  's1' 's2' 'e1' 'e2'});
% figure showing biological replicates of WT & each mutant
wtbaseline_datamat = [T.DupFrq( strcmp(T.chr,'SPCC1235.01.1')) T.DupFrq( strcmp(T.chr,'SPCC1235.01.2')) ] ; 
WTbaseline = mean( wtbaseline_datamat , 2 );  
%%
R = table();
R.allstrains = unique(T.chr(regexpcmp(T.chr,'SPCC123')));
for I = 1:height(R)
    dupfrq = T.DupFrq(strcmp(T.chr,R.allstrains{I}));
    log2r = log2( (dupfrq+0.1) ./ (WTbaseline+0.1)  );
    log2r_notzero = log2r(  WTbaseline>0 & dupfrq>0 ) ;
    R.mean_dupfreq(I) = mean(dupfrq);
    R.mean_log10_dupfreq(I) = mean(log10(dupfrq(dupfrq>0)));
    R.mean_dupfreq_nz(I) = mean(dupfrq(dupfrq>0));
    R.mean_not_zero(I) = mean(log2r_notzero);
    R.mean_all(I) = mean(log2r);
    [ ~,R.ttest_mean_not_zero_p(I) ]=ttest(log2r_notzero);
    [ ~,R.ttest_mean_all_p(I) ]=ttest(log2r);
    [ ~,R.ttest_paired_not_zero_p(I) ]=ttest(log2r_notzero , WTbaseline(WTbaseline>0 & dupfrq>0));
    [ ~,R.ttest_paired_all_p(I) ]=ttest(log2r  , WTbaseline );
    R.Pct_MHPs_with_Dup(I)  = 100*mean(T.DupFrq(strcmp(T.chr,R.allstrains{I}))>0) ; 
    R.log2r__Pct_MHPs_with_Dup_over_WT(I)  = log2( R.Pct_MHPs_with_Dup(I) ./ mean(mean(wtbaseline_datamat>0)*100))  ; 

end
R.Pct_MHPs_with_Dup = round(R.Pct_MHPs_with_Dup*10)/10 ; 
R    
R.genotype = regexprep( R.allstrains , 'SPCC1235.01.','');
R.genotype = regexprep( R.genotype , '\.[12]','');
R.genotype = regexprep( R.genotype , '^[12]$','WT');

G = grpstats( R , 'genotype' ,'mean' ,'datavars' , R.Properties.VariableNames(2:end-1));
G.genotype = regexprep( G.genotype , 'd$','\\Delta');

%writetable( R , '~/Downloads/AmpliconLib3__DuplicationFreqInMutants_RelativeToWT.xlsx') ;
%%

% Figure for publication
fh = figure('units','centimeters','position',[5 5  6 9]) ;
hold on ;
bar( G.mean_Pct_MHPs_with_Dup(1:3) ,'FaceColor',[.8 .8 .8]);
plot( [1 1 2 2 3 3 ] , R.Pct_MHPs_with_Dup(1:6) , 'ok' ,'MarkerSize',5)
xlim([0.5 3.5])
set(gca,'xtick',1:3)
set(gca,'xticklabel',G.genotype(1:3));
ylabel('% of MHPairs with an MTD');

fh = figure('units','centimeters','position',[5 5  6 9]) ;
hold on ;
bar( G.mean_mean_log10_dupfreq(1:3) ,'FaceColor',[.8 .8 .8]);
plot( [1 1 2 2 3 3 ] , R.mean_log10_dupfreq(1:6) , 'ok' ,'MarkerSize',5)
xlim([0.5 3.5])
set(gca,'xtick',1:3)
set(gca,'xticklabel',G.genotype(1:3));
ylabel('Mean MTD frequency (log_{10})');
ylim([0.7 1.3])


datamat = [T.DupFrq( strcmp(T.chr,'SPCC1235.01.1')) T.DupFrq( strcmp(T.chr,'SPCC1235.01.2')) ] ; 
WT = mean( datamat , 2 );  

datamat = [T.DupFrq( strcmp(T.chr,'SPCC1235.01.rad27d.1')) T.DupFrq( strcmp(T.chr,'SPCC1235.01.rad27d.2')) ] ; 
rad27 = mean( datamat , 2 );  

fh = figure('units','centimeters','position',[5 5  12 6]) ;
plot( WT+0.1 , rad27+0.1 , 'ok' ,'MarkerFaceColor',[0.8 0.8 0.8]);
xlabel('wild-type MTD frequency (per 10^6 reads)')
ylabel('rad27\Delta MTD frequency')
set(gca,'yscale','log')
set(gca,'xscale','log')
line(xlim,xlim)
set(gca,'xtick',[0.1 1 10 100 1e3 1e4 ]); set(gca,'xticklabel',[0 1 10 100 1e3 1e4])
set(gca,'ytick',[0.1 1 10 100 1e3 1e4 ]);set(gca,'yticklabel',[0  1 10 100 1e3 1e4])
xlim([0.1 1e4]);ylim(xlim)
legend({'one single MHPair in the gene SPCC1235.01'},'location','nw')
print('-dpng','~/Downloads/x.png','-r500')


fh = figure('units','centimeters','position',[5 5 4.5 4.5]) ;
plot( 0.1+T.DupFrq( strcmp(T.chr,'SPCC1235.01.1')) , 0.1+T.DupFrq( strcmp(T.chr,'SPCC1235.01.2'))  , 'ok' ,'MarkerFaceColor',[0.8 0.8 0.8]);
xlabel('WT MTD freq. rep 1')
ylabel('WT MTD freq. rep 2')
set(gca,'yscale','log')
set(gca,'xscale','log')
line(xlim,xlim,'LineWidth',2)
set(gca,'xtick',[0.1 1 10 100 1e3  ]); set(gca,'xticklabel',[0 1 10 100 1e3 ])
set(gca,'ytick',[0.1 1 10 100 1e3  ]);set(gca,'yticklabel',[0  1 10 100 1e3 ])
xlim([0.1 1e4]);ylim(xlim)
print('-dpng','~/Downloads/xWT.png','-r500')

%%
figure; bar(R.Pct_MHPs_with_Dup); set(gca,'xtick',1:height(R)); 
set(gca,'xticklabel',regexprep(R.allstrains,'.01',''))
xticklabel_rotate([],45)
ylabel('% of MHPs with duplication')
line( xlim , [R.Pct_MHPs_with_Dup(2) R.Pct_MHPs_with_Dup(2) ])
line( xlim , [R.Pct_MHPs_with_Dup(1) R.Pct_MHPs_with_Dup(1) ])

figure; barh(R.mean_log10_dupfreq); set(gca,'ytick',1:height(R)); 
set(gca,'yticklabel',regexprep(R.allstrains,'.01',''))
xlabel('MTD frequency (log_{10})')
line([R.mean_log10_dupfreq(2) R.mean_log10_dupfreq(2) ] , ylim )
line([R.mean_log10_dupfreq(1) R.mean_log10_dupfreq(1) ] , ylim)
xlim([0.3 max(xlim)])