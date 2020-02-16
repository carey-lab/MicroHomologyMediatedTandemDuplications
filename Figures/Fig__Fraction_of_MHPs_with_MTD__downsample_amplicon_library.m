% Fig__Fraction_of_MHPs_with_MTD__downsample_amplicon_library
% if we reduce the coverage, how many MHPairs would have an MTD? 
FIGDIR  = '~/Nutstore Files/Microhomology shared folder/Figures/Fig4/' ;
FN = '~/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/DataFromCluster/PombeAmpliconSeq_E2_alltsvs.txt';
T = readtable( FN ,'FileType','text','Format','%s%s%d%d%d%d%d%d%d%d%d','TreatAsEmpty','-');
T = T(ismember(T.Var2, {'ssp1' 'ssp2' 'SPCC1235.01' }) & ismember(T.Var1, {'lib_1' 'lib_2' })  ,:) ;

DIR = '~/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/DataFromCluster/';
G = readtable( [ DIR '10k_rm.sign.count.tsv'] ,'FileType','text');
GENOME_NDups = G.Var6 ; 
clear 'G' 

FigSize = [5 5 6 6] ; 
%% for testing
R = table();
R.xl = (logspace(0,5,1000) )' ; 


Q = T(ismember(T.Var2, {'ssp1' }) & ismember(T.Var1, {'lib_1'  })  ,:) ;
R.ndupreads_ssp1_lib1 = arrayfun( @(X)mean(round(Q.Var8 ./ X) > 0)*100, R.xl ) ; 
R.nreads_ssp1_lib1    = arrayfun( @(X)mean(round(Q.Var7 ./ X)), R.xl ) ; 

Q = T(ismember(T.Var2, {'ssp1' }) & ismember(T.Var1, {'lib_2'  })  ,:) ;
R.ndupreads_ssp1_lib2 = arrayfun( @(X)mean(round(Q.Var8 ./ X) > 0)*100, R.xl ) ; 
R.nreads_ssp1_lib2    = arrayfun( @(X)mean(round(Q.Var7 ./ X)), R.xl ) ; 

Q = T(ismember(T.Var2, { 'SPCC1235.01' }) & ismember(T.Var1, {'lib_1'  })  ,:) ;
R.ndupreads_SPCC1235_lib1 = arrayfun( @(X)mean(round(Q.Var8 ./ X) > 0)*100, R.xl ) ; 
R.nreads_SPCC1235_lib1    = arrayfun( @(X)mean(round(Q.Var7 ./ X)), R.xl ) ; 

Q = T(ismember(T.Var2, { 'SPCC1235.01' }) & ismember(T.Var1, {'lib_2'  })  ,:) ;
R.ndupreads_SPCC1235_lib2 = arrayfun( @(X)mean(round(Q.Var8 ./ X) > 0)*100, R.xl ) ; 
R.nreads_SPCC1235_lib2    = arrayfun( @(X)mean(round(Q.Var7 ./ X)), R.xl ) ; 

R.ndupreads_10k = arrayfun( @(X)mean(round(GENOME_NDups ./ X) > 0)*100, R.xl ) ; 
R.nreads_10k    = arrayfun( @(X)mean(round(1e4 ./ X)), R.xl ) ; 

%%
fh = figure('units','centimeters','position',FigSize);
hold on ; 
plot( R.nreads_SPCC1235_lib1 , R.ndupreads_SPCC1235_lib1 ,'-','LineWidth',3,'DisplayName' , 'SPCC1235 rep 1')
plot( R.nreads_SPCC1235_lib2 , R.ndupreads_SPCC1235_lib2 ,'-','LineWidth',3,'DisplayName' , 'SPCC1235 rep 2')
plot( R.nreads_ssp1_lib1 , R.ndupreads_ssp1_lib1 ,'-','LineWidth',3,'DisplayName' , 'ssp1 rep 1')
plot( R.nreads_ssp1_lib2 , R.ndupreads_ssp1_lib2 ,'-','LineWidth',3,'DisplayName' , 'ssp1 rep 2')
plot( R.nreads_10k , R.ndupreads_10k ,'-k','LineWidth',3,'DisplayName' , '10k genome')
set(gca,'xscale','log')
xlabel('Coverage (avg. # of reads)')
ylabel('% of MHPs with an MTD')
legend('location','nw','box','off')
axis tight; 
set(gca,'xtick',logspace(1,7,7))
xlim([99 max(xlim)])
print('-dpng',[ FIGDIR 'Fraction_of_MHPs_with_MTD__downsample_amplicon_library_lin'] ,'-r300');
set(gca,'yscale','log')
legend('off')
ylim([1e-3 max(ylim)])
set(gca,'ytick',logspace(-5,2,8))
print('-dpng',[ FIGDIR 'Fraction_of_MHPs_with_MTD__downsample_amplicon_library_log'] ,'-r300');
close  ; 

%%
X  = R.nreads_ssp1_lib1  ; 
Y  = R.ndupreads_ssp1_lib1 ;
Y2 = Y(X>1.5*10^5);
X2 = X(X>1.5*10^5);

%X  = R.nreads_SPCC1235_lib1 ; 
%Y  = R.ndupreads_SPCC1235_lib1 ;
[xData, yData] = prepareCurveData( X2, Y2 );
ft = fittype( 'power2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
[fitresult, gof] = fit( xData, yData, ft, opts );

xl = logspace(1,9,1000) ;
Y =  feval(fitresult,xl) ; 

fh = figure('units','centimeters','position',FigSize-[0 0 1 0]);
yyaxis left ; 

plot( xl , Y ,'-k','LineWidth',3,'DisplayName' , 'prediction : 10k genome')

ax = gca ; 
ax.YTick = 0:20:100 ;  
set(gca,'YColor','k') ; 
ylim([0 100])  ;  xlim([50 5e8])
ylabel('% of MHPs with an MTD')


yyaxis right
ax = gca ; 
plot( xl , Y ,'-k','LineWidth',3,'DisplayName' , 'prediction : 10k genome')
ylim([0 100])  ;  xlim([50 5e8])
set(gca,'YColor','k') ; 
ax.YTick = 0:20:100 ; yticklabels( {'0' '5e5' '10^6' '1.5e6' '2e6'  '2.5e6'}  ) ; 
%ax.YTick = 0:40:100 ; yticklabels( {'0'  '10^6'  '2e6'  }  ) ; 

ylabel('# of MTDs detected')
ax.XTick = (0:4) * 1e8  ;ax.XTickLabel = (0:4)   ;
xlabel('Coverage (×10^8)')
print('-dpng',[ FIGDIR 'Fraction_of_MHPs_with_MTD__downsample_amplicon_library_prediction'] ,'-r300');
close  ; 