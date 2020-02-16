DATADIR = '~/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/DataFromCluster/';
T0= readtable( [DATADIR '10k_rm.sign.count.tsv'],'FileType','text','Format','%s%d%d%d%d%d%d');
T1= readtable( [DATADIR 'trimmed_crop100.sign.count.tsv'],'FileType','text','Format','%s%d%d%d%d%d%d');
T2= readtable( [DATADIR 'trimmed_headcrop50.sign.count.tsv'],'FileType','text','Format','%s%d%d%d%d%d%d');


T0.Properties.VariableNames = {'chr' 's1' 'e1' 's2' 'e2' 'DupCounts' 'ColCounts'};
T1.Properties.VariableNames = {'chr' 's1' 'e1' 's2' 'e2' 'DupCounts' 'ColCounts'};
T2.Properties.VariableNames = {'chr' 's1' 'e1' 's2' 'e2' 'DupCounts' 'ColCounts'};

%%

T0.InterMHDist = T0.s2 - T0.e1 + 1;  
T1.InterMHDist = T1.s2 - T1.e1 + 1; 
T2.InterMHDist = T2.s2 - T2.e1 + 1; 
T0.MHLen = T0.e1 - T0.s1 + 1; 
T1.MHLen = T1.e1 - T1.s1 + 1; 
T2.MHLen = T2.e1 - T2.s1 + 1; 

T0.InterMHDistR = round(T0.InterMHDist./5)*5;
T1.InterMHDistR = round(T1.InterMHDist./5)*5;
T2.InterMHDistR = round(T2.InterMHDist./5)*5;

T0.MHLenR = T0.MHLen ; T0.MHLenR(T0.MHLenR>6) = 6 ; 
T1.MHLenR = T1.MHLen ; T1.MHLenR(T1.MHLenR>6) = 6 ; 
T2.MHLenR = T2.MHLen ; T2.MHLenR(T2.MHLenR>6) = 6 ; 

T0.HasDup = T0.DupCounts>0;
T1.HasDup = T1.DupCounts>0;
T2.HasDup = T2.DupCounts>0;

%%

G0 = grpstats(T0,{'InterMHDistR' 'MHLenR'},{'mean'},'DataVars','HasDup');
G1 = grpstats(T1,{'InterMHDistR' 'MHLenR'},{'mean'},'DataVars','HasDup');
G2 = grpstats(T2,{'InterMHDistR' 'MHLenR'},{'mean'},'DataVars','HasDup');
%%

fh = figure('units','centimeters','position',[5 5  8 25]) ; 
t = tiledlayout(3,1);
nexttile;
gscatter(G0.InterMHDistR,(100*G0.mean_HasDup)+1e-3,G0.MHLenR)
set(gca,'yscale','log')
title('2x150nt full length reads')
xlabel('nt between MHPair')
legend({'MHL=4' 'MHL=5' 'MHL>=6'} , 'location','NorthOutside')
set(gca,'ytick',[1e-4 1e-3 1e-2 1e-1 1])

nexttile;
gscatter(G1.InterMHDistR,(100*G1.mean_HasDup)+1e-3,G1.MHLenR)
set(gca,'yscale','log')
title('50nt removed from end of read')
xlabel('nt between MHPair')
legend('off')
set(gca,'ytick',[1e-4 1e-3 1e-2 1e-1 1])

nexttile;
gscatter(G2.InterMHDistR,(100*G2.mean_HasDup)+1e-3,G2.MHLenR)
set(gca,'yscale','log')
title('50nt removed from start of read')
xlabel('nt between MHPair')
legend('off')
set(gca,'ytick',[1e-4 1e-3 1e-2 1e-1 1])


ylabel(t,'% of MHPs with an MTD')
print('-dpng','~/Downloads/SupFig__ReadLenght_DoesNotDetermine_150nt_spacing','-r600')