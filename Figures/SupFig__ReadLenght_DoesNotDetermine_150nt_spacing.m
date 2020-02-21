% SupFig__ReadLenght_DoesNotDetermine_150nt_spacing
% determine if the 150nt inter-MH-Pair spacing from the diploid 10k sequencing is due to 
%   the 150nt readlength: trim 50nt from the front or back of each read, and recalculate
% 
%  also determine if shorter reads have more short inter-MHPs because insertions more likely to be softclipped
%
%
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

fh = figure('units','centimeters','position',[5 5  18 18]) ; 
t = tiledlayout(2,2);
nexttile;
gscatter(G0.InterMHDistR,(100*G0.mean_HasDup)+1e-3,G0.MHLenR)
set(gca,'yscale','log')
title('2x150nt full length reads')
xlabel('nt between MHPair','FontSize',15)
set(gca,'ytick',[1e-4 1e-3 1e-2 1e-1 1])
legend('off')

nexttile;
gscatter(G1.InterMHDistR,(100*G1.mean_HasDup)+1e-3,G1.MHLenR)
set(gca,'yscale','log')
title('50nt trimmed from end of read')
xlabel('nt between MHPair','FontSize',15)
legend('off')
set(gca,'ytick',[1e-4 1e-3 1e-2 1e-1 1])

nexttile;
gscatter(G2.InterMHDistR,(100*G2.mean_HasDup)+1e-3,G2.MHLenR)
set(gca,'yscale','log')
title('50nt trimmed from start of read')
xlabel('nt between MHPair','FontSize',15)
legend('off')
set(gca,'ytick',[1e-4 1e-3 1e-2 1e-1 1])
legend('off')

%xlabel(t,'nt between MHPair','FontSize',15)
ylabel(t,'% of MHPairs with an MTD','FontSize',15)
print('-dpng','~/Downloads/SupFig__ReadLenght_DoesNotDetermine_150nt_spacing','-r200')

legend({'MH Length = 4' 'MH Length = 5' 'MH Length >= 6'} , 'location','EastOutside')
print('-dpng','~/Downloads/SupFig__ReadLenght_DoesNotDetermine_150nt_spacing__lengend','-r400')
