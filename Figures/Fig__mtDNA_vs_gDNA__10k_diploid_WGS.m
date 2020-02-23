DATADIR = '~/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/DataFromCluster/';
T0= readtable( [DATADIR '10k_rm.sign.count.tsv'],'FileType','text','Format','%s%d%d%d%d%d%d');
%T0 = readtable( '~/Downloads/10k_rm_10pctA.sign.count.tsv' ,'FileType','text','Format','%s%d%d%d%d%d%d');
%T0= readtable( [DATADIR 'SRR7817502_rm.sign.count.tsv'],'FileType','text','Format','%s%d%d%d%d%d%d');

T0.Properties.VariableNames = {'chr' 's1' 'e1' 's2' 'e2' 'DupCounts' 'ColCounts'};
T0 = T0( ~strcmp(T0.chr,'mating_type_region'),:);
T0 = T0( ~strcmp(T0.chr,'chr_II_telomeric_gap'),:);
T0 = T0( ~strcmp(T0.chr,'MTR'),:);
T0 = T0( ~strcmp(T0.chr,'AB325691'),:);


T0.MHLen = T0.e1 - T0.s1 + 1;  
T0.InterMHDist = T0.s2 - T0.e1 + 1;  
T0.HasDup = T0.DupCounts>0;

T0.InterMHDistR = round(T0.InterMHDist./2)*2;
T0.InterMHDistR(T0.InterMHDistR>300)=300 ; 

T0.MHLenR = T0.MHLen ; T0.MHLenR(T0.MHLenR>7) = 7 ; 

%% %% %%  properly downsample
% Divide DupCounts by 10
%  calculate the sum of all DupCounts < 1
%  randomly assign these reads to all MHPs < 1
%   set the rest of the MHPs to 0
DNSMPL = 100 ; 
T = T0(strcmp(T0.chr,'MT') , :);
idx_has_dup = find(T.HasDup);
T_HasDup = T.HasDup ;
T_InterMHDistR = T.InterMHDistR ; 
T_MHLenR = T.MHLenR ; 
T_RS_HasDup = NaN(0);
T_RS_InterMHDistR = NaN(0);
T_RS_MHLenR = NaN(0);

% this works, work with vectors; 
%   G = grpstats( T.HasDup ,  [T.InterMHDistR T.MHLenR ] )
for N = 1:10000
    idx_to_keep = randsample( idx_has_dup , round(numel(idx_has_dup)/DNSMPL) ) ; 
    has_dup = zeros(numel(T_HasDup),1);
    has_dup(idx_to_keep) = 1; 
    T_RS_HasDup = vertcat(T_RS_HasDup , has_dup);
    T_RS_InterMHDistR = vertcat(T_RS_InterMHDistR , T_InterMHDistR);
    T_RS_MHLenR = vertcat(T_RS_MHLenR , T_MHLenR);
end
save('~/Downloads/downsample_MTdna.mat')
%%
R = table(); R.HasDup=T_RS_HasDup; R.InterMHDistR=T_RS_InterMHDistR;R.MHLenR=T_RS_MHLenR;
R = grpstats( R , {'InterMHDistR' 'MHLenR'} ,'mean' , 'DataVars' , 'HasDup');
save('~/Downloads/downsample_MTdnaR.mat' , 'R')
%%
figure; 
gscatter( R.InterMHDistR , R.mean_HasDup , R.MHLenR)
set(gca,'yscale','log')
ylabel('% of MHPs with an MTD')
xlabel('inter-MHP distance (nt)')
%% MT DNA has 800,000x coverage, vs 9,000x coverage for whole-gnome
%%
T0.InterMHDistR = round(T0.InterMHDist./10)*10;

T0.MHLenR = T0.MHLen ; T0.MHLenR(T0.MHLenR>7) = 7 ; 
G0 = grpstats(T0,{'InterMHDistR' 'MHLenR' 'chr'},{'mean'},'DataVars',{'HasDup' 'DupCounts'});

%G0.mean_HasDup = G0.mean_HasDup+1e-5 ;
%



fh = figure('units','centimeters','position',[5 5  18 18]) ; 
t = tiledlayout(2,2);

nexttile;
idx = G0.MHLenR==4; 
gscatter( G0.InterMHDistR(idx) , 100*G0.mean_HasDup(idx) , G0.chr(idx) )
set(gca,'yscale','log')
title('MH length = 4')
xlabel('nt between MHPair','FontSize',15)
set(gca,'ytick',[1e-4 1e-3 1e-2 1e-1 1])
legend('location','NorthOutside')
xlim([0 350])

nexttile;
idx = G0.MHLenR==5; 
gscatter( G0.InterMHDistR(idx) , 100*G0.mean_HasDup(idx) , G0.chr(idx) )
set(gca,'yscale','log')
title('MH length = 5')
xlabel('nt between MHPair','FontSize',15)
set(gca,'ytick',[1e-4 1e-3 1e-2 1e-1 1 10])
legend('location','NorthOutside')
xlim([0 350])

nexttile;
idx = G0.MHLenR==6; 
gscatter( G0.InterMHDistR(idx) , 100*G0.mean_HasDup(idx) , G0.chr(idx) )
set(gca,'yscale','log')
title('MH length = 6')
xlabel('nt between MHPair','FontSize',15)
set(gca,'ytick',[1e-4 1e-3 1e-2 1e-1 1 10])
legend('off')
xlim([0 350])


nexttile;
idx = G0.MHLenR==7; 
gscatter( G0.InterMHDistR(idx) , 100*G0.mean_HasDup(idx) , G0.chr(idx) );
legend('off');
set(gca,'yscale','log')
title('MH length >= 7')
xlabel('nt between MHPair','FontSize',15)
set(gca,'ytick',[1e-4 1e-3 1e-2 1e-1 1 10 ])
xlim([0 350])

ylabel( t  , '% of MHPs with an MTD' ,'FontSize',15)
print('-dpng','~/Downloads/mtDNA_10pct.png','-r200');
%close;