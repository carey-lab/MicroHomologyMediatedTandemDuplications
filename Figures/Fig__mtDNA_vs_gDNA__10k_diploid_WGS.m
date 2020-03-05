DATADIR = '~/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/DataFromCluster/';

FIGNAME = '~/Nutstore Files/Microhomology shared folder/Figures/Supplementary Figures/InterMHDist__mtDNA_vs_gDNA__10k' ; 

T0= readtable( [DATADIR '10k_rm.sign.count.tsv'],'FileType','text','Format','%s%d%d%d%d%f%d');
%T0 = readtable( '~/Downloads/10k_rm_10pctA.sign.count.tsv' ,'FileType','text','Format','%s%d%d%d%d%d%d');
HAP = readtable( [DATADIR 'SRR7817502_rm.sign.count.tsv'],'FileType','text','Format','%s%d%d%d%d%d%d');
HAP.Properties.VariableNames = {'chr' 's1' 'e1' 's2' 'e2' 'DupCounts' 'ColCounts'};
HAP = HAP( ~strcmp(HAP.chr,'MTR'),:);
HAP = HAP( ~strcmp(HAP.chr,'AB325691'),:);
HAP.InterMHDist = HAP.s2 - HAP.e1 + 1;  


T0.Properties.VariableNames = {'chr' 's1' 'e1' 's2' 'e2' 'DupCounts' 'ColCounts'};
T0 = T0( ~strcmp(T0.chr,'mating_type_region'),:);
T0 = T0( ~strcmp(T0.chr,'chr_II_telomeric_gap'),:);
T0 = T0( ~strcmp(T0.chr,'MTR'),:);
T0 = T0( ~strcmp(T0.chr,'AB325691'),:);


T0.MHLen = T0.e1 - T0.s1 + 1;  
T0.InterMHDist = T0.s2 - T0.e1 + 1;  
T0.HasDup = T0.DupCounts>0;

T0.InterMHDistR = round(T0.InterMHDist./5)*5;
T0.InterMHDistR(T0.InterMHDistR>300)=300 ; 

T0.MHLenR = T0.MHLen ; T0.MHLenR(T0.MHLenR>7) = 7 ; 

T = T0(strcmp(T0.chr,'MT') , :);

G0 = grpstats(T0(~strcmp(T0.chr,'MT'),:),{'InterMHDistR' 'MHLenR'},{'mean'},'DataVars',{'HasDup' 'DupCounts'});
G0a = grpstats(T0(~strcmp(T0.chr,'MT'),:),{'MHLenR'},{'mean'},'DataVars',{'HasDup' 'DupCounts'});
G0b = grpstats(T0(~strcmp(T0.chr,'MT'),:),{'InterMHDistR'},{'mean'},'DataVars',{'HasDup' 'DupCounts'});


%% %% %%  properly downsample
% Divide DupCounts by 10
%  calculate the sum of all DupCounts < 1
%  randomly assign these reads to all MHPs < 1
%   set the rest of the MHPs to 0
% MT DNA has 836678x coverage, vs 9139x coverage for whole-gnome 
tic ;
DNSMPL = 836678 / 9139 ; 
NSamples = 1000 ; 
idx_has_dup = find(T.HasDup);
T_RS_HasDup = cell(1,NSamples);
N_MHPs = height(T) ; 
N_MTDs_per_sample_max = round(sum(T.DupCounts./DNSMPL  )) ; % the max value
N_MTDs_per_sample_min = round(sum(T.HasDup)/DNSMPL )  ;  % the min value
% this works, work with vectors; 
%   G = grpstats( T.HasDup ,  [T.InterMHDistR T.MHLenR ] )
parfor N = 1:NSamples
    idx_with_MTD = randsample( idx_has_dup ,N_MTDs_per_sample_min ) ; 
    
    has_dup = zeros(1,N_MHPs);
    has_dup(idx_with_MTD) = 1; 
    
    T_RS_HasDup{N} = has_dup ; 

end
%save('~/Downloads/downsample_MTdna.mat')
%
R = table(); 
R.InterMHDistR = repmat( T.InterMHDistR , NSamples , 1) ;
R.MHLenR       = repmat( T.MHLenR , NSamples , 1) ;
R.HasDup       = cell2mat(T_RS_HasDup)' ;

RG = grpstats( R , {'InterMHDistR' 'MHLenR'} ,'mean' , 'DataVars' , 'HasDup');
RGa = grpstats( R , {'MHLenR'} ,'mean' , 'DataVars' , 'HasDup');
RGb = grpstats( R , {'InterMHDistR'} ,'mean' , 'DataVars' , 'HasDup');

%save('~/Downloads/downsample_MTdnaR.mat' , 'R')
fprintf('Running took %0.02f minutes\n' , toc / 60 ) ;

%% histogram of inter-MH distance for gDNA and for mtDNA
fh = figure('units','centimeters','position',[5 5  10 7]) ;
t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');

nexttile;
hold on ; 
mtDNA = T0.InterMHDist(T0.DupCounts>0 & strcmp(T0.chr,'MT')) ;
gDNA  = T0.InterMHDist(T0.DupCounts>0 & ~strcmp(T0.chr,'MT')) ; 
histogram(mtDNA, 10 ,'Normalization','Probability')
histogram(gDNA , 10  ,'Normalization','Probability') 
legend( {'mtDNA' 'gDNA'} , 'Location', 'nw')
xlabel('inter-MH distance')
ylabel('Fraction of MTDs')
title('10k diploid')

% and for haploid
nexttile;
hold on ; 
mtDNA = HAP.InterMHDist(HAP.DupCounts>0 & strcmp(HAP.chr,'MT')) ;
gDNA  = HAP.InterMHDist(HAP.DupCounts>0 & ~strcmp(HAP.chr,'MT')) ; 
histogram(mtDNA, 10 ,'Normalization','Probability')
histogram(gDNA , 10  ,'Normalization','Probability') 
legend( {'mtDNA' 'gDNA'} , 'Location', 'ne')
xlabel('inter-MH distance')
ylabel('Fraction of MTDs')
title('haploid')


%% ksdensity
fh = figure('units','centimeters','position',[5 5  10 7]) ;
t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');

nexttile;
hold on ; 
mtDNA = T0.InterMHDist(T0.DupCounts>0 & strcmp(T0.chr,'MT')) ;
gDNA  = T0.InterMHDist(T0.DupCounts>0 & ~strcmp(T0.chr,'MT')) ; 
ksdensity(mtDNA)
ksdensity(gDNA ) 
legend( {'mtDNA' 'gDNA'} , 'Location', 'nw')
xlabel('inter-MH distance')
ylabel('Fraction of MTDs')
title('10k diploid')

% and for haploid
nexttile;
hold on ; 
mtDNA = HAP.InterMHDist(HAP.DupCounts>0 & strcmp(HAP.chr,'MT')) ;
gDNA  = HAP.InterMHDist(HAP.DupCounts>0 & ~strcmp(HAP.chr,'MT')) ; 
ksdensity(mtDNA)
ksdensity(gDNA ) 
legend( {'mtDNA' 'gDNA'} , 'Location', 'ne')
xlabel('inter-MH distance')
ylabel('Fraction of MTDs')
title('haploid')

%% Figure 


clr = [0.8 0.3 0.3] ; 
fh = figure('units','centimeters','position',[5 5  10 14]) ;
t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');

for I  = 4:6
    nexttile;
    hold on ;
    idx = G0.MHLenR == I ;
    plot( G0.InterMHDistR(idx) , 100*G0.mean_HasDup(idx) ,'ok' ,'DisplayName' , 'gDNA','MarkerFaceColor',[.8 .8 .8])
    idx = RG.MHLenR == I ;
    plot( RG.InterMHDistR(idx) , 100*RG.mean_HasDup(idx) ,'o','Color',clr ,'DisplayName' , 'mtDNA','MarkerFaceColor',clr)
    set(gca,'yscale','log')
    title( sprintf('MH length = %d' , I ) )
    xlabel('nt between MHPair')
    set(gca,'ytick',[1e-4 1e-3 1e-2 1e-1 1])
    legend('location','best')
    ylim([2e-3 5e-1])
    xlim([20 300])
    
end

nexttile
hold on ;
plot( G0a.MHLenR , G0a.mean_HasDup*100 ,'o-k','DisplayName' , 'gDNA','MarkerFaceColor',[.8 .8 .8]) ; 
plot( RGa.MHLenR , RGa.mean_HasDup*100 ,'o-','Color',clr ,'DisplayName' , 'mtDNA','MarkerFaceColor',clr) ; 
set(gca,'xtick',4:7)
xlim([3.5 7.5])
xlabel('MH length (nt)')
set(gca,'yscale','lin')
legend('location','nw')
set(gca,'xticklabel',{'4' '5' '6' '>=7'});

ylabel(t,'% of MHPs with an MTD','FontSize',15)
print( '-dpng' , [ FIGNAME '' ] , '-r300' );
%close;