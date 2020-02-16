% plot the measured and predicted MTDs across the genome, and show a hot and cold spot


% load data and set output dirs
DATADIR = '/Users/lcarey/SynologyDrive/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/Sarah/MH_project/Manuscript-todo/ProcessedData/' ;
FIGDIR  = '~/Nutstore Files/Microhomology shared folder/Figures/Fig2 - cis-determinants of MTD through ultra-deep sequencing/' ;
FIGBASENAME = [FIGDIR 'MTD_observed_predict_across_chromosome_' ] ; 


% load whole-genome 1kb window data

T = readtable( [ DATADIR 'loc_obs_pre1k.txt' ] );
T.loc = T.loc ./ 1000  ;
locs_of_chr_breaks = find(diff(T.loc)<0) ; 
T.chr = NaN(height(T),1);
T.chr(1:locs_of_chr_breaks(1))=1;
T.chr(locs_of_chr_breaks(1)+1:locs_of_chr_breaks(2))=2;
T.chr(locs_of_chr_breaks(2)+1:end)=3;

T = T(1:locs_of_chr_breaks(1),:);
T.locR = round(T.loc./50)*50 ;


% for finding hot and cold spots
R = grpstats(T,{'chr' 'locR'},'mean');
R.obs_over_mhr = R.mean_sumObs ./ R.mean_sumMHR ; 
R = sortrows( R , 'obs_over_mhr' ,'ascend') ;

%% plot MTDs across chromosome I
clrs = cbrewer('qual','Set1',5);
fh = figure('units','centimeters','position',[5 5  30 8]) ;
t = tiledlayout(3,1,'TileSpacing','none','Padding','none') ; 
nexttile
plot( T.loc , T.sumMHR ,'Color',clrs(1,:) )
ylabel('# of MHPs')
axis tight; 
ylim([1501 3200])
set(gca,'YTick', 0:1e3:1e5 )
set(gca,'TickLength',[0 0])

nexttile
plot( T.loc , T.sumObs ,'Color',clrs(2,:) )
ylabel({'observed' '# of MTDs'})
axis tight; 
ylim([-0.5  6.5])
set(gca,'TickLength',[0 0])

nexttile
plot( T.loc , T.sumPre ,'Color',clrs(3,:) )
ylabel({'predicted' '# of MTDs'})
axis tight; 
ylim([-0.5  6.5])
set(gca,'TickLength',[0 0])
% 
% nexttile
% plot( T.loc , T.sumFloat ,'Color',clrs(4,:) )
% xlabel('Position on chromosome (kb)')
% ylabel('MTD prediction score')
% axis tight; 

title(t,'MHPairs and predicted & observed MTDs across chromosome I')
xlabel(t,'Position on chromosome (kb)')
print('-dpng',[FIGBASENAME 'A'],'-r150');
close ; 

%%
hot_window = [3950 4030]  ; 
Q = T(T.loc > hot_window(1) & T.loc < hot_window(2)  ,:);

R(R.obs_over_mhr>prctile(R.obs_over_mhr,95) &  R.mean_sumObs>prctile(R.mean_sumObs,95)  ,:)
Q = T(T.locR == 4000,:);

xt = 0:20:1e5 ; 

fh = figure('units','centimeters','position',[5 5  7 10]) ;
t = tiledlayout(3,1,'TileSpacing','none','Padding','none') ; 
nexttile
plot( Q.loc , Q.sumMHR ,'Color',clrs(1,:) )
ylabel('# of MHPs')
axis tight; 
set(gca,'xtick',xt);
%rectangle('Position',[gene_window(1) , min(ylim) , diff(gene_window) , max(ylim)])

nexttile
plot( Q.loc , Q.sumObs ,'Color',clrs(2,:) )
ylabel('observed MTDs')
axis tight; 
ylim([-0.5  6.5])
set(gca,'xtick',xt);

nexttile
plot( Q.loc , Q.sumPre ,'Color',clrs(3,:) )
ylabel('predicted MTDs')
axis tight; 
ylim([-0.5  6.5])
set(gca,'xtick',xt);

% 
% nexttile
% plot( T.loc , T.sumFloat ,'Color',clrs(4,:) )
% ylabel('MTD prediction score')
% axis tight; 

xlabel(t,'Position on chromosome (kb)')
title(t,'hot spot in chrI')
print('-dpng',[FIGBASENAME 'H'],'-r300');
close ; 


%% Find cold spot
R(R.obs_over_mhr<prctile(R.obs_over_mhr,10) &  R.mean_sumObs<prctile(R.mean_sumObs,5)  ,:)


Q = T(T.locR ==  1300 ,:);
xt = 0:20:1e5 ; 

fh = figure('units','centimeters','position',[5 5  7 10]) ;
t = tiledlayout(3,1,'TileSpacing','none','Padding','none') ; 
nexttile
plot( Q.loc , Q.sumMHR ,'Color',clrs(1,:) )
ylabel('# of MHPs')
axis tight; 
set(gca,'xtick',xt);
%rectangle('Position',[gene_window(1) , min(ylim) , diff(gene_window) , max(ylim)])

nexttile
plot( Q.loc , Q.sumObs ,'Color',clrs(2,:) )
ylabel('observed MTDs')
axis tight; 
ylim([-0.5  6.5])
set(gca,'xtick',xt);

nexttile
plot( Q.loc , Q.sumPre ,'Color',clrs(3,:) )
ylabel('predicted MTDs')
axis tight; 
ylim([-0.5  6.5])
set(gca,'xtick',xt);

% 
% nexttile
% plot( T.loc , T.sumFloat ,'Color',clrs(4,:) )
% ylabel('MTD prediction score')
% axis tight; 

xlabel(t,'Position on chromosome (kb)')
title(t,'cold spot in chrI')
print('-dpng',[FIGBASENAME 'C'],'-r300');
close ; 