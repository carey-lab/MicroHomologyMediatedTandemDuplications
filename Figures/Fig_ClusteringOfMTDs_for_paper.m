% a local clustering on the scale of 100nt
% MHPs that are close to an observed MTD are more likely to have an MTD
DIR = '/Users/lcarey/SynologyDrive/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/Sarah/MH_project/ProcessedData/';
%MTD = readtable([ DIR '10k.sign.count.tsv'] , 'FileType','text','Format','%s%d%d%d%d%d%d');
%MTD.Properties.VariableNames = {'chr' 's1' 'e1' 's2' 'e2' 'NDups' 'NCol'} ;

FIGDIR  = '~/Nutstore Files/Microhomology shared folder/Figures/Fig2 - cis-determinants of MTD through ultra-deep sequencing/' ;
FIGBASENAME = [FIGDIR 'MTD_Fig_ClusteringOfMHPs_With_MTDs_' ] ; 

D = readtable( [DIR 'distance_to_closest_MHR_with_dup.bed'] , 'FileType','text','Format','%s%d%d%d%d');
A = readtable( [DIR '10k.sign.count.tsv.MN_Spombe_972h-_Rep1.features.txt'], 'FileType','text','TreatAsEmpty','NA');
closest_MTD = D.Var5 ; 
HasDup = A.Var1 > 0 ; 
clear D A ; 

% access to very long tables is slow. use vectors instead
[closest_MTD , o] = sort(closest_MTD,'ascend');
HasDup = HasDup(o);
%% bar plot of % of MHPs with an MTD within 100nt, vs >100nt
%mypool = parpool();
%opt = statset('UseParallel',true);

within_100 = mean(HasDup(closest_MTD<=100))*100 ; 
not_within_100 = mean(HasDup(closest_MTD>100))*100 ; 
% too slow
%tic ; within_100 = 100 * bootstrp( 100 , @mean , within_100,'Options',opt) ; toc
%tic ; not_within_100 = 100 * bootstrp( 1e3 , @mean , not_within_100,'Options',opt) ; toc

%%
fh = figure('units','centimeters','position',[50 50 5 9]) ;
clrs = get(gca,'ColorOrder');
bar(  [0.9 2.1] , [within_100 not_within_100] , 'FaceColor',clrs(4,:)  );
ylabel('% of MHPs with an MTD')
set(gca,'xtick', [0.9 2.1] )
set(gca,'xticklabel',{'<100nt' '>100nt'});
xlabel({'Distance to nearest' 'MHP with an MTD'})
xlim([0.4 2.6])
box off ;
print('-dpng',[FIGBASENAME 'bars'] , '-r300');
close ;
%%  several different methods for binning, clustering, plotting
%  ksdensity looks the best

[fYkl,xYkl] = ksdensity(log10(double(closest_MTD(HasDup)))) ;
[fNkl,xNkl] = ksdensity(log10(double(closest_MTD(~HasDup))));


%[fY,xY] = ecdf(closest_MTD(HasDup)) ;
%[fN,xN] = ecdf(closest_MTD(~HasDup)) ;

%[fYkl,xYkl] = ksdensity(closest_MTD(HasDup),1:1e4);
%[fNkl,xNkl] = ksdensity(closest_MTD(~HasDup),1:1e4);


% window_size = 1e5;
% xl = 1:1e3:(numel(closest_MTD)-window_size) ; 
% y = NaN(numel(xl),1);
% x = NaN(numel(xl),1);
% parfor I = 1:numel(xl)
%     idx = xl(I):xl(I)+window_size ; 
%     y(I) = mean(HasDup( idx ) );
%     x(I) = mean(closest_MTD( idx ) );
% end
% xl = log10(x);


%% plot ksdensity histograms of distances for MHPs with or without an MTD

fh = figure('units','centimeters','position',[50 50 9 9]) ;
hold on ;
plot(xYkl,fYkl,'-','DisplayName','Has MTD','LineWidth',4)
plot(xNkl,fNkl,'-','DisplayName','no MTD','LineWidth',4)
legend('location','nw')
set(gca,'xtick',log10([3 10 25 50 100 300  1e3 3e3])) 
set(gca,'xticklabel',[3 10 25 50 100  300  1e3 3e3]) 
xlabel('nt to the nearest MHP with an observed MTD')
ylabel('Density of MHPs with or without an MTD')
xlim(log10([2  9e3]))

print('-dpng',[FIGBASENAME 'ksdensity'] , '-r300');
close ;
%% calculate the pct of MHPs with an MTD within each window
% discarded for ksdensity()

% %%
% %fh = figure('units','centimeters','position',[50 50 10 10]) ;
% fh = figure('units','centimeters','position',[50 50 15 6]) ;
% hold on; 
% plot(xl,y*100,'.','Color','k','LineWidth',2)
% ylabel('% of MHPs with an MTD')
% xlabel('nt to the nearest MHPair with an observed MTD')
% set(gca,'xtick',[1 2 3])
% set(gca,'xticklabel',[10 100 1000]) 
% axis tight; 
% xlim([ log10(5) log10(550)])
% set(gca,'xtick',[1 log10(25) log10(50) 2 log10(150) log10(250) log10(500) ])
% set(gca,'xticklabel',[10 25 50 100 150 250 500]) 
% line(xlim,[mean(HasDup) mean(HasDup)]*100,'Color',[.5 .5 .5],'LineStyle','--')
% print('-dpng',[FIGBASENAME 'A'] , '-r300');
% close ; 
% 
% fh = figure('units','centimeters','position',[50 50 20 8]) ;
% hold on; 
% plot(xl,log2(y./mean(HasDup)),'.k')
% line(xlim,[0 0],'Color',[.5 .5 .5],'LineStyle','--')
% ylabel('log2( this-window / whole-genome )')
% xlabel('nt to the nearest MHP with an MTD')
% set(gca,'xtick',[1 log10(25) log10(50) 2 log10(150) log10(250) log10(500) ])
% set(gca,'xticklabel',[10 25 50 100 150 250 500]) 
% axis tight; 
% grid on ;
% xlim([ log10(5) log10(550)])
% print('-dpng',[FIGBASENAME 'B'] , '-r300');
% close ; 