PROJDIR = '/Users/lcarey/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/' ;
FILE  = 'Sarah/MH_project/ProcessedData/10k.sign.count.tsv' ; 
T = readtable( [PROJDIR FILE ] , 'FileType', 'text' , 'ReadVariableNames',false,'Delimiter','\t');
T.Properties.VariableNames = {'chr' 's1' 'e1' 's2' 'e2' 'NDup' 'NColl'} ;
T.HasDup = T.NDup > 0 ; 

T.MHlen = T.e1 - T.s1 ; 
T.InterMHlen = T.s2 - T.e1 ; 

T.MHlen(T.MHlen>11)=11;
%%
G = grpstats( T , 'MHlen' , {'mean' 'sum'}, 'DataVars' , 'HasDup') ; 
%%
fh = figure('units','centimeters','position',[5 5  5 5]) ;
plot( G.MHlen , 100* G.mean_HasDup , '-ok')
xlabel('MH length (nt)')
ylabel('% of MHRs w/MTD')
set(gca,'ytick',0:.2:1)
set(gca,'xtick',0:11)
grid on ;

%%
G = grpstats( T , 'InterMHlen' , {'mean' 'sum'}, 'DataVars' , 'HasDup') ; 
%%
fh = figure('units','centimeters','position',[5 5  9 5]) ;
plot( G.InterMHlen , 100* G.mean_HasDup , '-.k')
xlabel('Inter-MH distance (nt)')
ylabel('% of MHRs w/MTD')
%set(gca,'ytick',0:.2:1)
set(gca,'xtick',0:50:1000)
grid on ;
xlim([0 400])

%% clustering
T = sortrows(T,'s1','ascend');
idx = find(T.HasDup);
idx_rnd = sort( randsample( height(T) , numel(idx)) , 'ascend') ;


real_distances = diff(T.s1(idx)) ;
rand_distances = diff(T.s1(idx_rnd)) ;
%%
Y1 = log10( double(rand_distances )) ; 
Y2 = log10( double(real_distances )) ; 
fh = figure('units','centimeters','position',[5 5  9 6]) ;
hold on ;
histogram( Y1 , 50 ,'FaceColor','k')
histogram( Y2 , 50)
axis tight; 
ylabel('# of inter-duplication events')
xlabel('nt between MHRs')
legend({'shuffled' 'real'},'location','nw')
set(gca,'xtick',0:5)
set(gca,'xticklabel',[1 10 100 1e3 1e4 1e5])

fh = figure('units','centimeters','position',[5 5  9 6]) ;
hold on ;
[f,x]=ecdf( Y1 ) ; plot(x,f,'-k','LineWidth',3)
[f,x]=ecdf( Y2 ); plot(x,f,'-','LineWidth',3)
axis tight; 
ylabel('# of inter-duplication events')
xlabel('nt between MHRs')
legend({'shuffled' 'real'},'location','nw')
set(gca,'xtick',0:5)
set(gca,'xticklabel',[1 10 100 1e3 1e4 1e5])




fh = figure('units','centimeters','position',[5 5  9 9]) ;
hold on ; 
bar( 1, 100 *  mean(Y1<log10(100)) , 'FaceColor'  , 'k')
bar( 2, 100 *  mean(Y2<log10(100)) , 'FaceColor'  , 'r')
set(gca,'xtick',1:2)
set(gca,'xticklabel',{'random' 'real'})
ylabel('% of distances < 100nt')
xlim([0.4 2.6])
