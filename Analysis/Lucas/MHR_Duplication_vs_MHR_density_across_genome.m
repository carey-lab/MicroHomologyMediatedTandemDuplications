T = readtable('~/Downloads/chrI_10k.sign.count.tsv','FileType','text','Format','%s%d%d%d%d%d%d');
T.Properties.VariableNames = {'chr' 's1' 'e1' 's2' 'e2' 'NDup' 'Ncol'} ;
T.Ncol = [] ;
T.HasDup = T.NDup > 0  ; 
T.MHlen = T.e1 - T.s1 ; 
%%
T.kb01 = round( mean([T.s1 T.e2],2) ./ 100 ) * 0.1 ; 
T.kb05 = round( mean([T.s1 T.e2],2) ./ 500 ) * 0.5 ; 
T.kb1  = round( mean([T.s1 T.e2],2) ./ 1000 ) ; 
T.kb5  = round( mean([T.s1 T.e2],2) ./ 5000 ) * 5 ; 
T.kb10 = round( mean([T.s1 T.e2],2) ./ 10000 ) * 10 ; 
T.kb20 = round( mean([T.s1 T.e2],2) ./ 20000 ) * 20 ; 

%%
G01 = grpstats( T , 'kb01' , {'mean' 'sum'}, 'DataVars' , 'HasDup') ; 
G05 = grpstats( T , 'kb05' , {'mean' 'sum'}, 'DataVars' , 'HasDup') ; 
G1 = grpstats( T , 'kb1' , {'mean' 'sum'}, 'DataVars' , 'HasDup') ; 
G5 = grpstats( T , 'kb5' , {'mean' 'sum'}, 'DataVars' , 'HasDup') ; 
G10 = grpstats( T , 'kb10' ,{ 'mean' 'sum'}, 'DataVars' , 'HasDup') ; 
G20 = grpstats( T , 'kb20' , {'mean' 'sum'}, 'DataVars' , 'HasDup') ; 

%%
fh = figure('units','centimeters','position',[5 5 40 40]) ;

subplot(2,1,1)
plotyy( G1.kb1 , G1.GroupCount , G1.kb1 , G1.sum_HasDup )
legend({'# of MHRs' '# of Duplications'}) 
xlabel('Position along chr I (kb)')
ylabel('# of MHRs or Duplications')
title('1 kb windows')

subplot(2,1,2)
[ax,h1,h2] = plotyy( G10.kb10 , G10.GroupCount , G10.kb10 , G10.sum_HasDup )
legend({'# of MHRs' '# of Duplications'}) 
xlabel('Position along chr I (kb)')
ylabel('# of MHRs or Duplications')
title('10 kb windows')
set(ax(2),'ytick',0:2:20)

%%
fh = figure('units','centimeters','position',[5 5 40 40]) ;
Y1 = 100*( G1.sum_HasDup ./ G1.GroupCount)  ; 
idx = Y1>(mean(Y1)+1*std(Y)) | Y1<(mean(Y1)-1*std(Y1)) ; 

subplot(2,1,1) ; hold on ;
plot( G1.kb1 , Y1, '.','Color',[.7 .7 .7]);
plot( G1.kb1(idx) , Y1(idx) , '.r');
xlabel('Position along chr I (kb)')
ylabel('% of MHRs w/a Duplication')
title('1 kb windows')
legend({'all' '>1std away from mean'})
set(gca,'ytick',0:0.05:1)
grid on ; 

%Y = log2( (1+G10.sum_HasDup) ./ (1+G10.GroupCount))  ; 
Y = 100*( G10.sum_HasDup ./ G10.GroupCount)  ; 
%Y  = Y ./ modefit(Y) ; 
idx = Y>(mean(Y)+1*std(Y)) | Y<(mean(Y)-1*std(Y)) ; 
subplot(2,1,2); hold on ;
plot( G10.kb10 , Y , '.','Color',[.7 .7 .7]);
plot( G10.kb10(idx) , Y(idx) , '.r');
xlabel('Position along chr I (kb)')
ylabel('% of MHRs w/a Duplication')
title('10 kb windows')
legend({'all' '>1std away from mean'})
set(gca,'ytick',0:0.05:1)
grid on ; 

figure; hold on; 
histogram( Y1 , 25 , 'Normalization','Probability')
histogram( Y , 25 , 'Normalization','Probability')
ylim([0 0.1])

%%
gs = {'G01' 'G05' 'G1' 'G5' 'G10' 'G20'};
gt = {'100bp' '500bp' '1kb' '5kb' '10kb' '20kb'};
figure; 
for I = 1:numel(gs)
subplot(2,3,I)
G = eval(gs{I});
boxplot(G.GroupCount, G.sum_HasDup,'symbol','','plotstyle','compact')
grid on ; 
title(gt{I})
xlabel('# of Duplications')
ylabel('# of MHRs')
axis tight; 
end


%%

fh = figure('units','centimeters','position',[5 5 10 10]) ;
X = G10.GroupCount  ;% X(X<prctile(X,1))=prctile(X,1) ;  X(X>prctile(X,99))=prctile(X,99) ;
Y = G10.sum_HasDup ;% Y(Y<prctile(Y,1))=prctile(Y,1) ;  Y(Y>prctile(Y,99))=prctile(Y,99) ; 
scatter(X,Y,'k','MarkerFaceColor',[.7 .7 .7],'MarkerFaceAlpha',0.5);
corr(X,Y)
title('10 kb')
xlabel('# of MHRs')
ylabel('# of duplications')

fh = figure('units','centimeters','position',[5 5 10 10]) ;
X = G1.GroupCount  ;% X(X<prctile(X,1))=prctile(X,1) ;  X(X>prctile(X,99))=prctile(X,99) ;
Y = G1.sum_HasDup ;% Y(Y<prctile(Y,1))=prctile(Y,1) ;  Y(Y>prctile(Y,99))=prctile(Y,99) ; 
scatter(X,Y,'k','MarkerFaceColor',[.7 .7 .7],'MarkerFaceAlpha',0.5);
corr(X,Y)
title('1 kb')
xlabel('# of MHRs')
ylabel('# of duplications')


%%
pL = prctile(G10.sum_HasDup,10) 
pH = prctile(G10.sum_HasDup,90) 

figure; hold on ;
histogram( G10.GroupCount(G10.sum_HasDup<=pL),25,'Normalization','Probability')
histogram( G10.GroupCount(G10.sum_HasDup>pH),25,'Normalization','Probability')
%%
pL = prctile(G1.sum_HasDup,10) 
pH = prctile(G1.sum_HasDup,90) 
pL = 0 ;  pH = 4 ; % visual

fh = figure('units','centimeters','position',[5 5 8 8]) ;
hold on ;
X = G1.GroupCount(G1.sum_HasDup<=pL); 
Y = G1.GroupCount(G1.sum_HasDup>=pH) ; 
histogram( X,25,'Normalization','Probability')
histogram( Y,25,'Normalization','Probability')
xlim([1500 3500])
xlabel('# of MHRs in the 1kb window')
ylabel('# of genomic loci')
legend({'0 duplications' '>= 4 duplications'})
[p,h,teststats] = ranksum(X,Y)
