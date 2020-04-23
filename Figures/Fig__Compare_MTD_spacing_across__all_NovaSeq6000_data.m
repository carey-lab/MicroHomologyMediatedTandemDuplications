
%% load data
DATADIR = '~/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/DataFromCluster/NovaSeq6000/';

FIGNAME = '~/Downloads/Fig__Compare_MTD_spacing_across__all_NovaSeq6000_data' ; 

T5= readtable( [DATADIR 'Spombe_NovaSeq6000_B5.sign.count.tsv'],'FileType','text','Format','%s%d%d%d%d%f%d');
T6= readtable( [DATADIR 'Spombe_NovaSeq6000_B6.sign.count.tsv'],'FileType','text','Format','%s%d%d%d%d%f%d');
T7= readtable( [DATADIR 'Spombe_NovaSeq6000_B7.sign.count.tsv'],'FileType','text','Format','%s%d%d%d%d%f%d');
T8= readtable( [DATADIR 'Spombe_NovaSeq6000_B8.sign.count.tsv'],'FileType','text','Format','%s%d%d%d%d%f%d');
%% set variable names
T5.Properties.VariableNames = {'chr' 's1' 'e1' 's2' 'e2' 'DupCounts' 'ColCounts'};
T6.Properties.VariableNames = {'chr' 's1' 'e1' 's2' 'e2' 'DupCounts' 'ColCounts'};
T5.InterMHDist = T5.s2 - T5.e1 + 1;  
T6.InterMHDist = T6.s2 - T6.e1 + 1;  
T5.HasDup = T5.DupCounts>0;
T6.HasDup = T6.DupCounts>0;
T5.MHLen = T5.e1 - T5.s1 +1 ;
T6.MHLen = T6.e1 - T6.s1 +1 ;


T7.Properties.VariableNames = {'chr' 's1' 'e1' 's2' 'e2' 'DupCounts' 'ColCounts'};
T8.Properties.VariableNames = {'chr' 's1' 'e1' 's2' 'e2' 'DupCounts' 'ColCounts'};
T7.InterMHDist = T7.s2 - T7.e1 + 1;  
T8.InterMHDist = T8.s2 - T8.e1 + 1;  
T7.HasDup = T7.DupCounts>0;
T8.HasDup = T8.DupCounts>0;
T7.MHLen = T7.e1 - T7.s1 +1 ;
T8.MHLen = T8.e1 - T8.s1 +1 ;
%% grp stats
N=20 ; 
T5.InterMHDistR = round(T5.InterMHDist./N)*N;
T5.InterMHDistR(T5.InterMHDistR>300)=300 ; 
T5.MHLenR = T5.MHLen ; T5.MHLenR(T5.MHLenR>7) = 7 ; 

T6.InterMHDistR = round(T6.InterMHDist./N)*N;
T6.InterMHDistR(T6.InterMHDistR>300)=300 ; 
T6.MHLenR = T6.MHLen ; T6.MHLenR(T6.MHLenR>7) = 7 ; 

T7.InterMHDistR = round(T7.InterMHDist./N)*N;
T7.InterMHDistR(T7.InterMHDistR>300)=300 ; 
T7.MHLenR = T7.MHLen ; T7.MHLenR(T7.MHLenR>7) = 7 ; 

T8.InterMHDistR = round(T8.InterMHDist./N)*N;
T8.InterMHDistR(T8.InterMHDistR>300)=300 ; 
T8.MHLenR = T8.MHLen ; T8.MHLenR(T8.MHLenR>7) = 7 ; 

G5 = grpstats(T5(~strcmp(T5.chr,'MT'),:),{'InterMHDistR' },{'mean' 'sum'},'DataVars',{'HasDup' });
G6 = grpstats(T6(~strcmp(T6.chr,'MT'),:),{'InterMHDistR'},{'mean' 'sum'},'DataVars',{'HasDup'});
G7 = grpstats(T7(~strcmp(T5.chr,'MT'),:),{'InterMHDistR' },{'mean' 'sum'},'DataVars',{'HasDup' });
G8 = grpstats(T8(~strcmp(T6.chr,'MT'),:),{'InterMHDistR'},{'mean' 'sum'},'DataVars',{'HasDup'});

G = G5;
G.PctDup_5 = G5.mean_HasDup * 100 ; 
G.PctDup_6 = G6.mean_HasDup * 100 ; 
G.PctDup_7 = G7.mean_HasDup * 100 ; 
G.PctDup_8 = G8.mean_HasDup * 100 ; 
G.GroupCount = []; 
G.mean_HasDup = [];
G.sum_HasDup = [] ; 

%%
figure; hold on ;
plot(G.InterMHDistR ,  G.PctDup_5 ,'-s','DisplayName' , 'Diploid LD388 midlog' ,'MarkerFaceColor',[.7 .7 .7]);
plot(G.InterMHDistR ,  G.PctDup_6 ,'-o','DisplayName' , 'Haploid DY1778 midlog');
plot(G.InterMHDistR ,  G.PctDup_7 ,'-s','DisplayName' , 'Diploid LD388 saturated','MarkerFaceColor',[.7 .7 .7]);
plot(G.InterMHDistR ,  G.PctDup_8 ,'-o','DisplayName' , 'Haploid DY1778 saturated');
legend('location','best')
xlabel('Inter-MH Distance')
ylabel('% of MHPs with an MTD')
%%
d5 = T5.InterMHDist(T5.HasDup & ~strcmp(T5.chr,'MT'));
d6 = T6.InterMHDist(T6.HasDup & ~strcmp(T6.chr,'MT'));
d7 = T7.InterMHDist(T7.HasDup & ~strcmp(T7.chr,'MT'));
d8 = T8.InterMHDist(T8.HasDup & ~strcmp(T8.chr,'MT'));
d = vertcat(d5,d6,d7,d8);

figure; hold on ;
[f,x]=ecdf(d5); plot(x,f,'-or','DisplayName' , 'Diploid LD388 midlog' );
[f,x]=ecdf(d6); plot(x,f,'-ob','DisplayName' , 'Haploid DY1778 midlog');
[f,x]=ecdf(d7); plot(x,f,'-sr','DisplayName' , 'Diploid LD388 saturated','MarkerFaceColor',[.7 .7 .7]);
[f,x]=ecdf(d8); plot(x,f,'-sb','DisplayName' , 'Haploid DY1778 saturated' ,'MarkerFaceColor',[.7 .7 .7]);
[f,x]=ecdf(T8.InterMHDist); plot(x,f,'-k','DisplayName' , 'all MHPs');
[f,x]=ecdf(d); plot(x,f,'-g','DisplayName' , 'all MTDs','LineWidth',2);

legend('location','best')
xlabel('Inter-MH Distance')
ylabel('Cummulative fraction of MHPs with an MTD')
grid on ;


figure; hold on ;
[f,x]=ksdensity(d5,0:2:500); plot(x,f,'-or','DisplayName' , 'Diploid LD388 midlog' );
[f,x]=ksdensity(d6,0:2:500); plot(x,f,'-ob','DisplayName' , 'Haploid DY1778 midlog');
[f,x]=ksdensity(d7,0:2:500); plot(x,f,'-sr','DisplayName' , 'Diploid LD388 saturated','MarkerFaceColor',[.7 .7 .7]);
[f,x]=ksdensity(d8,0:2:500); plot(x,f,'-sb','DisplayName' , 'Haploid DY1778 saturated' ,'MarkerFaceColor',[.7 .7 .7]);
[f,x]=ksdensity(T8.InterMHDist,0:500); plot(x,f,'-k','DisplayName' , 'all MHPs');
[f,x]=ksdensity(d,0:500); plot(x,f,'-g','DisplayName' , 'all MTDs','LineWidth',2);

legend('location','best')
xlabel('Inter-MH Distance')
ylabel('Fraction of MHPs with an MTD')
grid on ;


%%
xl = 0:50:526 ; 
figure; hold on ;
[f,x]=hist(d5,xl); plot(x,f,'-or','DisplayName' , 'Diploid LD388 midlog' );
[f,x]=hist(d6,xl); plot(x,f,'-ob','DisplayName' , 'Haploid DY1778 midlog');
[f,x]=hist(d7,xl); plot(x,f,'-sr','DisplayName' , 'Diploid LD388 saturated','MarkerFaceColor',[.7 .7 .7],'LineWidth',2);
[f,x]=hist(d8,xl); plot(x,f,'-sb','DisplayName' , 'Haploid DY1778 saturated' ,'MarkerFaceColor',[.7 .7 .7],'LineWidth',2);
%[f,x]=hist(T8.InterMHDist,xl); plot(x,f,'-k','DisplayName' , 'all MHPs');
%[f,x]=hist(d,xl); plot(x,f,'-g','DisplayName' , 'all MTDs','LineWidth',2);

legend('location','best')
xlabel('Inter-MH Distance')
ylabel('# of MHPs with an MTD')
grid on ;