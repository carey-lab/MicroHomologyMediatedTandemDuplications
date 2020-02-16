% Fig__Spombe_WholeGenome_10kDiploid_vs_SRR7817502haploid__EssentialGenes_SameRules
% we have two high-coverage whole-genome S. pombe datasets: 10k, & SRR7817502
% same MTDfrequency rules? 
% 
% SRR7817502 is liklely haploid. Fewer MTDs in essential genes? 
%
%

%% load data
F1 = '10k' ;
F2 = 'Ecoli' ;
DATADIR = '~/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/DataFromCluster/' ;
D = readtable( [DATADIR F1 '_rm.sign.count.tsv']  , 'FileType','text','Delimiter','\t','Format','%s%d%d%d%d%d%d');
H = readtable( [DATADIR F2 '.sign.count.tsv']  , 'FileType','text','Delimiter','\t','Format','%s%d%d%d%d%d%d');


D.Properties.VariableNames = {'chr' 's1' 'e1' 's2' 'e2' 'NDup' 'NCol'};
D.HasDup = D.NDup>0;
D.MHLen = D.e1 - D.s1 + 1 ; 
D.InterMHDist = D.s2 - D.e1 + 1 ; 


H.Properties.VariableNames = { 'chr' 's1' 'e1' 's2' 'e2' 'NDup' 'NCol'};
H.HasDup = H.NDup>0;
H.MHLen = H.e1 - H.s1 + 1 ; 
H.InterMHDist = H.s2 - H.e1 + 1 ; 

D.InterMHDistR = round(D.InterMHDist./10)*10 ;
H.InterMHDistR = round(H.InterMHDist./10)*10 ;

whos
%


figure; 

subplot(2,2,1)
hold on ;
GD = grpstats( D(D.MHLen==4,:) , 'InterMHDistR' ,'mean' ,'DataVars' ,'HasDup');
GH = grpstats( H(H.MHLen==4,:) , 'InterMHDistR' ,'mean' ,'DataVars' ,'HasDup');
plot( GD.InterMHDistR  , 100*GD.mean_HasDup ,'o-','DisplayName',F1)
plot( GH.InterMHDistR  , 100*GH.mean_HasDup ,'o-','DisplayName',F2)
legend('location','ne')
xlabel('Inter-MH distance (nt)')
ylabel('% of MHPs with duplication')
title('MHlength = 4')
xlim([0 375])

subplot(2,2,2)
hold on ;
GD = grpstats( D(D.MHLen==5,:) , 'InterMHDistR' ,'mean' ,'DataVars' ,'HasDup');
GH = grpstats( H(H.MHLen==5,:) , 'InterMHDistR' ,'mean' ,'DataVars' ,'HasDup');
plot( GD.InterMHDistR  , 100*GD.mean_HasDup ,'o-','DisplayName','10k diploid')
plot( GH.InterMHDistR  , 100*GH.mean_HasDup ,'o-','DisplayName','SRR7817502 haploid')
%legend('location','ne')
xlabel('Inter-MH distance (nt)')
ylabel('% of MHPs with duplication')
title('MHlength = 5')
xlim([0 375])

subplot(2,2,3)
hold on ;
GD = grpstats( D(D.MHLen==6,:) , 'InterMHDistR' ,'mean' ,'DataVars' ,'HasDup');
GH = grpstats( H(H.MHLen==6,:) , 'InterMHDistR' ,'mean' ,'DataVars' ,'HasDup');
plot( GD.InterMHDistR  , 100*GD.mean_HasDup ,'o-','DisplayName',F1)
plot( GH.InterMHDistR  , 100*GH.mean_HasDup ,'o-','DisplayName',F2)
legend('location','ne')
xlabel('Inter-MH distance (nt)')
ylabel('% of MHPs with duplication')
title('MHlength = 6')
xlim([0 375])

subplot(2,2,4)
hold on ;
GD = grpstats( D(D.MHLen>=7,:) , 'InterMHDistR' ,'mean' ,'DataVars' ,'HasDup');
GH = grpstats( H(H.MHLen>=7,:) , 'InterMHDistR' ,'mean' ,'DataVars' ,'HasDup');
plot( GD.InterMHDistR  , 100*GD.mean_HasDup ,'o-','DisplayName','10k diploid')
plot( GH.InterMHDistR  , 100*GH.mean_HasDup ,'o-','DisplayName','SRR7817502 haploid')
%legend('location','ne')
xlabel('Inter-MH distance (nt)')
ylabel('% of MHPs with duplication')
title('MHlength >= 7')
xlim([0 375])
