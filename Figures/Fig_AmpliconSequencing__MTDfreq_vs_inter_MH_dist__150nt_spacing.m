% Do the amplicon-sequencing results have the same 150nt spacing? 
%   they were done in haploids!!
% 
% 
%% load data
A = 'lib.sign.norm.tsv' ;
DATADIR = '~/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/DataFromCluster/' ;

FIGDIR = '~/Nutstore Files/Microhomology shared folder/Figures/Supplementary Figures/' ; 
%FIGNAME = [ FIGDIR 'SupFig__Haploid_vs_Diploid_vs_Ecoli_MTDfreq__InterMH_Distance'  ] ; 


A = readtable( [DATADIR A] ,'FileType','text','TreatAsEmpty','-');
A.Properties.VariableNames =  {'chr'      's1'      'e1'      's2'    'e2'   'Nreads' 'NDup' 'DupFrq' 'NCol' 'ColFrq'} ; 

A.HasDup = A.NDup>0;
A.MHLen = A.e1 - A.s1 + 1 ; 
A.InterMHDist = A.s2 - A.e1 + 1 ; 
%%
%%
n=5; 
A.MHLenR = A.MHLen;
A.MHLenR(A.MHLenR>6)=6 ;
A.InterMHDistR = round(double(A.InterMHDist)./n)*n ;


G = grpstats( A , {'MHLenR' 'InterMHDistR' } , {'mean' 'sum'} , 'DataVars' , {'HasDup' 'DupFrq'} );
G2 = grpstats( A(A.HasDup,:) , {'MHLenR' 'InterMHDistR' } , {'mean' 'sum'} , 'DataVars' , {'DupFrq'} );

figure; 
gscatter(G.InterMHDistR,100*G.mean_HasDup,G.MHLenR)
xlabel('Inter-MH distance')
ylabel('% of MHPs with an MTD')
title('Amplicon sequencing')