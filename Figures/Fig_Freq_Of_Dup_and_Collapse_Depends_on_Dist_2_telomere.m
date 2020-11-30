% Figure for rebuttal
%   showing that the frequency of duplications and collapse depends on the
%   distance to the telomere

figpath = '~/Downloads/'; 
figname = [ figpath 'Fig_Freq_Of_Dup_and_Collapse_Depends_on_Dist_2_telomere__' ];

Sp_Col_NonSub = 0.0077586 ;
Sp_Col_Sub  = 0.0875000 ;
Sp_Dup_NonSub = 0.0099137 ;
Sp_Dup_Sub    = 0.145258 ;

Sp_midlog_col =	11.15702479338843 ;
Sp_midlog_dup =	12.809917355371898; 
Sp_saturated_col =	10.78512396694215 ;
Sp_saturated_dup =	16.65289256198347 ;

%%
figure( 'Position',[ 50 50 150 250] ); 
hold on ; 
bar( 1.07 , Sp_Col_NonSub ,'FaceColor' ,  [0.3 0.9 0.9])
bar( 1.93 , Sp_Col_Sub ,'FaceColor' , [0.7 0.1 0.1])
bar( 3.07 , Sp_Dup_NonSub ,'FaceColor' , [0.3 0.9 0.9])
bar( 3.93 , Sp_Dup_Sub ,'FaceColor' , [0.7 0.1 0.1])
set(gca,'xtick',[1.5 3.5])
set(gca,'xticklabel',{'Collapse' 'Duplication'})
xlim([0.5 4.5])
ylabel('% of MHPs with a')
legend({'non-subtelomeric' 'subtelomeric'} , 'location','NorthOutside')
print('-dpng',[figname 'a.png'],'-r300')

%
figure( 'Position',[ 50 50 150 250] ); 
hold on ; 
bar( 1.07 , Sp_midlog_col ,'FaceColor' ,  [0.9 0.3 0.9])
bar( 1.93 , Sp_midlog_dup ,'FaceColor' , [0.1 0.1 0.1])
bar( 3.07 , Sp_saturated_col ,'FaceColor' , [0.9 0.3 0.9])
bar( 3.93 , Sp_saturated_dup ,'FaceColor' , [0.1 0.1 0.1])
set(gca,'xtick',[1.5 3.5])
set(gca,'xticklabel',{'mid-log' , 'saturated'})
xlabel('growth state')
xlim([0.5 4.5])
ylabel('subtelomeric freq. / non-subtel. frequency')
legend({'Collapse' 'Duplication'} , 'location','NorthOutside')
axis tight; 
print('-dpng',[figname 'b.png'],'-r300')
