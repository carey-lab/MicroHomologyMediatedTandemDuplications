% Fig_NaturalIsolates_TandemDuplications_WithWithout_MH
%   what % of >10bp insertions in wild S. pombe isolates have MH? 
FIGDIR  = '~/Nutstore Files/Microhomology shared folder/Figures/Fig4/' ;
FN = [ FIGDIR '57 nature isolates_sequencing summary.xlsx' ] ;

T = readtable( FN , 'Sheet' , 2);
T = T( : , 1:13) ; 
T.Properties.VariableNames{6} = 'REF_N' ; 
T.Properties.VariableNames{7} = 'ALT_N' ; 
T = T(2:end,:);
T.REF_N = str2double(T.REF_N) ; 
T.ALT_N = str2double(T.ALT_N) ; 
T.HasMH = ~cellfun(@isempty,T.MH);

T.REF_len = cellfun(@length,T.REF);
T.ALT_len = cellfun(@length,T.ALT);
T.IsInsertion = T.ALT_len > T.REF_len ; 

T.Location = categorical(T.Location);

T = T( T.IsInsertion , :); 

%%
fh = figure('units','centimeters','position',[5 5  12 5]) ;
FS = 12 ; 
clrs = get(gca,'ColorOrder');

%clr1 = [118 214 255] ./ 255 ; 
%clr2 = [.8 .8 .8] ; 
clr1 = ( clrs(4,:) ).^0.05;
clr2 = ( clrs(5,:) ).^0.05;

idx_1x_2x = T.REF_N==1 & T.ALT_N==2 ; 
idx_1x_Xx = T.REF_N==1 & T.ALT_N>2 ; 
idx_Xx_Xx = T.REF_N>1 & T.ALT_N>T.REF_N ; 


data0X = [ sum( T.REF_N==0 & T.HasMH ) ; sum( T.REF_N==0 & ~T.HasMH ) ] ;
data12 = [ sum( idx_1x_2x & T.HasMH ) ; sum( idx_1x_2x & ~T.HasMH ) ] ;
data1X = [ sum( idx_1x_Xx & T.HasMH ) ; sum( idx_1x_Xx & ~T.HasMH ) ] ;
dataXX = [ sum( idx_Xx_Xx & T.HasMH ) ; sum( idx_Xx_Xx & ~T.HasMH ) ] ;

t = tiledlayout(1,4,'TileSpacing','Compact','Padding','Compact');

nexttile ; 
bh = bar([ data0X [0;0]]' , 'stacked' ,'FaceColor','flat') ;
axis tight; xlim([0.5 1.5])
set(gca,'xticklabel','simple insertion')
box off
bh(1).FaceColor = clr1 ; bh(2).FaceColor = clr2 ;
text( 0.8 , max(ylim)*0.15 , {'with' 'MH'},'FontSize',FS)
text( 0.8 , max(ylim)*0.75 , {' no' 'MH'},'FontSize',FS)


nexttile ; 
bh = bar([ data12 [0;0]]' , 'stacked' ) ;
axis tight; xlim([0.5 1.5])
set(gca,'xticklabel','1x ? 2x')
box off
bh(1).FaceColor = clr1.^4 ; bh(2).FaceColor = clr2.^4 ;
text( 0.8 , max(ylim)*0.15 , {'with' 'MH'},'FontSize',FS)
text( 0.8 , max(ylim)*0.75 , {' no' 'MH'},'FontSize',FS)

nexttile ; 
bh = bar([ data1X [0;0]]' , 'stacked' ) ;
axis tight; xlim([0.5 1.5])
set(gca,'xticklabel','1x ? >2x')
box off
bh(1).FaceColor = (clr1) ; bh(2).FaceColor = (clr2) ;
text( 0.8 , max(ylim)*0.15 , {'with' 'MH'},'FontSize',FS)
text( 0.8 , max(ylim)*0.75 , {' no' 'MH'},'FontSize',FS)

nexttile ; 
bh = bar([ dataXX [0;0]]' , 'stacked' ) ;
axis tight; xlim([0.5 1.5])
set(gca,'xticklabel','>1x ? >Nx')
box off
bh(1).FaceColor = (clr1) ; bh(2).FaceColor = (clr2) ;
text( 0.8 , max(ylim)*0.15 , {'with' 'MH'},'FontSize',FS)
text( 0.8 , max(ylim)*0.75 , {' no' 'MH'},'FontSize',FS)


ylabel(t,'Number of insertions')

print('-dpng' , [ FIGDIR 'NaturalIsolates_TandemDuplications_WithWithout_MH_barplots' ] , '-r300');
close ; 
