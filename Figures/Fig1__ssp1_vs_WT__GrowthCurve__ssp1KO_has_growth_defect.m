PROJDIR = '/Users/lcarey/Nutstore Files/Microhomology shared folder/';
FN = [ PROJDIR 'Figures/Fig1/growth for ssp1 MTD-AGGCA.xlsx' ];

FIGNAME = [ PROJDIR  'Figures/Fig1/Fig1__ssp1_vs_WT__GrowthCurve__ssp1KO_has_growth_defect' ] ;
T = readtable( FN , 'Sheet',1);
data = table2array(T(:,2:end));
hours = [0 3 6 9 12] ; 

FigSize = [5 5 7 5] ; 
%%
fh = figure('units','centimeters','position', FigSize) ;
hold on; 

plot( hours , data(3,:),'.-','DisplayName',T.Var1{3},'Color','b','LineWidth',1)
plot( hours , data(4,:),'.-','DisplayName',T.Var1{4},'Color','b','LineWidth',1)
plot( hours , data(5,:),'.-','DisplayName',T.Var1{5},'Color','b','LineWidth',1)
plot( hours , data(6,:),'.-','DisplayName',T.Var1{6},'Color','b','LineWidth',1)
plot( hours , data(1,:),'.-','DisplayName',T.Var1{1},'Color','r','LineWidth',2)
plot( hours , data(2,:),'.-','DisplayName',T.Var1{2},'Color','r','LineWidth',2)
xlabel('Time (hrs)')
ylabel('Cell Density (OD_{600})')
set(gca,'yscale','lin')
axis tight ;
set(gca,'ytick',0:0.2:2);

l2r = log2(T.x12./T.x0)
[~,p] = ttest2(l2r(1:2),l2r(3:end))

print('-dpng',FIGNAME,'-r600');
close ; 

%% I don't understand why this doesn't give the expected results!!
% Set up fittype and options.
ft = fittype( 'exp1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Robust = 'Bisquare';
exponentials = NaN(size(data,1),1);

% Set up fittype and options.
ft = fittype( {'x', '1'}, 'independent', 'x', 'dependent', 'y', 'coefficients', {'m', 'b'} );


for I = 1:numel(exponentials)
[fitresult, gof] = fit( log(data(I,:)') , hours', ft);
exponentials(I) = fitresult.m  ; 
end
exponentials
[~,p] = ttest2(exponentials(1:2),exponentials(3:end))