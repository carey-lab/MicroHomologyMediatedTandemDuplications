% Figure S4. Characterization of the logistic regression model for predicting MTDs from cis MHP features. 
%

FIGDIR = '~/Nutstore Files/Microhomology shared folder/Figures/Supplementary Figures/';
FIGNAME = [FIGDIR 'model feature importance' ] ;

T = readtable('~/Develop/MicroHomologyMediatedIndels/Data/diploid_model_parameters.txt' ,'Delimiter','\t','ReadVariableNames',true);
fh = figure('units','centimeters','position',[5 5 7.5 6]) ;
T = sortrows(T,'coeff','descend');
clrs = parula(height(T)+3);
gh = gscatter(  abs(log10(T.p)) , T.coeff , T.feature , clrs , 'vps^o' , 18,'filled');
for I = 1:numel(gh)
    gh(I).MarkerFaceColor = clrs(I,:);
    gh(I).DisplayName = T.feature{I};
end
xlim([6 16.5])

ylabel({'Feature importance' '(model coefficient)'})
xlabel('p-value   (abs(log_{10}))')
legend('location','nw')
print('-dpng' , FIGNAME , '-r300');
%close ; 