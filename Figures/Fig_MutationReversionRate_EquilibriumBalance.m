% simulate mutation and reversion at varying reversion rates
%   show that MTDs will, by nature, always be subclonal

% High reversion rates do not strongly affect measured mutant frequencies during short-term outgrowth from a single cell 
% into a population of ~108 cells (these experiments) 
% but does reduce the mutant genotype frequency within the population during long-term evolution. 

%% https://www.apsnet.org/edcenter/disimpactmngmnt/topc/PopGenetics/Pages/Mutation.aspx
% Assume that u = 1 x 10-5 and v = 1 x 10-6 per generation. 
% These are typical forward and backward mutation rates. 
% Let the frequency of A1 (call it p) = 0.99, and the frequency of A2 (q) = 0.01.
% What is the new frequency of A2 after one generation of mutation?

FIGDIR  = '~/Nutstore Files/Microhomology shared folder/Figures/Fig4/' ;
addpath_recurse('~/Develop/MicroHomologyMediatedIndels')

% = 0.01 + [(10-5)(0.99) - (10-6)(0.01)] = 0.01 + 9.89 x 10-6 approximately = 0.01001.
LW = 3 ; 
WTfreq = 1 ; 
mutation_rate = 10^-7 ;
number_of_generations = 30 ; 

fh = figure('units','centimeters','position',[5 5 12 6]) ;hold on ;
t = tiledlayout(1,2);
nexttile ; hold on ;
plot( simulate_mutation_reversion(number_of_generations,mutation_rate,1e-7)  , '-','Color',[.5 .5 .5] ,'DisplayName','10^{-7}','LineWidth',LW)
plot( simulate_mutation_reversion(number_of_generations,mutation_rate,0.01)  , '-','DisplayName','0.01','LineWidth',LW)
plot( simulate_mutation_reversion(number_of_generations,mutation_rate,0.001)  , '-','DisplayName','0.001','LineWidth',LW)
xlabel('Number of generations')
ylabel('Mutant frequency (10^{-6})')
%title( sprintf('             mutation rate = 10^{%d}' , log10(mutation_rate) ) );
title('       short-term')
grid off ;
lh = legend('location','best','box','off');
%text( 1 , max(ylim)*0.97 , 'reversion rate')
set(gca,'ytick',(0:3).*1e-6 ) ;  set(gca,'yticklabel',0:3) ; ylim([0 max(get(gca,'ytick'))])
text(1 , max(ylim)*0.96 , '\underline{reversion rate}', 'FontSize', 12, 'Interpreter', 'latex','FontName','Helvetica')

nexttile ; hold on ;
number_of_generations = 3e3 ; mutation_rate = 10^-7 ;
plot( simulate_mutation_reversion(number_of_generations,mutation_rate,1e-7)  , '-','Color',[.5 .5 .5] ,'DisplayName','10^{-7}','LineWidth',LW)
plot( simulate_mutation_reversion(number_of_generations,mutation_rate,0.01)  , '-','DisplayName','0.01','LineWidth',LW)
plot( simulate_mutation_reversion(number_of_generations,mutation_rate,0.001)  , '-','DisplayName','0.001','LineWidth',LW)
xlabel('Number of generations')
ylabel('Mutant frequency (10^{-4})')
%title( sprintf('             mutation rate = 10^{%d}' , log10(mutation_rate) ) );
grid off ;
legend('location','best','box','off')
title('       long-term')
set(gca,'ytick',(0:3).*1e-4 ) ;  set(gca,'yticklabel',0:3) ; ylim([0 max(get(gca,'ytick'))])
text(100 , max(ylim)*0.96 , '\underline{reversion rate}', 'FontSize', 12, 'Interpreter', 'latex','FontName','Helvetica')
print('-dpng',[FIGDIR 'MTDs_are_always_subclonal__simulation_mutation_reversion__MutationFrequencies'] ,'-r300') ;
close all ; 
%% %% old steady-state equilibrium simulation
% 
% % reading
% % https://www.apsnet.org/edcenter/disimpactmngmnt/topc/PopGenetics/Pages/Mutation.aspx
% % http://www.nyu.edu/projects/fitch/courses/evolution/html/genetic_drift.html#Exercises
% 
% % A mutation (A*) arises (repeatedly) at a slow rate (u) from an allele (A) that already exists in the population 
% %     at a frequency f(A) = p.  
% % But there is also a rate of reversion (v) back to A.  
% % The change in q is up-vq = u(1-q) - vq.  
% % At equilibrium, this change is 0.  
% % What is the predicted equilibrium frequency of A* (q) if u = 10-5 and v = 10-7?
% 
% 
% % change in q is up-vq = u(1-q) - vq.   ; at equilibrium, = 0 
% % up-vq  = 0
% % u(1-q) - vq = 0
% 
% 
% %u = 10^-5 ;
% %v = 10^-7 ;
% syms u v  p q ;
% eqn1 =  u*(1-q) - v*q == 0 ;
% eqn2 =  u*p-v*q  == 0 ;
% eqn2 =  u*(1-q)-v*q  == 0 ;
% 
% sol1 = solve( eqn1  , q  , 'ReturnConditions',true) ;
% sol2 = solve( eqn2  , q  , 'ReturnConditions',true) ;
% 
% syms f1(u,v) f2(u,v) ;
% f1(u,v) = sol1.q ;
% f2(u,v) = sol2.q ;
% 
% 
% xV = -9:0.1:-2;
% resultsmat = NaN(numel(xV));
% for I = 1:numel(xV)
%     resultsmat(I,:) = double( f1( 10^xV(I) , 10.^xV) ) ; 
% end
% 
% %%
% fh = figure('units','centimeters','position',[5 5 8 7]) ;
% hold on ; 
% %imagesc(resultsmat);
% colormap(vertcat(  [.7 .7 .7], hot(100)) )
% X = resultsmat ; 
% % X =  tril(resultsmat')' ; % remove all values where mutation < reversion 
% % X(X==0)=NaN ;
% for I = 1:size(X,1)
%     for J = 1:size(X,1)
%         if I>J
%             X(I,J) = NaN ; 
%         end
%     end
% end
% 
% X = log10(X);
% imagesc( X  ); 
% 
% ch = colorbar ;
% ch.TickLabels =  cellfun(@(X) sprintf('10^{%s}' , X)  , ch.TickLabels  , 'UniformOutput'  , false) ; 
% xlabel('Reversion rate')
% ylabel('MTD rate');
% ch.Label.String = 'MTD frequency at equilibruium' ;
% tickPositions = linspace( min(xticks()) , max(xticks) , 8) ; 
% tickLabels = linspace( min(xV()) , max(xV) , 8) ;
% tickLabels = arrayfun(@(X) sprintf('10^{%d}' , X)  , tickLabels  , 'UniformOutput'  , false) ; 
% set(gca,'xtick',tickPositions);
% set(gca,'ytick',tickPositions);
% set(gca,'xticklabels',tickLabels);
% set(gca,'yticklabels',tickLabels);
% axis tight ;
% 
% 
% [~,idx] = min(abs(10.^xV - 10^-3) ) ; 
% xTickLoc_Closest_to_ReversionRate = idx ; 
% 
% [~,idx] = min(abs(10.^xV - 10^-8) ) ; 
% yTickLoc_Closest_to_ReversionRate = idx ; 
% 
% plot( xTickLoc_Closest_to_ReversionRate , yTickLoc_Closest_to_ReversionRate , 'ok','MarkerFaceColor',[.7 .7 .7]);
% 
% %%
% 
% fh = figure('units','centimeters','position',[5 5 8 7]) ;
% hold on ; 
% imagesc( tril(resultsmat')');
% ch = colorbar ;
% xlabel('Reversion rate')
% ylabel('MTD rate');
% axis tight;
% %a1 = f1(10^-5,10.^xV) ;
% 
% % figure; 
% % hold on ;
% % plot(xV , double(a1) ,'ok') ;
% 
% % xlabel('v (rate of reversion)') ;
% % ylabel('q (A* freq, mutation freq)') ;