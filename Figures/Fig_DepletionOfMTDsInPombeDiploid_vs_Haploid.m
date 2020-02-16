% enrichment or depletion for MTDs in essential genes
% make essential    in   10kCoverageSequencingData
% % of MTDs in essential genes
d1 = [ 1262 4604 ; 335 1457] ;

[~,p1,~] = fishertest(d1)

d1(:,1)./d1(:,2)
