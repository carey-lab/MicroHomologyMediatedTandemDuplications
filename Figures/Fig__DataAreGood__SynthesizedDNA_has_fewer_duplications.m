% compare the Duplication frequency for chemically synthesized vs plasmid
% vs gDNA for ssp1
% LBC December 2019
%
% for a figure showing that our data are not a sequencing artifact

%% load data
DATADIR = '~/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/DataFromCluster/' ;
FIGUREDIR = '~/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/FIGURES/' ;

T = readtable( [ DATADIR  'PombeAmpliconSeq_E4_alltsvs.txt'] ,'TreatAsEmpty','-');
T.Properties.VariableNames = {'lib' 'chr' 's1' 'e1' 's2' 'e2' 'ReadDepth' 'DupCounts' 'DupFreq' 'CollapseCounts' 'CollapseFreq' };
T.DupCounts(isnan(T.DupCounts))=0 ;
T.DupFreq(isnan(T.DupFreq))=0 ;
T.HasDup = T.DupCounts>0 ;
T.MHlen = T.e1-T.s1+1 ;

% old version used to takes as input, result of
% Develop/MicroHomologyMediatedIndels/Analysis/Lucas/MH_Event_Frequency_as_a_function_of_MHlen_and_distance('~/Downloads/alltsvs.txt' )

%% make barplot
FIGURE_NAME = [ FIGUREDIR 'Fig__DataAreGood__SynthesizedDNA_has_fewer_duplications___PombeAmpliconSeq_E4_alltsvs_'] ;

idx_chem = strcmp(T.chr , 'ssp1.short.1') | regexpcmp(T.chr , 'ssp1.short...PCR') | strcmp(T.chr , 'ssp1.long.PCR') ;
idx_not_chem = regexpcmp(T.chr , 'ssp1.dup.')  ;
idx_plasmid = regexpcmp(T.chr , 'ssp1.fragment.')  ;

T.class = zeros(height(T),1);
T.class(idx_not_chem) = 1 ; % yeast
T.class(idx_chem) = 2 ;
T.class(idx_plasmid) = 3 ;

T = T(T.class>0,:);
T = T( regexpcmp(T.chr,'^ssp'),:);

G = grpstats(T,'chr',{'median' 'max'},'datavars',{'e2' 'ReadDepth'}); G = sortrows(G,'max_e2','ascend');
T = T( ismember(T.chr , G.chr(G.median_ReadDepth>20000)),:);

%for all MHLs
for MHL = 4:7
    
    T.idx_mhl_short   = T.MHlen== MHL ;
    idx = T.idx_mhl_short ;
    
    b1 = bootstrp( 1000 , @mean , T.DupFreq(idx & T.class==1 & T.DupFreq>0)) ;
    b2 = bootstrp( 1000 , @mean , T.DupFreq(idx & T.class==2 & T.DupFreq>0)) ;
    b3 = bootstrp( 1000 , @mean , T.DupFreq(idx & T.class==3 & T.DupFreq>0)) ;
    
    [~,p2] = ttest2( log10(T.DupFreq(idx & T.class==1 & T.DupFreq>0)) , log10(T.DupFreq(idx & T.class==2 & T.DupFreq>0)));
    [~,p3] = ttest2( log10(T.DupFreq(idx & T.class==1 & T.DupFreq>0)) , log10(T.DupFreq(idx & T.class==3 & T.DupFreq>0)));
    
    titlestr = sprintf('MH Length = %d\nN=%d N=%d N=%d' , MHL , sum(idx & T.class==1) , sum(idx & T.class==2), sum(idx & T.class==3) ) ; 
    
    G = grpstats( T( idx ,:) , 'class' , {'mean'} , 'DataVars' , {'DupFreq' 'HasDup' 'ReadDepth'}) ;
    fh = figure('units','centimeters','position',[5 5 6 7]) ;
    hold on ;
    bh = bar( [ mean(b1) mean(b2) mean(b3) ] ,'FaceColor',[.8 .8 .8]) ;
    errorbar( 1:3 , mean([ b1 b2 b3] ,1) , std([ b1 b2 b3] ,1) ,'.k', 'LineWidth',2);
    if p2<0.001
        text( 2 , max(ylim)*0.9 , '***')
    elseif p2<0.01
        text( 2 , max(ylim)*0.9 , '**')
     elseif p2<0.05
         text( 2 , max(ylim)*0.9 , '*')
    end
    if p3<0.001
        text( 3 , max(ylim)*0.9 , '***')
    elseif p3<0.01
        text( 3 , max(ylim)*0.9 , '**')
    elseif p3<0.05
         text( 3 , max(ylim)*0.9 , '*')
    end
    
    xlim([0.4 3.6])
    set(gca,'xtick',1:3)
    set(gca,'xticklabel',{'gDNA' 'chem' 'plasmid'})
    xlabel('DNA source')
    ylabel('Duplication frequency (10^{-6})')
    title(titlestr)
    print('-dpng' , [ FIGURE_NAME num2str(MHL)] , '-r300');
    close ;
    
    
    % and for the % of MHRs in which we observe a duplication event
    b1 = bootstrp( 1000 , @mean , T.HasDup(idx & T.class==1)) ;
    b2 = bootstrp( 1000 , @mean , T.HasDup(idx & T.class==2)) ;
    b3 = bootstrp( 1000 , @mean , T.HasDup(idx & T.class==3)) ;

    [~,p2] = fishertest( [ sum( T.HasDup(idx & T.class==1) ) sum( ~T.HasDup(idx & T.class==1) ) ; ...
        sum( T.HasDup(idx & T.class==2) ) sum( ~T.HasDup(idx & T.class==2) ) ] ) ; 
    [~,p2] = fishertest( [ sum( T.HasDup(idx & T.class==1) ) sum( ~T.HasDup(idx & T.class==1) ) ; ...
        sum( T.HasDup(idx & T.class==3) ) sum( ~T.HasDup(idx & T.class==3) ) ] )    ; 
    G = grpstats( T( idx ,:) , 'class' , {'mean'} , 'DataVars' , {'DupFreq' 'HasDup' 'ReadDepth'}) ;
    fh = figure('units','centimeters','position',[5 5 6 7]) ;
    hold on ;
    bar( 100*[G.mean_HasDup] ,'FaceColor',[.8 .8 .8])
    errorbar( 1:3 , mean(100*[ b1 b2 b3] ,1) , std(100*[ b1 b2 b3] ,1) ,'.k', 'LineWidth',2);
    if p2<0.001
        text( 2 , max(ylim)*0.9 , '***')
    elseif p2<0.01
        text( 2 , max(ylim)*0.9 , '**')
     elseif p2<0.05
         text( 2 , max(ylim)*0.9 , '*')
    end
    if p3<0.001
        text( 3 , max(ylim)*0.9 , '***')
    elseif p3<0.01
        text( 3 , max(ylim)*0.9 , '**')
    elseif p3<0.05
         text( 3 , max(ylim)*0.9 , '*')
    end
    xlim([0.4 3.6])
    set(gca,'xtick',1:3)
    set(gca,'xticklabel',{'gDNA' 'chem' 'plasmid'})
    xlabel('DNA source')
    ylabel('% of MHRs w/duplication')
    title(titlestr)
    print('-dpng' , [ FIGURE_NAME num2str(MHL) '_2' ] , '-r300');
    close ;
    
end
%% using all MHRs for main figure
idx = T.MHlen >2 ; 

b1 = bootstrp( 1000 , @mean , T.DupFreq(idx & T.class==1 & T.DupFreq>0)) ;
b2 = bootstrp( 1000 , @mean , T.DupFreq(idx & T.class==2 & T.DupFreq>0)) ;
b3 = bootstrp( 1000 , @mean , T.DupFreq(idx & T.class==3 & T.DupFreq>0)) ;

[~,p2] = ttest2( log10(T.DupFreq(idx & T.class==1 & T.DupFreq>0)) , log10(T.DupFreq(idx & T.class==2 & T.DupFreq>0)));
[~,p3] = ttest2( log10(T.DupFreq(idx & T.class==1 & T.DupFreq>0)) , log10(T.DupFreq(idx & T.class==3 & T.DupFreq>0)));

titlestr = sprintf('N=%d N=%d N=%d' , sum(idx & T.class==1) , sum(idx & T.class==2), sum(idx & T.class==3) ) ;

G = grpstats( T( idx ,:) , 'class' , {'mean'} , 'DataVars' , {'DupFreq' 'HasDup' 'ReadDepth'}) ;
fh = figure('units','centimeters','position',[5 5 6 7]) ;
hold on ;
bh = bar( [ mean(b1) mean(b2) mean(b3) ] ,'FaceColor',[.8 .8 .8]) ;
errorbar( 1:3 , mean([ b1 b2 b3] ,1) , std([ b1 b2 b3] ,1) ,'.k', 'LineWidth',2);
    if p2<0.001
        text( 2 , max(ylim)*0.9 , '***')
    elseif p2<0.01
        text( 2 , max(ylim)*0.9 , '**')
     elseif p2<0.05
         text( 2 , max(ylim)*0.9 , '*')
    end
    if p3<0.001
        text( 3 , max(ylim)*0.9 , '***')
    elseif p3<0.01
        text( 3 , max(ylim)*0.9 , '**')
    elseif p3<0.05
         text( 3 , max(ylim)*0.9 , '*')
    end

xlim([0.4 3.6])
set(gca,'xtick',1:3)
set(gca,'xticklabel',{'gDNA' 'chem' 'plasmid'})
xlabel('DNA source')
ylabel('Duplication frequency (10^{-6})')
title(titlestr)
print('-dpng' , [ FIGURE_NAME '__MainFig'] , '-r300');
close ;
%%
TotalDupEventCounts_Native = nansum(T.DupCounts(strcmp(T.chr,'ssp1.dup.2')))   % 12,555
TotalDupEventCounts_SYN = nansum(T.DupCounts(strcmp(T.chr,'ssp1.short.1.PCR')))   % 2

% $ samtools view -c lib_14.bam ssp1.dup.2    ->  43420306
% $ samtools view -c lib_1.bam ssp1.short.1   ->     73303
%12,555 / 43420306 = 0.000289150

% p = 1.2478e-09

% $ samtools view -c lib_2.bam ssp1.short.1.PCR  ->  2405009

% X = [#_reads_Native #_dup_events_native
%     #_reads_synth  #_dup_events_synth ]

X = [ 43420306 12555 ;...
    2405009 2] ;

[~,p,festats] = fishertest( X )

% ../sequencing_raw_data/E4/lib_2.fa:1:>ssp1.short.1.PCR

p =  4.5889e-289

%
%
% %% SSP1 FLASH results
% %$ cut -f 2 lib_15.readsAndPairs.tab |grep GAGAACGTGAAGGAAGTTCGTTAACTCACTCATGGACTTTTCAACCTGGTAAGCATAACCAGCGTCTTTATTCTGATAATTTTCAAGAGGCTCAGCGCCAGTGGAAGCGCCTGCAAGAATGGGGCGAGGTGAAGGAAACAAAA |wc -l
% % 237712 events
% fid=fopen('~/Downloads/ssp1_89_dup');
% C = textscan(fid,'%s');
% fclose(fid);
% C=C{:};
% C = C( ~regexpcmp( C , 'AGAGGTAGACGATAGCGAAGTTCCTCCTTCTGTCTTTCCTGAATATCCCGTCCACAAGGCCATCCAGAGAGGTAGACGATAGCGAAGTTCCTCCTTCTGTCTTTCCTGAATATCCCGTCCACAAGGCCATCCAGAAAACGTCCGATTCATTTCGTAAACGGAACTACAGCGCGGGAGATTATGTA'));
%
%
% m=multialign(C) ; showalignment(m)