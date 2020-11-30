pn = {'MHlen' 'interMH'  'GCcon'  'interGCcon' 'dis_to_telomere_log10' 'POLII' 'DAVID'  'bintclosestMHR' 'replication_time' };
pn = [ pn T.Properties.VariableNames(regexpcmp(T.Properties.VariableNames,'^K')) ] ; 
pn = {'MHlen' 'InGeneALL' 'InGeneANY' 'GCcon' }
%pn = {'ntclosestMHR' 'logntclosestMHR' 'bintclosestMHR'   } ;

KFoldCV = 2 ;


sc = NaN(numel(pn),1);
R = table;
R.nP = NaN(1e4,1);
R.predNams  = cell(1e4,1);
R.predNams2 = cell(1e4,1);
R.vA = NaN(1e4,1);
c = 0 ;
for K = [1 2]
    pn2 = nchoosek( pn , K);
    for I = 1:nrows(pn2)
        predNams = pn2(I,:);
        c = c+1 ;
        trainingData = Q ;
        inputTable = trainingData;
        %predictorNames = {'MHlen'  'GCcon'   'dis_to_telomere_log10' 'POLII' 'interMH' 'interGCcon' 'bintclosestMHR'  }; %   'dis_to_telomere_log10'  'replication_time'  };
        predictorNames = predNams ;
        predictors = inputTable(:, predictorNames);
        response = inputTable.dup;
        isCategoricalPredictor = [false, false, false, false];
        
        % Train a classifier
        % This code specifies all the classifier options and trains the classifier.
        template = templateTree(...
            'MaxNumSplits', 20);
        classificationEnsemble = fitcensemble(...
            predictors, ...
            response, ...
            'Method', 'AdaBoostM1', ...
            'NumLearningCycles', 30, ...
            'Learners', template, ...
            'LearnRate', 0.1, ...
            'ClassNames', [false; true]);
        
        % Create the result struct with predict function
        predictorExtractionFcn = @(t) t(:, predictorNames);
        ensemblePredictFcn = @(x) predict(classificationEnsemble, x);
        trainedClassifier.predictFcn = @(x) ensemblePredictFcn(predictorExtractionFcn(x));
        
        % Add additional fields to the result struct
        trainedClassifier.RequiredVariables = predictorNames ;
        trainedClassifier.ClassificationEnsemble = classificationEnsemble;
        trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2020a.';
        trainedClassifier.HowToPredict = sprintf('To make predictions on a new table, T, use: \n  yfit = c.predictFcn(T) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nThe table, T, must contain the variables returned by: \n  c.RequiredVariables \nVariable formats (e.g. matrix/vector, datatype) must match the original training data. \nAdditional variables are ignored. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');
        
        % Extract predictors and response
        % This code processes the data into the right shape for training the
        % model.
        inputTable = trainingData;
        predictors = inputTable(:, predictorNames);
        response = inputTable.dup;
        isCategoricalPredictor = [false, false, false, false];
        
        % Perform cross-validation
        partitionedModel = crossval(trainedClassifier.ClassificationEnsemble, 'KFold', KFoldCV );
        
        % Compute validation predictions
        [validationPredictions, validationScores] = kfoldPredict(partitionedModel);
        
        % Compute validation accuracy
        validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');
        
        
        
        R.predNams{c} = predNams ;
        R.predNams2{c} = strjoin(predNams,',') ;
        R.vA(c) = validationAccuracy ;
        R.nP(c) = K ;
        
        fprintf('%d\t%d\t%0.02f\t%s\n' , c , K , validationAccuracy, strjoin(predNams,',') ) ;
        
    end
end
R = R( ~isnan(R.nP) , : );
R = sortrows(R,'vA','descend');
%% %% %% %% GRAPH best features %% %% %%

Q = R( max(R.vA) - R.vA < 0.01  , :);
vn = horzcat(Q.predNams{:}) ; 
[a,b] = count_unique(vn);
b = 100* (b ./ height(Q)) ;
Q = table();
Q.vn = a ; 
Q.pct = b ; 
Q = sortrows(Q,'pct','descend');

figure ; 
bar( Q.pct )
set(gca,'xtick',1:height(Q))
set(gca,'xticklabels',Q.vn)
xtickangle(45);
ylabel('% of models with this feature')