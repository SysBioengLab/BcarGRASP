function [conditionLimits,chosenNormalization] = processConditionData(genesDeltaCqs,referenceCondition,targetCondition,conditionNames,conditionNamesFull,targetNames,normalizationNames,imagesFolder)
% This function estimates the mean and 95% CI of each gene in the
% normalized conditions compared to a reference condition. It also saves
% which set of reference genes where used for the normalization

% INPUTS
% - genesDeltaCqs: DeltaCqs of target genes for every reference gene within
% a same condition
% - referenceCondition: number of the reference condition
% - targetCondition: number of the target condition to be processed
% - conditionNames: names of the conditions
% - conditionNamesFull: full names of the conditions
% - targetNames: names of the target genes
% - normalizationNames: names of the reference genes
% - imagesFolder: folder to store images

% OUTPUTS
% - conditionLimits: 9%CI and mean relative transcript quantities
% - chosenNormalization: boolean stating reference genes used

%% Transform into cells with conditions instead of genes
genesN = length(targetNames);
allDeltaCqs = cell(genesN,6);
for i = 1:6
    for j = 1:genesN
        allDeltaCqs{j,i} = genesDeltaCqs{j,1}{i,1};
    end
end


referenceDeltaCqs = allDeltaCqs(:,referenceCondition);

conditionDeltaCqs = allDeltaCqs(:,targetCondition);

%%  Create combinatorial vectors

normalizationN = length(normalizationNames);

logicalVectors = cell(1,normalizationN);
for i = 1:normalizationN
    logicalVectors{1,i} = [0,1];
end
combinatorialVectors = combvec(logicalVectors{:,:});
combinatorialVectors = boolean(combinatorialVectors(:,2:end));
combinatorialN = size(combinatorialVectors,2);

%% Attempt all possible combinations for a condition

cis = cell(genesN,1);
confidenceIntervals = cell(genesN,1);
stats = cell(genesN,1);

conditionDeltaCqsOutliersAll = cell(genesN,1);
referenceDeltaCqsOutliersAll = cell(genesN,1);
conditionDeltaDeltaCqsAll = cell(genesN,1);
relativeVariability = cell(genesN,1);

conditionDeltaCqsOutliers = boolean(zeros(genesN,3));
referenceDeltaCqsOutliers = boolean(zeros(genesN,3));
conditionDeltaDeltaCqs = zeros(genesN,1);
conditionDeltaDeltaCqsLowerCI = zeros(genesN,1);
conditionDeltaDeltaCqsUpperCI = zeros(genesN,1);

for i = 1:genesN
    for j = 1:combinatorialN
        conditionDeltaCqsPrefiltered = mean(conditionDeltaCqs{i,1}(combinatorialVectors(:,j),:),1,'omitnan');
        conditionDeltaCqsOutliersAll{i,1}(j,:) = isoutlier(conditionDeltaCqsPrefiltered);
        conditionDeltaCqsSelected = conditionDeltaCqsPrefiltered;
        conditionDeltaCqsSelected(conditionDeltaCqsOutliersAll{i,1}(j,:)) = nan;

        referenceDeltaCqsPrefiltered = mean(referenceDeltaCqs{i,1}(combinatorialVectors(:,j),:),1,'omitnan');
        referenceDeltaCqsOutliersAll{i,1}(j,:) = isoutlier(referenceDeltaCqsPrefiltered);
        referenceDeltaCqsSelected = referenceDeltaCqsPrefiltered;
        referenceDeltaCqsSelected(referenceDeltaCqsOutliersAll{i,1}(j,:)) = nan;

        conditionDeltaDeltaCqsAll{i,1}(j,1) = mean(conditionDeltaCqsSelected,2,'omitnan') - mean(referenceDeltaCqsSelected,2,'omitnan');
        [~,~,cis{i,1}{j,1},stats{i,1}{j,1}] = ttest2(conditionDeltaCqsSelected',referenceDeltaCqsSelected,'Vartype','unequal','Alpha',0.05);
        confidenceIntervals{i,1}(j,1) = cis{i,1}{j,1}(1,:);
        confidenceIntervals{i,1}(j,2) = cis{i,1}{j,1}(2,:);
    end
    relativeVariability{i,1} = abs((10.^-confidenceIntervals{i,1}(:,1) - 10.^-confidenceIntervals{i,1}(:,2))./(10.^-conditionDeltaDeltaCqsAll{i,1}));
end

[minRelativeVariability,bestCombination] = min(max([relativeVariability{:,:}],[],2));
chosenNormalization = combinatorialVectors(:,bestCombination)';

for i = 1:genesN
    conditionDeltaCqsOutliers(i,:) = conditionDeltaCqsOutliersAll{i,1}(bestCombination,:);
    referenceDeltaCqsOutliers(i,:) = referenceDeltaCqsOutliersAll{i,1}(bestCombination,:);
    conditionDeltaDeltaCqs(i,1) = conditionDeltaDeltaCqsAll{i,1}(bestCombination,1);
    conditionDeltaDeltaCqsLowerCI(i,1) = confidenceIntervals{i,1}(bestCombination,1);
    conditionDeltaDeltaCqsUpperCI(i,1) = confidenceIntervals{i,1}(bestCombination,2);
end

%% Transform relations in linear scale

conditionRelativeMean = 10.^-conditionDeltaDeltaCqs;
conditionRelativeUpperCI = 10.^-conditionDeltaDeltaCqsLowerCI;
conditionRelativeLowerCI = 10.^-conditionDeltaDeltaCqsUpperCI;

conditionLowerError = conditionRelativeMean - conditionRelativeLowerCI;
conditionUpperError = conditionRelativeUpperCI - conditionRelativeMean;

%% Plot relative Concentrations
figureRelativeConcentrations = figure();
categories = categorical(targetNames);
categories = reordercats(categories,targetNames);
bar(categories,conditionRelativeMean,...
    'FaceColor','w','EdgeColor','k','LineWidth',0.5)
hold on
errorbar(categories,conditionRelativeMean',conditionLowerError',conditionUpperError',...
    'LineWidth',0.5,'Color','k','LineStyle','none')
title(['Relative Expression ',conditionNamesFull{targetCondition},' relative to condition ',conditionNamesFull{referenceCondition}])
ylabel('Relative Expression');
xlabel('Genes')
hold off
outputFileName = fullfile(imagesFolder,'relative_expression',[conditionNames{referenceCondition},'_',conditionNames{targetCondition}]);
saveas(figureRelativeConcentrations,outputFileName,'jpeg')

%% Save outliers data for the condition
% referenceOutliersTable = table(targetNames,referenceDeltaCqsOutliers);
% referenceOutliersTable = splitvars(referenceOutliersTable,'referenceDeltaCqsOutliers','newVariableNames',["Replicate 1","Replicate 2","Replicate 3"]);
% referenceOutliersTable = renamevars(referenceOutliersTable,"targetNames","Gene");
% writetable(referenceOutliersTable,'output\ReferenceOutliers.xlsx','Sheet',[conditionNames{referenceCondition},'_',conditionNames{targetCondition}]);
% 
% conditionOutliersTable = table(targetNames,conditionDeltaCqsOutliers);
% conditionOutliersTable = splitvars(conditionOutliersTable,'conditionDeltaCqsOutliers','newVariableNames',["Replicate 1","Replicate 2","Replicate 3"]);
% conditionOutliersTable = renamevars(conditionOutliersTable,"targetNames","Gene");
% writetable(conditionOutliersTable,'output\ConditionsOutliers.xlsx','Sheet',[conditionNames{referenceCondition},'_',conditionNames{targetCondition}]);

%% Save Relative Variability Table
% 
% relativeVariabilities = [relativeVariability{:,:}];
% combinatorials = (1:combinatorialN)';
% relativeVariabilitiesTable = table(combinatorials,relativeVariabilities);
% relativeVariabilitiesTable = splitvars(relativeVariabilitiesTable,'relativeVariabilities','NewVariableNames',targetNames(:,:)');
% relativeVariabilitiesTable = renamevars(relativeVariabilitiesTable,"combinatorials","Combinatorial ID");
% writetable(relativeVariabilitiesTable,'output\RelativeVariabilites.xlsx','Sheet',[conditionNames{referenceCondition},'_',conditionNames{targetCondition}]);

%% Save Combinatorial Table

% combinatorialCombinations = combinatorialVectors';
% combinatorialsTable = table(combinatorials,combinatorialCombinations);
% combinatorialsTable = splitvars(combinatorialsTable,'combinatorialCombinations','NewVariableNames',normalizationNames(:,:)');
% combinatorialsTable = renamevars(combinatorialsTable,"combinatorials","Combinatorial ID");
% writetable(combinatorialsTable,'output\RelativeVariabilites.xlsx','Sheet',"Normalizations Combinations")

%% Define Condition Limits

conditionLimits = [conditionRelativeLowerCI,conditionRelativeMean,conditionRelativeUpperCI];

end
