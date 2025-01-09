function conditionsDeltaCqs = processGeneData(targetFilename,targetName,normalizationFilenames,normalizationNames,referenceCondition,dataFolder,imagesFolder)
% This function returns estimated detla Cqs of a gene for all conditions,
% normalizing using every reference gene included in the analysis

% INPUTS
% - targetFilename: Aria report filename for a specific gene
% - targetName: name of the respective gene
% - normalizationFilenames: Aria report filenames for reference gene
% - normalizationNames: names of the respective reference genes
% - referenceCondition: number of the reference condition
% - dataFolder: folder with the Aria reports
% - imagesFolder: folder to store images

% OUTPUTS
% - conditionLimits: 9%CI and mean relative transcript quantities
% - chosenNormalization: boolean stating reference genes used


%% Define Conditions that will be relativized
relativeConditions = (1:6)';
relativeConditions = relativeConditions(relativeConditions~=referenceCondition);

%% Retrieve Tables Information

targetTableCqs = transformTable(targetFilename,dataFolder);
normalizationTablesCqs = transformTables(normalizationFilenames,dataFolder);

%% Calculate Average of Technical Replicates
[targetSamplesAdjustedCqsMean,targetSamplesOutliers] = calculateSamplesAdjustedCqsMean(targetTableCqs,targetName,imagesFolder);
[normalizationMultipleSamplesAdjustedCqsMean,normalizationMultipleSamplesOutliers] = calculateMultipleSamplesAdjustedCqsMean(normalizationTablesCqs,normalizationNames,imagesFolder);


%% Normalize Biological Replicates by Reference Genes and obtain mean replicated

normalizationN = length(normalizationMultipleSamplesAdjustedCqsMean);
targetConditionsAdjustedCqsMean = cell(6,1);
normalizationConditionsAdjustedCqsMean = cell(6,1);

conditionsDeltaCqs = cell(6,1);
conditionsDeltaCqsRelStd = ones(6,normalizationN);

for i = 1:6
    targetConditionsAdjustedCqsMean{i,1} = repmat(targetSamplesAdjustedCqsMean(i*3-2:i*3,:)',normalizationN,1);
    for j = 1:normalizationN 
        normalizationConditionsAdjustedCqsMean{i,1}(j,:) = normalizationMultipleSamplesAdjustedCqsMean{j,1}(i*3-2:i*3,:)';
    end
    conditionsDeltaCqs{i,1} = targetConditionsAdjustedCqsMean{i,1} - normalizationConditionsAdjustedCqsMean{i,1};
    conditionsDeltaCqsRelStd(i,:) = (abs(std(conditionsDeltaCqs{i,1},0,2,'omitnan')./mean(conditionsDeltaCqs{i,1},2,'omitnan')))';
end

end