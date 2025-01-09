function [referenceMultipleSamplesAdjustedCqsMean,referenceMultipleSamplesOutliers] = calculateMultipleSamplesAdjustedCqsMean(referenceTablesCqs,referenceNames,imagesFolder)
% This script performs the function calculateSamplesAdjustedCqsMean, but
% for multiple target genes

% INPUTS
% - referenceTablesCqs: tables Cqs of multiple reference genes
% - referenceNames: names of reference genes
% - imagesFolder: folder to stores images of linear regressions

% OUTPUTS
% - referenceMultipleSamplesAdjustedCqsMean: average adjusted Cqs of multiple
% reference genes
% - referenceMultipleSamplesOutliers: outliers of the Cqs discarded to
% calculate their averages

referenceMultipleSamplesAdjustedCqsMean = cell(length(referenceTablesCqs),1);
referenceMultipleSamplesOutliers = cell(length(referenceTablesCqs),1);

for i = 1:length(referenceTablesCqs)
    [referenceMultipleSamplesAdjustedCqsMean{i,:},referenceMultipleSamplesOutliers{i,:}] = calculateSamplesAdjustedCqsMean(referenceTablesCqs{i,:},referenceNames{i,1},imagesFolder);
end

end