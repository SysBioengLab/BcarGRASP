function [samplesAdjustedCqsMean,samplesOutliers] = calculateSamplesAdjustedCqsMean(tableCqs,targetName,imagesFolder)
% This script estimates the adjusted Cqs using the Common Base method for a
% specific gene but for multiple conditions using various options of
% reference genes

% INPUTS
% - tableCqs: Cqs obtained for a PCR plate in matrix format
% - targetName: name of the gene target
% - imagesFolder: folder to store images of linear regressions

% OUTPUTS
% - samlesAdjustedCqsMean: average adjusted Cqs to reference genes by
% separate
% - samplesOutliers: states which samples where discarded to estimate the
% mean values of adjusted Cqs


%% Process data of Cqs

standardsCqs = tableCqs(1:7,10:12);

% standardsOutliers = isoutlier(standardsCqs,2);
% standardsCqs(standardsOutliers) = nan;

standardsCqsVector = reshape(standardsCqs,[],1);
standardsLogQuantity = repmat(linspace(-5,-11,7),1,3);

standardsLM = fitlm(standardsCqsVector,standardsLogQuantity);

samplesCqs = tableCqs(1:6,1:9);
samplesCqsVector = reshape(samplesCqs',[],1);
samplesCqsMatrix = reshape(samplesCqsVector,3,[])';

samplesOutliers = isoutlier(samplesCqsMatrix,2);
samplesCqsMatrix(samplesOutliers) = nan;
samplesCqsVector = reshape(samplesCqsMatrix',[],1);

slope = standardsLM.Coefficients{2,1};
efficiency = 10^-slope;

samplesCqsMean = mean(samplesCqsMatrix,2,'omitnan');
samplesCqsN = sum(~isnan(samplesCqsMatrix),2);
samplesCqsSE = std(samplesCqsMatrix,0,2,'omitnan')./(samplesCqsN).^(1/2);
samplesCqsCI = samplesCqsSE.*tinv(0.975,samplesCqsN-1);
samplesCqsLCI = samplesCqsMean - samplesCqsSE.*tinv(0.975,samplesCqsN-1);
samplesCqsUCI = samplesCqsMean + samplesCqsSE.*tinv(0.975,samplesCqsN-1);
samplesAdjustedCqsMean = log10(efficiency).*samplesCqsMean;

gene = {targetName};
intercept = standardsLM.Coefficients{1,1};
slopeSE = standardsLM.Coefficients{2,2};
freedomDegrees = standardsLM.DFE;
adjustedRsquared = standardsLM.Rsquared.Adjusted;
tScore = tinv(0.975,freedomDegrees);
efficiencyPercentage = ((10^-(slope))-1)*100;
efficiencyPercentageLowerCI = ((10^-(slope+tScore*slopeSE))-1)*100;
efficiencyPercentageUpperCI = ((10^-(slope-tScore*slopeSE))-1)*100;
standardsApplied = standardsLM.Variables(~isnan(standardsLM.Variables.x1),:);
standardsCqMin = min(standardsApplied{:,1});
standardsCqMax = max(standardsApplied{:,1});
standardsDilutionMin = max(-standardsApplied{:,2});
standardsDilutionMax = min(-standardsApplied{:,2});
ntcCq = min(tableCqs(8,10:12),[],'omitnan');

%% Export table

linearRegressionTable = table(gene,slope,intercept,slopeSE,freedomDegrees,adjustedRsquared,efficiencyPercentage,efficiencyPercentageLowerCI,efficiencyPercentageUpperCI, ...
    standardsCqMin,standardsCqMax,standardsDilutionMin,standardsDilutionMax,ntcCq);
writetable(linearRegressionTable,'output\LinearRegressionStatistics.xlsx','Sheet',targetName);


%% Plot linear regressions

% % [samplesQuantityVector,~] = predict(standardsLM,samplesCqsVector,'Alpha',0.05,'Prediction','observation','Simultaneous',false);
% [samplesMeanQuantity,~] = predict(standardsLM,samplesCqsMean,'Alpha',0.05,'Prediction','observation','Simultaneous',false);
% figureLM = figure('visible','on');
% plotAdded(standardsLM)
% hold on
% % plot(samplesCqsVector,samplesQuantityVector,'.k')
% plot(samplesCqsMean,samplesMeanQuantity,'.k')
% errorbar(samplesCqsMean,samplesMeanQuantity,samplesCqsCI,'horizontal','LineStyle','none','Color',[0 0 0])
% xlabel('Detection Cq')
% ylabel('log10 Relative Standard Quantity')
% title(['Linear Regression of Gene: ',targetName])
% legend('Fitted Standards',['Fit: y=',num2str(slope),'x'],'Regression 95% CIs','Samples mean','Samples 95% CIs')
% grid on
% hold off
% outputFileName = fullfile(imagesFolder,'linear_regressions',targetName);
% saveas(figureLM,outputFileName,'jpeg')
% saveas(figureLM,outputFileName,'svg')

end


