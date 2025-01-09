%% Process transcript data
% This script process the Cq data obtained from the Aria qPCR reports. The
% files where previously generated using the software of the equipment to
% generate Cq values at a defined fluorescence threshold.

% It requires all the report files for genes in the network and the
% reference genes for normalization. The code only processed one reference
% condition at a time, but does not overwrite the results of other
% reference conditions. To change the reference condition to be used, the
% variable 'referenceCondition' be modified to match the inputs

% The main output of the code are the mean and confidence intervals of the
% relativized non-reference conditions. They are stored in
% 'output/ConditionsRelativeLimits'. Other outputs are also stored in the
% folder 'output'. 'AppliedNormalizations' includes which set of genes
% were used as normalizers for each non-reference condition. 'Linear 
% Regression Statistics' summarizes the results for the linear regressions
% performed for qPCR data. 'images/GeneRelativeExpression' stores the
% boxplots of relative expression for one gene. 'images/RelativeExpression'
% includes the boxplots of relative expression for a
% non-reference condition. 'Images/Linear Regressions' has the plots of the linear
% regressions obtained for the qPCR fluorescence data.

%% Clean variables and close windows
clear, close all
rng('default');                 % for reproducibility
format longE
addpath('additional_code')


%% Define Inputs
conditionNames = {'2A';'2B';'3A';'3B';'4A';'4B'};                                       % names of the conditions (2,3,4 stand for the strain; A and B for low and high growth rate respectively)
strains = {'2';'2';'3';'3';'4';'4'};                                                    % strain names
growthRates = {'0.101';'0.254';'0.101';'0.254';'0.101';'0.254'};                        % growth rates
conditionNamesFull = {'\beta-car2 \mu=0.101 h^{-1}';'\beta-car2 \mu=0.254 h^{-1}';...
    '\beta-car3 \mu=0.101 h^{-1}';'\beta-car3 \mu=0.254 h^{-1}';
    '\beta-car4 \mu=0.101 h^{-1}';'\beta-car4 \mu=0.254 h^{-1}'};                       % condition full names

dataFolder = '..\..\rt-qpcr_data\aria_reports';                                         % folder with the Aria Reports
imagesFolder = 'output\images';                                                         % folder to store the images

referenceCondition = 5;                                                                 % reference condition on which to perform relativization

targetFilenames = {'bcar-s-erg10-24-03-13-report.xlsx';
    'bcar-r-erg13-24-03-12-report.xlsx';
    'bcar-g-hmg1-24-03-08-report.xlsx';
    'bcar-q-hmg2-24-03-12-report.xlsx';
    'bcar-p-erg12-24-03-11-report.xlsx';
    'bcar-o-erg8-24-03-11-report.xlsx';
    'bcar-f-mvd1-24-03-07-report.xlsx';
    'bcar-n-idi1-24-03-11-report.xlsx';
    'bcar-l-erg20-24-03-09-report.xlsx';
    'bcar-m-bts1-24-03-10-report.xlsx';
    'bcar-d-crte-24-03-06-report.xlsx';
    'bcar-k-erg9-24-03-09-report.xlsx';
    'bcar-e-crti-24-03-07-report.xlsx';
    'bcar-c-crtyb-24-03-06-report.xlsx'};                                               % Names of the files created with Aria Reports for the target genes
targetNames = {'ERG10';'ERG13';'HMG1';'HMG2';'ERG12';'ERG8';'MVD1';...
    'IDI1';'ERG20';'BTS1';'CrtE';'ERG9';'CrtI';'CrtYB'};                                % Name of target genes (in the pathway)

normalizationFilenames = {'bcar-a-taf10-24-03-13-report.xlsx';
    'bcar-b-ubc6-24-03-14-report.xlsx';
    'bcar-i-alg9-24-03-08-report.xlsx';
    'bcar-h-act1-24-03-13-report.xlsx'};                                                % Name of the files created with Aria Reports for the reference genes
normalizationNames = {'TAF10';
    'UBC6';
    'ALG9';
    'ACT1'};                                                                            % Name of reference genes (normalizers)

%% Process Gene Information
genesN = length(targetFilenames);
genesDeltaCqs = cell(genesN,1);
for i = 1:genesN
    genesDeltaCqs(i,1) = {processGeneData(targetFilenames{i,1},targetNames{i,1},normalizationFilenames,normalizationNames,referenceCondition,dataFolder,imagesFolder)};
end


%% Process Condition Information and save them
relativeConditions = (1:6)';
relativeConditions = relativeConditions(relativeConditions ~= referenceCondition,:);
conditionsLimits = cell(6,1);
chosenNormalizations = zeros(6,length(normalizationFilenames));

for i = relativeConditions'
    targetCondition = i;
    [conditionsLimits{i,1},chosenNormalizations(i,:)] = processConditionData(genesDeltaCqs,referenceCondition,targetCondition,conditionNames,conditionNamesFull,targetNames,normalizationNames,imagesFolder);
end

%% Save selection of normalizations

chosenNormalizationsTable = table(strains,growthRates,chosenNormalizations);
chosenNormalizationsTable = chosenNormalizationsTable(relativeConditions,:);
chosenNormalizationsTable = splitvars(chosenNormalizationsTable,'chosenNormalizations','NewVariableNames',normalizationNames');
chosenNormalizationsTable = renamevars(chosenNormalizationsTable,["strains","growthRates"],["Strain","Growth Rate (1/h)"]);
writetable(chosenNormalizationsTable,'output\AppliedNormalizations.xlsx','Sheet',conditionNames{referenceCondition});

%% Save Limits of conditions

conditionsLimits{referenceCondition,1} = ones(genesN,3);
conditionsLimitsTable = table(targetNames,[conditionsLimits{:,:}]);
writetable(conditionsLimitsTable,'output\ConditionsRelativeLimits.xlsx','Sheet',conditionNames{referenceCondition});

%% Plot Relative Gene Expression

conditionsRelativeLimits = conditionsLimits(relativeConditions,:);

genesRelativeMean = cell(genesN,1);
genesRelativeLowerCI = cell(genesN,1);
genesRelativeUpperCI = cell(genesN,1);

for i = 1:genesN
    for j = 1:5
        genesRelativeLowerCI{i,1}(j,1) = conditionsRelativeLimits{j,1}(i,1);
        genesRelativeMean{i,1}(j,1) = conditionsRelativeLimits{j,1}(i,2);
        genesRelativeUpperCI{i,1}(j,1) = conditionsRelativeLimits{j,1}(i,3);
    end
end

for i = 1:genesN
    relativeMean = genesRelativeMean{i,1};
    lowerError = genesRelativeMean{i,1} - genesRelativeLowerCI{i,1};
    upperError = genesRelativeUpperCI{i,1} - genesRelativeMean{i,1};

    figureRelativeConcentrations = figure();
    bar(categorical(conditionNamesFull(relativeConditions)),relativeMean,...
        'FaceColor','w','EdgeColor','k','LineWidth',0.5)
    hold on
    errorbar((1:5)',relativeMean,lowerError,upperError,...
        'LineWidth',0.5,'Color','k','LineStyle','none')
    title(['Relative Expression ',targetNames{i,:},' relative to condition ',conditionNamesFull{referenceCondition}])
    ylabel('Relative Expression');
    xlabel('Conditions')
    hold off
    outputFileName = fullfile(imagesFolder,'gene_relative_expression',[conditionNames{referenceCondition},'_',targetNames{i,1}]);
    saveas(figureRelativeConcentrations,outputFileName,'jpeg')
end


