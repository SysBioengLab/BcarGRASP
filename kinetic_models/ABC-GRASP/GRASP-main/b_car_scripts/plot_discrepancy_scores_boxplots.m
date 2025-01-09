%% Creates boxplots with the discrepancy scores for each ensemble
% It creates  the boxplots of the discrepancy scores for each ensemble,
% distiguishing between the experimental conditions obtained for beta-car2
% and beta-car3. It also plot the average score for each model particle

% It requires the backup files generated using the codes stored in the
% 'ensemble_sampling' folder.

% The output boxplots are stored in the
% 'GRASP-main/io/figure/scores_boxplot' folder

%% Clean variables and state matlab standards

clear, close all
rng('default');                 % for reproducibility
format longE
digits(256)                     % to prevent errors in calculations

%% Input

rejectionFileDirs = {'../../ensemble_sampling/Simple_D010/backup_b_car_rejection_model_simple_D010';
    '../../ensemble_sampling/Regulated_D010/backup_b_car_rejection_model_regulated_D010';
    '../../ensemble_sampling/Detailed_D010/backup_b_car_rejection_model_detailed_D010';
    '../../ensemble_sampling/Simple_D025/backup_b_car_rejection_model_simple_D025';
    '../../ensemble_sampling/Regulated_D025/backup_b_car_rejection_model_regulated_D025';
    '../../ensemble_sampling/Detailed_D025/backup_b_car_rejection_model_detailed_D025'};                    % direction of the files

conditionNames = {'4D010';'4D010';'4D010';'4D025';'4D025';'4D025'};                                         % condition names
structureTypeNames = {'simple';' regulated';
    'detailed';'simple';
    'regulated';'detailed'};                                                                                % structure type of each ensemble model
scoresNames = {'\beta-car2';'\beta-car3';'Average'};                                                        % names of the scores, 2 is for scores of beta-car2, 3 is for scores of beta-car3, average for the mean of both
outputFilenamesUnfiltered = {'2_unfiltered';'3_unfiltered';'average_unfiltered'};                           % filenames of the ouput figures for prior dicrepancy scores, 2 is for scores of beta-car2, 3 is for scores of beta-car3, average for the mean of both
outputFilenamesFiltered = {'2_filtered';'3_filtered';'average_filtered'};                                   % filenames of the ouput figures for posterior dicrepancy scores, 2 is for scores of beta-car2, 3 is for scores of beta-car3, average for the mean of both

modelsN = 10000; %nValidModels if all models want to be checked                                             % number of models to be considered from the prior (only first N models)
thresholdProportion = 0.01;  % 1 if all the generated models want to be checked                             % proportion of models to be selected from the prior to constitute the posterior

%% Extract scores

discrepancyScores = cell(3,1);
discrepancyScoresFiltered = cell(3,1);
discrepancyScoreThreshold = zeros(1,length(rejectionFileDirs));

for ix = 1:length(rejectionFileDirs)
    rejectionFile = load(rejectionFileDirs{ix,1});
    tolScores = rejectionFile.tolScoreBackup;
    tolScoresValid = rejectionFile.isModelValidBackup;
    discrepancyScoresTemp = tolScores(tolScoresValid,:);
    discrepancyScoresTemp = discrepancyScoresTemp(1:modelsN,:);
    discrepancyScores{1,1}(:,ix) = discrepancyScoresTemp(:,2);
    discrepancyScores{2,1}(:,ix) = discrepancyScoresTemp(:,1);
    discrepancyScores{3,1}(:,ix) = mean(discrepancyScoresTemp,2);
    discrepancyScoreThreshold(1,ix) = prctile(mean(discrepancyScoresTemp,2),thresholdProportion*100);
    filterTemp = (mean(discrepancyScoresTemp,2) <= discrepancyScoreThreshold(1,ix));
    discrepancyScoresFiltered{1,1}(:,ix) = discrepancyScoresTemp(filterTemp,2);
    discrepancyScoresFiltered{2,1}(:,ix) = discrepancyScoresTemp(filterTemp,1);
    discrepancyScoresFiltered{3,1}(:,ix) = mean(discrepancyScoresTemp(filterTemp,:),2);
end

%% Create boxplots

for ix = 1:length(discrepancyScores)
    boxplotFigure = figure();
    boxplotFigure.Position = [500 200 300 200];
    boxplot(discrepancyScores{ix,1},'Symbol','k.', 'Colors','k','PlotStyle','traditional','MedianStyle','line', ...
        'Labels',structureTypeNames,'Jitter',0.5,'OutlierSize',1,'Positions',[(1:3),(5:7)]);
    hold on
    if ix == 3
        for jx = [(1:3),(5:7)]
            plot([jx-0.2;jx+0.2],[discrepancyScoreThreshold(1,ix);discrepancyScoreThreshold(1,ix)],'r-','LineWidth',0.5)
        end
    end
    title({'Boxplot of discrepancy scores'})
    xticks([(1:3),(5:7)])
    ylabel([scoresNames{ix,1}],'Interpreter','tex')
    xlabel('Condition and Structure Type')
    %ylim([0 3])
    grid('on')
    set(gca,'units','points','position',[30,50,180,90])
    set(gca,'TickDir','out')
    set(gca,'LineWidth',1)
    set(gca,'FontSize',8,'FontName','arial')
    set(gca,'TickLabelInterpreter','tex')
    hold off
    saveas(boxplotFigure,fullfile('..','io','figures','scores_boxplot','png',[outputFilenamesUnfiltered{ix,1},'_boxplot.png']))
    saveas(boxplotFigure,fullfile('..','io','figures','scores_boxplot','jpg',[outputFilenamesUnfiltered{ix,1},'_boxplot.jpg']))
    saveas(boxplotFigure,fullfile('..','io','figures','scores_boxplot','svg',[outputFilenamesUnfiltered{ix,1},'_boxplot.svg']))
    saveas(boxplotFigure,fullfile('..','io','figures','scores_boxplot','mat',[outputFilenamesUnfiltered{ix,1},'_boxplot.fig']))
end

for ix = 1:length(discrepancyScoresFiltered)
    boxplotFigure = figure();
    boxplotFigure.Position = [500 200 300 200];
    boxplot(discrepancyScoresFiltered{ix,1},'Symbol','k.', 'Colors','k','PlotStyle','traditional','MedianStyle','line', ...
        'Labels',structureTypeNames,'Jitter',0.5,'OutlierSize',5,'Positions',[(1:3),(5:7)]);
    hold on
    title({'Boxplot of filtered discrepancy scores'})
    xticks([(1:3),(5:7)])
    ylabel([scoresNames{ix,1}],'Interpreter','tex')
    xlabel('Condition and Structure Type')
    %ylim([0 3])
    grid('on')
    set(gca,'units','points','position',[30,50,180,90])
    set(gca,'TickDir','out')
    set(gca,'LineWidth',1)
    set(gca,'FontSize',8,'FontName','arial')
    set(gca,'TickLabelInterpreter','tex')
    hold off
    saveas(boxplotFigure,fullfile('..','io','figures','scores_boxplot','png',[outputFilenamesFiltered{ix,1},'_boxplot.png']))
    saveas(boxplotFigure,fullfile('..','io','figures','scores_boxplot','jpg',[outputFilenamesFiltered{ix,1},'_boxplot.jpg']))
    saveas(boxplotFigure,fullfile('..','io','figures','scores_boxplot','svg',[outputFilenamesFiltered{ix,1},'_boxplot.svg']))
    saveas(boxplotFigure,fullfile('..','io','figures','scores_boxplot','mat',[outputFilenamesFiltered{ix,1},'_boxplot.fig']))
end


