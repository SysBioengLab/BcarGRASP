%% Plot heatmap of mean relative transcript concentrations at high growth rate
% This script plots the relative transcript concentrations of strains
% beta-car3 and beta-car2 at high growth rate, using the strain beta-car4 as
% reference

% It uses as input the data generated from 'processData', and stored in
% 'Output/Conditions Relative Limits.xlsx'

% The output is a figure with the heatmap for when beta-car4 is used as
% reference, another figure with the values inversed to consider beta-car2
% as reference for illustrative purposes, and a figure to plot not only the
% means, but also the confidence intervals when beta-car4 is used as
% reference. All images are stored in 'Output/Images/Heatmaps'. The
% filename indicates the reference condition and the growth rate.

%% Clean variables and state matlab standards

clear, close all
rng('default');                 % for reproducibility
format longE
digits(256)

%% Create colormap for heatmap

cmap_cm = [linspace(1,0.02,50)',zeros(50,1),linspace(1,0.02,50)';
    0,0,0;
    zeros(50,1),linspace(0.02,1,50)',linspace(0.02,1,50)'];

conditionFullNames = {'\beta-car2 \mu = 0.254h^{-1}';
    '\beta-car3 \mu = 0.254h^{-1}';
    '\beta-car4 \mu = 0.254h^{-1}'};

conditionNames = {'\beta-car2';
    '\beta-car3';
    '\beta-car4'};

geneFullNames = {'ERG10';'ERG13';'HMG1';'HMG2';'ERG12';'ERG8';'MVD1';'IDI1';
    'ERG20';'BTS1';'CrtE';'ERG9';'CrtI';'CrtYB'};

labelsFullNames = {'95% Lower CI';'Mean';'95% Upper CI';
    '95% Lower CI';'Mean';'95% Upper CI';
    '95% Lower CI';'Mean';'95% Upper CI'};

cLimitsExpression = [-2 2];
cValuesExpression = linspace(cLimitsExpression(:,1),cLimitsExpression(:,2),9);
cCellExpression = num2cell(cValuesExpression);
%cCellExpression{1} = ['<',num2str(cLimitsExpression(:,1))];
%cCellExpression{end} = ['>',num2str(cLimitsExpression(:,2))];

cLimitsExpressionComplete = [-3 3];
cValuesExpressionComplete = linspace(cLimitsExpressionComplete(:,1),cLimitsExpressionComplete(:,2),9);
cCellExpressionComplete = num2cell(cValuesExpressionComplete);
%cCellExpressionComplete{1} = ['<',num2str(cLimitsExpressionComplete(:,1))];
%cCellExpressionComplete{end} = ['>',num2str(cLimitsExpressionComplete(:,2))];

%% Get Data

expressionData = readtable('output\ConditionsRelativeLimits.xlsx','Sheet','4B');

meanExpressionColumns = [5,11,17];
meanExpression = expressionData{:,meanExpressionColumns+1};
meanLogExpression = log2(meanExpression);

minMeanLogExpression = min(meanLogExpression,[],'all');
maxMeanLogExpression = max(meanLogExpression,[],'all');


%% Plot Heatmap relatvized to condition 4 D01

nx = length(conditionFullNames);
ny = length(geneFullNames);
mx = linspace(0.5,nx+0.5,nx+1);
my = linspace(0.5,ny+0.5,ny+1);

expressionFigure = figure();
expressionFigure.Position = [400 250 250 250];
imagesc(meanLogExpression); hold on
set(gca,'DefaultAxesTickLabelInterpreter','tex')
set(gca,'TickLabelInterpreter','tex')
set(gca,'xticklabel',[],'yticklabel',[],'xtick',1:length(conditionFullNames),'ytick',1:length(geneFullNames))
set(gca,'xticklabel',conditionNames,'yticklabel',geneFullNames)
set(gca,'units','points','position',[60,25,125,125])
set(gca,'TickDir','out')
set(gca,'LineWidth',1)
set(gca,'XTickLabelRotation',0)
xlabel('Strains')
ylabel({'Transcripts log2-fold change';'in average relative expression'})
title({'Gene Relative Expression to ',conditionFullNames{3,1}})
set(gca,'FontSize',8,'FontName','arial')
caxis(cLimitsExpression)
colormap(expressionFigure,cmap_cm)
cb = colorbar;
set(cb,'LineWidth',1)
set(cb,'TickDir','out')
set(cb,'FontSize',8,'FontName','arial')
set(cb,'Ticks',cValuesExpression)
set(cb,'TickLabels',cCellExpression)
mg = mesh(mx,my,zeros([ny,nx]+1));
mg.FaceColor = 'none';
mg.EdgeColor = 'k';
hold off

saveas(expressionFigure,fullfile('output\images','heatmaps','svg','heatmap_D025_4D025.svg'))
saveas(expressionFigure,fullfile('output\images','heatmaps','jpg','heatmap_D025_4D025.jpg'))


%% Plot Heatmap relativized to condition 2 D01

meanExpression2D01 = meanExpression./repmat(meanExpression(:,1),1,3);
meanLogExpression2D01 = log2(meanExpression2D01);

minMeanLogExpression2D01 = min(meanLogExpression2D01,[],'all');
maxMeanLogExpression2D01 = max(meanLogExpression2D01,[],'all');

expressionFigure2 = figure();
expressionFigure2.Position = [400 250 250 250];
imagesc(meanLogExpression2D01); hold on
set(gca,'DefaultAxesTickLabelInterpreter','tex')
set(gca,'TickLabelInterpreter','tex')
set(gca,'xticklabel',[],'yticklabel',[],'xtick',1:length(conditionFullNames),'ytick',1:length(geneFullNames))
set(gca,'xticklabel',conditionNames,'yticklabel',geneFullNames)
set(gca,'units','points','position',[60,25,125,125])
set(gca,'TickDir','out')
set(gca,'LineWidth',1)
set(gca,'XTickLabelRotation',0)
xlabel('Strains')
ylabel({'Transcripts log2-fold change';' in average relative expression'})
title({'Gene Relative Expression to ',conditionFullNames{1,1}})
set(gca,'FontSize',8,'FontName','arial')
caxis(cLimitsExpression)
colormap(expressionFigure2,cmap_cm)
cb = colorbar;
set(cb,'LineWidth',1)
set(cb,'TickDir','out')
set(cb,'FontSize',8,'FontName','arial')
set(cb,'Ticks',cValuesExpression)
set(cb,'TickLabels',cCellExpression)
mg = mesh(mx,my,zeros([ny,nx]+1));
mg.FaceColor = 'none';
mg.EdgeColor = 'k';
hold off

saveas(expressionFigure2,fullfile('output\images','heatmaps','svg','heatmap_D025_2D025.svg'))
saveas(expressionFigure2,fullfile('output\images','heatmaps','jpg','heatmap_D025_2D025.jpg'))

%% Plot Complete Heatmap with confidence Intervals

meanExpressionComplete = expressionData{:,[2:4,8:10,14:16]+3};
meanLogExpressionComplete = log2(meanExpressionComplete);

minMeanLogExpressionComplete = min(meanLogExpressionComplete,[],'all');
maxMeanLogExpressionComplete = max(meanLogExpressionComplete,[],'all');

nx = length(labelsFullNames);
ny = length(geneFullNames);
mx = linspace(0.5,nx+0.5,nx+1);
my = linspace(0.5,ny+0.5,ny+1);

nx2 = length(labelsFullNames);
ny2 = length(geneFullNames);
mx2 = linspace(0.5,nx2+0.5,nx2/3+1);
my2 = linspace(0.5,ny2+0.5,ny2+1);

expressionFigure3 = figure();
expressionFigure3.Position = [400 250 700 450];
imagesc(meanLogExpressionComplete); hold on
set(gca,'DefaultAxesTickLabelInterpreter','tex')
set(gca,'TickLabelInterpreter','tex')
set(gca,'xticklabel',[],'yticklabel',[],'xtick',1:length(labelsFullNames),'ytick',1:length(geneFullNames))
set(gca,'xticklabel',labelsFullNames,'yticklabel',geneFullNames)
set(gca,'units','points','position',[50,75,450,225])
set(gca,'TickDir','out')
set(gca,'LineWidth',1)
set(gca,'XTickLabelRotation',45)
xlabel('Conditions')
ylabel("Transcripts log2-fold change in average relative expression")
title(['Gene Relative Expression to ',conditionFullNames{3,1}])
set(gca,'FontSize',8,'FontName','arial')
caxis(cLimitsExpressionComplete)
colormap(expressionFigure3,cmap_cm)
cb = colorbar;
set(cb,'LineWidth',1)
set(cb,'TickDir','out')
set(cb,'FontSize',8,'FontName','arial')
set(cb,'Ticks',cValuesExpressionComplete)
set(cb,'TickLabels',cCellExpressionComplete)
% mg = mesh(mx,my,zeros([ny,nx]+1));
% mg.FaceColor = 'none';
% mg.EdgeColor = 'k';
mg2 = mesh(mx2,my2,zeros([ny2,nx2/3]+1));
mg2.FaceColor = 'none';
mg2.EdgeColor = 'k';
mg2.LineWidth = 2;
hold off

saveas(expressionFigure3,fullfile('output\images','heatmaps','svg','heatmap_D025_4D025_Complete.svg'))
saveas(expressionFigure3,fullfile('output\images','heatmaps','jpg','heatmap_D025_4D025_Complete.jpg'))
