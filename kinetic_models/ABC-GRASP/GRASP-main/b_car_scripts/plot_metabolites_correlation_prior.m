%% Plot heatmaps of metabolite distribution correlations, prior models
% This script generates a heatmap of the metabolite distributions Spearman
% rank correlation coeffcients. It only creates one heatmap at a time, The
% variable 'i' in the 'Load data' section can be modified to define the
% experimental condition and structure type to be used. This script uses
% the prior ensembles.

% It requires as input a prior ensemble model located in the folder
% 'GRASP-main/io/ulfiltered_output'.Other variables are defined to include
% nomenclatures in their figures and filenames

% The output are correlation heatmaps. They are stored in the folder
% 'GRASP-main/io/figures/metabolites_correlation_prior', and their
% filenames indicate the condition and the type of model structure.


%% Clean variables and state matlab standards

clear, close all
rng('default');                 % for reproducibility
format longE
digits(256)                     % to prevent errors in calculations
%% Input

modelIDs = {'b_car_rejection_model_detailed_D010';'b_car_rejection_model_detailed_D025';...
    'b_car_rejection_model_regulated_D010';'b_car_rejection_model_regulated_D025';...
    'b_car_rejection_model_simple_D010';'b_car_rejection_model_simple_D025'};                   % names of the models in the specified folder

fileNames = {'detailed_4D010';'detailed_4D025';...
    'regulated_4D010';'regulated_4D025';...
    'simple_4D010';'simple_4D025'};                                                             % filenames for the figures

modelTypes = {'detailed_D010';'detailed_D025';...
    'regulated_D010';'regulated_D025';...
    'simple_D010';'simple_D025'};                                                               % structure types and growth rate

conditionNames = {'4D010';'4D025';'4D010';'4D025';'4D010';'4D025'};                             % experimental conditions (strain and growth rate)

conditionFullNames = {'\beta-car4 \mu = 0.101h^{-1}';'\beta-car4 \mu = 0.254h^{-1}';...
    '\beta-car4 \mu = 0.101h^{-1}';'\beta-car4 \mu = 0.254h^{-1}';...
    '\beta-car4 \mu = 0.101h^{-1}';'\beta-car4 \mu = 0.254h^{-1}'};                             % full name of conditions


cmap_bp = [ones(50,1),linspace(0.5,0.99,50)',linspace(0,0.98,50)';
    1,1,1;
    linspace(0.98,0,50)',linspace(0.98,0,50)',ones(50,1)];                                      % heatmaps colors

% cmap_bp = [ones(50,1),linspace(0.5,0.99,50)',linspace(0,0.98,50)';
%     1,1,1;
%     linspace(0.99,0.5,50)',linspace(0.98,0,50)',ones(50,1)];

% cmap_bp = [linspace(0.5,0.99,50)',linspace(0,0.98,50)',ones(50,1);
%     1,1,1;
%     linspace(0.98,00,50)',linspace(0.98,00,50)',linspace(0.98,0,50)'];

%% Load Data
% Modfify variable 'i' to choose a different ensemble model

i = 6;
modelID = modelIDs{i,1};
modelType = modelTypes{i,1};
conditionName = conditionNames{i,1};
conditionFullName = conditionFullNames{i,1};
fileName = fileNames{i,1};
load(fullfile('..', 'io', 'unfiltered_output', [modelID, '.mat']));
ensemble = ensembleUnfiltered;

%% Create Main Variables

% Optimization & simulation parameters

metNames = ensemble.mets;%(ensemble.metsActive);
modelsN = length(ensemble.populations.strucIdx);

metN = length(metNames);

for i= 1:length(metNames)
    metTemp = strsplit(metNames{i},'m_');
    metNames{i} = metTemp{2};
end

metaboliteConcentrations = zeros(modelsN,metN);

for i = 1:modelsN
    metaboliteConcentrations(i,:) = log10(ensemble.populations.models(i).metConcRef');
end

[corrMatrix,pvalueMatrix] = corr(metaboliteConcentrations,'Type','Spearman');

%% Plot Correlation

nx = metN;
ny = metN;
mx = linspace(0.5,nx+0.5,nx+1);
my = linspace(0.5,ny+0.5,ny+1);

corrFigure = figure('Name',conditionName);
corrFigure.Position = [500 300 300 300];
imagesc(rot90(corrMatrix)); hold on
set(gca,'DefaultAxesTickLabelInterpreter','none')
set(gca,'TickLabelInterpreter','none')
set(gca,'xticklabel',[],'yticklabel',[],'xtick',1:metN,'ytick',1:metN)
set(gca,'xticklabel',metNames,'yticklabel',flip(metNames))
set(gca,'units','points','position',[40,40,150,150])
set(gca,'TickDir','out')
set(gca,'LineWidth',0.5)
set(gca,'XTickLabelRotation',90)
xlabel('Metabolites','FontSize',7)
ylabel('Metabolites','FontSize',7)
title({'Metabolites Concentration correlation coefficients: ',conditionFullName})
set(gca,'FontSize',6,'FontName','arial')
caxis([-1 1])
colormap(corrFigure,cmap_bp)
cb = colorbar;
set(cb,'LineWidth',0.5)
set(cb,'TickDir','out')
set(cb,'FontSize',7,'FontName','arial')
cbparameters = cb.Position;
set(cb,'Position',[0.875,cbparameters(2),0.025,cbparameters(4)])
mg = mesh(mx,my,zeros([nx,ny]+1));
mg.FaceColor = 'none';
mg.EdgeColor = 'k';
hold off

%% Export Results

outputFig = fullfile('..','io','figures','metabolites_correlation_prior','mat',[fileName,'_metCorrelation.fig']);
savefig(corrFigure,outputFig)
outputJpg = fullfile('..','io','figures','metabolites_correlation_prior','jpg',[fileName,'_metCorrelation.jpg']);
saveas(corrFigure,outputJpg)
outputSvg = fullfile('..','io','figures','metabolites_correlation_prior','svg',[fileName,'_metCorrelation.svg']);
saveas(corrFigure,outputSvg)
outputPng = fullfile('..','io','figures','metabolites_correlation_prior','png',[fileName,'_metCorrelation.png']);
saveas(corrFigure,outputPng)

corrTable1 = array2table(corrMatrix,"VariableNames",metNames);
corrTable2 = cell2table(metNames,"VariableNames","Metabolite");
corrTable = [corrTable2,corrTable1];
OutputXls = fullfile('..','io','figures','metabolites_correlation_unfiltered','xls','metCorrelation.xlsx');
writetable(corrTable,OutputXls,'Sheet',fileName);


