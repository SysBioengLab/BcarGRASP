%% Plot independent fluxes of posterior ensemble models
% This scripts creates boxplots with for the three model structures used
% within one figure. It generates plots for the fluxes EX_b_car, EX_lyc,
% and ERG9b. Each plot represents a different reference conditions (low or
% high growth rate). Each plot contains the fluxes simulated for each 
% experimental condition.

% It requires as input the filtered ensembles generated with
% 'filter_ensemble' stored at 'GRASP-main/io/output'

% The output figures are stored in the folder
% 'GRASP-main/io/figures/fluxes_posterior_boxplot'

%% Clean variables and state matlab standards

clear, close all
rng('default');                 % for reproducibility
format longE
digits(256)                     % to prevent error in calculations

%% Path

addpath(fullfile('..', 'matlab_code', 'analysisFxns'), ...
        fullfile('..', 'matlab_code', 'ensembleFxns'), ...
        fullfile('..', 'matlab_code', 'patternFxns'));

%% Input

% Generate plots for high growth rate
ensembleFilteredDirs = {'../io/output/b_car_rejection_model_simple_D025.mat';
    '../io/output/b_car_rejection_model_regulated_D025.mat';
    '../io/output/b_car_rejection_model_detailed_D025.mat'};

referenceConditionName = '4D025';
adjustedConditionsNames = {'3D025';'2D025'};
referenceconditionFullName = '\beta-car4 \mu = 0.254h^{-1}';

% % Generate plots for low growth rate
% ensembleFilteredDirs = {'../io/output/b_car_rejection_model_simple_D010.mat';
%     '../io/output/b_car_rejection_model_regulated_D010.mat';
%     '../io/output/b_car_rejection_model_detailed_D010.mat'};
% 
% referenceConditionName = '4D010';
% adjustedConditionsNames = {'3D010';'2D010'};
% referenceconditionFullName = '\beta-car4 \mu = 0.101h^{-1}';


modelTypes = {'simple';'regulated';'detailed'};

% boxplotNames = {'\beta-car3 simple';'\beta-car3 regulated';'\beta-car3 detailed';
%     '\beta-car2 simple';'\beta-car2 regulated';'\beta-car2 detailed'};
boxplotNames = {'simple';'regulated';'detailed';
    'simple';'regulated';'detailed'};

adjustedConditionsFullNames = {'\beta-car3';'\beta-car2'};

targetFluxes = {'EX_lyc';'EX_b_car';'ERG9b'};                                       % names of fluxes to be plotted
colorPalette = [0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0 0.4470 0.7410];         % colors for the fluxes


%% Load Data

ensembleFiltered = cell(length(ensembleFilteredDirs),1);

for ix = 1:length(ensembleFilteredDirs)
    ensembleFilteredFile = load(ensembleFilteredDirs{ix,1});
    ensembleFiltered{ix,1} = ensembleFilteredFile.ensembleFiltered;
end

%% Extract data

structuresN = length(ensembleFiltered);
conditionsN = size(ensembleFiltered{1,1}.populations.simFluxes{1,1},2);

metNames = ensembleFiltered{ix,1}.mets;
fluxNames = ensembleFiltered{ix,1}.rxns;

for ix = 1:numel(metNames)
    metTemp = strsplit(metNames{ix},'m_');
    metNames{ix} = metTemp{2};
end

for ix = 1:numel(fluxNames)
    fluxTemp = strsplit(fluxNames{ix},'r_');
    fluxNames{ix} = fluxTemp{2};
end

indexERG9b = find(strcmp(fluxNames,'ERG9b'));
indexEX_lyc = find(strcmp(fluxNames,'EX_lyc'));
indexEX_b_car = find(strcmp(fluxNames,'EX_b_car'));

modelsN = zeros(length(ensembleFiltered),1);
for ix = 1:structuresN
    modelsN(ix,1) = length(ensembleFiltered{ix,1}.populations.simFluxes);
end

simulatedFluxes = cell(conditionsN,3);

for ix = 1:structuresN
    for jx = 1:conditionsN
        for hx = 1:modelsN(ix,1)
            simulatedFluxes{jx,1}(hx,ix) = ensembleFiltered{ix,1}.populations.simFluxes{hx,1}(indexEX_lyc,jx);
            simulatedFluxes{jx,2}(hx,ix) = ensembleFiltered{ix,1}.populations.simFluxes{hx,1}(indexEX_b_car,jx);
            simulatedFluxes{jx,3}(hx,ix) = ensembleFiltered{ix,1}.populations.simFluxes{hx,1}(indexERG9b,jx);
        end
    end
end

 
expMeanFlux(:,1) = ensembleFiltered{1,1}.expFluxes(indexEX_lyc,:)';
expMeanFlux(:,2) = ensembleFiltered{1,1}.expFluxes(indexEX_b_car,:);
expMeanFlux(:,3) = ensembleFiltered{1,1}.expFluxes(indexERG9b,:);

expStdFlux(:,1) = ensembleFiltered{1,1}.expFluxesStd(indexEX_lyc,:)';
expStdFlux(:,2) = ensembleFiltered{1,1}.expFluxesStd(indexEX_b_car,:);
expStdFlux(:,3) = ensembleFiltered{1,1}.expFluxesStd(indexERG9b,:);


for ix = 1:size(simulatedFluxes,2)
    
    boxplotFigure = figure();
    boxplotFigure.Position = [500 200 300 200];
    rectanglePosition = [0.75,expMeanFlux(1,ix)-2*expStdFlux(1,ix),2.5,4*expStdFlux(1,ix)];
    rectangle('Position',rectanglePosition,'FaceColor',colorPalette(ix,:),'EdgeColor',colorPalette(ix,:))
    hold on
    rectanglePosition = [4.75,expMeanFlux(2,ix)-2*expStdFlux(2,ix),2.5,4*expStdFlux(2,ix)];
    rectangle('Position',rectanglePosition,'FaceColor',colorPalette(ix,:),'EdgeColor',colorPalette(ix,:))
    boxplot([simulatedFluxes{1,ix},simulatedFluxes{2,ix}],'Symbol','k.', 'Colors','k','PlotStyle','traditional','MedianStyle','line','Labels',boxplotNames, ...
        'OutlierSize',5,'Positions',[(1:length(modelTypes)),(1:length(modelTypes))+4])
    title({'Boxplot of the adjusted fluxes for reference ',referenceconditionFullName})
    ylabel([targetFluxes{ix,1}],'Interpreter','none')
    xlabel('Adjusted Condition and Structure Type')
    ylimActual = get(gca,'yLim');
    ylim([0 ylimActual(1,2)])
    xticks([(1:length(modelTypes)),(1:length(modelTypes))+4])
    grid('on')
    set(gca,'units','points','position',[30,50,180,60])
    set(gca,'TickDir','out')
    set(gca,'LineWidth',1)
    set(gca,'FontSize',7,'FontName','arial')
    set(gca,'TickLabelInterpreter','tex')
    hold off
    saveas(boxplotFigure,fullfile('..','io','figures','fluxes_posterior_boxplot','png',[referenceConditionName,'_',targetFluxes{ix,1},'_boxplot.png']))
    saveas(boxplotFigure,fullfile('..','io','figures','fluxes_posterior_boxplot','jpg',[referenceConditionName,'_',targetFluxes{ix,1},'_boxplot.jpg']))
    saveas(boxplotFigure,fullfile('..','io','figures','fluxes_posterior_boxplot','svg',[referenceConditionName,'_',targetFluxes{ix,1},'_boxplot.svg']))
    saveas(boxplotFigure,fullfile('..','io','figures','fluxes_posterior_boxplot','mat',[referenceConditionName,'_',targetFluxes{ix,1},'_boxplot.fig']))
    
end
