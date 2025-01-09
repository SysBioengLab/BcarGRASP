%% Plot fluxes obtained from the first set of simulations, selected simulations
% This script creates the boxplots and barplots of the fluxes EX_b_car,
% EX_lyc, ERG9b, and ERG10 for each of the simulations performed. It also
% plots the log-ratio between EX_b_car and EX_lyc. Each boxplot includes
% a selected set of the perturbations placed on an ensemble model. Besides
% the fluxes for each simulated model particle, they also include the
% confidence intervals of the reference condition from where the particles
% were originally sampled. Only one ensemble model is plotted at a time. To
% select the ensemble, the variable 'ix' in the 'Load data' can be modified.

% It requires the files of the simulations, both for step changes in
% enzymes and metabolites. The location of the filefolder can be specified.
% Other variables are defined to include nomenclatures in their figures and
% filenames.

% The output are the flux distribution boxplots of the first set of
% simulations. They are stored in
% 'GRASP-main/io/figure/simulations_boxplots' and 
% 'GRASP-main/io/figure/simulations_barplots'. Their filenames indicate the
% structure type, reference condition and flux/ratio.

%% Clean variables and state matlab standards

clear, close all
rng('default');                 % for reproducibility
format longE
digits(256)                     % to prevent errors in calculations

%% Path
addpath(fullfile('..', 'matlab_code', 'analysisFxns'), ...
        fullfile('..', 'matlab_code', 'ensembleFxns'), ...
        fullfile('..', 'matlab_code', 'patternFxns'));



%% Input

enzymesFluxesDirs = {'../../ensemble_simulation/enzymes_simple_D010/backup_b_car_simulation_enzymes_model_simple_D010.mat';
    '../../ensemble_simulation/enzymes_simple_D025/backup_b_car_simulation_enzymes_model_simple_D025.mat';
    '../../ensemble_simulation/enzymes_regulated_D010/backup_b_car_simulation_enzymes_model_regulated_D010.mat';
    '../../ensemble_simulation/enzymes_regulated_D025/backup_b_car_simulation_enzymes_model_regulated_D025.mat';
    '../../ensemble_simulation/enzymes_detailed_D010/backup_b_car_simulation_enzymes_model_detailed_D010.mat';
    '../../ensemble_simulation/enzymes_detailed_D025/backup_b_car_simulation_enzymes_model_detailed_D025.mat'};                             % filenames of enzyme modification simulations

metabolitesFluxesDirs = {'../../ensemble_simulation/metabolites_simple_D010/backup_b_car_simulation_metabolites_model_simple_D010.mat';
    '../../ensemble_simulation/metabolites_simple_D025/backup_b_car_simulation_metabolites_model_simple_D025.mat';
    '../../ensemble_simulation/metabolites_regulated_D010/backup_b_car_simulation_metabolites_model_regulated_D010.mat';
    '../../ensemble_simulation/metabolites_regulated_D025/backup_b_car_simulation_metabolites_model_regulated_D025.mat';
    '../../ensemble_simulation/metabolites_detailed_D010/backup_b_car_simulation_metabolites_model_detailed_D010.mat';
    '../../ensemble_simulation/metabolites_detailed_D025/backup_b_car_simulation_metabolites_model_detailed_D025.mat'};                     % filenames of metabolite modification simulations

referenceConditionNames = {'4D010';'4D025';'4D010';'4D025';'4D010';'4D025'};                                                                % reference conditions
referenceConditionFullNames = {'\beta-car4 \mu = 0.101h^{-1}';'\beta-car4 \mu = 0.254h^{-1}';
    '\beta-car4 \mu = 0.101h^{-1}';'\beta-car4 \mu = 0.254h^{-1}';
    '\beta-car4 \mu = 0.101h^{-1}';'\beta-car4 \mu = 0.254h^{-1}'};                                                                         % reference condition full names

modelTypes = {'simple';'simple';'regulated';'regulated';'detailed';'detailed'};                                                             % structure type

perturbationNames = {'-50%';'-30%';'+40%';'+100%';'-67%';'-33%';'+33%';'+67%';'-33%';'+33%';'-50%';'+50%';'-30%';'+40%';
    '-50%';'-30%';'+40%';'+100%';'-30%';'+40%';'-30%';'+40%';'-30%';'+40%'};                                                                % x label of perturbations

enzymesPerturbations = 1:14;                                                                                                                % enzyme perturbations to be included
metabolitesPerturbations = 1:10;                                                                                                            % metabolite perturbations to be included

targetFluxes = {'EX_lyc';'EX_b_car';'ERG9b';'ERG10'};                                                                                       % fluxes to be ploted
colorPaletteReference = [0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0 0.4470 0.7410;0.4660 0.6740 0.1880;0.4940 0.1840 0.5560];              % first color palette
% colorPalettePerturbations = [0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840;0.4940 0.1840 0.5560];   % second color palette


%% Load Data
% Modify ix to choose the ensemble to be plotted

ix = 6;
enzymesFluxesDir = enzymesFluxesDirs{ix,1};
metabolitesFluxesDir = metabolitesFluxesDirs{ix,1};
referenceConditionName = referenceConditionNames{ix,1};
referenceConditionFullName = referenceConditionFullNames{ix,1};
modelType = modelTypes{ix,1};

enzymesFluxesFile = load(enzymesFluxesDir);
metabolitesFluxesFile = load(metabolitesFluxesDir);
ensemble = enzymesFluxesFile.ensemble;
simFluxesEnzymes = enzymesFluxesFile.simFluxes;
simFluxesMetabolites = metabolitesFluxesFile.simFluxes;


%% Extract data
simFluxes = {simFluxesEnzymes{enzymesPerturbations},simFluxesMetabolites{metabolitesPerturbations}}';
metNames = ensemble.mets;
fluxNames = ensemble.rxns;

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
indexERG10 = find(strcmp(fluxNames,'ERG10'));

perturbationsN = length(perturbationNames);%size(simFluxes,1);
modelsN = zeros(perturbationsN,1);
for ix = 1:perturbationsN
    modelsN(ix,1) = size(simFluxes{ix,1},1);
end

simulatedFluxes = cell(4,1);

for ix = 1:perturbationsN
    for jx = 1:modelsN(ix,1)
        simulatedFluxes{1,1}(jx,ix) = simFluxes{ix,1}(jx,indexEX_lyc);
        simulatedFluxes{2,1}(jx,ix) = simFluxes{ix,1}(jx,indexEX_b_car);
        simulatedFluxes{3,1}(jx,ix) = simFluxes{ix,1}(jx,indexERG9b);
        simulatedFluxes{4,1}(jx,ix) = simFluxes{ix,1}(jx,indexERG10);
    end
end


expMeanFlux(:,1) = ensemble.fluxRef(indexEX_lyc,:)';
expMeanFlux(:,2) = ensemble.fluxRef(indexEX_b_car,:);
expMeanFlux(:,3) = ensemble.fluxRef(indexERG9b,:);
expMeanFlux(:,4) = ensemble.fluxRef(indexERG10,:);

expStdFlux(:,1) = ensemble.fluxRefStd(indexEX_lyc,:)';
expStdFlux(:,2) = ensemble.fluxRefStd(indexEX_b_car,:);
expStdFlux(:,3) = ensemble.fluxRefStd(indexERG9b,:);
expStdFlux(:,4) = ensemble.fluxRefStd(indexERG10,:);

%% Boxplot of simulation fluxes

for ix = 1:size(simulatedFluxes,1)
    
    boxplotFigure = figure();
    boxplotFigure.Position = [500 200 300 200];
    rectanglePosition = [0.6,expMeanFlux(1,ix)-2*expStdFlux(1,ix),perturbationsN-0.2,4*expStdFlux(1,ix)];
    rectangle('Position',rectanglePosition,'FaceColor',colorPaletteReference(ix,:),'EdgeColor','none')
    hold on
%     plot([0.6;perturbationsN+0.4], [expMeanFlux(1,ix) expMeanFlux(1,ix)],'k-')
    boxplot([simulatedFluxes{ix,1}],'Symbol','k.', 'Colors','k','PlotStyle','traditional','MedianStyle','line','Labels',perturbationNames, ...
        'OutlierSize',5,'Positions',1:perturbationsN)
    title({'Boxplot of the simulated perturbations for reference ',referenceConditionFullName})
    ylabel([targetFluxes{ix,1}],'Interpreter','none')
    xlabel('Performed Perturbation')
    ylimActual = get(gca,'yLim');
    ylim([0 ylimActual(1,2)])
    xlim([0 perturbationsN+1])
    xticks(1:perturbationsN)
    grid('on')
    set(findobj(gca,'type','line'),'linew',0.5)
    set(gca,'units','points','position',[30,50,180,60])
    set(gca,'TickDir','out')
    set(gca,'LineWidth',1)
    set(gca,'FontSize',7,'FontName','arial')
    set(gca,'TickLabelInterpreter','tex')
%     h = findobj(gca,'Tag','Box');
%     for jx=1:length(h)
%         patch(get(h(jx),'XData'),get(h(jx),'YData'),colorPalettePerturbations(jx,:),'FaceAlpha',0.9);
%     end
    hold off
    saveas(boxplotFigure,fullfile('..','io','figures','simulations_boxplot','png',[modelType,'_',referenceConditionName,'_',targetFluxes{ix,1},'_boxplot.png']))
    saveas(boxplotFigure,fullfile('..','io','figures','simulations_boxplot','jpg',[modelType,'_',referenceConditionName,'_',targetFluxes{ix,1},'_boxplot.jpg']))
    saveas(boxplotFigure,fullfile('..','io','figures','simulations_boxplot','svg',[modelType,'_',referenceConditionName,'_',targetFluxes{ix,1},'_boxplot.svg']))
    saveas(boxplotFigure,fullfile('..','io','figures','simulations_boxplot','mat',[modelType,'_',referenceConditionName,'_',targetFluxes{ix,1},'_boxplot.fig']))   
end

%% Boxplot of simulation ratios

simulatedRatioBcarLyc = simulatedFluxes{2,1}./simulatedFluxes{1,1};
newMeanRatioBcarLyc = mean(log10(simulatedRatioBcarLyc),1);
meanRatioBcarLyc = expMeanFlux(:,2)/expMeanFlux(:,1);
oldMeanRatioBcarLyc = log10(meanRatioBcarLyc);

boxplotFigure = figure();
boxplotFigure.Position = [500 200 300 200];
lineplot = plot([0.6;perturbationsN+0.4],[meanRatioBcarLyc;meanRatioBcarLyc],'-','LineWidth',1,'Color',colorPaletteReference(5,:));
hold on
boxplot(simulatedRatioBcarLyc,'Symbol','k.', 'Colors','k','PlotStyle','traditional','MedianStyle','line','Labels',perturbationNames, ...
    'OutlierSize',5,'Positions',1:perturbationsN)
title({'Boxplot of the simulated perturbations for reference ',referenceConditionFullName})
ylabel('Log-ratio EX_b_car/EX_lyc','Interpreter','none')
xlabel('Performed Perturbation')
ylimActual = get(gca,'yLim');
ylim([0 ylimActual(1,2)])
%ylim([0 max(simulatedRatioBcarLyc(:,8))+0.5])
%ylim([0 150])
xlim([0 perturbationsN+1])
xticks(1:perturbationsN)
xtickangle(90)
grid('on')
set(findobj(gca,'type','line'),'linew',0.5)
set(gca,'units','points','position',[30,50,180,60])
set(gca,'TickDir','out')
set(gca,'LineWidth',1)
set(gca,'FontSize',7,'FontName','arial')
set(gca,'TickLabelInterpreter','tex')
ax = gca;
ax.YAxis.Scale = "log";
lineplot.LineWidth = 2;
hold off
saveas(boxplotFigure,fullfile('..','io','figures','simulations_boxplot','png',[modelType,'_',referenceConditionName,'_','ratio','_boxplot.png']))
saveas(boxplotFigure,fullfile('..','io','figures','simulations_boxplot','jpg',[modelType,'_',referenceConditionName,'_','ratio','_boxplot.jpg']))
saveas(boxplotFigure,fullfile('..','io','figures','simulations_boxplot','svg',[modelType,'_',referenceConditionName,'_','ratio','_boxplot.svg']))
saveas(boxplotFigure,fullfile('..','io','figures','simulations_boxplot','mat',[modelType,'_',referenceConditionName,'_','ratio','_boxplot.fig']))   


%% Barplot of simulation fluxes

meanFluxes = cell(4,1);
stdFluxes = cell(4,1);
relativeIncrease = cell(4,1);

for ix = 1:length(meanFluxes)
    meanFluxes{ix,1} = mean(simulatedFluxes{ix,1},1);
    stdFluxes{ix,1} = std(simulatedFluxes{ix,1},0,1);
    relativeIncrease{ix,1} = meanFluxes{ix,1}./expMeanFlux(:,ix);
end

meanRatio = mean(log10(simulatedRatioBcarLyc),1);
stdRatio = std(log10(simulatedRatioBcarLyc),0,1);

for ix = 1:length(meanFluxes)
    barplotFigure = figure();
    barplotFigure.Position = [500 200 300 200];
    %rectanglePosition = [0.6,expMeanFlux(1,ix)-2*expStdFlux(1,ix),perturbationsN-0.2,4*expStdFlux(1,ix)];
    %rectangle('Position',rectanglePosition,'FaceColor',colorPaletteReference(ix,:),'EdgeColor',[1 0 0],'LineWidth',1,'LineStyle','-')
    hold on
    lineplot1 = plot([0;perturbationsN+1],[expMeanFlux(1,ix);expMeanFlux(1,ix)],'-','Color',[1 0 0]);
    lineplot2 = plot([0;perturbationsN+1],[expMeanFlux(1,ix);expMeanFlux(1,ix)]-2*expStdFlux(1,ix),'--','Color',[1 0 0]);
    lineplot3 = plot([0;perturbationsN+1],[expMeanFlux(1,ix);expMeanFlux(1,ix)]+2*expStdFlux(1,ix),'--','Color',[1 0 0]);
    bar(meanFluxes{ix,1},0.6,'FaceColor',colorPaletteReference(ix,:),'EdgeColor',[0 0 0],'LineWidth',0.5)
    errorbar(1:perturbationsN,meanFluxes{ix,1},stdFluxes{ix,1},'LineStyle','none','LineWidth',0.5,'CapSize',0.5,'Color',[0 0 0])
    title({'Barplot of the simulated perturbations for reference ',referenceConditionFullName})
    ylabel([targetFluxes{ix,1}],'Interpreter','none')
    xlabel('Performed Perturbation')
    ylimActual = get(gca,'yLim');
    ylim([0 ylimActual(1,2)])
    xlim([0 perturbationsN+1])
    xticks(1:perturbationsN)
    xticklabels(perturbationNames)
    grid('on')
    set(findobj(gca,'type','line'),'linew',0.5)
    set(gca,'Box','on')
    set(gca,'units','points','position',[30,50,180,60])
    set(gca,'TickDir','out')
    set(gca,'LineWidth',1)
    set(gca,'FontSize',7,'FontName','arial')
    set(gca,'TickLabelInterpreter','tex')
    lineplot1.LineWidth = 1;
    lineplot2.LineWidth = 1;
    lineplot3.LineWidth = 1;
    saveas(barplotFigure,fullfile('..','io','figures','simulations_barplot','png',[modelType,'_',referenceConditionName,'_',targetFluxes{ix,1},'_boxplot.png']))
    saveas(barplotFigure,fullfile('..','io','figures','simulations_barplot','jpg',[modelType,'_',referenceConditionName,'_',targetFluxes{ix,1},'_boxplot.jpg']))
    saveas(barplotFigure,fullfile('..','io','figures','simulations_barplot','svg',[modelType,'_',referenceConditionName,'_',targetFluxes{ix,1},'_boxplot.svg']))
    saveas(barplotFigure,fullfile('..','io','figures','simulations_barplot','mat',[modelType,'_',referenceConditionName,'_',targetFluxes{ix,1},'_boxplot.fig'])) 
end

%% Barplot of simulation ratios

barplotFigure = figure();
barplotFigure.Position = [500 200 300 200];
hold on
lineplot = plot([0.6;perturbationsN+0.4],log10([meanRatioBcarLyc;meanRatioBcarLyc]),'-','LineWidth',1,'Color',[1 0 0]);
bar(meanRatio,0.6,'FaceColor',colorPaletteReference(5,:),'EdgeColor',[0 0 0],'LineWidth',0.5)
errorbar(1:perturbationsN,meanRatio,stdRatio,'LineStyle','none','LineWidth',0.5,'CapSize',0.5,'Color',[0 0 0])
title({'Barplot of the simulated perturbations for reference ',referenceConditionFullName})
ylabel('Log-ratio EX_b_car/EX_lyc','Interpreter','none')
xlabel('Performed Perturbation')
ylimActual = get(gca,'yLim');
%ylim([0 ylimActual(1,2)])
xlim([0 perturbationsN+1])
xticks(1:perturbationsN)
xticklabels(perturbationNames)
grid('on')
set(findobj(gca,'type','line'),'linew',0.5)
set(gca,'Box','on')
set(gca,'units','points','position',[30,50,180,60])
set(gca,'TickDir','out')
set(gca,'LineWidth',1)
set(gca,'FontSize',7,'FontName','arial')
set(gca,'TickLabelInterpreter','tex')
lineplot.LineWidth = 1;
saveas(barplotFigure,fullfile('..','io','figures','simulations_barplot','png',[modelType,'_',referenceConditionName,'_',targetFluxes{ix,1},'_boxplot.png']))
saveas(barplotFigure,fullfile('..','io','figures','simulations_barplot','jpg',[modelType,'_',referenceConditionName,'_',targetFluxes{ix,1},'_boxplot.jpg']))
saveas(barplotFigure,fullfile('..','io','figures','simulations_barplot','svg',[modelType,'_',referenceConditionName,'_',targetFluxes{ix,1},'_boxplot.svg']))
saveas(barplotFigure,fullfile('..','io','figures','simulations_barplot','mat',[modelType,'_',referenceConditionName,'_',targetFluxes{ix,1},'_boxplot.fig'])) 



