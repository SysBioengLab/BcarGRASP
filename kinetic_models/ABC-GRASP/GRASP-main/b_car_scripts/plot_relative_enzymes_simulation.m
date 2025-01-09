%% Plot relative enzyme abundances distributions, simulation results
% This script creates the boxplots of relative enzyme abundances
% distributions. Each figure contains the results of a simulation
% performed. The 'Input' section must be modified according to the
% simulation that wants to be used for the plot.

% It requires a file with the performed simulation. The location of the
% filefolder can be specified. Other variables are defined to include
% nomenclatures in their figures and filenames.

% The output are the relative enzyme abundance boxplots. They are stored in
% 'GRASP-main/io/figures/relative_enzymes_simulation_boxplots', and their
% filenames indicate the condition and the type of model structure.

%% Clean variables and state matlab standards

clear, close all
rng('default');                 % for reproducibility
format longE
digits(256)                     % to prevent errors in calculations

%% Input

simulationFileDirs = {'../../ensemble_simulation/enzymes_detailed_D010/backup_b_car_simulation_enzymes_model_detailed_D010';
    '../../ensemble_simulation/enzymes_detailed_D025/backup_b_car_simulation_enzymes_model_detailed_D025';};                    % filenames of the simulations

modelTypes = {'detailed';'detailed'};                                                                                           % structure types
conditionNames = {'D010';'D025'};                                                                                               % growth rates
conditionFullNames = {'detailed \mu = 0.101h^{-1}';'detailed \mu = 0.254h^{-1}'};                                               % full name of conditions
referenceName = {'4'};                                                                                                          % strain used as reference

simulationN = 8;                                                                                                                % number of the simulation to be plotted

%% Extract Relative Enzyme Abundances

enzymeNames = cell(6,1);
enzymesRelativeSimulation = cell(length(simulationFileDirs),1);

for ix = 1:length(simulationFileDirs)
    simulationFile = load(simulationFileDirs{ix,1});
    enzymeNames{ix,1} = simulationFile.ensemble.rxns;
    enzymeNames{ix,1}([21,22],:) = [];
    keepEnzymes = 1:length(enzymeNames{ix,1});
    metsN = length(simulationFile.ensemble.mets);
    for jx = 1:length(enzymeNames{ix,1})
        rxnTemp = strsplit(enzymeNames{ix,1}{jx},'r_');
        enzymeNames{ix,1}{jx} = rxnTemp{2};
    end
    if strcmp(modelTypes{ix,1},'detailed')
        keepEnzymes = [1,2,3,4,5,6,7,8,9,11,12,13,15,17];
        enzymeNames{ix,1}([10,14,16,18,19,20],:) = [];
        enzymeNames{ix,1}(9,:) = {'ERG20'};
        enzymeNames{ix,1}(12,:) = {'ERG9'};
        enzymeNames{ix,1}(13,:) = {'CrtI'};
        enzymeNames{ix,1}(14,:) = {'CrtYB'};
    end
    enzymesRelativeSimulation{ix,1} = simulationFile.xopts{simulationN,1}(:,keepEnzymes+metsN);
end

%% Plot Boxplots

for ix = 1:length(enzymesRelativeSimulation)
    boxplotFigure = figure();
    boxplotFigure.Position = [500 200 300 250];
    boxplot(log2(enzymesRelativeSimulation{ix,1}),'Symbol','k.', 'Colors','k','PlotStyle','traditional','MedianStyle','line', ...
        'Labels',enzymeNames{ix,1}(1:end),'Widths',0.5,'Jitter',0.5,'OutlierSize',2,'Positions',(1:length(enzymeNames{ix,1})));
    title({'Boxplot relative enzymes abundances of best simulation',conditionFullNames{ix,1}})
    ylabel('Log_{2} relative abundance')
    %ylabel('Relative Abundance')
    xlabel('Enzymes')
    xticks((1:length(enzymeNames{ix,1})))
    grid('on')
    set(gca,'units','points','position',[30,50,180,90])
    set(gca,'TickDir','out')
    set(gca,'LineWidth',1)
    set(gca,'FontSize',8,'FontName','arial')
    set(gca,'TickLabelInterpreter','none')
    hold off
    saveas(boxplotFigure,fullfile('..','io','figures','relative_enzymes_simulation_boxplot','png',[conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.png']))
    saveas(boxplotFigure,fullfile('..','io','figures','relative_enzymes_simulation_boxplot','jpg',[conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.jpg']))
    saveas(boxplotFigure,fullfile('..','io','figures','relative_enzymes_simulation_boxplot','svg',[conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.svg']))
    saveas(boxplotFigure,fullfile('..','io','figures','relative_enzymes_simulation_boxplot','mat',[conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.fig']))
end






