%% Plot relative metabolite concetration distributions, simulation results
% This script creates the boxplots of relative metabolite concentration
% distributions. Each figure contains the results of a simulation
% performed. The 'Input' section must be modified according to the
% simulation that wants to be used for the plot.

% It requires a file with the performed simulation. The location of the
% filefolder can be specified. Other variables are defined to include
% nomenclatures in their figures and filenames.

% The output are the relative metabolite concentration boxplots. They are stored in
% 'GRASP-main/io/figures/relative_metabolites_simulation_boxplots', and their
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

simulationN = 8;  

%% Extract Metabolite Concentrations

metaboliteNames = cell(length(simulationFileDirs),1);
metabolitesRelativeSimulation = cell(length(simulationFileDirs),1);

for ix = 1:length(simulationFileDirs)
    simulationFile = load(simulationFileDirs{ix,1});
    metaboliteNames{ix,1} = simulationFile.ensemble.mets;
    for jx = 1:length(metaboliteNames{ix,1})
        rxnTemp = strsplit(metaboliteNames{ix,1}{jx},'m_');
        metaboliteNames{ix,1}{jx} = rxnTemp{2};
    end
    metabolitesRelativeSimulation{ix,1} = simulationFile.xopts{simulationN,1}(:,1:length(metaboliteNames{ix,1}));
end

%% Plot Boxplots


for ix = 1:length(metabolitesRelativeSimulation)
    boxplotFigure = figure();
    boxplotFigure.Position = [500 200 300 250];
    boxplot(log2(metabolitesRelativeSimulation{ix,1}),'Symbol','k.', 'Colors','k','PlotStyle','traditional','MedianStyle','line', ...
        'Labels',metaboliteNames{ix,1}(1:end),'Widths',0.5,'Jitter',0.5,'OutlierSize',2,'Positions',(1:length(metaboliteNames{ix,1})));
    title({'Boxplot relative metabolite concentration of best simulation',conditionFullNames{ix,1}})
    ylabel('Log_{2} relative concentrations')
    xlabel('Metabolites')
    xticks((1:length(metaboliteNames{ix,1})))
    grid('on')
    set(gca,'units','points','position',[30,50,180,90])
    set(gca,'TickDir','out')
    set(gca,'LineWidth',1)
    set(gca,'FontSize',8,'FontName','arial')
    set(gca,'TickLabelInterpreter','none')
    hold off
    saveas(boxplotFigure,fullfile('..','io','figures','relative_metabolites_simulation_boxplot','png',[conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.png']))
    saveas(boxplotFigure,fullfile('..','io','figures','relative_metabolites_simulation_boxplot','jpg',[conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.jpg']))
    saveas(boxplotFigure,fullfile('..','io','figures','relative_metabolites_simulation_boxplot','svg',[conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.svg']))
    saveas(boxplotFigure,fullfile('..','io','figures','relative_metabolites_simulation_boxplot','mat',[conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.fig']))
end






