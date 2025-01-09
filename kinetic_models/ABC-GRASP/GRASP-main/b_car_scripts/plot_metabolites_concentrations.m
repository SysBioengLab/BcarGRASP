%% Plot metabolite concentration distributions, prior and posterior
% This scripts creates the boxplots of metabolite concentration
% distributions. Each figure contains an experimental condition, it can be
% a reference condition or a non-reference condition. Every plot includes
% the prior and posterior distributions

% It requires the files with the posterior models located in
% 'GRASP-main/io/output', and the prior models located in
% 'GRASP-main/io/unfiltered_output'. Other variables are defined to include
% nomenclatures in their figures and filenames

% The output are the metabolite concentrations boxplots. They are stored in
% 'GRASP-main/io/figures/metabolites_concentrations_boxplots', and their
% filenames indicate the condition and the type of model structure. The
% black boxplots represent the prior, while blue boxplots correspond to the
% posterior.

%% Clean variables and state matlab standards

clear, close all
rng('default');                 % for reproducibility
format longE
digits(256)                     % to prevent errors in calculations

%% Input

filteredFileDirs = {'../io/output/b_car_rejection_model_simple_D010';
    '../io/output/b_car_rejection_model_regulated_D010';
    '../io/output/b_car_rejection_model_detailed_D010';
    '../io/output/b_car_rejection_model_simple_D025';
    '../io/output/b_car_rejection_model_regulated_D025';
    '../io/output/b_car_rejection_model_detailed_D025'};                            % filenames of the posterior models

unfilteredFileDirs = {'../io/unfiltered_output/b_car_rejection_model_simple_D010';
    '../io/unfiltered_output/b_car_rejection_model_regulated_D010';
    '../io/unfiltered_output/b_car_rejection_model_detailed_D010';
    '../io/unfiltered_output/b_car_rejection_model_simple_D025';
    '../io/unfiltered_output/b_car_rejection_model_regulated_D025';
    '../io/unfiltered_output/b_car_rejection_model_detailed_D025'};|                % filenames of the prior models

modelTypes = {'simple';'regulated';'detailed';'simple';'regulated';'detailed'};     % structure types
conditionNames = {'D010';'D010';'D010';'D025';'D025';'D025'};                       % growth rates
conditionFullNames = {'simple \mu = 0.101h^{-1}';'regulated \mu = 0.101h^{-1}';
    'detailed \mu = 0.101h^{-1}';'simple \mu = 0.254h^{-1}';
    'regulated \mu = 0.254h^{-1}';'detailed \mu = 0.254h^{-1}'};                    % full names of the conditions
referenceName = {'4'};                                                              % reference strain
adjustedNames = {'3';'2'};                                                          % non-reference strains

balancedN = 17;                                                                     % number of balanced metabolites

%% Extract Metabolite Concentrations

metabolitesFiltered = cell(length(filteredFileDirs),1);
metabolitesUnfiltered = cell(length(unfilteredFileDirs),1);
metaboliteNames = cell(6,1);
metabolitesRelativeFiltered = cell(length(filteredFileDirs),length(adjustedNames));
metabolitesRelativeUnfiltered = cell(length(unfilteredFileDirs),length(adjustedNames));

for ix = 1:length(filteredFileDirs)
    filteredFile = load(filteredFileDirs{ix,1});
    unfilteredFile = load(unfilteredFileDirs{ix,1});
    filteredN = length(filteredFile.ensembleFiltered.populations.models);
    unfilteredN = length(unfilteredFile.ensembleUnfiltered.populations.models);
    for jx = 1:filteredN
        metabolitesFiltered{ix,1}(jx,:) =  filteredFile.ensembleFiltered.populations.models(jx).metConcRef';
    end
    for jx = 1:unfilteredN
        metabolitesUnfiltered{ix,1}(jx,:) =  unfilteredFile.ensembleUnfiltered.populations.models(jx).metConcRef';
    end
    metaboliteNames{ix,1} = filteredFile.ensembleFiltered.mets;
    for jx = 1:length(metaboliteNames{ix,1})
        rxnTemp = strsplit(metaboliteNames{ix,1}{jx},'m_');
        metaboliteNames{ix,1}{jx} = rxnTemp{2};
    end
    for jx = 1:filteredN
        for kx = 1:length(adjustedNames)
            metabolitesRelativeFiltered{ix,kx}(jx,:) = filteredFile.ensembleFiltered.populations.xopt{jx,1}(1:length(metaboliteNames{ix,1}),kx)';
        end
    end
    for jx = 1:unfilteredN
        for kx = 1:length(adjustedNames)
            metabolitesRelativeUnfiltered{ix,kx}(jx,:) = unfilteredFile.ensembleUnfiltered.populations.xopt{jx,1}(1:length(metaboliteNames{ix,1}),kx)';
        end
    end
end

%% Plot Boxplots

for ix = 1:length(metabolitesUnfiltered)
    boxplotFigure = figure();
    boxplotFigure.Position = [500 200 540 300];
    boxplot(log10(metabolitesUnfiltered{ix,1}),'Symbol','k.', 'Colors','k','PlotStyle','traditional','MedianStyle','line', ...
        'Labels',metaboliteNames{ix,1}(1:end),'Widths',0.5,'Jitter',0.25,'OutlierSize',0.5,'Positions',(1:length(metaboliteNames{ix,1}))*3-2);
    hold on
    boxplot(log10(metabolitesFiltered{ix,1}),'Symbol','b.', 'Colors','b','PlotStyle','traditional','MedianStyle','line', ...
        'Labels',metaboliteNames{ix,1}(1:end),'Widths',0.5,'Jitter',0.5,'OutlierSize',0.5,'Positions',(1:length(metaboliteNames{ix,1}))*3-1);
    title({'Boxplot Metabolite Concentrations',['\beta-car',referenceName{1,1},' ',conditionFullNames{ix,1}]})
    ylabel('Log_{10} Concentrations [M]')
    xlabel('Metabolites')
    xticks((1:length(metaboliteNames{ix,1}))*3-1.5)
    grid('on')
    set(gca,'units','points','position',[50,50,350,125])
    set(gca,'TickDir','out')
    set(gca,'LineWidth',1)
    set(gca,'FontSize',8,'FontName','arial')
    set(gca,'TickLabelInterpreter','none')
    hold off
    saveas(boxplotFigure,fullfile('..','io','figures','metabolites_concentrations_boxplot','png',[referenceName{1,1},conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.png']))
    saveas(boxplotFigure,fullfile('..','io','figures','metabolites_concentrations_boxplot','jpg',[referenceName{1,1},conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.jpg']))
    saveas(boxplotFigure,fullfile('..','io','figures','metabolites_concentrations_boxplot','svg',[referenceName{1,1},conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.svg']))
    saveas(boxplotFigure,fullfile('..','io','figures','metabolites_concentrations_boxplot','mat',[referenceName{1,1},conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.fig']))
end

for jx = 1:length(adjustedNames)
    for ix = 1:length(metabolitesUnfiltered)
        boxplotFigure = figure();
        boxplotFigure.Position = [500 200 540 300];
        boxplot(log10(metabolitesUnfiltered{ix,1}.*metabolitesRelativeUnfiltered{ix,jx}),'Symbol','k.', 'Colors','k','PlotStyle','traditional','MedianStyle','line', ...
            'Labels',metaboliteNames{ix,1}(1:end),'Widths',0.5,'Jitter',0.25,'OutlierSize',0.5,'Positions',(1:length(metaboliteNames{ix,1}))*3-2);
        hold on
        boxplot(log10(metabolitesFiltered{ix,1}.*metabolitesRelativeFiltered{ix,jx}),'Symbol','b.', 'Colors','b','PlotStyle','traditional','MedianStyle','line', ...
            'Labels',metaboliteNames{ix,1}(1:end),'Widths',0.5,'Jitter',0.5,'OutlierSize',0.5,'Positions',(1:length(metaboliteNames{ix,1}))*3-1);
        title({'Boxplot Metabolite Concentrations',['\beta-car',adjustedNames{jx,1},' ',conditionFullNames{ix,1}]})
        ylabel('Log_{10} Concentrations [M]')
        xlabel('Metabolites')
        xticks((1:length(metaboliteNames{ix,1}))*3-1.5)
        grid('on')
        set(gca,'units','points','position',[50,50,350,125])
        set(gca,'TickDir','out')
        set(gca,'LineWidth',1)
        set(gca,'FontSize',8,'FontName','arial')
        set(gca,'TickLabelInterpreter','none')
        hold off
        saveas(boxplotFigure,fullfile('..','io','figures','metabolites_concentrations_boxplot','png',[adjustedNames{jx,1},conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.png']))
        saveas(boxplotFigure,fullfile('..','io','figures','metabolites_concentrations_boxplot','jpg',[adjustedNames{jx,1},conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.jpg']))
        saveas(boxplotFigure,fullfile('..','io','figures','metabolites_concentrations_boxplot','svg',[adjustedNames{jx,1},conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.svg']))
        saveas(boxplotFigure,fullfile('..','io','figures','metabolites_concentrations_boxplot','mat',[adjustedNames{jx,1},conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.fig']))
    end
end





