%% Plot relative enzyme abundances distributions, posterior models
% This script creates the boxplots of relative enzyme abundances
% distributions. Each figure contains a non-reference experimental 
% condition, since reference conditions have a relative enzyme value of 1
% Every plot includes only the posterior distributions.

% It requires the files with the posterior models located in
% 'GRASP-main/io/output'. Other variables are defined to include
% nomenclatures in their figures and filenames

% The output are the relative enzyme abundance boxplots. They are stored in
% 'GRASP-main/io/figures/relative_enzymes_posterior_boxplots', and their
% filenames indicate the condition and the type of model structure.

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
    '../io/output/b_car_rejection_model_detailed_D025'};                            % filenames of posterior models

modelTypes = {'simple';'regulated';'detailed';'simple';'regulated';'detailed'};     % structure type
conditionNames = {'D010';'D010';'D010';'D025';'D025';'D025'};                       % growth rates
conditionFullNames = {'simple \mu = 0.101h^{-1}';'regulated \mu = 0.101h^{-1}';
    'detailed \mu = 0.101h^{-1}';'simple \mu = 0.254h^{-1}';
    'regulated \mu = 0.254h^{-1}';'detailed \mu = 0.254h^{-1}'};                    % full names of conditions
adjustedNames = {'3';'2'};                                                          % non-reference strains


%% Extract Relative Enzyme Abundances

enzymeNames = cell(6,1);
enzymesRelativeFiltered = cell(length(filteredFileDirs),length(adjustedNames));

for ix = 1:length(filteredFileDirs)
    filteredFile = load(filteredFileDirs{ix,1});
    filteredN = length(filteredFile.ensembleFiltered.populations.models);
    enzymeNames{ix,1} = filteredFile.ensembleFiltered.rxns;
    enzymeNames{ix,1}([21,22],:) = [];
    keepEnzymes = 1:length(enzymeNames{ix,1});
    metsN = length(filteredFile.ensembleFiltered.mets);
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
    for jx = 1:filteredN
        for kx = 1:length(adjustedNames)
            enzymesRelativeFiltered{ix,kx}(jx,:) = filteredFile.ensembleFiltered.populations.xopt{jx,1}(keepEnzymes+metsN,kx)';
        end
    end
end

%% Plot Boxplots

for jx = 1:length(adjustedNames)
    for ix = 1:length(enzymesRelativeFiltered)
        boxplotFigure = figure();
        boxplotFigure.Position = [500 200 300 250];
        boxplot(log2(enzymesRelativeFiltered{ix,jx}),'Symbol','k.', 'Colors','k','PlotStyle','traditional','MedianStyle','line', ...
            'Labels',enzymeNames{ix,1}(1:end),'Widths',0.5,'Jitter',0.5,'OutlierSize',2,'Positions',(1:length(enzymeNames{ix,1})));
        title({'Boxplot relative enzymes abundances',['\beta-car',adjustedNames{jx,1},' ',conditionFullNames{ix,1}]})
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
        saveas(boxplotFigure,fullfile('..','io','figures','relative_enzymes_posterior_boxplot','png',[adjustedNames{jx,1},conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.png']))
        saveas(boxplotFigure,fullfile('..','io','figures','relative_enzymes_posterior_boxplot','jpg',[adjustedNames{jx,1},conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.jpg']))
        saveas(boxplotFigure,fullfile('..','io','figures','relative_enzymes_posterior_boxplot','svg',[adjustedNames{jx,1},conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.svg']))
        saveas(boxplotFigure,fullfile('..','io','figures','relative_enzymes_posterior_boxplot','mat',[adjustedNames{jx,1},conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.fig']))
    end
end





