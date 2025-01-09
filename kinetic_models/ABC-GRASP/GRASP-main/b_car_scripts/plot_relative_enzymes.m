%% Plot relative enzyme abundances distributions, prior and posterior
% This script creates the boxplots of relative enzyme abundances
% distributions. Each figure contains a non-reference experimental 
% condition, since reference conditions have a relative enzyme value of 1
% Every plot includes the prior and posterior distributions.

% It requires the files with the posterior models located in
% 'GRASP-main/io/output', and the prior models located in
% 'GRASP-main/io/unfiltered_output'. Other variables are defined to include
% nomenclatures in their figures and filenames

% The output are the relative enzyme abundance boxplots. They are stored in
% 'GRASP-main/io/figures/relative_enzymes_boxplots', and their
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
    '../io/output/b_car_rejection_model_detailed_D025'};                                % filenames of posterior models

unfilteredFileDirs = {'../io/unfiltered_output/b_car_rejection_model_simple_D010';
    '../io/unfiltered_output/b_car_rejection_model_regulated_D010';
    '../io/unfiltered_output/b_car_rejection_model_detailed_D010';
    '../io/unfiltered_output/b_car_rejection_model_simple_D025';
    '../io/unfiltered_output/b_car_rejection_model_regulated_D025';
    '../io/unfiltered_output/b_car_rejection_model_detailed_D025'};                     % filenames of prior models

modelTypes = {'simple';'regulated';'detailed';'simple';'regulated';'detailed'};         % structure types
conditionNames = {'D010';'D010';'D010';'D025';'D025';'D025'};                           % growth rates
conditionFullNames = {'simple \mu = 0.101h^{-1}';'regulated \mu = 0.101h^{-1}';
    'detailed \mu = 0.101h^{-1}';'simple \mu = 0.254h^{-1}';
    'regulated \mu = 0.254h^{-1}';'detailed \mu = 0.254h^{-1}'};                        % full condition names
adjustedNames = {'3';'2'};                                                              % non-reference strains


%% Extract Relative Enzyme Abundances

enzymeNames = cell(6,1);
enzymesRelativeFiltered = cell(length(filteredFileDirs),length(adjustedNames));
enzymesRelativeUnfiltered = cell(length(unfilteredFileDirs),length(adjustedNames));

for ix = 1:length(filteredFileDirs)
    filteredFile = load(filteredFileDirs{ix,1});
    unfilteredFile = load(unfilteredFileDirs{ix,1});
    filteredN = length(filteredFile.ensembleFiltered.populations.models);
    unfilteredN = length(unfilteredFile.ensembleUnfiltered.populations.models);
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
    for jx = 1:unfilteredN
        for kx = 1:length(adjustedNames)
            enzymesRelativeUnfiltered{ix,kx}(jx,:) = unfilteredFile.ensembleUnfiltered.populations.xopt{jx,1}(keepEnzymes+metsN,kx)';
        end
    end
end

%% Plot Boxplots

for jx = 1:length(adjustedNames)
    for ix = 1:length(enzymesRelativeUnfiltered)
        boxplotFigure = figure();
        boxplotFigure.Position = [500 200 425 300];
        boxplot(log2(enzymesRelativeUnfiltered{ix,jx}),'Symbol','k.', 'Colors','k','PlotStyle','traditional','MedianStyle','line', ...
            'Labels',enzymeNames{ix,1}(1:end),'Widths',0.5,'Jitter',0.25,'OutlierSize',0.5,'Positions',(1:length(enzymeNames{ix,1}))*3-2);
        hold on
        boxplot(log2(enzymesRelativeFiltered{ix,jx}),'Symbol','b.', 'Colors','b','PlotStyle','traditional','MedianStyle','line', ...
            'Labels',enzymeNames{ix,1}(1:end),'Widths',0.5,'Jitter',0.5,'OutlierSize',0.5,'Positions',(1:length(enzymeNames{ix,1}))*3-1);
        title({'Boxplot Relative Enzymes Abundances',['\beta-car',adjustedNames{jx,1},' ',conditionFullNames{ix,1}]})
        ylabel('Log_{2} Relative Abundance')
        %ylabel('Relative Abundance')
        xlabel('Enzymes')
        xticks((1:length(enzymeNames{ix,1}))*3-1.5)
        grid('on')
        set(gca,'units','points','position',[50,50,250,125])
        set(gca,'TickDir','out')
        set(gca,'LineWidth',1)
        set(gca,'FontSize',8,'FontName','arial')
        set(gca,'TickLabelInterpreter','tex')
        hold off
        saveas(boxplotFigure,fullfile('..','io','figures','relative_enzymes_boxplot','png',[adjustedNames{jx,1},conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.png']))
        saveas(boxplotFigure,fullfile('..','io','figures','relative_enzymes_boxplot','jpg',[adjustedNames{jx,1},conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.jpg']))
        saveas(boxplotFigure,fullfile('..','io','figures','relative_enzymes_boxplot','svg',[adjustedNames{jx,1},conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.svg']))
        saveas(boxplotFigure,fullfile('..','io','figures','relative_enzymes_boxplot','mat',[adjustedNames{jx,1},conditionNames{ix,1},'_',modelTypes{ix,1},'_boxplot.fig']))
    end
end





