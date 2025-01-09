%% This script creates the mean control heatmaps of the different conditions
% It only creates the heatmaps for one condition at a time, the condition
% can be selected by changing the 'i' variable at the 'Load data' section
% to match the inputs of the desired condition

% It requires that the CR analysis of the ensemble has been performed
% previously, ant that the data is stored in GRASP-main/io/cr with the name
% of the model ID with '_crResults' added to the end of the filename

% The heatmaps are stored in GRASP-main/io/figures/control_mean
% The name of of the figures is based on the input variable 'fileName'. The
% ending of the output figures and their meaning are:
%   _crMetabolites: Response coefficients of metabolite concentrations
%   _crFluxes: Response coefficients of fluxes
%   _mcaMetabolites: Concentration control coefficients
%   _mcaFluces: Flux control coefficients

% The mean coefficients of each figure are stored in an excel file, where
% each sheet is the value of the coefficients for a specific condition

% Unlike 'plot_control_mean', in the heatmaps only the fluxes and
% concentrations associated with lycopene and beta-carotene are displayed.
% The names of the enzymes, reactions, fluxes and metabolites are included

%% Clean variables and state matlab standards

clear, close all
rng('default');                 % for reproducibility
format longE
digits(256)

%% Path
addpath(fullfile('..', 'matlab_code', 'analysisFxns'), ...
        fullfile('..', 'matlab_code', 'ensembleFxns'), ...
        fullfile('..', 'matlab_code', 'patternFxns'));

%% Input

modelIDs = {'b_car_rejection_model_detailed_D010';'b_car_rejection_model_detailed_D025';...
    'b_car_rejection_model_regulated_D010';'b_car_rejection_model_regulated_D025';...
    'b_car_rejection_model_simple_D010';'b_car_rejection_model_simple_D025'};                   % IDs of models

fileNames = {'detailed_4D010';'detailed_4D025';...
    'regulated_4D010';'regulated_4D025';...
    'simple_4D010';'simple_4D025'};                                                             % filenames of the outputs

modelTypes = {'detailed';'detailed';...
    'regulated';'regulated';...
    'simple';'simple'};                                                                         % model structures types

conditionNames = {'4D010';'4D025';'4D010';'4D025';'4D010';'4D025'};                             % names of the conditions

conditionFullNames = {'\beta-car4 \mu = 0.101h^{-1}';'\beta-car4 \mu = 0.254h^{-1}';...
    '\beta-car4 \mu = 0.101h^{-1}';'\beta-car4 \mu = 0.254h^{-1}';...
    '\beta-car4 \mu = 0.101h^{-1}';'\beta-car4 \mu = 0.254h^{-1}'};                             % full names of the conditions

fluxSelections = {[21,22];[21,22];[21,22];[21,22];[21,22];[21,22];[21,22]};                     % fluxes to be displayed for each condition
metsSelections = {[15,17];[15,17];[15,17];[15,17];[15,17];[15,17];[15,17]};                     % concentrations to be displayed for each condition

%% Load Data and Create Heatmaps
% Modify variable 'i' to choose the ensemble to be analyzed

load('cmap_rgb.mat')
cmap = cmap;

i = 1;
modelID = modelIDs{i,1};
modelType = modelTypes{i,1};
conditionName = conditionNames{i,1};
conditionFullName = conditionFullNames{i,1};
fileName = fileNames{i,1};
fluxSelection = fluxSelections{i,1};
metsSelection = metsSelections{i,1};
load(fullfile('..', 'io', 'output', [modelID, '.mat']));
load(fullfile('..','io','cr',[modelID,'_crResults.mat']));
ensemble = ensembleFiltered;

%% Create Main Variables

% Optimization & simulation parameters
freeVarsN = numel(ensemble.freeVars);

metNames = ensemble.mets(ensemble.metsBalanced);
fluxNames = ensemble.rxns;
enzNames = mcaResults.enzNames;

fluxN = numel(ensemble.fluxRef);
metN = length(metNames);
enzN = length(enzNames);

metIndexes = 1:numel(ensemble.metsBalanced);

% Plot parameters
cLimitsFluxes = [-2 2];
cLimitsMetabolites = [-2 2];
cValuesFluxes = linspace(cLimitsFluxes(:,1),cLimitsFluxes(:,2),5);
cValuesMetabolites = linspace(cLimitsMetabolites(:,1),cLimitsMetabolites(:,2),5);
cCellFluxes = num2cell(cValuesFluxes);
cCellMetabolites = num2cell(cValuesMetabolites);
% cCellFluxes{1} = ['<',num2str(cLimitsFluxes(:,1))];
% cCellFluxes{end} = ['>',num2str(cLimitsFluxes(:,2))];
% cCellMetabolites{1} = ['<',num2str(cLimitsMetabolites(:,1))];
% cCellMetabolites{end} = ['>',num2str(cLimitsMetabolites(:,2))];

% Remove prefixes
for i= 1:length(metNames)
    metTemp = strsplit(metNames{i},'m_');
    metNames{i} = metTemp{2};
end
for i = 1:length(fluxNames)
    rxnTemp = strsplit(fluxNames{i},'r_');
    fluxNames{i} = rxnTemp{2};
end
for i = 1:length(enzNames)
    rxnTemp = strsplit(enzNames{i},'r_');
    enzNames{i} = rxnTemp{2};
end

%% Plot Fluxes vs Fluxes

mcavAll = mcaResults.vControl{1,1};
modelsN = size(mcavAll,1)/fluxN;

mcav = zeros(fluxN,fluxN);
mcavFluxes = cell(fluxN,1);

for i = 1:fluxN
    for j = 1:modelsN
        mcavFluxes{i,1}(j,:) = mcavAll((j-1)*fluxN+i,:);
    end
    mcav(i,:) =  mean(mcavFluxes{i,1},1);
end

CrtBacol = mcav(:,20);
CrtBbcol = mcav(:,19);
CrtYacol = mcav(:,17);
CrtYbcol = mcav(:,18);
mcav(:,17) = CrtBacol;
mcav(:,18) = CrtBbcol;
mcav(:,19) = CrtYacol;
mcav(:,20) = CrtYbcol;

CrtBarow = mcav(20,:);
CrtBbrow = mcav(19,:);
CrtYarow = mcav(17,:);
CrtYbrow = mcav(18,:);
mcav(17,:) = CrtBarow;
mcav(18,:) = CrtBbrow;
mcav(19,:) = CrtYarow;
mcav(20,:) = CrtYbrow;

fluxNames{17,:} = 'CrtBa';
fluxNames{18,:} = 'CrtBb';
fluxNames{19,:} = 'CrtYa';
fluxNames{20,:} = 'CrtYb';

nx = fluxN-2;
ny = length(fluxSelection);
mx = linspace(0.5,nx+0.5,nx+1);
my = linspace(0.5,ny+0.5,ny+1);

fluxFigure = figure('Name',conditionName);
fluxFigure.Position = [500 300 250 110];
imagesc(mcav(fluxSelection,1:(end-2))); hold on
set(gca,'DefaultAxesTickLabelInterpreter','none')
set(gca,'TickLabelInterpreter','none')
set(gca,'xticklabel',[],'yticklabel',[],'xtick',1:(fluxN-2),'ytick',1:(length(fluxSelection)))
set(gca,'xticklabel',fluxNames(1:(end-2)),'yticklabel',fluxNames(fluxSelection))
set(gca,'units','points','position',[50,40,125,17.857])
set(gca,'TickDir','out')
set(gca,'LineWidth',0.5)
set(gca,'XTickLabelRotation',90)
xlabel('Fluxes')
ylabel('Fluxes')
%title(['Flux control coefficients: ',conditionFullName])
set(gca,'FontSize',7,'FontName','arial')
caxis(cLimitsFluxes)
colormap(fluxFigure,cmap)
cb = colorbar;
set(cb,'Location','northoutside')
set(cb,'LineWidth',0.5)
set(cb,'TickDir','out')
set(cb,'FontSize',7,'FontName','arial')
set(cb,'Ticks',cValuesFluxes)
set(cb,'TickLabels',cCellFluxes)
cbparameters = cb.Position;
set(cb,'Position',[cbparameters(1),0.75,cbparameters(3),0.05])
mg = mesh(mx,my,zeros([ny,nx]+1));
mg.FaceColor = 'none';
mg.EdgeColor = 'k';
hold off

outputFluxFig = fullfile('..','io','figures','control_mean_selection','mat',[fileName,'_mcaFluxes_selection.fig']);
savefig(fluxFigure,outputFluxFig)
outputFluxJpg = fullfile('..','io','figures','control_mean_selection','jpg',[fileName,'_mcaFluxes_selection.jpg']);
saveas(fluxFigure,outputFluxJpg)
outputFluxSvg = fullfile('..','io','figures','control_mean_selection','svg',[fileName,'_mcaFluxes_selection.svg']);
saveas(fluxFigure,outputFluxSvg)
outputFluxPng = fullfile('..','io','figures','control_mean_selection','png',[fileName,'_mcaFluxes_selection.png']);
saveas(fluxFigure,outputFluxPng)

%% Plot Metabolites vs. Fluxes

mcaxAll = mcaResults.xControl{1,1};
modelsN = size(mcaxAll,1)/metN;

mcax = zeros(metN,fluxN);
mcaxMetabolites = cell(metN,1);

for i = 1:metN
    for j = 1:modelsN
        mcaxMetabolites{i,1}(j,:) = mcaxAll((j-1)*metN+i,:);
    end
    mcax(i,:) =  mean(mcaxMetabolites{i,1},1);
end

CrtBamet = mcax(:,20);
CrtBbmet = mcax(:,19);
CrtYamet = mcax(:,17);
CrtYbmet = mcax(:,18);
mcax(:,17) = CrtBamet;
mcax(:,18) = CrtBbmet;
mcax(:,19) = CrtYamet;
mcax(:,20) = CrtYbmet;

nx = fluxN-2;
ny = length(metsSelection);
mx = linspace(0.5,nx+0.5,nx+1);
my = linspace(0.5,ny+0.5,ny+1);

metFigure = figure('Name',conditionName);
metFigure.Position = [500 300 250 110];
imagesc(mcax(metsSelection,1:(end-2))); hold on
set(gca,'DefaultAxesTickLabelInterpreter', 'none')
set(gca,'TickLabelInterpreter','none')
set(gca,'xticklabel',[],'yticklabel',[],'ytick',1:length(metsSelection),'xtick',1:(fluxN-2))
set(gca,'yticklabel',metNames(metsSelection),'xticklabel',fluxNames(1:(end-2)))
set(gca,'units','points','position',[50,40,125,17.857])
set(gca,'TickDir','out')
set(gca,'LineWidth',0.5)
set(gca,'XTickLabelRotation',90)
ylabel('Metabolites')
xlabel('Fluxes')
%title(['Concentration control coefficients: ',conditionFullName])
set(gca,'FontSize',7,'FontName','arial')
caxis(cLimitsMetabolites)
colormap(metFigure,cmap)
cb = colorbar;
set(cb,'Location','northoutside')
set(cb,'LineWidth',0.5)
set(cb,'TickDir','out')
set(cb,'FontSize',7,'FontName','arial')
set(cb,'Ticks',cValuesMetabolites)
set(cb,'TickLabels',cCellMetabolites)
cbparameters = cb.Position;
set(cb,'Position',[cbparameters(1),0.75,cbparameters(3),0.05])
mg = mesh(mx,my,zeros([ny,nx]+1));
mg.FaceColor = 'none';
mg.EdgeColor = 'k';
hold off

outputMetFig = fullfile('..','io','figures','control_mean_selection','mat',[fileName,'_mcaMetabolites_selection.fig']);
savefig(metFigure,outputMetFig)
outputMetJpg = fullfile('..','io','figures','control_mean_selection','jpg',[fileName,'_mcaMetabolites_selection.jpg']);
saveas(metFigure,outputMetJpg)
outputMetSvg = fullfile('..','io','figures','control_mean_selection','svg',[fileName,'_mcaMetabolites_selection.svg']);
saveas(metFigure,outputMetSvg)
outputMetPng = fullfile('..','io','figures','control_mean_selection','png',[fileName,'_mcaMetabolites_selection.png']);
saveas(metFigure,outputMetPng)


%% Plot Fluxes vs Enzymes

mcaevAll = mcaResults.eResponse{1,1};
modelsN = size(mcaevAll,1)/fluxN;

mcaev = zeros(fluxN,enzN);
mcaevFluxes = cell(fluxN,1);

for i = 1:fluxN
    for j = 1:modelsN
        mcaevFluxes{i,1}(j,:) = mcaevAll((j-1)*fluxN+i,:);
    end
    mcaev(i,:) =  mean(mcaevFluxes{i,1},1);
end

CrtBaerow = mcaev(20,:);
CrtBberow = mcaev(19,:);
CrtYaerow = mcaev(17,:);
CrtYberow = mcaev(18,:);
mcaev(17,:) = CrtBaerow;
mcaev(18,:) = CrtBberow;
mcaev(19,:) = CrtYaerow;
mcaev(20,:) = CrtYberow;


if strcmp('simple',modelType) || strcmp('regulated',modelType)

    CrtBaecol = mcaev(:,20);
    CrtBbecol = mcaev(:,19);
    CrtYaecol = mcaev(:,17);
    CrtYbecol = mcaev(:,18);
    mcaev(:,17) = CrtBaecol;
    mcaev(:,18) = CrtBbecol;
    mcaev(:,19) = CrtYaecol;
    mcaev(:,20) = CrtYbecol;

    enzNames{17,:} = 'CrtBa';
    enzNames{18,:} = 'CrtBb';
    enzNames{19,:} = 'CrtYa';
    enzNames{20,:} = 'CrtYb';

end

if strcmp('detailed',modelType)

   enzNames{9,:} = 'ERG20';
   enzNames{12,:} = 'ERG9';
   enzNames{13,:} = 'CrtI';
   enzNames{14,:} = 'CrtYB';

end

nx = size(mcaev,2)-2;
ny = length(fluxSelection);
mx = linspace(0.5,nx+0.5,nx+1);
my = linspace(0.5,ny+0.5,ny+1);

efluxFigure = figure('Name',conditionName);
efluxFigure.Position = [500 300 250 110];
imagesc(mcaev(fluxSelection,1:(end-2))); hold on
set(gca,'DefaultAxesTickLabelInterpreter', 'none')
set(gca,'TickLabelInterpreter','none')
set(gca,'xticklabel',[],'yticklabel',[],'xtick',1:(enzN-2),'ytick',1:length(fluxSelection))
set(gca,'xticklabel',enzNames(1:(end-2)),'yticklabel',fluxNames(fluxSelection))
set(gca,'units','points','position',[50,40,125,17.857])
set(gca,'TickDir','out')
set(gca,'LineWidth',0.5)
set(gca,'XTickLabelRotation',90)
ylabel('Fluxes')
xlabel('Enzymes')
%title(['Enzyme flux control coefficients: ',conditionFullName])
set(gca,'FontSize',7,'FontName','arial')
caxis(cLimitsFluxes)
colormap(efluxFigure,cmap)
cb = colorbar;
set(cb,'Location','northoutside')
set(cb,'LineWidth',0.5)
set(cb,'TickDir','out')
set(cb,'FontSize',7,'FontName','arial')
set(cb,'Ticks',cValuesFluxes)
set(cb,'TickLabels',cCellFluxes)
cbparameters = cb.Position;
set(cb,'Position',[cbparameters(1),0.75,cbparameters(3),0.05])
mg = mesh(mx,my,zeros([ny,nx]+1));
mg.FaceColor = 'none';
mg.EdgeColor = 'k';
hold off

outputeFluxFig = fullfile('..','io','figures','control_mean_selection','mat',[fileName,'_crFluxes_selection.fig']);
savefig(efluxFigure,outputeFluxFig)
outputeFluxJpg = fullfile('..','io','figures','control_mean_selection','jpg',[fileName,'_crFluxes_selection.jpg']);
saveas(efluxFigure,outputeFluxJpg)
outputeFluxSvg = fullfile('..','io','figures','control_mean_selection','svg',[fileName,'_crFluxes_selection.svg']);
saveas(efluxFigure,outputeFluxSvg)
outputeFluxPng = fullfile('..','io','figures','control_mean_selection','png',[fileName,'_crFluxes_selection.png']);
saveas(efluxFigure,outputeFluxPng)

%% Plot Metabolites vs Enzymes

mcaexAll = mcaResults.xResponse{1,1};
modelsN = size(mcaexAll,1)/metN;

mcaex = zeros(metN,enzN);
mcaexMetabolites = cell(metN,1);

for i = 1:metN
    for j = 1:modelsN
        mcaexMetabolites{i,1}(j,:) = mcaexAll((j-1)*metN+i,:);
    end
    mcaex(i,:) =  mean(mcaexMetabolites{i,1},1);
end

if strcmp('simple',modelType) || strcmp('regulated',modelType)
    CrtBaemet = mcaex(:,20);
    CrtBbemet = mcaex(:,19);
    CrtYaemet = mcaex(:,17);
    CrtYbemet = mcaex(:,18);
    mcaex(:,17) = CrtBaemet;
    mcaex(:,18) = CrtBbemet;
    mcaex(:,19) = CrtYaemet;
    mcaex(:,20) = CrtYbemet;
end 

nx = size(mcaex,2)-2;
ny = length(metsSelection);
mx = linspace(0.5,nx+0.5,nx+1);
my = linspace(0.5,ny+0.5,ny+1);

emetFigure = figure('Name',conditionName);
emetFigure.Position = [500 300 250 110];
imagesc(mcaex(metsSelection,1:(enzN-2))); hold on
set(gca, 'DefaultAxesTickLabelInterpreter', 'none')
set(gca,'TickLabelInterpreter','none')
set(gca,'xticklabel',[],'yticklabel',[],'ytick',1:length(metsSelection),'xtick',1:(enzN-2))
set(gca,'yticklabel',metNames(metsSelection),'xticklabel',enzNames(1:(end-2)))
set(gca,'units','points','position',[50,40,125,17.857])
set(gca,'TickDir','out')
set(gca,'LineWidth',0.5)
set(gca,'XTickLabelRotation',90)
ylabel('Metabolites')
xlabel('Enzymes')
%title(['Enzyme concentration control coefficients: ',conditionFullName])
set(gca,'FontSize',7,'FontName','arial')
caxis(cLimitsMetabolites)
colormap(emetFigure,cmap)
cb = colorbar;
set(cb,'Location','northoutside')
set(cb,'LineWidth',0.5)
set(cb,'TickDir','out')
set(cb,'FontSize',7,'FontName','arial')
set(cb,'Ticks',cValuesMetabolites)
set(cb,'TickLabels',cCellMetabolites)
cbparameters = cb.Position;
set(cb,'Position',[cbparameters(1),0.75,cbparameters(3),0.05])
mg = mesh(mx,my,zeros([ny,nx]+1));
mg.FaceColor = 'none';
mg.EdgeColor = 'k';
hold off

outputeMetFig = fullfile('..','io','figures','control_mean_selection','mat',[fileName,'_crMetabolites_selection.fig']);
savefig(emetFigure,outputeMetFig)
outputeMetJpg = fullfile('..','io','figures','control_mean_selection','jpg',[fileName,'_crMetabolites_selection.jpg']);
saveas(emetFigure,outputeMetJpg)
outputeMetSvg = fullfile('..','io','figures','control_mean_selection','svg',[fileName,'_crMetabolites_selection.svg']);
saveas(emetFigure,outputeMetSvg)
outputeMetPng = fullfile('..','io','figures','control_mean_selection','png',[fileName,'_crMetabolites_selection.png']);
saveas(emetFigure,outputeMetPng)

%% Save table with values

% eFluxTable1 = array2table(mcaev,"VariableNames",enzNames(:));
% eFluxTable2 = cell2table(fluxNames,"VariableNames","Flux");
% eFluxTable = [eFluxTable2,eFluxTable1];
% OutputeFluxXls = fullfile('..','io','figures','control_mean_selection','xls','crFluxes.xlsx');
% writetable(eFluxTable,OutputeFluxXls,'Sheet',fileName);
% 
% eMetTable1 = array2table(mcaex,"VariableNames",enzNames(:));
% eMetTable2 = cell2table(metNames,"VariableNames","Flux");
% eMetTable = [eMetTable2,eMetTable1];
% OutputeMetXls = fullfile('..','io','figures','control_mean_selection','xls','crMetabolites.xlsx');
% writetable(eMetTable,OutputeMetXls,'Sheet',fileName);




