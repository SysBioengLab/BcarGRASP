%% Test of significance for the mean concentration response coefficients
% This script tests the significance of the mean concentration response
% coefficients. Each plot only includes one ensemble model. The Bonferroni
% correction is applied within the code. To choose the ensemble to be
% tested, the variable 'i' in 'Load data' must be modified. 

% It requires as inputs the folders for the posterior ensembles and the
% control and response outputs. It also needs other variables to define the
% nomenclatures of the plots and files generated. It is necessary to state
% the numbers of tests applied in total per ensemble model, also the
% significance to test. Although this last variable is not used for the
% plots.

% The output are heatmaps for the significance of the concentration response
% coefficients. These are stored in
% 'GRASP-main/io/figure/metabolite_control_significance_selection'. Their
% filenames indicate the type of sampling, reference condition, and
% structure type. An excel file is generated in the same folder with the
% values of the test

%% Clean variables and state matlab standards

clear, close all
rng('default');                 % for reproducibility
format longE
digits(256)                     % to prevent errors in calculation

%% Input


ensembleDirs = {'../io/output/b_car_rejection_model_simple_D010';
    '../io/output/b_car_rejection_model_simple_D025';
    '../io/output/b_car_rejection_model_regulated_D010';
    '../io/output/b_car_rejection_model_regulated_D025';
    '../io/output/b_car_rejection_model_detailed_D010';
    '../io/output/b_car_rejection_model_detailed_D025'};
controlDirs = {'../io/cr/b_car_rejection_model_simple_D010_crResults';                          % filenames of the ensembles
    '../io/cr/b_car_rejection_model_simple_D025_crResults';
    '../io/cr/b_car_rejection_model_regulated_D010_crResults';
    '../io/cr/b_car_rejection_model_regulated_D025_crResults';
    '../io/cr/b_car_rejection_model_detailed_D010_crResults';
    '../io/cr/b_car_rejection_model_detailed_D025_crResults'};                                  % filenames of the control and response results
modelTypes = {'simple';'simple';'regulated';'regulated';'detailed';'detailed'};                 % structure types
conditionNames = {'4D010';'4D025';'4D010';'4D025';'4D010';'4D025'};                             % reference conditions
conditionFullNames = {'\beta-car4 \mu = 0.101h^{-1}';'\beta-car4 \mu = 0.254h^{-1}';
    '\beta-car4 \mu = 0.101h^{-1}';'\beta-car4 \mu = 0.254h^{-1}';
    '\beta-car4 \mu = 0.101h^{-1}';'\beta-car4 \mu = 0.254h^{-1}'};                             % full names of conditions
samplingTypes = {'rejection';'rejection';'rejection';'rejection';'rejection';'rejection'};      % sampling type used
metsSelections = {[15,17];[15,17];[15,17];[15,17];[15,17];[15,17];[15,17]};                     % metabolites to be tested
testsNs = {80;80;80;80;56;56};                                                                  % number of tests applied per ensemble

significance = 0.05;                                                                            % significance for the single test in excel
nbootstrap = 1000000;                                                                           % number of boostrap resamplings

%% Load Data

i = 6;
modelType = modelTypes{i,1};
conditionName = conditionNames{i,1};
conditionFullName = conditionFullNames{i,1};
samplingType = samplingTypes{i,1};
metsSelection = metsSelections{i,1};

ensembleDir = ensembleDirs{i,1};
ensembleFile = load(ensembleDir);
ensemble = ensembleFile.ensembleFiltered;
testsN = testsNs{i,1};


controlDir = controlDirs{i,1};
controlFile = load(controlDir);
responseUnprocessed = controlFile.mcaResults.xResponse{1,1};

%% Create colormap for heatmaps

cmap = [[linspace(0,0.98,50)',linspace(0,0.98,50)',ones(50,1)];
    ones(1,3);
    [ones(50,1),linspace(0.98,0,50)',linspace(0.98,0,50)']];

%% Extract reaction names and arrange metabolites coefficients

metNames = ensemble.mets(ensemble.metsBalanced);
metN = length(metNames);
enzNames = controlFile.mcaResults.enzNames;
enzN = length(enzNames);

obsN = size(responseUnprocessed,1)/metN;

for i= 1:metN
    metTemp = strsplit(metNames{i},'m_');
    metNames{i} = metTemp{2};
end
for i= 1:enzN
    enzTemp = strsplit(enzNames{i},'r_');
    enzNames{i} = enzTemp{2};
end

if strcmp('simple',modelType) || strcmp('regulated',modelType)

    CrtBaecol = responseUnprocessed(:,20);
    CrtBbecol = responseUnprocessed(:,19);
    CrtYaecol = responseUnprocessed(:,17);
    CrtYbecol = responseUnprocessed(:,18);
    responseUnprocessed(:,17) = CrtBaecol;
    responseUnprocessed(:,18) = CrtBbecol;
    responseUnprocessed(:,19) = CrtYaecol;
    responseUnprocessed(:,20) = CrtYbecol;

    enzNames{17,:} = 'CrtBa';
    enzNames{18,:} = 'CrtBb';
    enzNames{19,:} = 'CrtYa';
    enzNames{20,:} = 'CrtYb';

end

if contains('detailed',modelType)

   enzNames{9,:} = 'ERG20';
   enzNames{12,:} = 'ERG9';
   enzNames{13,:} = 'CrtI';
   enzNames{14,:} = 'CrtYB';

end

response = reshape(reshape(responseUnprocessed(:,1:(end-2))',1,[]),(enzN-2)*metN,[])';
selectionVector = reshape(((repmat(1:1:(enzN-2),length(metsSelection),1)+repmat((enzN-2)*(metsSelection-1)',1,(enzN-2))))',1,[]);
response = response(:,selectionVector);
minResponse = sum(abs(min(response,[],1))>100);
maxResponse = sum(abs(max(response,[],1))>100);

enzN = enzN-2;
metN = length(metsSelection);
enzNames = enzNames(1:end-2);
metNames = metNames(metsSelection');

%% Adjust siginificance to the number of tests

% testsN = testsN+metN*enzN;
correctedSignificance = significance/testsN;


%% Plot heatmap of the mean for visualization

responseMean = reshape(mean(response,1),enzN,[])';

figure()
imagesc(responseMean)
set(gca,'xticklabel',[],'yticklabel',[],'xtick',1:enzN,'ytick',1:metN)
set(gca,'xticklabel',enzNames,'yticklabel',metNames)
colormap(cmap)
caxis([-max(abs(responseMean),[],"all"),max(abs(responseMean),[],"all")])

%% Plot histograms of values

% if ~exist(fullfile('..','io','figures','metabolite_control_significance_selection','histograms_values'), 'dir')
%    mkdir(fullfile('..','io','figures','metabolite_control_significance_selection','histograms_values'))
% end
% for ix = 1:size(response,2)
%     metPlot = fix((ix-1)/enzN)+1;
%     enzPlot = mod((ix-1),enzN)+1;
%     figure('visible','off')
%     histogram(response(:,ix),100,'Normalization','probability')
%     title(['Concentration Response Coefficient of ',metNames{metPlot,1},' by ',enzNames{enzPlot,1}])
%     saveas(gcf,fullfile('..','io','figures','metabolite_control_significance_selection','histograms_values',[samplingType,'_',modelType,'_',conditionName,'_',metNames{metPlot,1},'_',enzNames{enzPlot,1},'.png']))
% end

%% Perform Bootstrapping of the control coefficients

meanResponses = bootstrp(nbootstrap,@(x) mean(x),response);

% if ~exist(fullfile('..','io','figures','metabolite_control_significance_selection','histograms_bootstrap'), 'dir')
%    mkdir(fullfile('..','io','figures','metabolite_control_significance_selection','histograms_bootstrap'))
% end
% for ix = 1:size(response,2)
%     metPlot = fix((ix-1)/enzN)+1;
%     enzPlot = mod((ix-1),enzN)+1;
%     figure('visible','off')
%     histogram(meanResponses(:,ix),100,'Normalization','probability')
%     title(['Concentration Response Coefficient of ',metNames{metPlot,1},' by ',enzNames{enzPlot,1}])
%     saveas(gcf,fullfile('..','io','figures','metabolite_control_significance_selection','histograms_bootstrap',[samplingType,'_',modelType,'_',conditionName,'_',metNames{metPlot,1},'_',enzNames{enzPlot,1},'.png']))
% end

meanResponse = mean(meanResponses,1)';

meanResponseLCI = prctile(meanResponses,100*(correctedSignificance/2),1)';
meanResponseUCI = prctile(meanResponses,100*(1-correctedSignificance/2),1)';
resultsSignificance = (meanResponseUCI.*meanResponseLCI)>0;

%% Calculate p-value

[~,pValuesIndexes] = min(abs(sort(meanResponses,1)),[],1);
pValues = min(pValuesIndexes./(size(meanResponses,1)).*2,abs(1-pValuesIndexes./(size(meanResponses,1))).*2)';

%% Plot Heatmap of the P-values

pValuesMatrix = reshape(pValues,enzN,[])';
pValuesSignificance = zeros(metN,enzN);
pValuesSignificance(pValuesMatrix < 0.1/testsN) = 1;
pValuesSignificance(pValuesMatrix < 0.05/testsN) = 2;
pValuesSignificance(pValuesMatrix < 0.01/testsN) = 3;
pValuesSignificance(pValuesMatrix < 0.001/testsN) = 4;

cmap = repmat(((5:-1:1)/5)',1,3);

nx = enzN;
ny = metN;
mx = linspace(0.5,nx+0.5,nx+1);
my = linspace(0.5,ny+0.5,ny+1);

figure('Position',[500 300 250 110])
imagesc(pValuesSignificance); hold on
set(gca,'DefaultAxesTickLabelInterpreter','none')
set(gca,'TickLabelInterpreter','none')
set(gca,'xticklabel',[],'yticklabel',[],'xtick',1:enzN,'ytick',1:metN)
set(gca,'xticklabel',enzNames,'yticklabel',metNames)
set(gca,'units','points','position',[50,40,125,17.857])
set(gca,'TickDir','out')
set(gca,'LineWidth',0.5)
set(gca,'XTickLabelRotation',90)
colormap(cmap)
xlabel('Enzymes')
ylabel('Metabolites')
%title({'Escaled Significance of Mean Concentration Control Coefficient',[modelType,' ',conditionFullName]})
set(gca,'FontSize',7,'FontName','arial')
caxis([0 4])
colormap(cmap)
cb = colorbar;
set(cb,'Location','northoutside')
set(cb,'LineWidth',0.5)
set(cb,'TickDir','out')
set(cb,'FontSize',7,'FontName','arial')
set(cb,'Ticks',[0.4 1.2 2 2.8 3.6])
set(cb,'TickLabels',{'> 0.1';'0.05-0.1';'0.01-0.05';'0.001-0.01';'< 0.001'})
cbparameters = cb.Position;
set(cb,'Position',[cbparameters(1),0.75,cbparameters(3),0.05])
mg = mesh(mx,my,zeros([ny,nx]+1));
mg.FaceColor = 'none';
mg.EdgeColor = 'k';
hold off
saveas(gcf,fullfile('..','io','figures','metabolite_control_significance_selection','png',[samplingType,'_',modelType,'_',conditionName,'.png']))
saveas(gcf,fullfile('..','io','figures','metabolite_control_significance_selection','jpg',[samplingType,'_',modelType,'_',conditionName,'.jpg']))
saveas(gcf,fullfile('..','io','figures','metabolite_control_significance_selection','svg',[samplingType,'_',modelType,'_',conditionName,'.svg']))
saveas(gcf,fullfile('..','io','figures','metabolite_control_significance_selection','mat',[samplingType,'_',modelType,'_',conditionName,'.fig']))


%% Create vectors for output

originalSignificances = repmat(significance,metN*enzN,1);
testsNVector = repmat(testsN,metN*enzN,1);
correctedSignificances = repmat(correctedSignificance,metN*enzN,1);

metNamesColumn = reshape(repmat(metNames',enzN,1),[],1);
enzNamesColumn = repmat(enzNames,metN,1);

%% Create Table with the results

outputTable = table(metNamesColumn,enzNamesColumn,originalSignificances,testsNVector,correctedSignificances,meanResponseLCI,meanResponse,meanResponseUCI,resultsSignificance,pValues);
outputTable = renamevars(outputTable,'metNamesColumn','Concentration Control Coefficient');
outputTable = renamevars(outputTable,'enzNamesColumn','Enzyme');
outputTable = renamevars(outputTable,'originalSignificances','Original Significance');
outputTable = renamevars(outputTable,'testsNVector','N of Tests Correction');
outputTable = renamevars(outputTable,'correctedSignificances','Corrected Significance');
outputTable = renamevars(outputTable,'meanResponseLCI','Mean CCC Lower CI');
outputTable = renamevars(outputTable,'meanResponse','Mean CCC Estimated');
outputTable = renamevars(outputTable,'meanResponseUCI','Mean CCC Upper CI');
outputTable = renamevars(outputTable,'resultsSignificance','Statistical Significance');
outputTable = renamevars(outputTable,'pValues','Mean p-Values');

outputFile = fullfile('..','io','figures',"metabolite_control_significance_selection","xls","metabolite_control_significance.xlsx");
writetable(outputTable,outputFile,'Sheet',[samplingType,'_',modelType,'_',conditionName])





